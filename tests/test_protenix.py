"""Protenix: what the wrapper writes for `protenix pred`, reads back from it, and declares.

Nothing here runs the model. The upstream layout and flags are pinned from its source:
`runner/batch_inference.py` for the command line, `runner/dumper.py` for the files a prediction
leaves, `runner/inference.py` for what happens to a query that fails (a text file in `ERR/`, exit
status 0), and `docs/infer_json_format.md` for the input JSON.

Not covered, because it needs the model: the shape and token order of a real `token_pair_pae`,
and whether upstream accepts the entities this builds on a real run.
"""

import csv
import json
import os
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPTS = os.path.join(REPO_ROOT, "pipe_scripts")
sys.path.insert(0, SCRIPTS)

from biopipelines.protenix import MODEL_NAMES, Protenix  # noqa: E402
import pipe_protenix_config as config_script  # noqa: E402
import pipe_protenix_weights  # noqa: E402

POSTPROCESS = os.path.join(SCRIPTS, "pipe_protenix_postprocess.py")
CONFIG = os.path.join(SCRIPTS, "pipe_protenix_config.py")
PROTEIN = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ"


@pytest.fixture
def tool(local_config, isolated_cwd):
    """A Protenix with only the attributes validate_params and the section builders read."""
    defaults = dict(proteins=object(), ssDNA=None, dsDNA=None, ssRNA=None, dsRNA=None,
                    ligands=None, msas=None, use_msa_server=True, model_name="protenix-v2",
                    num_diffusion_samples=None, seeds=None, top_only=True, covalent_bonds=None,
                    include_pae=True, no_kernels=False)

    def _make(**kwargs):
        made = Protenix.__new__(Protenix)
        for name, value in {**defaults, **kwargs}.items():
            setattr(made, name, value)
        return made
    return _make


def _predict_section(made):
    made.queries_folder = "/q"
    made.prediction_folder = "/out"
    made.container_prefix = lambda: ""
    return Protenix._generate_predict_section(made)


class TestItRefusesWhatItCannotDo:
    """Each of these would otherwise run and return something quietly wrong."""

    def test_no_input_axis_is_refused(self, tool):
        with pytest.raises(ValueError, match="at least one input axis"):
            Protenix.validate_params(tool(proteins=None))

    def test_an_unknown_model_is_refused_with_the_accepted_names(self, tool):
        with pytest.raises(ValueError, match="protenix-v2"):
            Protenix.validate_params(tool(model_name="protenix-v3"))

    @pytest.mark.parametrize("name", ["protenix_base_constraint_v0.5.0", "protenix_mini_esm_v0.5.0"])
    def test_the_checkpoints_that_need_more_than_sequences_are_not_offered(self, tool, name):
        """The constraint model needs a constraint input and the ESM ones precompute embeddings."""
        assert name not in MODEL_NAMES
        with pytest.raises(ValueError, match="model_name"):
            Protenix.validate_params(tool(model_name=name))

    def test_every_offered_model_is_accepted(self, tool):
        for name in MODEL_NAMES:
            Protenix.validate_params(tool(model_name=name))

    def test_supplied_alignments_and_the_server_are_mutually_exclusive(self, tool):
        """Upstream searches every protein chain of a query once any chain lacks a path, and overwrites the rest."""
        with pytest.raises(ValueError, match="use_msa_server must be False"):
            Protenix.validate_params(tool(msas=object(), use_msa_server=True))

    def test_supplied_alignments_with_the_server_off_are_accepted(self, tool):
        Protenix.validate_params(tool(msas=object(), use_msa_server=False))

    def test_an_alignment_stream_that_is_not_a3m_is_refused_at_configuration(self, tool):
        made = tool(msas=object(), use_msa_server=False)
        made.msas_stream = type("S", (), {"format": "csv"})()
        with pytest.raises(ValueError, match="a3m"):
            Protenix.validate_params(made)

    @pytest.mark.parametrize("seeds", [[], [1.5], [True], "101"])
    def test_bad_seeds_are_refused(self, tool, seeds):
        with pytest.raises(ValueError, match="seeds"):
            Protenix.validate_params(tool(seeds=seeds))

    @pytest.mark.parametrize("samples", [0, -1, 2.5, True])
    def test_bad_sample_counts_are_refused(self, tool, samples):
        with pytest.raises(ValueError, match="num_diffusion_samples"):
            Protenix.validate_params(tool(num_diffusion_samples=samples))


class TestCovalentBondValidation:

    BOND = {"entity1": 1, "position1": 5, "atom1": "SG", "entity2": 2, "position2": 1, "atom2": "C1"}

    def test_a_complete_bond_is_accepted(self, tool):
        Protenix.validate_params(tool(covalent_bonds=[self.BOND]))

    def test_upstream_documents_ids_as_strings_and_they_are_accepted(self, tool):
        bond = {**self.BOND, "entity1": "1", "position1": "5", "copy1": 1, "copy2": "1"}
        Protenix.validate_params(tool(covalent_bonds=[bond]))

    @pytest.mark.parametrize("missing", ["entity1", "position2", "atom1"])
    def test_a_missing_key_is_named(self, tool, missing):
        bond = {k: v for k, v in self.BOND.items() if k != missing}
        with pytest.raises(ValueError, match=missing):
            Protenix.validate_params(tool(covalent_bonds=[bond]))

    def test_an_unknown_key_is_refused(self, tool):
        """The deprecated left_/right_ spellings would run, upstream says it may drop them."""
        with pytest.raises(ValueError, match="left_entity"):
            Protenix.validate_params(tool(covalent_bonds=[{**self.BOND, "left_entity": 1}]))

    def test_one_copy_without_the_other_is_refused(self, tool):
        with pytest.raises(ValueError, match="copy1 and copy2"):
            Protenix.validate_params(tool(covalent_bonds=[{**self.BOND, "copy1": 1}]))

    @pytest.mark.parametrize("value", [0, -2, "x", True])
    def test_entities_and_positions_count_from_one(self, tool, value):
        with pytest.raises(ValueError, match="entity1"):
            Protenix.validate_params(tool(covalent_bonds=[{**self.BOND, "entity1": value}]))

    def test_an_empty_list_is_refused(self, tool):
        with pytest.raises(ValueError, match="non-empty"):
            Protenix.validate_params(tool(covalent_bonds=[]))


class TestThePredictCommand:
    """Pinned because a wrong flag is not a soft failure: click rejects it or the setting is ignored."""

    def test_the_command_and_the_flags_upstream_spells(self, tool):
        text = _predict_section(tool())
        for fragment in ("protenix pred ", '-i "/q"', '-o "/out"', "-n protenix-v2",
                         "--use_default_params true", "--use_msa true",
                         "--need_atom_confidence true"):
            assert fragment in text, fragment

    def test_the_default_model_is_v2(self, tool):
        assert "-n protenix-v2" in _predict_section(tool())
        import inspect
        assert inspect.signature(Protenix.__init__).parameters["model_name"].default == "protenix-v2"

    def test_another_model_is_passed_by_name(self, tool):
        assert "-n protenix_mini_default_v0.5.0" in _predict_section(
            tool(model_name="protenix_mini_default_v0.5.0"))

    def test_default_params_are_always_requested(self, tool):
        """Without it the CLI's own 10 cycles and 200 steps override the 4 and 5 the mini and tiny models need."""
        assert "--use_default_params true" in _predict_section(tool(model_name="protenix_tiny_default_v0.5.0"))

    def test_seeds_and_samples_are_passed_only_when_set(self, tool):
        plain = _predict_section(tool())
        assert " -s " not in plain and " -e " not in plain
        text = _predict_section(tool(seeds=[1, 2, 3], num_diffusion_samples=4))
        assert " -s 1,2,3" in text and " -e 4" in text

    def test_the_msa_service_is_switched_off_explicitly(self, tool):
        """Upstream defaults --use_msa to true, so omitting the flag does not disable it."""
        assert "--use_msa false" in _predict_section(tool(use_msa_server=False))

    def test_supplied_alignments_keep_msa_features_on_without_a_server_flag(self, tool):
        text = _predict_section(tool(msas=object(), use_msa_server=False))
        assert "--use_msa true" in text
        assert "msa_server" not in text

    def test_pae_asks_upstream_for_the_full_data_file_and_nothing_else_does(self, tool):
        assert "--need_atom_confidence true" in _predict_section(tool(include_pae=True))
        assert "need_atom_confidence" not in _predict_section(tool(include_pae=False))

    def test_no_kernels_selects_the_torch_implementations(self, tool):
        text = _predict_section(tool(no_kernels=True))
        assert "--trimul_kernel torch --triatt_kernel torch" in text
        assert "kernel" not in _predict_section(tool(no_kernels=False))

    def test_a_container_prefix_is_applied_to_the_command_only(self, tool):
        made = tool()
        made.queries_folder, made.prediction_folder = "/q", "/out"
        made.container_prefix = lambda: "apptainer exec x.sif "
        assert "apptainer exec x.sif protenix pred" in Protenix._generate_predict_section(made)


class TestTheDeclaredOutputs:

    @staticmethod
    def _declare(made, monkeypatch, ids=("q",)):
        import biopipelines.protenix as module
        monkeypatch.setattr(module, "predict_output_ids_with_provenance", lambda **k: (list(ids), {}))
        made._axis_kwargs = lambda: {}
        made.stream_path = lambda name, f: f"/x/{name}/{f}"
        made.structures_map_csv = made.pae_map_csv = made.confidence_csv = made.missing_csv = "/x/m.csv"
        made.output_folder = "/x"
        made.missing_table_info = lambda path: None
        return made.get_output_files()

    def test_top_only_declares_the_bare_ids(self, tool, monkeypatch):
        out = self._declare(tool(), monkeypatch, ids=("a", "b"))
        assert out["structures"].ids == ["a", "b"]
        assert out["structures"].format == "cif"

    def test_the_default_declares_upstreams_five_samples_for_one_seed(self, tool, monkeypatch):
        out = self._declare(tool(top_only=False), monkeypatch)
        assert out["structures"].ids == ["q_<1..5>"]

    def test_every_seed_is_counted_in_the_declared_sample_ids(self, tool, monkeypatch):
        """seeds=[100, 101] at 3 samples writes six structures; declaring three wired downstream to half."""
        out = self._declare(tool(top_only=False, seeds=[100, 101], num_diffusion_samples=3), monkeypatch)
        assert out["structures"].ids == ["q_<1..6>"]

    def test_pae_mirrors_the_structures_ids_and_is_an_npz(self, tool, monkeypatch):
        out = self._declare(tool(top_only=False, num_diffusion_samples=2), monkeypatch)
        assert out["pae"].ids == out["structures"].ids == ["q_<1..2>"]
        assert out["pae"].format == "npz"
        assert out["pae"].files[0].endswith("pae/<id>.npz")

    def test_pae_is_not_declared_when_it_is_not_written(self, tool, monkeypatch):
        """Declaring a file nothing writes would fail a run that worked."""
        out = self._declare(tool(include_pae=False), monkeypatch)
        assert len(out["pae"]) == 0

    def test_the_confidence_table_has_the_documented_columns(self, tool, monkeypatch):
        out = self._declare(tool(), monkeypatch)
        assert out["tables"]["confidence"].info.columns == [
            "id", "plddt", "ptm", "iptm", "gpde", "has_clash", "disorder", "ranking_score",
            "seed", "sample"]


def _axis(entity_type, records, mode="each"):
    return {"entity_type": entity_type, "mode": mode,
            "iterated": records if mode == "each" else [],
            "static": [] if mode == "each" else records, "static_first": False}


class TestTheQueryJson:

    def test_each_axis_gives_one_query_per_pair_with_the_framework_ids(self):
        elements = {
            "proteins": _axis("protein", [{"id": "p1", "sequence": PROTEIN}, {"id": "p2", "sequence": PROTEIN}]),
            "ligands": _axis("ligand", [{"id": "l1", "smiles": "CCO"}]),
        }
        jobs, skipped = config_script.build_queries(elements, None, None)
        assert sorted(jobs) == ["p1+l1", "p2+l1"] and not skipped
        job = jobs["p1+l1"]
        assert job["name"] == "p1+l1"
        assert job["sequences"] == [{"proteinChain": {"sequence": PROTEIN, "count": 1}},
                                    {"ligand": {"ligand": "CCO", "count": 1}}]

    def test_entities_follow_the_axis_order_whatever_iterates(self):
        """covalent_bonds names entities by position, so the position cannot depend on Each versus Bundle."""
        elements = {
            "ligands": _axis("ligand", [{"id": "l1", "smiles": "CCO"}]),
            "proteins": _axis("protein", [{"id": "p1", "sequence": PROTEIN}, {"id": "p2", "sequence": PROTEIN}],
                              mode="bundle"),
            "dsDNA": _axis("dsdna", [{"id": "d1", "sequence": "ACGTAC"}]),
        }
        jobs, _ = config_script.build_queries(elements, None, None)
        (job,) = jobs.values()
        assert [next(iter(e)) for e in job["sequences"]] == [
            "proteinChain", "proteinChain", "dnaSequence", "dnaSequence", "ligand"]

    def test_a_double_stranded_axis_adds_the_reverse_complement_strand(self):
        elements = {"dsDNA": _axis("dsdna", [{"id": "d1", "sequence": "AAAG"}])}
        (job,) = config_script.build_queries(elements, None, None)[0].values()
        assert [e["dnaSequence"]["sequence"] for e in job["sequences"]] == ["AAAG", "CTTT"]

    def test_rna_uses_the_rna_entity_and_single_strands_stay_single(self):
        elements = {"ssRNA": _axis("ssrna", [{"id": "r1", "sequence": "GGCAUC"}])}
        (job,) = config_script.build_queries(elements, None, None)[0].values()
        assert job["sequences"] == [{"rnaSequence": {"sequence": "GGCAUC", "count": 1}}]

    def test_a_ccd_ligand_is_prefixed_and_a_smiles_ligand_is_not(self):
        ccd = config_script.build_entity(
            "ligand", {"id": "l1", "ccd": "NAG", "smiles": "CC", "source": "rcsb"}, None)
        assert ccd == {"ligand": {"ligand": "CCD_NAG", "count": 1}}
        smiles = config_script.build_entity("ligand", {"id": "l2", "smiles": "CCO"}, None)
        assert smiles == {"ligand": {"ligand": "CCO", "count": 1}}

    def test_an_ion_is_a_ccd_ligand(self):
        assert config_script.build_entity("ligand", {"id": "i", "ccd": "MG"}, None) == {
            "ligand": {"ligand": "CCD_MG", "count": 1}}

    @pytest.mark.parametrize("record, kind, fragment", [
        ({"id": "p_empty", "sequence": float("nan")}, "protein", "p_empty"),
        ({"id": "p_u", "sequence": "MKTU"}, "protein", "U"),
        ({"id": "d_u", "sequence": "ACGU"}, "ssdna", "U"),
        ({"id": "r_t", "sequence": "ACGT"}, "ssrna", "T"),
        ({"id": "d_low", "sequence": "acgt"}, "ssdna", "a"),
    ])
    def test_a_sequence_upstream_would_reject_is_named_by_its_id(self, record, kind, fragment):
        with pytest.raises(config_script.QueryError, match=fragment):
            config_script.build_entity(kind, record, None)

    def test_a_bad_query_is_skipped_and_the_rest_are_still_written(self):
        elements = {"proteins": _axis("protein", [{"id": "good", "sequence": PROTEIN},
                                                  {"id": "bad", "sequence": "MKT*"}])}
        jobs, skipped = config_script.build_queries(elements, None, None)
        assert list(jobs) == ["good"]
        assert "bad" in skipped and "*" in skipped["bad"]

    def test_a_ligand_without_chemistry_skips_its_query(self):
        elements = {"proteins": _axis("protein", [{"id": "p", "sequence": PROTEIN}]),
                    "ligands": _axis("ligand", [{"id": "hetatm", "code": "LIG"}])}
        jobs, skipped = config_script.build_queries(elements, None, None)
        assert not jobs and "p+hetatm" in skipped

    def test_an_iterating_axis_with_no_records_builds_no_query_instead_of_crashing(self):
        assert config_script.build_queries({"proteins": _axis("protein", [])}, None, None) == ({}, {})


class TestAlignments:

    @staticmethod
    def _elements():
        return {"proteins": _axis("protein", [{"id": "p1", "sequence": PROTEIN},
                                              {"id": "p2", "sequence": PROTEIN}])}

    def test_a_supplied_alignment_becomes_an_absolute_unpaired_path(self, tmp_path):
        a3m = tmp_path / "p1.a3m"
        a3m.write_text(">q\n" + PROTEIN + "\n")
        jobs, skipped = config_script.build_queries(self._elements(), {"p1": str(a3m)}, None)
        assert jobs["p1"]["sequences"][0]["proteinChain"]["unpairedMsaPath"] == str(a3m)
        assert "pairedMsaPath" not in jobs["p1"]["sequences"][0]["proteinChain"]

    def test_a_chain_without_an_alignment_skips_its_query_rather_than_reaching_the_server(self, tmp_path):
        """A path that is absent makes upstream search for every chain, the opposite of supplying them."""
        a3m = tmp_path / "p1.a3m"
        a3m.write_text(">q\n" + PROTEIN + "\n")
        jobs, skipped = config_script.build_queries(self._elements(), {"p1": str(a3m)}, None)
        assert list(jobs) == ["p1"]
        assert "no alignment" in skipped["p2"]

    def test_an_alignment_file_that_does_not_exist_skips_its_query(self, tmp_path):
        jobs, skipped = config_script.build_queries(
            self._elements(), {"p1": str(tmp_path / "gone.a3m"), "p2": str(tmp_path / "gone2.a3m")}, None)
        assert not jobs and "does not exist" in skipped["p1"]

    def test_without_supplied_alignments_no_path_is_written(self):
        jobs, _ = config_script.build_queries(self._elements(), None, None)
        assert "unpairedMsaPath" not in jobs["p1"]["sequences"][0]["proteinChain"]

    def test_an_alignment_the_model_cannot_read_is_named_not_ignored(self, tmp_path):
        table = tmp_path / "msas_map.csv"
        pd.DataFrame([{"id": "p1", "file": str(tmp_path / "p1.csv")}]).to_csv(table, index=False)
        stream = tmp_path / "msas.json"
        stream.write_text(json.dumps({"map_table": str(table)}))
        with pytest.raises(ValueError, match="a3m"):
            config_script.load_msa_lookup(str(stream))

    def test_the_lookup_maps_ids_to_absolute_paths(self, tmp_path):
        table = tmp_path / "msas_map.csv"
        pd.DataFrame([{"id": "p1", "file": "rel/p1.a3m"}]).to_csv(table, index=False)
        stream = tmp_path / "msas.json"
        stream.write_text(json.dumps({"map_table": str(table)}))
        assert config_script.load_msa_lookup(str(stream)) == {"p1": os.path.abspath("rel/p1.a3m")}


class TestCovalentBondsInTheQuery:

    BOND = {"entity1": 1, "position1": 5, "atom1": "SG", "entity2": 2, "position2": 1, "atom2": "C1"}

    @classmethod
    def _elements(cls):
        return {"proteins": _axis("protein", [{"id": "p", "sequence": PROTEIN}]),
                "ligands": _axis("ligand", [{"id": "l", "ccd": "NAG"}])}

    def test_the_bonds_are_copied_into_every_query_verbatim(self):
        jobs, skipped = config_script.build_queries(self._elements(), None, [self.BOND])
        assert jobs["p+l"]["covalent_bonds"] == [self.BOND] and not skipped

    def test_no_bonds_leaves_the_key_out(self):
        jobs, _ = config_script.build_queries(self._elements(), None, None)
        assert "covalent_bonds" not in jobs["p+l"]

    def test_a_bond_naming_a_third_entity_skips_the_query(self):
        """Upstream would index past the end of its sequences list."""
        jobs, skipped = config_script.build_queries(
            self._elements(), None, [{**self.BOND, "entity2": 3}])
        assert not jobs and "3" in skipped["p+l"] and "2 entities" in skipped["p+l"]


class TestTheConfigScriptEndToEnd:
    """The script as the generated bash calls it, against a combinatorics config on disk."""

    @staticmethod
    def _run(tmp_path, proteins, extra=()):
        from biopipelines.combinatorics import AxisConfig, CombinatoricsConfig
        table = tmp_path / "proteins.csv"
        pd.DataFrame(proteins).to_csv(table, index=False)
        config = CombinatoricsConfig(axes={"proteins": AxisConfig(
            name="proteins", mode="each", sources=[{"path": str(table), "iterate": True}],
            entity_type="protein")})
        config.save(str(tmp_path / "combinatorics.json"))
        out = tmp_path / "cfg"
        done = subprocess.run(
            [sys.executable, CONFIG, "--combinatorics-config", str(tmp_path / "combinatorics.json"),
             "--queries-dir", str(out / "queries"), "--query-ids-csv", str(out / "ids.csv"),
             "--skipped-csv", str(out / "skipped.csv"), *extra],
            capture_output=True, text=True)
        return out, done

    def test_one_json_per_query_and_an_id_list_that_includes_the_skipped(self, tmp_path):
        out, done = self._run(tmp_path, [{"id": "a", "sequence": PROTEIN}, {"id": "b", "sequence": "MK*"}])
        assert done.returncode == 0, done.stdout + done.stderr
        assert sorted(os.listdir(out / "queries")) == ["a.json"]
        payload = json.loads((out / "queries" / "a.json").read_text())
        assert isinstance(payload, list) and payload[0]["name"] == "a"
        assert pd.read_csv(out / "ids.csv")["id"].tolist() == ["a", "b"]
        assert pd.read_csv(out / "skipped.csv")["id"].tolist() == ["b"]

    def test_a_stale_queries_folder_is_cleared(self, tmp_path):
        stale = tmp_path / "cfg" / "queries"
        stale.mkdir(parents=True)
        (stale / "old.json").write_text("[]")
        (stale / "a-update-msa.json").write_text("[]")
        out, done = self._run(tmp_path, [{"id": "a", "sequence": PROTEIN}])
        assert done.returncode == 0, done.stdout + done.stderr
        assert os.listdir(out / "queries") == ["a.json"]

    def test_nothing_left_to_predict_is_a_step_error(self, tmp_path):
        out, done = self._run(tmp_path, [{"id": "b", "sequence": "MK*"}])
        assert done.returncode == 1
        assert "every query was skipped" in done.stdout
        assert pd.read_csv(out / "skipped.csv")["id"].tolist() == ["b"]


# -- postprocess -----------------------------------------------------------------

def _sample(pred, query, seed, sample, ranking, pae=True, n_token=3, layout="current"):
    """One sample the way runner/dumper.py leaves it: cif, summary JSON and (optionally) full-data JSON."""
    if layout == "current":
        folder = pred / query / f"seed_{seed}" / "predictions"
        prefix = query
    else:
        folder = pred / query / str(seed)
        prefix = f"{query}_{seed}"
    folder.mkdir(parents=True, exist_ok=True)
    (folder / f"{prefix}_sample_{sample}.cif").write_text(f"data_{query}_{seed}_{sample}\n")
    (folder / f"{prefix}_summary_confidence_sample_{sample}.json").write_text(json.dumps({
        "plddt": 80.0 + sample, "gpde": 1.5, "ptm": 0.8, "iptm": 0.7, "has_clash": False,
        "disorder": 0.0, "ranking_score": ranking, "chain_ptm": [0.8]}))
    if pae is True:
        matrix = np.full((n_token, n_token), seed * 10 + sample, dtype=float).tolist()
        (folder / f"{prefix}_full_data_sample_{sample}.json").write_text(json.dumps({
            "atom_plddt": [0.9, 0.9], "token_pair_pae": matrix, "token_pair_pde": matrix}))
    elif pae is not False:
        (folder / f"{prefix}_full_data_sample_{sample}.json").write_text(json.dumps(pae))


def _postprocess(tmp_path, queries, *flags, skipped=None, with_pae=True):
    out = tmp_path / "out"
    (tmp_path / "ids.csv").write_text("id\n" + "".join(f"{q}\n" for q in queries))
    args = [sys.executable, POSTPROCESS, "--prediction-folder", str(tmp_path / "pred"),
            "--query-ids-csv", str(tmp_path / "ids.csv"),
            "--structures-dir", str(out / "structures"),
            "--structures-map-csv", str(out / "structures" / "m.csv"),
            "--confidence-csv", str(out / "c.csv"), "--local-missing-csv", str(out / "x.csv")]
    if skipped is not None:
        (tmp_path / "skipped.csv").write_text("id,cause\n" + "".join(f"{k},{v}\n" for k, v in skipped.items()))
        args += ["--skipped-csv", str(tmp_path / "skipped.csv")]
    if with_pae:
        args += ["--pae-dir", str(out / "pae"), "--pae-map-csv", str(out / "pae" / "pae_map.csv")]
    return out, subprocess.run([*args, *flags], capture_output=True, text=True)


def _matrix(path):
    with np.load(path) as data:
        assert list(data.files) == ["pae"], "the array key is the one boltz writes"
        return data["pae"]


class TestThePostprocess:

    def test_top_only_keeps_the_best_structure_across_seeds_not_the_first(self, tmp_path):
        """Upstream ranks per seed, so sample 0 of the first seed is only the best of that seed."""
        pred = tmp_path / "pred"
        _sample(pred, "q", 101, 0, ranking=0.5)
        _sample(pred, "q", 101, 1, ranking=0.1)
        _sample(pred, "q", 102, 0, ranking=0.9)
        out, done = _postprocess(tmp_path, ["q"])
        assert done.returncode == 0, done.stdout + done.stderr
        assert (out / "structures" / "q.cif").read_text() == "data_q_102_0\n"
        row = pd.read_csv(out / "c.csv").iloc[0]
        assert (row["id"], row["seed"], row["sample"], row["ranking_score"]) == ("q", 102, 0, 0.9)
        assert row["plddt"] == 80.0 and row["ptm"] == 0.8 and row["iptm"] == 0.7 and row["gpde"] == 1.5

    def test_every_structure_is_numbered_by_seed_then_sample_when_all_are_kept(self, tmp_path):
        pred = tmp_path / "pred"
        for seed in (102, 101):
            for sample in (1, 0):
                _sample(pred, "q", seed, sample, ranking=0.5 - 0.1 * sample)
        out, done = _postprocess(tmp_path, ["q"], "--top-only", "False")
        assert done.returncode == 0, done.stdout + done.stderr
        table = pd.read_csv(out / "c.csv")
        assert table["id"].tolist() == ["q_1", "q_2", "q_3", "q_4"]
        assert list(zip(table["seed"], table["sample"])) == [(101, 0), (101, 1), (102, 0), (102, 1)]
        assert (out / "structures" / "q_3.cif").read_text() == "data_q_102_0\n"
        assert pd.read_csv(out / "structures" / "m.csv")["id"].tolist() == table["id"].tolist()

    def test_the_pae_belongs_to_the_structure_on_disk_and_is_an_npz_with_the_boltz_key(self, tmp_path):
        pred = tmp_path / "pred"
        _sample(pred, "q", 101, 0, ranking=0.2, n_token=4)
        _sample(pred, "q", 101, 1, ranking=0.9, n_token=4)
        out, done = _postprocess(tmp_path, ["q"])
        assert done.returncode == 0, done.stdout + done.stderr
        matrix = _matrix(out / "pae" / "q.npz")
        assert matrix.shape == (4, 4) and matrix.dtype == np.float32
        assert float(matrix[0, 0]) == 1011.0, "the matrix of sample 1, the best one, not of sample 0"
        assert (out / "structures" / "q.cif").read_text() == "data_q_101_1\n"

    def test_pae_ids_mirror_the_structures_ids_when_all_are_kept(self, tmp_path):
        pred = tmp_path / "pred"
        _sample(pred, "q", 101, 0, ranking=0.9)
        _sample(pred, "q", 101, 1, ranking=0.1)
        out, done = _postprocess(tmp_path, ["q"], "--top-only", "False")
        assert done.returncode == 0, done.stdout + done.stderr
        structures = pd.read_csv(out / "structures" / "m.csv")
        pae = pd.read_csv(out / "pae" / "pae_map.csv")
        assert structures["id"].tolist() == pae["id"].tolist() == ["q_1", "q_2"]
        assert all(os.path.exists(f) for f in pae["file"])

    def test_the_bulky_full_data_files_are_deleted_once_the_matrix_is_out(self, tmp_path):
        pred = tmp_path / "pred"
        _sample(pred, "q", 101, 0, ranking=0.9)
        _sample(pred, "q", 101, 1, ranking=0.1)
        out, done = _postprocess(tmp_path, ["q"])
        assert done.returncode == 0, done.stdout + done.stderr
        assert not list(pred.rglob("*_full_data_sample_*.json"))
        assert list(pred.rglob("*_summary_confidence_sample_*.json")), "the scores stay"

    def test_without_the_pae_flags_nothing_extra_is_written_or_deleted(self, tmp_path):
        pred = tmp_path / "pred"
        _sample(pred, "q", 101, 0, ranking=0.9, pae=False)
        out, done = _postprocess(tmp_path, ["q"], with_pae=False)
        assert done.returncode == 0, done.stdout + done.stderr
        assert not (out / "pae").exists()
        assert sorted(p.name for p in out.iterdir()) == ["c.csv", "structures", "x.csv"]

    def test_a_structure_with_no_usable_pae_is_kept_and_only_left_out_of_the_pae_map(self, tmp_path):
        pred = tmp_path / "pred"
        _sample(pred, "ok", 101, 0, ranking=0.5)
        _sample(pred, "absent", 101, 0, ranking=0.5, pae=False)
        _sample(pred, "ragged", 101, 0, ranking=0.5, pae={"token_pair_pae": [[1.0, 2.0, 3.0], [1.0, 2.0, 3.0]]})
        _sample(pred, "nokey", 101, 0, ranking=0.5, pae={"atom_plddt": [0.5]})
        out, done = _postprocess(tmp_path, ["ok", "absent", "ragged", "nokey"])
        assert done.returncode == 0, done.stdout + done.stderr
        assert pd.read_csv(out / "pae" / "pae_map.csv")["id"].tolist() == ["ok"]
        assert pd.read_csv(out / "structures" / "m.csv")["id"].tolist() == ["ok", "absent", "ragged", "nokey"]
        assert done.stdout.count("WARNING: no PAE for") == 3
        assert pd.read_csv(out / "x.csv").empty, "a missing matrix is not a missing structure"

    def test_a_leading_singleton_dimension_is_squeezed(self, tmp_path):
        pred = tmp_path / "pred"
        _sample(pred, "q", 101, 0, ranking=0.5, pae={"token_pair_pae": [[[0.5, 1.0], [1.0, 0.5]]]})
        out, done = _postprocess(tmp_path, ["q"])
        assert done.returncode == 0, done.stdout + done.stderr
        assert _matrix(out / "pae" / "q.npz").shape == (2, 2)

    def test_a_failed_query_keeps_its_row_and_carries_upstreams_error_as_the_cause(self, tmp_path):
        """Upstream writes ERR/<query>.txt and exits 0, so the files are the only evidence."""
        pred = tmp_path / "pred"
        _sample(pred, "ok", 101, 0, ranking=0.5)
        (pred / "ERR").mkdir(parents=True)
        (pred / "ERR" / "big.txt").write_text(
            "[Rank 0] big failed: protenix-v2 model does not support n_token > 2560.\nTraceback ...\n")
        out, done = _postprocess(tmp_path, ["ok", "big", "silent"])
        assert done.returncode == 0, done.stdout + done.stderr
        assert pd.read_csv(out / "c.csv")["id"].tolist() == ["ok", "big", "silent"]
        missing = pd.read_csv(out / "x.csv")
        assert missing["id"].tolist() == ["big", "silent"]
        assert set(missing["kind"]) == {"failure"} and set(missing["removed_by"]) == {"Protenix"}
        assert "n_token > 2560" in missing.loc[0, "cause"]
        assert "no structure" in missing.loc[1, "cause"]
        assert pd.read_csv(out / "structures" / "m.csv")["id"].tolist() == ["ok"]

    def test_a_skipped_query_is_a_failure_row_with_the_reason_it_was_skipped(self, tmp_path):
        pred = tmp_path / "pred"
        _sample(pred, "ok", 101, 0, ranking=0.5)
        out, done = _postprocess(tmp_path, ["ok", "bad"], skipped={"bad": "bad: letter U not allowed"})
        assert done.returncode == 0, done.stdout + done.stderr
        missing = pd.read_csv(out / "x.csv")
        assert missing["id"].tolist() == ["bad"] and "letter U" in missing.loc[0, "cause"]
        assert pd.read_csv(out / "c.csv")["id"].tolist() == ["ok", "bad"]

    def test_a_missing_summary_leaves_nan_scores_but_keeps_the_structure(self, tmp_path):
        pred = tmp_path / "pred"
        _sample(pred, "q", 101, 0, ranking=0.5)
        os.remove(pred / "q" / "seed_101" / "predictions" / "q_summary_confidence_sample_0.json")
        out, done = _postprocess(tmp_path, ["q"])
        assert done.returncode == 0, done.stdout + done.stderr
        row = pd.read_csv(out / "c.csv").iloc[0]
        assert row["id"] == "q" and pd.isna(row["plddt"]) and pd.isna(row["ranking_score"])
        assert (out / "structures" / "q.cif").exists()

    def test_a_query_whose_name_starts_another_does_not_take_its_files(self, tmp_path):
        pred = tmp_path / "pred"
        _sample(pred, "p1", 101, 0, ranking=0.5)
        _sample(pred, "p1_2", 101, 0, ranking=0.5)
        out, done = _postprocess(tmp_path, ["p1", "p1_2"], "--top-only", "False")
        assert done.returncode == 0, done.stdout + done.stderr
        assert pd.read_csv(out / "c.csv")["id"].tolist() == ["p1_1", "p1_2_1"]

    def test_ids_with_a_plus_sign_are_found(self, tmp_path):
        """The framework's combined ids look like `p1+l1`."""
        pred = tmp_path / "pred"
        _sample(pred, "p1+l1", 101, 0, ranking=0.5)
        out, done = _postprocess(tmp_path, ["p1+l1"])
        assert done.returncode == 0, done.stdout + done.stderr
        assert (out / "structures" / "p1+l1.cif").exists()

    def test_the_layout_upstreams_docs_describe_is_read_too(self, tmp_path):
        pred = tmp_path / "pred"
        _sample(pred, "q", 7, 0, ranking=0.4, layout="docs")
        out, done = _postprocess(tmp_path, ["q"])
        assert done.returncode == 0, done.stdout + done.stderr
        row = pd.read_csv(out / "c.csv").iloc[0]
        assert (row["seed"], row["sample"]) == (7, 0)
        assert (out / "pae" / "q.npz").exists()

    def test_nothing_produced_is_a_step_error(self, tmp_path):
        out, done = _postprocess(tmp_path, ["a", "b"])
        assert done.returncode == 1
        assert "no Protenix prediction produced a structure" in done.stdout
        assert pd.read_csv(out / "x.csv")["id"].tolist() == ["a", "b"]

    def test_an_empty_id_list_is_a_step_error(self, tmp_path):
        out, done = _postprocess(tmp_path, [])
        assert done.returncode == 1


# -- install ---------------------------------------------------------------------

class TestTheInstall:

    @pytest.fixture
    def script(self, local_config, monkeypatch):
        monkeypatch.setattr(Protenix, "_install_env", classmethod(lambda cls, mgr: "protenix"))
        return lambda **folders: Protenix._install_script(
            {"biopipelines": "/repo", **folders}, env_manager="mamba")

    def test_the_shortcut_for_an_existing_env_still_ensures_the_weights(self, script):
        shortcut = script().split("already installed, skipping")[1].split("exit 0")[0]
        assert "pipe_protenix_weights.py" in shortcut and "--model-name" in shortcut

    def test_it_verifies_the_cli_not_just_an_import(self, script):
        assert "protenix --version" in script()

    def test_the_fetch_script_runs_inside_the_env(self, script):
        fetch = os.path.join("/repo", "pipe_scripts", "pipe_protenix_weights.py")
        assert f'mamba run -n protenix python "{fetch}"' in script()

    def test_the_cache_folder_becomes_the_root_upstream_writes_to(self, script):
        text = script(ProtenixCache="/cache/Protenix")
        assert 'export PROTENIX_ROOT_DIR="/cache/Protenix"' in text
        assert text.index("PROTENIX_ROOT_DIR") < text.index("pipe_protenix_weights.py")

    def test_without_a_cache_folder_the_upstream_default_applies(self, script):
        assert "PROTENIX_ROOT_DIR" not in script()

    def test_the_default_checkpoint_is_the_one_fetched(self, script):
        assert '--model-name "protenix-v2"' in script()

    def test_a_checkpoint_upstream_refuses_warns_but_does_not_fail_the_install(self, script):
        """Status 3 from the fetch script; any other non-zero status still stops it."""
        text = script()
        assert 'if [ "$FETCH_RC" -eq 3 ]; then' in text
        assert "WARNING: Protenix is installed without the protenix-v2 checkpoint" in text
        assert 'elif [ "$FETCH_RC" -ne 0 ]; then' in text
        assert 'ERROR: the CCD files could not be downloaded"\n    exit 1' in text

    def test_the_status_the_script_exits_with_is_the_one_the_install_tests_for(self):
        from biopipelines.protenix import CHECKPOINT_UNAVAILABLE
        assert CHECKPOINT_UNAVAILABLE == pipe_protenix_weights.CHECKPOINT_UNAVAILABLE == 3

    def test_the_run_script_points_upstream_at_the_same_cache(self, tool):
        made = tool()
        made.folders = {"ProtenixCache": "/cache/Protenix"}
        assert Protenix._cache_root(made) == "/cache/Protenix"
        made.folders = {}
        assert Protenix._cache_root(made) is None


class TestTheWeightsScript:

    COMMON = pipe_protenix_weights.COMMON_FILES

    @classmethod
    def _targets(cls, tmp_path):
        targets = {name: str(tmp_path / "common" / name) for name in cls.COMMON}
        targets["protenix-v2"] = str(tmp_path / "checkpoint" / "protenix-v2.pt")
        return targets, {name: f"https://example.invalid/{name}" for name in targets}

    @staticmethod
    def _download(refuse=(), calls=None):
        """upstream's download_from_url, writing the file it was asked for, or refusing like an HTTP 403."""
        import urllib.error

        def download(url, path, check_weight):
            if calls is not None:
                calls.append(os.path.basename(url))
            with open(path, "w") as handle:
                handle.write("partial")
            if os.path.basename(url) in refuse:
                raise urllib.error.HTTPError(url, 403, "Forbidden", {}, None)
        return download

    def test_the_targets_are_the_ccd_files_and_the_named_checkpoint(self):
        targets = pipe_protenix_weights.targets_for(
            "protenix-v2",
            {name: f"/root/common/{name}" for name in self.COMMON},
            {"load_checkpoint_dir": "/root/checkpoint"})
        assert set(targets) == {*self.COMMON, "protenix-v2"}
        assert targets["protenix-v2"] == os.path.join("/root/checkpoint", "protenix-v2.pt")

    def test_everything_absent_is_downloaded_to_its_final_path_with_no_partial_left(self, tmp_path):
        targets, urls = self._targets(tmp_path)
        assert pipe_protenix_weights.fetch_missing("protenix-v2", targets, urls, self._download()) is None
        assert all(os.path.exists(path) for path in targets.values())
        assert not list(tmp_path.rglob("*.part"))

    def test_what_is_already_there_is_not_downloaded_again(self, tmp_path):
        targets, urls = self._targets(tmp_path)
        for path in targets.values():
            os.makedirs(os.path.dirname(path), exist_ok=True)
            open(path, "w").close()
        calls = []
        assert pipe_protenix_weights.fetch_missing(
            "protenix-v2", targets, urls, self._download(calls=calls)) is None
        assert calls == []

    def test_a_refused_checkpoint_is_reported_and_leaves_the_ccd_files_in_place(self, tmp_path):
        """protenix-v2 answered HTTP 403 where the other checkpoints downloaded: the install must carry on."""
        targets, urls = self._targets(tmp_path)
        reason = pipe_protenix_weights.fetch_missing(
            "protenix-v2", targets, urls, self._download(refuse=("protenix-v2",)))
        assert reason and "403" in reason and "protenix-v2" in reason
        assert not os.path.exists(targets["protenix-v2"])
        assert not list(tmp_path.rglob("*.part")), "a partial file would be mistaken for the checkpoint next time"
        assert all(os.path.exists(targets[name]) for name in self.COMMON)

    def test_a_ccd_file_that_cannot_be_fetched_is_an_error(self, tmp_path):
        targets, urls = self._targets(tmp_path)
        with pytest.raises(RuntimeError, match="ccd_components_file"):
            pipe_protenix_weights.fetch_missing(
                "protenix-v2", targets, urls, self._download(refuse=("ccd_components_file",)))

    def test_a_target_upstream_lists_no_url_for_is_an_error(self, tmp_path):
        targets, urls = self._targets(tmp_path)
        del urls["pdb_cluster_file"]
        with pytest.raises(RuntimeError, match="pdb_cluster_file"):
            pipe_protenix_weights.fetch_missing("protenix-v2", targets, urls, self._download())


# -- through a pipeline ---------------------------------------------------------

class TestThroughAPipeline:

    @staticmethod
    def _built(new_pipeline, **kwargs):
        """(the Protenix instance, the output its constructor returned, the generated step script)."""
        from pathlib import Path
        from biopipelines.sequence import Sequence
        pipeline = new_pipeline("protenix")
        with pipeline:
            output = Protenix(proteins=Sequence(seq=PROTEIN, type="protein", ids="p1"), **kwargs)
            script_path = pipeline.save()
        steps = sorted(Path(script_path).parent.glob("*_Protenix.sh"))
        script = "\n".join(step.read_text(encoding="utf-8", errors="replace") for step in steps)
        return pipeline.tools[-1], output, script

    def test_the_generated_step_runs_config_then_prediction_then_collection(
            self, local_config, isolated_cwd, new_pipeline):
        tool, output, script = self._built(new_pipeline)
        config = script.index("pipe_protenix_config.py")
        predict = script.index("protenix pred")
        collect = script.index("pipe_protenix_postprocess.py")
        assert config < predict < collect
        assert f'--queries-dir "{tool.queries_folder}"' in script
        assert f'-i "{tool.queries_folder}"' in script
        assert f'--pae-map-csv "{output.streams.pae.map_table}"' in script

    def test_the_streams_and_tables_are_declared(self, local_config, isolated_cwd, new_pipeline):
        _tool, output, _script = self._built(new_pipeline)
        assert output.streams.structures.ids == output.streams.pae.ids == ["p1"]
        assert output.streams.pae.format == "npz"
        assert output.tables.confidence.info.columns[0] == "id"
        assert output.tables.missing is not None

    def test_pae_off_leaves_the_flags_and_the_stream_out(self, local_config, isolated_cwd, new_pipeline):
        _tool, output, script = self._built(new_pipeline, include_pae=False)
        assert len(output.streams.pae) == 0
        assert "--pae-dir" not in script and "need_atom_confidence" not in script

    def test_covalent_bonds_are_written_at_configuration_and_handed_to_the_config_script(
            self, local_config, isolated_cwd, new_pipeline):
        bond = {"entity1": 1, "position1": 5, "atom1": "SG", "entity2": 2, "position2": 1, "atom2": "C1"}
        tool, _output, script = self._built(new_pipeline, covalent_bonds=[bond])
        assert json.loads(open(tool.covalent_bonds_file).read()) == [bond]
        assert f'--covalent-bonds-json "{tool.covalent_bonds_file}"' in script

    def test_a_misspelled_parameter_is_reported_rather_than_forwarded(
            self, local_config, isolated_cwd, new_pipeline, capsys):
        """Protenix forwards no unknown keys, so a typo must be flagged and must not reach the command."""
        _tool, _output, script = self._built(new_pipeline, num_samples=3)
        assert "unknown constructor parameter(s): 'num_samples'" in capsys.readouterr().err
        assert "num_samples" not in script and " -e " not in script

    def test_it_advertises_no_affinity_parameter(self):
        """There is no affinity head; a parameter named for one would be accepted and ignored."""
        import inspect
        assert "affinity" not in set(inspect.signature(Protenix.__init__).parameters)
