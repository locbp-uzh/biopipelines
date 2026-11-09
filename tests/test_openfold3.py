"""OpenFold3: the second co-folder, and what it deliberately is not.

The tool exists so a prediction can be cross-checked against a model that is not Boltz2, which
only works if the two take the same inputs. So most of what is worth pinning here is agreement
with Boltz2 — the same axes, the same ligand chemistry decision, the same double-stranded
expansion — plus the two places they must NOT agree, because OpenFold3 has no affinity head and
no covalent linkage and pretending otherwise would hand a user a silent wrong answer.

The upstream flag spellings are pinned too. The docs' auto-generated option table renders them
with dashes while both worked examples use underscores; a wrong flag is not a soft failure,
`run_openfold` either rejects the run or ignores the setting.
"""

import json
import os
import sys

import pytest

from biopipelines.openfold3 import OpenFold3
from tests._helpers import sys_path_prepended

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
with sys_path_prepended(os.path.join(REPO_ROOT, "pipe_scripts")):
    import pipe_openfold3_config

pytestmark = pytest.mark.usefixtures("pipe_scripts_on_path")


@pytest.fixture
def tool(local_config, isolated_cwd):
    def _make(**kwargs):
        made = OpenFold3.__new__(OpenFold3)
        made.proteins = kwargs.pop("proteins", object())
        for axis in ("ssDNA", "dsDNA", "ssRNA", "dsRNA", "ligands", "msas"):
            setattr(made, axis, kwargs.pop(axis, None))
        made.use_msa_server = kwargs.pop("use_msa_server", True)
        made.output_format = kwargs.pop("output_format", "cif")
        made.top_only = kwargs.pop("top_only", True)
        made.template = kwargs.pop("template", None)
        made.template_chain_ids = kwargs.pop("template_chain_ids", None)
        made.low_mem = kwargs.pop("low_mem", False)
        for name in ("num_diffusion_samples", "num_model_seeds", "seeds",
                     "inference_ckpt_name", "inference_ckpt_path", "devices", "runner_yaml"):
            setattr(made, name, kwargs.pop(name, None))
        assert not kwargs, f"unexpected: {sorted(kwargs)}"
        return made
    return _make


class TestItRefusesWhatItCannotDo:
    """Each of these would otherwise run and return something quietly wrong."""

    def test_precomputed_msas_and_the_server_are_mutually_exclusive(self, tool):
        with pytest.raises(ValueError, match="use_msa_server must be False"):
            OpenFold3.validate_params(tool(msas=object(), use_msa_server=True))

    def test_precomputed_msas_with_the_server_off_are_accepted(self, tool):
        OpenFold3.validate_params(tool(msas=object(), use_msa_server=False))

    def test_contradicting_seed_specifications_are_refused(self, tool):
        """The count says how many to draw, the list says which. Picking one silently would make
        a run that was asked to be reproducible by seed reproducible by a different seed."""
        with pytest.raises(ValueError, match="contradicts"):
            OpenFold3.validate_params(tool(seeds=[1, 2, 3], num_model_seeds=5))

    def test_agreeing_seed_specifications_are_accepted(self, tool):
        OpenFold3.validate_params(tool(seeds=[1, 2, 3], num_model_seeds=3))

    def test_two_checkpoints_are_refused(self, tool):
        with pytest.raises(ValueError, match="not both"):
            OpenFold3.validate_params(tool(inference_ckpt_name="a", inference_ckpt_path="/b.pt"))

    def test_no_input_axis_is_refused(self, tool):
        with pytest.raises(ValueError, match="at least one input axis"):
            OpenFold3.validate_params(tool(proteins=None))

    def test_an_unknown_output_format_is_refused(self, tool):
        with pytest.raises(ValueError, match="output_format"):
            OpenFold3.validate_params(tool(output_format="mmtf"))

    def test_template_chains_without_a_template_are_refused(self, tool):
        """Silently ignoring it would leave the user believing a template was applied."""
        with pytest.raises(ValueError, match="template_chain_ids"):
            OpenFold3.validate_params(tool(template_chain_ids=["A"]))


class TestTheFlagSpellings:
    """Pinned because the docs disagree with themselves and a wrong flag is not a soft failure."""

    def test_the_predict_command_uses_underscore_flags(self, tool):
        made = tool()
        made.queries_json = "/q.json"
        made.prediction_folder = "/out"
        made.runner_yaml_file = "/r.yaml"
        made.extra_args_bash = lambda: ""
        made.extra_args_echo = lambda: ""
        made.container_prefix = lambda: ""
        text = OpenFold3._generate_predict_section(made)
        for flag in ("--query_json=", "--output_dir=", "--runner_yaml=", "--use_msa_server="):
            assert flag in text, f"{flag} is how run_openfold spells it"
        assert "--query-json" not in text

    def test_the_msa_server_is_switched_off_explicitly(self, tool):
        """Upstream defaults it to True, so omitting the flag does not disable it."""
        made = tool(use_msa_server=False, msas=object())
        made.queries_json = "/q.json"
        made.prediction_folder = "/out"
        made.runner_yaml_file = "/r.yaml"
        made.extra_args_bash = lambda: ""
        made.extra_args_echo = lambda: ""
        made.container_prefix = lambda: ""
        assert "--use_msa_server=False" in OpenFold3._generate_predict_section(made)


class TestTheQueryJson:

    def test_a_double_stranded_axis_becomes_two_chains(self):
        from pipe_openfold3_config import build_chain
        forward = build_chain("dsdna", {"id": "d1", "sequence": "ACGT"}, "A", None, None, None)
        reverse = build_chain("dsdna", {"id": "d1", "sequence": "ACGT"}, "B", None, None, None,
                              rev_comp=True)
        assert forward["sequence"] == "ACGT"
        assert reverse["sequence"] == "ACGT", "reverse complement of ACGT is itself"
        second = build_chain("dsdna", {"id": "d2", "sequence": "AAAG"}, "B", None, None, None,
                             rev_comp=True)
        assert second["sequence"] == "CTTT"

    def test_dna_and_rna_map_onto_the_upstream_molecule_types(self):
        from pipe_openfold3_config import build_chain
        assert build_chain("ssdna", {"id": "d", "sequence": "ACGT"}, "A", None, None,
                           None)["molecule_type"] == "dna"
        assert build_chain("ssrna", {"id": "r", "sequence": "ACGU"}, "A", None, None,
                           None)["molecule_type"] == "rna"

    def test_a_ccd_ligand_uses_ccd_codes_and_a_smiles_ligand_uses_smiles(self):
        from pipe_openfold3_config import build_chain
        ccd = build_chain("ligand", {"id": "l1", "ccd": "NAG", "smiles": "CC", "source": "rcsb"},
                          "A", None, None, None)
        assert ccd["ccd_codes"] == ["NAG"] and "smiles" not in ccd
        smiles = build_chain("ligand", {"id": "l2", "smiles": "CCO"}, "A", None, None, None)
        assert smiles["smiles"] == "CCO" and "ccd_codes" not in smiles

    def test_an_empty_polymer_is_refused_by_the_id_that_caused_it(self):
        """An empty FASTA arrives as float nan; unchecked it becomes a query with no residues
        and the failure surfaces inside the model with nothing naming the input."""
        from pipe_openfold3_config import build_chain
        with pytest.raises(ValueError, match="p_bad"):
            build_chain("protein", {"id": "p_bad", "sequence": float("nan")}, "A", None, None, None)

    def test_chain_labels_continue_past_z(self):
        from pipe_openfold3_config import chain_ids
        assert chain_ids(1, 0) == ["A"]
        assert chain_ids(1, 25) == ["Z"]
        assert chain_ids(1, 26) == ["AA"]

    def test_a_supplied_a3m_is_staged_under_a_name_upstream_reads(self, tmp_path):
        """Upstream parses an a3m only when its stem is a database it names (colabfold_main, uniref90_hits, ...): `p1.a3m` was skipped and the chain failed with an IndexError."""
        source = tmp_path / "p1.a3m"
        source.write_text(">p1\nMKV\n")
        staged = pipe_openfold3_config.stage_alignment(str(source), "p1", str(tmp_path / "msas"))
        assert os.path.basename(staged) == "colabfold_main.a3m"
        assert open(staged).read() == ">p1\nMKV\n"

    def test_each_chain_gets_a_folder_of_its_own(self, tmp_path):
        staged = {}
        for name in ("p1", "p2"):
            (tmp_path / f"{name}.a3m").write_text(f">{name}\nMKV\n")
            staged[name] = pipe_openfold3_config.stage_alignment(
                str(tmp_path / f"{name}.a3m"), name, str(tmp_path / "msas"))
        assert staged["p1"] != staged["p2"]
        assert open(staged["p2"]).read() == ">p2\nMKV\n"

    def test_a_preparsed_npz_is_passed_through_unchanged(self, tmp_path):
        source = tmp_path / "p1.npz"
        source.write_bytes(b"x")
        assert pipe_openfold3_config.stage_alignment(str(source), "p1", str(tmp_path / "msas")) == str(source)

    def test_an_msa_the_model_cannot_read_is_named_not_ignored(self, tmp_path):
        """Our msas streams are often CSV (Boltz's format); OpenFold3 reads a3m/sto/npz. Passing
        one through would fold every chain single-sequence and report nothing amiss."""
        with pytest.raises(ValueError, match="a3m"):
            pipe_openfold3_config.stage_alignment(str(tmp_path / "p1.csv"), "p1", str(tmp_path / "msas"))


class TestTheDeclaredOutputs:

    def test_top_only_off_declares_a_sample_suffix(self, tool, monkeypatch):
        """The postprocess writes <id>_1..N then; declaring the bare id would make the
        completion check demand files that were never meant to exist."""
        made = tool(top_only=False, num_diffusion_samples=3)
        made._axis_kwargs = lambda: {}
        made.stream_path = lambda *a: "/x/<id>.cif"
        made.stream_map_path = lambda n: "/x/map.csv"
        made.structures_map_csv = "/x/map.csv"
        made.confidence_csv = "/x/c.csv"
        made.missing_csv = "/x/m.csv"
        made.output_folder = "/x"
        made.missing_table_info = lambda p: None
        import biopipelines.openfold3 as module
        monkeypatch.setattr(module, "predict_output_ids_with_provenance", lambda **k: (["q1"], {}))
        out = OpenFold3.get_output_files(made)
        assert out["structures"].ids == ["q1_<1..3>"]


def test_it_advertises_no_affinity_and_no_covalent_parameters():
    """OpenFold3 has neither. A parameter named for one would be accepted and ignored, which is
    how a user ends up believing they ran a covalent prediction."""
    import inspect
    names = set(inspect.signature(OpenFold3.__init__).parameters)
    assert "affinity" not in names
    assert "covalent_linkage" not in names


class TestFindingsFromThe151Review:

    def test_every_seed_is_counted_in_the_declared_sample_ids(self, tool, monkeypatch):
        """seeds=[100, 101] with 3 samples writes _1.._6; declaring _1.._3 wired downstream to half."""
        import biopipelines.openfold3 as module
        monkeypatch.setattr(module, "predict_output_ids_with_provenance", lambda **k: (["q"], {}))
        made = tool(top_only=False, seeds=[100, 101], num_diffusion_samples=3)
        made._axis_kwargs = lambda: {}
        made.stream_path = lambda name, f: f
        made.structures_map_csv = made.confidence_csv = made.missing_csv = "x.csv"
        made.output_folder = "out"
        made.missing_table_info = lambda path: None
        assert made.get_output_files()["structures"].ids == ["q_<1..6>"]

    def test_a_sample_count_hidden_in_runner_yaml_is_refused_for_per_sample_ids(self, tool):
        with pytest.raises(ValueError, match="needs num_diffusion_samples"):
            OpenFold3.validate_params(tool(top_only=False, runner_yaml="r.yaml"))

    def test_a_csv_msa_stream_is_refused_at_configuration(self, tool):
        made = tool(msas=object(), use_msa_server=False)
        made.msas_stream = type("S", (), {"format": "csv"})()
        with pytest.raises(ValueError, match="a3m"):
            OpenFold3.validate_params(made)

    def test_the_install_shortcut_still_ensures_the_weights(self, local_config, monkeypatch):
        monkeypatch.setattr(OpenFold3, "_install_env", classmethod(lambda cls, mgr: "OpenFold3Env"))
        script = OpenFold3._install_script({"biopipelines": "/repo"}, env_manager="mamba")
        shortcut = script.split("already installed, skipping")[1].split("exit 0")[0]
        assert "setup_openfold --non-interactive" in shortcut


def test_a_query_with_no_folder_is_reported_missing(tmp_path):
    """run_openfold makes no folder for a query that failed before inference; listing folders dropped it."""
    import subprocess
    pred = tmp_path / "pred"
    (pred / "ok" / "seed_1").mkdir(parents=True)
    (pred / "ok" / "seed_1" / "ok_seed_1_sample_1_model.cif").write_text("data_x\n")
    queries = tmp_path / "queries.json"
    queries.write_text(json.dumps({"queries": {"ok": {}, "lost": {}}}))
    out = tmp_path / "out"
    script = os.path.join(REPO_ROOT, "pipe_scripts", "pipe_openfold3_postprocess.py")
    done = subprocess.run([sys.executable, script, "--prediction-folder", str(pred),
                           "--combinatorics-config", str(tmp_path / "none.json"),
                           "--queries-json", str(queries),
                           "--structures-dir", str(out / "s"), "--structures-map-csv", str(out / "m.csv"),
                           "--confidence-csv", str(out / "c.csv"), "--local-missing-csv", str(out / "x.csv")],
                          capture_output=True, text=True)
    assert done.returncode == 0, done.stderr
    assert "lost" in (out / "x.csv").read_text()


class TestEveryDeclaredSampleSlotIsNamed:
    """With top_only off the declared ids are <query>_1..N, so a slot upstream never filled
    has to appear as a failure: the bare query id names none of them, and a silent gap
    leaves the confidence table short of the rows the contract promises."""

    POSTPROCESS = os.path.join(REPO_ROOT, "pipe_scripts", "pipe_openfold3_postprocess.py")

    @staticmethod
    def _sample(pred, query, sample):
        folder = pred / query / "seed_1"
        folder.mkdir(parents=True, exist_ok=True)
        (folder / f"{query}_seed_1_sample_{sample}_model.cif").write_text(f"data_{query}_{sample}\n")

    def _run(self, tmp_path, queries, per_query):
        import subprocess
        pred = tmp_path / "pred"
        pred.mkdir(exist_ok=True)
        queries_json = tmp_path / "queries.json"
        queries_json.write_text(json.dumps({"queries": {q: {} for q in queries}}))
        out = tmp_path / "out"
        done = subprocess.run(
            [sys.executable, self.POSTPROCESS, "--prediction-folder", str(pred),
             "--combinatorics-config", str(tmp_path / "none.json"),
             "--queries-json", str(queries_json),
             "--structures-dir", str(out / "s"), "--structures-map-csv", str(out / "m.csv"),
             "--confidence-csv", str(out / "c.csv"), "--local-missing-csv", str(out / "x.csv"),
             "--top-only", "False", "--structures-per-query", str(per_query)],
            capture_output=True, text=True)
        return out, done

    def test_a_query_that_produced_fewer_samples_than_declared_names_each_id_it_owes(self, tmp_path):
        import pandas as pd
        self._sample(tmp_path / "pred", "q", 1)
        out, done = self._run(tmp_path, ["q"], per_query=3)
        assert done.returncode == 0, done.stdout + done.stderr
        assert pd.read_csv(out / "m.csv")["id"].tolist() == ["q_1"]
        missing = pd.read_csv(out / "x.csv")
        assert missing["id"].tolist() == ["q_2", "q_3"]
        assert set(missing["kind"]) == {"failure"}
        assert "1 of 3" in missing.loc[0, "cause"]
        assert pd.read_csv(out / "c.csv")["id"].tolist() == ["q_1", "q_2", "q_3"]

    def test_a_query_that_produced_nothing_names_every_declared_sample_id(self, tmp_path):
        import pandas as pd
        out, done = self._run(tmp_path, ["q"], per_query=2)
        assert done.returncode == 1
        missing = pd.read_csv(out / "x.csv")
        assert missing["id"].tolist() == ["q_1", "q_2"]
        assert set(missing["kind"]) == {"failure"}

    def test_the_declared_count_and_the_flag_come_from_one_helper(self, local_config, isolated_cwd):
        """The ids the wrapper declares and the number the postprocess is told must not drift apart."""
        from biopipelines.mock import Mock
        from biopipelines.pipeline import Pipeline
        with Pipeline(project="TestSuite", job="of3", description="x", on_the_fly=False,
                      local_output=True, config="local") as pipeline:
            seqs = Mock(ids=["p"], streams={"sequences": {"format": "csv", "file": "<id>.csv"}})
            OpenFold3(proteins=seqs.streams.sequences, top_only=False,
                      num_diffusion_samples=3, seeds=[1, 2])
            tool = pipeline.tools[-1]
        assert tool._structures_per_query() == 6
        assert "--structures-per-query 6" in tool._generate_postprocess_section()


def test_an_iterating_axis_with_no_records_builds_no_query_instead_of_crashing():
    elements = {"proteins": {"mode": "each", "entity_type": "protein", "iterated": [],
                             "static": [], "static_first": False}}
    assert pipe_openfold3_config.build_queries(elements, None, None, None) == ({"queries": {}}, {})



def test_msas_alone_turn_the_server_off(local_config, isolated_cwd):
    """The docstring promised it; the default use_msa_server=True made msas= raise instead."""
    from biopipelines.mock import Mock
    from biopipelines.pipeline import Pipeline
    with Pipeline(project="TestSuite", job="of3", description="x", on_the_fly=False,
                  local_output=True, config="local") as pipeline:
        seqs = Mock(ids=["p"], streams={"sequences": {"format": "csv", "file": "<id>.csv"}})
        msas = Mock(ids=["p"], streams={"msas": {"format": "a3m", "file": "<id>.a3m"}})
        OpenFold3(proteins=seqs.streams.sequences, msas=msas.streams.msas)
        assert pipeline.tools[-1].use_msa_server is False
        with pytest.raises(ValueError, match="use_msa_server must be False"):
            OpenFold3(proteins=seqs.streams.sequences, msas=msas.streams.msas, use_msa_server=True)


class TestThePaeStream:
    """The full confidence JSON, copied per kept structure, for a downstream reader of `pae`.

    Layout pinned from upstream's writer: `<prefix>_confidences.json` holds `plddt` (per atom),
    `pde` and `pae` (token x token). The real files and their token order are not covered here.
    """

    @staticmethod
    def _sample(pred, query, sample, ranking, with_json=True):
        """One sample's structure, aggregated scores and (optionally) full confidence JSON."""
        folder = pred / query / "seed_1"
        folder.mkdir(parents=True, exist_ok=True)
        prefix = f"{query}_seed_1_sample_{sample}"
        (folder / f"{prefix}_model.cif").write_text(f"data_{query}_{sample}\n")
        (folder / f"{prefix}_confidences_aggregated.json").write_text(
            json.dumps({"sample_ranking_score": ranking}))
        if with_json:
            (folder / f"{prefix}_confidences.json").write_text(json.dumps({
                "plddt": [90.0, 91.0], "pde": [[0.0, 1.0], [1.0, 0.0]],
                "pae": [[0.5, float(sample)], [float(sample), 0.5]]}, indent=4))

    @staticmethod
    def _collect(tmp_path, queries, *flags):
        import subprocess
        out = tmp_path / "out"
        (tmp_path / "queries.json").write_text(json.dumps({"queries": {q: {} for q in queries}}))
        script = os.path.join(REPO_ROOT, "pipe_scripts", "pipe_openfold3_postprocess.py")
        return out, subprocess.run(
            [sys.executable, script, "--prediction-folder", str(tmp_path / "pred"),
             "--combinatorics-config", str(tmp_path / "none.json"),
             "--queries-json", str(tmp_path / "queries.json"),
             "--structures-dir", str(out / "structures"),
             "--structures-map-csv", str(out / "structures" / "m.csv"),
             "--confidence-csv", str(out / "c.csv"), "--local-missing-csv", str(out / "x.csv"),
             *flags], capture_output=True, text=True)

    @staticmethod
    def _pae_flags(tmp_path):
        return ["--pae-dir", str(tmp_path / "out" / "pae"),
                "--pae-map-csv", str(tmp_path / "out" / "pae" / "pae_map.csv")]

    def test_the_kept_sample_json_is_copied_unchanged_under_the_structure_id(self, tmp_path):
        """The best sample by ranking score is kept, so its JSON, not the first sample's, is the one copied."""
        pred = tmp_path / "pred"
        self._sample(pred, "q1", 1, ranking=0.2)
        self._sample(pred, "q1", 2, ranking=0.9)
        out, done = self._collect(tmp_path, ["q1"], *self._pae_flags(tmp_path))
        assert done.returncode == 0, done.stderr
        source = pred / "q1" / "seed_1" / "q1_seed_1_sample_2_confidences.json"
        assert (out / "pae" / "q1.json").read_bytes() == source.read_bytes()
        assert json.loads((out / "pae" / "q1.json").read_text())["pae"] == [[0.5, 2.0], [2.0, 0.5]]
        assert (out / "structures" / "q1.cif").read_text() == "data_q1_2\n"

    def test_every_sample_is_copied_when_all_are_kept(self, tmp_path):
        pred = tmp_path / "pred"
        self._sample(pred, "q1", 1, ranking=0.2)
        self._sample(pred, "q1", 2, ranking=0.9)
        out, done = self._collect(tmp_path, ["q1"], "--top-only", "False", *self._pae_flags(tmp_path))
        assert done.returncode == 0, done.stderr
        structures = sorted(p.stem for p in (out / "structures").glob("*.cif"))
        assert structures == sorted(p.stem for p in (out / "pae").glob("*.json")) == ["q1_1", "q1_2"]

    def test_the_map_lists_only_ids_that_have_a_json(self, tmp_path):
        import pandas as pd
        pred = tmp_path / "pred"
        self._sample(pred, "q1", 1, ranking=0.5)
        self._sample(pred, "q2", 1, ranking=0.5, with_json=False)
        out, done = self._collect(tmp_path, ["q1", "q2"], *self._pae_flags(tmp_path))
        assert done.returncode == 0, done.stderr
        assert pd.read_csv(out / "pae" / "pae_map.csv")["id"].tolist() == ["q1"]
        assert "q2" in done.stdout
        assert pd.read_csv(out / "structures" / "m.csv")["id"].tolist() == ["q1", "q2"]

    def test_without_the_flags_nothing_extra_is_written(self, tmp_path):
        """A run that does not ask for the stream (a runner_yaml that disabled the file) is unchanged."""
        self._sample(tmp_path / "pred", "q1", 1, ranking=0.5)
        out, done = self._collect(tmp_path, ["q1"])
        assert done.returncode == 0, done.stderr
        assert not (out / "pae").exists()
        assert sorted(p.name for p in out.iterdir()) == ["c.csv", "structures", "x.csv"]

    @staticmethod
    def _declared(new_pipeline, tmp_path, runner=None, **kwargs):
        from pathlib import Path
        from biopipelines.sequence import Sequence
        if runner is not None:
            (tmp_path / "runner.yaml").write_text(runner)
            kwargs["runner_yaml"] = str(tmp_path / "runner.yaml")
        pipeline = new_pipeline("of3_pae")
        with pipeline:
            tool = OpenFold3(proteins=Sequence(seq="MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ",
                                               type="protein", ids="p1"), **kwargs)
            script_path = pipeline.save()
        steps = sorted(Path(script_path).parent.glob("*_OpenFold3.sh"))
        return tool, "\n".join(step.read_text(encoding="utf-8", errors="replace") for step in steps)

    def test_the_stream_is_declared_by_default_with_the_structures_ids(
            self, local_config, isolated_cwd, new_pipeline):
        tool, script = self._declared(new_pipeline, isolated_cwd)
        assert tool.streams.pae.ids == tool.streams.structures.ids == ["p1"]
        assert tool.streams.pae.format == "json"
        assert tool.streams.pae.files[0].endswith(os.path.join("pae", "<id>.json"))
        assert tool.streams.pae.map_table.endswith(os.path.join("pae", "pae_map.csv"))
        assert f'--pae-dir "{os.path.dirname(tool.streams.pae.map_table)}"' in script
        assert f'--pae-map-csv "{tool.streams.pae.map_table}"' in script

    def test_every_sample_is_declared_when_all_are_kept(self, local_config, isolated_cwd, new_pipeline):
        tool, _script = self._declared(new_pipeline, isolated_cwd, top_only=False,
                                       num_diffusion_samples=2, num_model_seeds=3)
        assert tool.streams.pae.ids == tool.streams.structures.ids == ["p1_<1..6>"]

    @pytest.mark.parametrize("runner", [
        "output_writer_settings:\n  write_full_confidence_scores: False\n",
        "output_writer_settings:\n  full_confidence_output_format: npz\n",
    ])
    def test_a_runner_yaml_that_removes_the_json_removes_the_stream(
            self, local_config, isolated_cwd, new_pipeline, runner):
        """Declaring a file run_openfold will not write would fail a run that worked."""
        tool, script = self._declared(new_pipeline, isolated_cwd, runner=runner)
        assert len(tool.streams.pae) == 0
        assert "--pae-dir" not in script

    def test_a_runner_yaml_that_leaves_the_json_alone_keeps_the_stream(
            self, local_config, isolated_cwd, new_pipeline):
        tool, _script = self._declared(
            new_pipeline, isolated_cwd,
            runner="output_writer_settings:\n  write_full_confidence_scores: True\n"
                   "  full_confidence_output_format: json\n")
        assert tool.streams.pae.ids == ["p1"]


class TestFindingsFromThe160Review:

    @staticmethod
    def _declare(made, monkeypatch, runner_text=None, tmp_path=None):
        import biopipelines.openfold3 as module
        monkeypatch.setattr(module, "predict_output_ids_with_provenance", lambda **k: (["q"], {}))
        if runner_text is not None:
            (tmp_path / "runner.yaml").write_text(runner_text)
            made.runner_yaml = str(tmp_path / "runner.yaml")
        made._axis_kwargs = lambda: {}
        made.stream_path = lambda name, f: f
        made.structures_map_csv = made.confidence_csv = made.missing_csv = "x.csv"
        made.output_folder = "out"
        made.missing_table_info = lambda path: None
        return made.get_output_files()

    def test_the_seeds_of_runner_yaml_are_counted_in_the_declared_sample_ids(self, tool, monkeypatch, tmp_path):
        """The runner YAML's seeds replace `seeds`, so three of them at two samples write six structures per query."""
        made = tool(top_only=False, num_diffusion_samples=2)
        out = self._declare(made, monkeypatch, "experiment_settings:\n  seeds: [1, 2, 3]\n", tmp_path)
        assert out["structures"].ids == ["q_<1..6>"]

    def test_the_seeds_of_runner_yaml_win_over_the_seeds_argument(self, tool, monkeypatch, tmp_path):
        made = tool(top_only=False, num_diffusion_samples=2, seeds=[9])
        out = self._declare(made, monkeypatch, "experiment_settings:\n  seeds: [1, 2, 3]\n", tmp_path)
        assert out["structures"].ids == ["q_<1..6>"]

    def test_a_runner_yaml_without_seeds_leaves_the_count_to_the_arguments(self, tool, monkeypatch, tmp_path):
        made = tool(top_only=False, num_diffusion_samples=2, num_model_seeds=3)
        out = self._declare(made, monkeypatch, "model_update:\n  presets: [predict]\n", tmp_path)
        assert out["structures"].ids == ["q_<1..6>"]

    @pytest.mark.parametrize("seeds", ["3", "[]"])
    def test_seeds_in_runner_yaml_that_are_not_a_list_are_refused_for_per_sample_ids(self, tool, tmp_path, seeds):
        (tmp_path / "runner.yaml").write_text(f"experiment_settings:\n  seeds: {seeds}\n")
        made = tool(top_only=False, num_diffusion_samples=2, runner_yaml=str(tmp_path / "runner.yaml"))
        with pytest.raises(ValueError, match="non-empty list"):
            OpenFold3.validate_params(made)

    def test_a_seed_count_that_contradicts_the_seeds_of_runner_yaml_is_refused(self, tool, tmp_path):
        (tmp_path / "runner.yaml").write_text("experiment_settings:\n  seeds: [1, 2]\n")
        made = tool(top_only=False, num_diffusion_samples=2, num_model_seeds=3,
                    runner_yaml=str(tmp_path / "runner.yaml"))
        with pytest.raises(ValueError, match="contradicts"):
            OpenFold3.validate_params(made)

    def test_the_seeds_of_runner_yaml_do_not_matter_while_only_the_best_sample_is_kept(self, tool, tmp_path):
        (tmp_path / "runner.yaml").write_text("experiment_settings:\n  seeds: 3\n")
        OpenFold3.validate_params(tool(top_only=True, runner_yaml=str(tmp_path / "runner.yaml")))

    def test_an_alignment_the_model_cannot_read_skips_only_its_query(self, tmp_path):
        import biopipelines.id_map_utils as id_map_utils
        (tmp_path / "p1.a3m").write_text(">p1\nMKT\n")
        (tmp_path / "p2.csv").write_text("key,sequence\n-1,MKT\n")
        msas = id_map_utils.MsaLookup([{"id": "p1", "file": str(tmp_path / "p1.a3m")},
                                       {"id": "p2", "file": str(tmp_path / "p2.csv")}])
        elements = {"proteins": {"entity_type": "protein", "mode": "each", "static_first": False, "static": [],
                                 "iterated": [{"id": "p1", "sequence": "MKT"}, {"id": "p2", "sequence": "MAA"}]}}
        payload, skipped = pipe_openfold3_config.build_queries(
            elements, pipe_openfold3_config.Alignments(msas, str(tmp_path / "msas")), None, None)
        assert list(payload["queries"]) == ["p1"]
        assert list(skipped) == ["p2"] and "a3m" in skipped["p2"]

    def test_the_sample_ids_are_composed_by_the_shared_helper(self, tool, monkeypatch):
        import biopipelines.openfold3 as module
        calls = []
        monkeypatch.setattr(module, "generate_multiplied_ids_pattern",
                            lambda ids, suffix: calls.append((list(ids), suffix)) or [f"{i}#{suffix}" for i in ids])
        out = self._declare(tool(top_only=False, num_diffusion_samples=3), monkeypatch)
        assert calls == [(["q"], "<1..3>")] and out["structures"].ids == ["q#<1..3>"]

    def test_a_failure_row_names_its_step_so_a_later_step_of_the_same_tool_excuses_it(self, tmp_path):
        import importlib.util
        import subprocess
        spec = importlib.util.spec_from_file_location(
            "pipe_check_completion", os.path.join(REPO_ROOT, "pipe_scripts", "pipe_check_completion.py"))
        checker = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(checker)

        (tmp_path / "pred").mkdir()
        queries = tmp_path / "queries.json"
        queries.write_text(json.dumps({"queries": {"lost": {}}}))
        out = tmp_path / "out"
        missing = out / "001_OpenFold3" / "execution" / "local_missing.csv"
        script = os.path.join(REPO_ROOT, "pipe_scripts", "pipe_openfold3_postprocess.py")
        done = subprocess.run([sys.executable, script, "--prediction-folder", str(tmp_path / "pred"),
                               "--combinatorics-config", str(tmp_path / "none.json"),
                               "--queries-json", str(queries), "--structures-dir", str(out / "s"),
                               "--structures-map-csv", str(out / "m.csv"), "--confidence-csv", str(out / "c.csv"),
                               "--local-missing-csv", str(missing)], capture_output=True, text=True)
        assert done.returncode == 1 and missing.exists(), done.stdout + done.stderr
        expected = {"tables": {"missing": {"path": str(missing)}}}
        assert checker._load_expected_missing_ids(expected, "001_OpenFold3", "OpenFold3") == []
        assert checker._load_expected_missing_ids(expected, "003_OpenFold3", "OpenFold3") == ["lost"]


class TestProvenanceColumns:
    """`<axis>.id` columns on the maps and the confidence table, read by exact query id and never by splitting it."""

    PROTEIN = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ"

    @staticmethod
    def _axis(entity_type, records):
        return {"entity_type": entity_type, "mode": "each", "iterated": records, "static": [],
                "static_first": False}

    def _elements(self):
        return {"proteins": self._axis("protein", [{"id": "a_1", "sequence": self.PROTEIN},
                                                   {"id": "b+x", "sequence": self.PROTEIN}]),
                "ligands": self._axis("ligand", [{"id": "l_2", "smiles": "CCO"}])}

    def test_the_builder_records_what_each_axis_contributed_to_every_query(self):
        provenance = {}
        pipe_openfold3_config.build_queries(self._elements(), None, None, None, provenance)
        assert provenance == {"a_1+l_2": {"proteins": "a_1", "ligands": "l_2"},
                              "b+x+l_2": {"proteins": "b+x", "ligands": "l_2"}}

    def test_a_skipped_query_is_recorded_too(self):
        import biopipelines.id_map_utils as id_map_utils
        provenance = {}
        _payload, skipped = pipe_openfold3_config.build_queries(
            self._elements(), pipe_openfold3_config.Alignments(
                id_map_utils.MsaLookup([{"id": "a_1", "file": "/nowhere/a_1.a3m"}]), "/x"),
            None, None, provenance)
        assert set(skipped) == set(provenance) == {"a_1+l_2", "b+x+l_2"}

    @staticmethod
    def _collect(tmp_path, ids, provenance_text, *flags):
        import subprocess
        pred = tmp_path / "pred"
        for query in ids:
            (pred / query / "seed_1").mkdir(parents=True, exist_ok=True)
            (pred / query / "seed_1" / f"{query}_seed_1_sample_1_model.cif").write_text("data_x\n")
            (pred / query / "seed_1" / f"{query}_seed_1_sample_1_confidences.json").write_text('{"pae": [[0.0]]}')
        (tmp_path / "queries.json").write_text(json.dumps({"queries": {q: {} for q in ids}}))
        (tmp_path / "provenance.csv").write_text(provenance_text)
        out = tmp_path / "out"
        script = os.path.join(REPO_ROOT, "pipe_scripts", "pipe_openfold3_postprocess.py")
        return out, subprocess.run(
            [sys.executable, script, "--prediction-folder", str(pred),
             "--combinatorics-config", str(tmp_path / "none.json"),
             "--queries-json", str(tmp_path / "queries.json"),
             "--provenance-csv", str(tmp_path / "provenance.csv"),
             "--structures-dir", str(out / "structures"), "--structures-map-csv", str(out / "m.csv"),
             "--confidence-csv", str(out / "c.csv"), "--local-missing-csv", str(out / "x.csv"), *flags],
            capture_output=True, text=True)

    def test_every_row_names_the_input_of_each_axis_it_was_built_from(self, tmp_path):
        import pandas as pd
        text = "id,proteins.id,ligands.id\na_1+l_2,a_1,l_2\nb+x+l_2,b+x,l_2\n"
        out, done = self._collect(tmp_path, ["a_1+l_2", "b+x+l_2"], text,
                                  "--pae-dir", str(tmp_path / "out" / "pae"),
                                  "--pae-map-csv", str(tmp_path / "out" / "pae" / "pae_map.csv"))
        assert done.returncode == 0, done.stdout + done.stderr
        frame = pd.read_csv(out / "m.csv").set_index("id")
        assert frame.loc["b+x+l_2", "proteins.id"] == "b+x", "a + inside an input id is not a separator"
        assert frame.loc["a_1+l_2", "ligands.id"] == "l_2"
        assert list(pd.read_csv(out / "c.csv").columns[-2:]) == ["proteins.id", "ligands.id"]
        pae = pd.read_csv(out / "pae" / "pae_map.csv").set_index("id")
        assert list(pae.columns) == ["file", "proteins.id", "ligands.id"]
        assert pae.loc["b+x+l_2", "proteins.id"] == "b+x"

    def test_a_failed_query_keeps_its_provenance_in_the_confidence_table(self, tmp_path):
        import pandas as pd
        (tmp_path / "skipped.csv").write_text("id,cause\nlost+l_2,no alignment\n")
        text = "id,proteins.id,ligands.id\na_1+l_2,a_1,l_2\nlost+l_2,lost,l_2\n"
        out, done = self._collect(tmp_path, ["a_1+l_2"], text, "--skipped-csv", str(tmp_path / "skipped.csv"))
        assert done.returncode == 0, done.stdout + done.stderr
        confidence = pd.read_csv(out / "c.csv").set_index("id")
        assert confidence.loc["lost+l_2", "proteins.id"] == "lost"

    def test_the_pipeline_hands_the_provenance_file_from_the_config_step_to_the_postprocess(
            self, local_config, isolated_cwd, new_pipeline):
        import re
        _tool, script = TestThePaeStream._declared(new_pipeline, isolated_cwd)
        paths = re.findall(r'--provenance-csv "([^"]+)"', script)
        assert len(paths) == 2 and paths[0] == paths[1]
        # The pipe writes it at execution time, so it belongs beside local_missing.csv, not in _configuration.
        assert paths[0].endswith(os.path.join("_execution", "query_provenance.csv"))

    def test_the_confidence_table_declares_the_columns(self, tool, monkeypatch):
        import biopipelines.openfold3 as module
        monkeypatch.setattr(module, "predict_output_ids_with_provenance",
                            lambda **k: (["a+l"], {"proteins": ["a"], "ligands": ["l"]}))
        made = tool()
        made._axis_kwargs = lambda: {}
        made.stream_path = lambda name, f: f
        made.structures_map_csv = made.confidence_csv = made.missing_csv = "x.csv"
        made.output_folder = "out"
        made.missing_table_info = lambda path: None
        assert made.get_output_files()["tables"]["confidence"].info.columns[-2:] == ["proteins.id", "ligands.id"]
