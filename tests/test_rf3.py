"""RF3: RoseTTAFold3 co-folding through Foundry's `rf3 fold`.

The input axes are those of Boltz2 and OpenFold3, so most of what is worth pinning is where RF3 differs: the
Hydra override spelling of its command line, the explicit chain type on nucleic acids (a string such as "ACGT"
is also a protein), the on-disk layout of its outputs (read from upstream's source, not from a real run), and the
early stop, which is a deliberate drop and not a failure.
"""

import gzip
import json
import os
import subprocess
import sys

import pytest

from biopipelines.id_map_utils import MsaLookup
from biopipelines.rf3 import RF3
from tests._helpers import sys_path_prepended

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
with sys_path_prepended(os.path.join(REPO_ROOT, "pipe_scripts")):
    import pipe_rf3_config
    import pipe_rf3_postprocess as post


@pytest.fixture
def tool(local_config, isolated_cwd):
    def _make(**kwargs):
        made = RF3.__new__(RF3)
        made.proteins = kwargs.pop("proteins", object())
        for axis in ("ssDNA", "dsDNA", "ssRNA", "dsRNA", "ligands", "msas"):
            setattr(made, axis, kwargs.pop(axis, None))
        made.msas_stream = kwargs.pop("msas_stream", None)
        made.n_recycles = kwargs.pop("n_recycles", None)
        made.diffusion_batch_size = kwargs.pop("diffusion_batch_size", 5)
        made.num_steps = kwargs.pop("num_steps", None)
        made.seed = kwargs.pop("seed", None)
        made.early_stopping_plddt_threshold = kwargs.pop("early_stopping_plddt_threshold", 0.5)
        made.top_only = kwargs.pop("top_only", True)
        made.checkpoint_path = kwargs.pop("checkpoint_path", None)
        assert not kwargs, f"unexpected: {sorted(kwargs)}"
        return made
    return _make


class TestItRefusesWhatItCannotDo:

    def test_no_input_axis_is_refused(self, tool):
        with pytest.raises(ValueError, match="at least one input axis"):
            RF3.validate_params(tool(proteins=None))

    @pytest.mark.parametrize("name", ["n_recycles", "diffusion_batch_size", "num_steps"])
    @pytest.mark.parametrize("bad", [0, -1, True, 2.5])
    def test_a_count_must_be_a_positive_integer(self, tool, name, bad):
        with pytest.raises(ValueError, match=name):
            RF3.validate_params(tool(**{name: bad}))

    @pytest.mark.parametrize("bad", [-0.1, 1.5, True, "0.5"])
    def test_the_early_stop_threshold_is_a_fraction(self, tool, bad):
        with pytest.raises(ValueError, match="early_stopping_plddt_threshold"):
            RF3.validate_params(tool(early_stopping_plddt_threshold=bad))

    @pytest.mark.parametrize("ok", [0, 0.0, 0.5, 1])
    def test_zero_turns_early_stopping_off_and_is_accepted(self, tool, ok):
        RF3.validate_params(tool(early_stopping_plddt_threshold=ok))

    def test_one_recycle_is_refused_while_early_stopping_is_on(self, tool):
        """Upstream takes the first recycle for the early-stopping check and then needs another one to fold with; with one recycle every input fails with `Recycling generator produced no outputs`."""
        with pytest.raises(ValueError, match="n_recycles"):
            RF3.validate_params(tool(n_recycles=1))

    def test_one_recycle_is_accepted_once_early_stopping_is_off(self, tool):
        RF3.validate_params(tool(n_recycles=1, early_stopping_plddt_threshold=0))

    def test_two_recycles_are_accepted_with_early_stopping_on(self, tool):
        RF3.validate_params(tool(n_recycles=2))

    def test_msas_without_proteins_are_refused(self, tool):
        """An alignment has nothing to attach to; accepting it would leave the user believing one was used."""
        with pytest.raises(ValueError, match="proteins axis is required"):
            RF3.validate_params(tool(proteins=None, ligands=object(), msas=object()))

    def test_a_csv_msa_stream_is_refused(self, tool):
        from types import SimpleNamespace
        with pytest.raises(ValueError, match="a3m or fasta"):
            RF3.validate_params(tool(msas=object(), msas_stream=SimpleNamespace(format="csv")))

    def test_a_checkpoint_path_cannot_carry_shell_expansion(self, tool):
        with pytest.raises(ValueError, match="checkpoint_path"):
            RF3.validate_params(tool(checkpoint_path="$(rm -rf x)/rf3.ckpt"))

    def test_an_integer_seed_is_required(self, tool):
        with pytest.raises(ValueError, match="seed"):
            RF3.validate_params(tool(seed="7"))
        RF3.validate_params(tool(seed=7))


class TestTheCommandLine:
    """Hydra spells overrides `key=value`; a `--flag` would be rejected by `rf3 fold`."""

    @staticmethod
    def _section(made):
        made.queries_json = "/q.json"
        made.prediction_folder = "/out"
        made.folders = {"RF3": "/weights"}
        made.extra_args_bash_tokens = lambda: []
        made.extra_args_echo = lambda: ""
        made.container_prefix = lambda: ""
        return made

    def test_the_defaults_are_always_stated_and_nothing_else_is(self, tool):
        made = self._section(tool())
        assert RF3._predict_overrides(made) == [
            'inputs="/q.json"', 'out_dir="/out"', "diffusion_batch_size=5",
            "early_stopping_plddt_threshold=0.5"]

    def test_optional_settings_appear_only_when_given(self, tool):
        made = self._section(tool(n_recycles=3, num_steps=50, seed=7, checkpoint_path="/m.ckpt",
                                  diffusion_batch_size=2, early_stopping_plddt_threshold=0))
        overrides = RF3._predict_overrides(made)
        for expected in ("n_recycles=3", "num_steps=50", "seed=7", 'ckpt_path="/m.ckpt"',
                         "diffusion_batch_size=2", "early_stopping_plddt_threshold=0"):
            assert expected in overrides
        assert not any(o.startswith("--") for o in overrides)

    def test_forwarded_overrides_ride_along(self, tool):
        made = self._section(tool())
        made.extra_args_bash_tokens = lambda: ['"template_selection=[A]"']
        assert RF3._predict_overrides(made)[-1] == '"template_selection=[A]"'

    def test_the_weights_are_found_through_the_search_path_and_checked_first(self, tool):
        text = RF3._generate_predict_section(self._section(tool()))
        assert 'export FOUNDRY_CHECKPOINT_DIRS="/weights"' in text
        assert "rf3 fold" in text
        assert text.index("*.ckpt") < text.index("rf3 fold"), "a missing checkpoint must fail before the model loads"


class TestTheInputsJson:

    def test_a_protein_carries_its_msa_and_its_chain_label(self, tmp_path):
        alignment = tmp_path / "p1.a3m"
        alignment.write_text(">p1\nMKT\n")
        built = pipe_rf3_config.build_component(
            "protein", {"id": "p1", "sequence": "MKT"}, "A", MsaLookup([{"id": "p1", "file": str(alignment)}]))
        assert built == {"seq": "MKT", "chain_id": "A", "msa_path": str(alignment)}

    def test_a_protein_without_an_alignment_has_no_msa_key(self):
        built = pipe_rf3_config.build_component("protein", {"id": "p1", "sequence": "MKT"}, "A", None)
        assert "msa_path" not in built

    def test_nucleic_acids_name_their_chain_type(self):
        """Upstream infers a type from the letters, and ACGT reads as a protein as well as DNA."""
        dna = pipe_rf3_config.build_component("ssdna", {"id": "d", "sequence": "ACGT"}, "A", None)
        rna = pipe_rf3_config.build_component("ssrna", {"id": "r", "sequence": "ACGU"}, "A", None)
        assert dna["chain_type"] == "polydeoxyribonucleotide"
        assert rna["chain_type"] == "polyribonucleotide"

    def test_the_second_strand_of_a_double_stranded_axis_is_the_reverse_complement(self):
        elements = {"dsDNA": {"mode": "each", "entity_type": "dsdna", "static_first": False, "static": [],
                              "iterated": [{"id": "d", "sequence": "AAAG"}]}}
        (example,), _ = pipe_rf3_config.build_examples(elements, None)
        assert [c["seq"] for c in example["components"]] == ["AAAG", "CTTT"]
        assert [c["chain_id"] for c in example["components"]] == ["A", "B"]
        assert all(c["chain_type"] == "polydeoxyribonucleotide" for c in example["components"])

    def test_a_ccd_ligand_uses_ccd_code_and_a_smiles_ligand_uses_smiles(self):
        ccd = pipe_rf3_config.build_component(
            "ligand", {"id": "l1", "ccd": "NAG", "smiles": "CC", "source": "rcsb"}, "B", None)
        assert ccd == {"ccd_code": "NAG", "chain_id": "B"}
        smiles = pipe_rf3_config.build_component("ligand", {"id": "l2", "smiles": "CCO"}, "B", None)
        assert smiles == {"smiles": "CCO", "chain_id": "B"}

    def test_an_example_is_named_by_the_id_the_wrapper_predicted(self):
        elements = {
            "proteins": {"mode": "each", "entity_type": "protein", "static_first": False, "static": [],
                         "iterated": [{"id": "p1", "sequence": "MKT"}, {"id": "p2", "sequence": "MAA"}]},
            "ligands": {"mode": "each", "entity_type": "ligand", "static_first": False, "static": [],
                        "iterated": [{"id": "l1", "smiles": "CCO"}]},
        }
        examples, _ = pipe_rf3_config.build_examples(elements, None)
        assert [e["name"] for e in examples] == ["p1+l1", "p2+l1"]
        assert [c["chain_id"] for c in examples[0]["components"]] == ["A", "B"]

    def test_a_bundle_folds_every_member_in_one_example(self):
        elements = {"proteins": {"mode": "bundle", "entity_type": "protein", "static_first": False,
                                 "iterated": [], "static": [{"id": "a", "sequence": "MKT"},
                                                            {"id": "b", "sequence": "MAA"}]}}
        (example,), _ = pipe_rf3_config.build_examples(elements, None)
        assert [c["seq"] for c in example["components"]] == ["MKT", "MAA"]

    def test_an_empty_sequence_is_refused_by_the_id_that_caused_it(self):
        with pytest.raises(ValueError, match="'p9'"):
            pipe_rf3_config.build_component("protein", {"id": "p9", "sequence": float("nan")}, "A", None)

    def test_an_iterating_axis_with_no_records_builds_no_example(self):
        elements = {"proteins": {"mode": "each", "entity_type": "protein", "iterated": [],
                                 "static": [], "static_first": False}}
        assert pipe_rf3_config.build_examples(elements, None) == ([], {})

    @pytest.mark.parametrize("name", ["a.a3m", "a.a3m.gz", "a.fasta"])
    def test_the_alignment_formats_upstream_reads_are_accepted(self, name, tmp_path):
        (tmp_path / name).write_text(">p1\nMKT\n")
        msas = MsaLookup([{"id": "p1", "file": str(tmp_path / name)}])
        built = pipe_rf3_config.build_component("protein", {"id": "p1", "sequence": "MKT"}, "A", msas)
        assert built["msa_path"] == str(tmp_path / name)

    def test_an_alignment_upstream_cannot_read_is_named_not_ignored(self):
        msas = MsaLookup([{"id": "p1", "file": "/x/p1.sto"}])
        with pytest.raises(ValueError, match=r"\.sto"):
            pipe_rf3_config.build_component("protein", {"id": "p1", "sequence": "MKT"}, "A", msas)


class TestAnAlignmentThatCannotBeUsedSkipsOnlyItsExample:
    """A row whose file was purged, or whose format upstream cannot read, must not take the other examples with it."""

    @staticmethod
    def _elements():
        records = [{"id": "p1", "sequence": "MKT"}, {"id": "p2", "sequence": "MAA"}, {"id": "p3", "sequence": "MGG"}]
        return {"proteins": {"mode": "each", "entity_type": "protein", "static_first": False, "static": [],
                             "iterated": records}}

    def test_a_row_whose_file_is_gone_is_skipped_and_the_examples_after_it_are_built(self, tmp_path):
        for name in ("p1", "p3"):
            (tmp_path / f"{name}.a3m").write_text(f">{name}\nMKT\n")
        msas = MsaLookup([{"id": name, "file": str(tmp_path / f"{name}.a3m")} for name in ("p1", "p2", "p3")])
        examples, skipped = pipe_rf3_config.build_examples(self._elements(), msas)
        assert [e["name"] for e in examples] == ["p1", "p3"]
        assert list(skipped) == ["p2"]
        assert "does not exist" in skipped["p2"] and str(tmp_path / "p2.a3m") in skipped["p2"]

    def test_a_row_in_a_format_upstream_cannot_read_is_skipped_not_fatal(self, tmp_path):
        for name in ("p1.a3m", "p2.sto", "p3.a3m"):
            (tmp_path / name).write_text(">x\nMKT\n")
        msas = MsaLookup([{"id": f"p{i}", "file": str(tmp_path / name)}
                          for i, name in enumerate(("p1.a3m", "p2.sto", "p3.a3m"), start=1)])
        examples, skipped = pipe_rf3_config.build_examples(self._elements(), msas)
        assert [e["name"] for e in examples] == ["p1", "p3"]
        assert list(skipped) == ["p2"] and ".sto" in skipped["p2"]

    def test_the_skipped_example_is_a_failure_row_of_the_local_missing_table(self, tmp_path):
        import pandas as pd
        from biopipelines.biopipelines_io import write_skipped_queries
        msas = MsaLookup([{"id": "p1", "file": str(tmp_path / "purged.a3m")}])
        elements = {"proteins": {"mode": "each", "entity_type": "protein", "static_first": False, "static": [],
                                 "iterated": [{"id": "p1", "sequence": "MKT"}]}}
        _examples, skipped = pipe_rf3_config.build_examples(elements, msas)
        write_skipped_queries(str(tmp_path / "skipped.csv"), skipped)
        (tmp_path / "queries.json").write_text("[]")
        out = tmp_path / "out"
        done = subprocess.run(
            [sys.executable, os.path.join(REPO_ROOT, "pipe_scripts", "pipe_rf3_postprocess.py"),
             "--prediction-folder", str(tmp_path / "pred"), "--queries-json", str(tmp_path / "queries.json"),
             "--skipped-csv", str(tmp_path / "skipped.csv"), "--structures-dir", str(out / "structures"),
             "--structures-map-csv", str(out / "m.csv"), "--confidence-csv", str(out / "c.csv"),
             "--pae-dir", str(out / "pae"), "--pae-map-csv", str(out / "pae_map.csv"),
             "--local-missing-csv", str(out / "001_RF3" / "execution" / "local_missing.csv"),
             "--diffusion-batch-size", "2"], capture_output=True, text=True)
        assert done.returncode == 1
        missing = pd.read_csv(out / "001_RF3" / "execution" / "local_missing.csv")
        assert missing.loc[0, "kind"] == "failure" and "does not exist" in missing.loc[0, "cause"]
        assert pd.read_csv(out / "c.csv")["id"].tolist() == ["p1"]


def _write_example(pred, name, ranking, *, gz=False, with_confidences=True):
    """The tree `rf3 fold` writes for one example: top-level files plus one folder per sample."""
    folder = pred / name
    folder.mkdir(parents=True, exist_ok=True)

    def write_set(directory, stem, score, mark):
        extension = "cif.gz" if gz else "cif"
        model = directory / f"{stem}_model.{extension}"
        if gz:
            with gzip.open(model, "wt") as handle:
                handle.write(f"data_{mark}\n")
        else:
            model.write_text(f"data_{mark}\n")
        (directory / f"{stem}_summary_confidences.json").write_text(json.dumps({
            "overall_plddt": 0.8 + score / 100, "ptm": 0.7, "iptm": 0.6, "overall_pae": 4.0,
            "overall_pde": 3.0, "has_clash": False, "ranking_score": score}))
        if with_confidences:
            (directory / f"{stem}_confidences.json").write_text(json.dumps({
                "atom_plddts": [0.9], "pae": [[0.5, score], [score, 0.5]], "token_chain_ids": ["A", "A"]}))

    scores = list(ranking)
    best = max(range(len(scores)), key=lambda i: scores[i])
    write_set(folder, name, scores[best], f"{name}_sample{best}")
    (folder / f"{name}_ranking_scores.csv").write_text(
        "seed,sample,ranking_score\n" + "".join(f"0,{i},{s}\n" for i, s in enumerate(scores)))
    for index, score in enumerate(scores):
        sample_dir = folder / f"seed-0_sample-{index}"
        sample_dir.mkdir()
        write_set(sample_dir, f"{name}_seed-0_sample-{index}", score, f"{name}_sample{index}")


class TestCollectingAnExample:

    def test_the_top_ranked_sample_is_kept_and_labeled(self, tmp_path):
        _write_example(tmp_path, "q1", [0.2, 0.9, 0.5])
        result = post.collect_example(str(tmp_path), "q1", True, 3)
        (output,) = result["outputs"]
        assert output["id"] == "q1" and result["dropped"] == []
        assert (output["seed"], output["sample"]) == (0, 1)
        assert output["summary"]["ranking_score"] == 0.9
        assert output["summary"]["has_clash"] is False

    def test_every_sample_is_surfaced_when_all_are_kept(self, tmp_path):
        _write_example(tmp_path, "q1", [0.2, 0.9])
        result = post.collect_example(str(tmp_path), "q1", False, 2)
        assert [o["id"] for o in result["outputs"]] == ["q1_1", "q1_2"]
        assert [o["summary"]["ranking_score"] for o in result["outputs"]] == [0.2, 0.9]

    def test_samples_upstream_did_not_write_are_declared_missing_not_forgotten(self, tmp_path):
        _write_example(tmp_path, "q1", [0.2, 0.9])
        result = post.collect_example(str(tmp_path), "q1", False, 5)
        assert [d["id"] for d in result["dropped"]] == ["q1_3", "q1_4", "q1_5"]
        assert {d["kind"] for d in result["dropped"]} == {"failure"}

    def test_a_gzipped_model_is_found(self, tmp_path):
        _write_example(tmp_path, "q1", [0.5], gz=True)
        (output,) = post.collect_example(str(tmp_path), "q1", True, 1)["outputs"]
        assert output["structure"].endswith("_model.cif.gz")

    def test_an_example_with_no_folder_is_a_failure_that_names_it(self, tmp_path):
        result = post.collect_example(str(tmp_path), "lost", True, 5)
        assert result["outputs"] == []
        (dropped,) = result["dropped"]
        assert (dropped["id"], dropped["kind"]) == ("lost", "failure")
        assert "no output folder" in dropped["cause"]

    def test_an_early_stop_recorded_in_the_ranking_csv_is_a_filter(self, tmp_path):
        folder = tmp_path / "q1"
        folder.mkdir()
        (folder / "q1_ranking_scores.csv").write_text("early_stopped,mean_plddt\nTrue,0.31\n")
        (dropped,) = post.collect_example(str(tmp_path), "q1", True, 5, 0.5)["dropped"]
        assert dropped["kind"] == "filter" and dropped["early_stopped"] is True
        assert "0.310" in dropped["cause"] and "0.5" in dropped["cause"]

    def test_an_early_stop_recorded_in_score_csv_is_a_filter_too(self, tmp_path):
        """The 0.2.0 release writes score.csv rather than a ranking row."""
        folder = tmp_path / "q1"
        folder.mkdir()
        (folder / "score.csv").write_text("early_stopped,mean_plddt\nTrue,0.2\n")
        (dropped,) = post.collect_example(str(tmp_path), "q1", True, 5)["dropped"]
        assert dropped["kind"] == "filter"

    def test_a_folder_with_nothing_in_it_is_a_failure_not_an_early_stop(self, tmp_path):
        (tmp_path / "q1").mkdir()
        (dropped,) = post.collect_example(str(tmp_path), "q1", True, 5)["dropped"]
        assert dropped["kind"] == "failure" and dropped["early_stopped"] is False

    def test_a_missing_summary_leaves_blanks_not_a_dropped_structure(self, tmp_path):
        _write_example(tmp_path, "q1", [0.5])
        os.remove(tmp_path / "q1" / "q1_summary_confidences.json")
        (output,) = post.collect_example(str(tmp_path), "q1", True, 1)["outputs"]
        assert set(output["summary"].values()) == {None}


class TestThePostprocessStep:

    @staticmethod
    def _run(tmp_path, names, *flags):
        queries = tmp_path / "queries.json"
        queries.write_text(json.dumps([{"name": n, "components": []} for n in names]))
        out = tmp_path / "out"
        script = os.path.join(REPO_ROOT, "pipe_scripts", "pipe_rf3_postprocess.py")
        done = subprocess.run(
            [sys.executable, script, "--prediction-folder", str(tmp_path / "pred"),
             "--queries-json", str(queries), "--structures-dir", str(out / "structures"),
             "--structures-map-csv", str(out / "structures" / "m.csv"),
             "--confidence-csv", str(out / "tables" / "c.csv"),
             "--pae-dir", str(out / "pae"), "--pae-map-csv", str(out / "pae" / "pae_map.csv"),
             "--local-missing-csv", str(out / "001_RF3" / "execution" / "local_missing.csv"),
             "--diffusion-batch-size", "2", *flags],
            capture_output=True, text=True)
        return out, done

    def test_each_artifact_answers_its_own_question(self, tmp_path):
        import pandas as pd
        pred = tmp_path / "pred"
        _write_example(pred, "ok", [0.2, 0.9])
        (pred / "stopped").mkdir()
        (pred / "stopped" / "stopped_ranking_scores.csv").write_text("early_stopped,mean_plddt\nTrue,0.3\n")
        out, done = self._run(tmp_path, ["ok", "stopped", "lost"], "--early-stopping-plddt-threshold", "0.5")
        assert done.returncode == 0, done.stderr

        confidence = pd.read_csv(out / "tables" / "c.csv").set_index("id")
        assert list(confidence.index) == ["ok", "stopped", "lost"], "one row per id that entered"
        assert confidence.loc["ok", "ranking_score"] == 0.9 and confidence.loc["ok", "sample"] == 1
        assert confidence.loc["ok", "pae"] == 4.0 and confidence.loc["ok", "plddt"] == pytest.approx(0.809)
        assert bool(confidence.loc["stopped", "early_stopped"]) is True
        assert pd.isna(confidence.loc["stopped", "plddt"])
        assert bool(confidence.loc["lost", "early_stopped"]) is False

        assert pd.read_csv(out / "structures" / "m.csv")["id"].tolist() == ["ok"]
        assert pd.read_csv(out / "pae" / "pae_map.csv")["id"].tolist() == ["ok"]
        missing = pd.read_csv(out / "001_RF3" / "execution" / "local_missing.csv").set_index("id")
        assert missing.loc["stopped", "kind"] == "filter" and missing.loc["lost", "kind"] == "failure"
        assert missing.loc["lost", "removed_by"] == "001_RF3"

    def test_the_pae_file_is_the_confidences_json_unchanged(self, tmp_path):
        _write_example(tmp_path / "pred", "ok", [0.2, 0.9])
        out, done = self._run(tmp_path, ["ok"])
        assert done.returncode == 0, done.stderr
        source = tmp_path / "pred" / "ok" / "ok_confidences.json"
        assert (out / "pae" / "ok.json").read_bytes() == source.read_bytes()
        assert json.loads((out / "pae" / "ok.json").read_text())["pae"] == [[0.5, 0.9], [0.9, 0.5]]

    def test_a_gzipped_model_is_staged_decompressed(self, tmp_path):
        _write_example(tmp_path / "pred", "ok", [0.5], gz=True)
        out, done = self._run(tmp_path, ["ok"])
        assert done.returncode == 0, done.stderr
        assert (out / "structures" / "ok.cif").read_text() == "data_ok_sample0\n"

    def test_top_only_off_numbers_the_samples(self, tmp_path):
        _write_example(tmp_path / "pred", "ok", [0.2, 0.9])
        out, done = self._run(tmp_path, ["ok"], "--top-only", "False")
        assert done.returncode == 0, done.stderr
        assert sorted(p.stem for p in (out / "structures").glob("*.cif")) == ["ok_1", "ok_2"]
        assert sorted(p.stem for p in (out / "pae").glob("*.json")) == ["ok_1", "ok_2"]

    def test_every_row_names_the_input_of_each_axis_it_was_built_from(self, tmp_path):
        import pandas as pd
        _write_example(tmp_path / "pred", "a_1+l_2", [0.2, 0.9])
        _write_example(tmp_path / "pred", "b+x+l_2", [0.5])
        (tmp_path / "provenance.csv").write_text(
            "id,proteins.id,ligands.id\na_1+l_2,a_1,l_2\nb+x+l_2,b+x,l_2\nlost+l_2,lost,l_2\n")
        out, done = self._run(tmp_path, ["a_1+l_2", "b+x+l_2", "lost+l_2"],
                              "--provenance-csv", str(tmp_path / "provenance.csv"), "--top-only", "False")
        assert done.returncode == 0, done.stderr
        structures = pd.read_csv(out / "structures" / "m.csv").set_index("id")
        assert structures.loc["a_1+l_2_2", "proteins.id"] == "a_1" and structures.loc["a_1+l_2_2", "ligands.id"] == "l_2"
        assert structures.loc["b+x+l_2_1", "proteins.id"] == "b+x", "a + inside an input id is not a separator"
        pae = pd.read_csv(out / "pae" / "pae_map.csv").set_index("id")
        assert pae.loc["a_1+l_2_1", "proteins.id"] == "a_1"
        confidence = pd.read_csv(out / "tables" / "c.csv").set_index("id")
        assert list(confidence.columns[-2:]) == ["proteins.id", "ligands.id"]
        assert confidence.loc["lost+l_2_1", "proteins.id"] == "lost", "an id that failed keeps its provenance"

    def test_without_a_provenance_file_no_column_is_added(self, tmp_path):
        import pandas as pd
        _write_example(tmp_path / "pred", "ok", [0.5])
        out, done = self._run(tmp_path, ["ok"])
        assert done.returncode == 0, done.stderr
        assert list(pd.read_csv(out / "structures" / "m.csv").columns) == ["id", "file"]

    def test_the_config_step_records_the_provenance_of_every_example_including_the_skipped(self, tmp_path):
        provenance = {}
        elements = {
            "proteins": {"mode": "each", "entity_type": "protein", "static_first": False, "static": [],
                         "iterated": [{"id": "a_1", "sequence": "MKT"}, {"id": "b+x", "sequence": "MAA"}]},
            "ligands": {"mode": "each", "entity_type": "ligand", "static_first": False, "static": [],
                        "iterated": [{"id": "l_2", "smiles": "CCO"}]}}
        msas = MsaLookup([{"id": "a_1", "file": str(tmp_path / "gone.a3m")}])
        examples, skipped = pipe_rf3_config.build_examples(elements, msas, provenance)
        assert not examples and set(skipped) == {"a_1+l_2", "b+x+l_2"}
        assert provenance == {"a_1+l_2": {"proteins": "a_1", "ligands": "l_2"},
                              "b+x+l_2": {"proteins": "b+x", "ligands": "l_2"}}

    def test_every_example_failing_exits_non_zero(self, tmp_path):
        out, done = self._run(tmp_path, ["lost"])
        assert done.returncode == 1
        assert "no RF3 prediction produced a structure" in done.stdout

    def test_every_example_early_stopping_is_not_a_failed_step(self, tmp_path):
        """The tool worked; it selected nothing."""
        (tmp_path / "pred" / "q").mkdir(parents=True)
        (tmp_path / "pred" / "q" / "q_ranking_scores.csv").write_text("early_stopped,mean_plddt\nTrue,0.1\n")
        out, done = self._run(tmp_path, ["q"])
        assert done.returncode == 0, done.stderr


class TestTheDeclaredOutputs:

    @staticmethod
    def _declared(new_pipeline, **kwargs):
        from pathlib import Path
        from biopipelines.sequence import Sequence
        pipeline = new_pipeline("rf3_outputs")
        with pipeline:
            tool = RF3(proteins=Sequence(seq="MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ",
                                         type="protein", ids="p1"), **kwargs)
            script_path = pipeline.save()
        steps = sorted(Path(script_path).parent.glob("*_RF3.sh"))
        return tool, "\n".join(step.read_text(encoding="utf-8", errors="replace") for step in steps)

    def test_structures_and_pae_share_ids_and_formats(self, local_config, isolated_cwd, new_pipeline):
        tool, _script = self._declared(new_pipeline)
        assert tool.streams.structures.ids == tool.streams.pae.ids == ["p1"]
        assert tool.streams.structures.format == "cif" and tool.streams.pae.format == "json"
        assert tool.streams.pae.files[0].endswith(os.path.join("pae", "<id>.json"))
        assert tool.streams.pae.map_table.endswith(os.path.join("pae", "pae_map.csv"))

    def test_every_sample_is_declared_when_all_are_kept(self, local_config, isolated_cwd, new_pipeline):
        """The completion check demands exactly the declared ids, so they must say how many samples there are."""
        tool, _script = self._declared(new_pipeline, top_only=False, diffusion_batch_size=3)
        assert tool.streams.structures.ids == tool.streams.pae.ids == ["p1_<1..3>"]

    def test_the_tables_are_declared_with_their_columns(self, local_config, isolated_cwd):
        from biopipelines.mock import Mock
        from biopipelines.pipeline import Pipeline
        with Pipeline(project="TestSuite", job="rf3_tables", description="x", on_the_fly=False,
                      local_output=True, config="local") as pipeline:
            seqs = Mock(ids=["p"], streams={"sequences": {"format": "csv", "file": "<id>.csv"}})
            RF3(proteins=seqs.streams.sequences)
            tables = pipeline.tools[-1].get_output_files()["tables"]
        assert tables["confidence"].info.columns == [
            "id", "plddt", "ptm", "iptm", "pae", "pde", "has_clash", "ranking_score",
            "early_stopped", "seed", "sample", "proteins.id"]
        assert tables["missing"].info.columns == ["id", "removed_by", "kind", "cause"]

    def test_the_script_runs_config_then_rf3_then_postprocess(self, local_config, isolated_cwd, new_pipeline):
        _tool, script = self._declared(new_pipeline, seed=3)
        assert script.index("pipe_rf3_config.py") < script.index("rf3 fold") < script.index("pipe_rf3_postprocess.py")
        assert "seed=3" in script
        assert "--early-stopping-plddt-threshold 0.5" in script
        assert "--diffusion-batch-size 5" in script

    def test_the_postprocess_flags_are_all_known_to_the_pipe_script(self, local_config, isolated_cwd, new_pipeline):
        """A flag the wrapper sends and the pipe does not read would fail the step on the cluster."""
        import re
        _tool, script = self._declared(new_pipeline)
        block = script.split("pipe_rf3_postprocess.py")[1].split("\n\n")[0]
        sent = set(re.findall(r"--[a-z-]+", block))
        with open(os.path.join(REPO_ROOT, "pipe_scripts", "pipe_rf3_postprocess.py")) as handle:
            known = set(re.findall(r'add_argument\("(--[a-z-]+)"', handle.read()))
        assert sent <= known

    def test_the_config_step_hands_its_provenance_file_to_the_postprocess(self, local_config, isolated_cwd, new_pipeline):
        import re
        _tool, script = self._declared(new_pipeline)
        paths = re.findall(r'--provenance-csv "([^"]+)"', script)
        assert len(paths) == 2 and paths[0] == paths[1]
        assert paths[0].endswith(os.path.join("_configuration", "query_provenance.csv"))

    def test_the_config_flags_are_all_known_to_the_pipe_script(self, local_config, isolated_cwd, new_pipeline):
        import re
        _tool, script = self._declared(new_pipeline)
        block = script.split("pipe_rf3_config.py")[1].split("\n\n")[0]
        sent = set(re.findall(r"--[a-z-]+", block))
        with open(os.path.join(REPO_ROOT, "pipe_scripts", "pipe_rf3_config.py")) as handle:
            known = set(re.findall(r'add_argument\("(--[a-z-]+)"', handle.read()))
        assert sent <= known

    def test_it_advertises_no_affinity_parameter(self):
        import inspect
        assert "affinity" not in inspect.signature(RF3.__init__).parameters


class TestTheInstall:
    """The foundry environment is shared with RFdiffusion3; the weights and the shortcut are not."""

    @pytest.fixture(autouse=True)
    def _environment(self, local_config, monkeypatch):
        monkeypatch.setattr(RF3, "_install_env", classmethod(lambda cls, manager: "foundry"))

    def test_it_downloads_its_own_weights_not_rfd3s(self):
        script = RF3._install_script({"biopipelines": "/repo", "RF3": "/data/rf3"}, env_manager="mamba")
        assert 'foundry install rf3 --checkpoint-dir "/data/rf3"' in script
        assert "rfd3" not in script

    def test_the_shortcut_requires_an_rf3_checkpoint_not_just_the_environment(self):
        script = RF3._install_script({"biopipelines": "/repo", "RF3": "/data/rf3"}, env_manager="mamba")
        shortcut = script.split("already installed, skipping")[0]
        assert '/data/rf3"/*.ckpt' in shortcut

    def test_the_environment_is_created_only_when_absent(self):
        script = RF3._install_script({"biopipelines": "/repo", "RF3": "/data/rf3"}, env_manager="mamba")
        assert "environment already exists, skipping creation" in script

    def test_a_config_without_the_folder_is_refused_by_name(self):
        with pytest.raises(ValueError, match="repositories.RF3"):
            RF3._install_script({"biopipelines": "/repo"}, env_manager="mamba")
