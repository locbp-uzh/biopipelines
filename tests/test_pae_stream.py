# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""The `pae` stream of Boltz2 and AlphaFold: one full predicted aligned error matrix per structure id.

A downstream tool (IPSAE) reads the whole matrix, which neither wrapper used to keep: Boltz2 recorded
only scalar confidences, and AlphaFold left the ColabFold scores JSON in the raw folder. The contract
pinned here is that the ids mirror the `structures` stream exactly (including `top_only=False` and
several diffusion samples), that the map lists only ids whose file exists, and that nothing else the
two tools write changes.
"""

import json
import os
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPTS = os.path.join(REPO_ROOT, "pipe_scripts")
BOLTZ_POST = os.path.join(SCRIPTS, "pipe_boltz_postprocessing.py")
AF_CONFIDENCE = os.path.join(SCRIPTS, "pipe_alphafold_confidence.py")
AF_UNSANITIZE = os.path.join(SCRIPTS, "pipe_unsanitize_ids.py")

sys.path.insert(0, SCRIPTS)
import pipe_check_completion  # noqa: E402


def _run(script, *args):
    done = subprocess.run([sys.executable, script, *[str(a) for a in args]],
                          capture_output=True, text=True)
    assert done.returncode == 0, done.stdout + done.stderr
    return done.stdout


# ── Boltz2 post-processing on a fake prediction folder ───────────────────────

BOLTZ_CONFIDENCE_COLUMNS = ["id", "input_file", "confidence_score", "ptm", "iptm",
                            "complex_plddt", "complex_iplddt"]


def _boltz_prediction(root, config_id, samples=1, pae_for=None):
    """boltz_results_<id>/predictions/<id>/ with `samples` models; model k carries a PAE matrix filled with k."""
    folder = root / "execution" / f"boltz_results_{config_id}" / "predictions" / config_id
    folder.mkdir(parents=True)
    for k in range(samples):
        (folder / f"{config_id}_model_{k}.pdb").write_text(f"MODEL {k}\nEND\n")
        scores = {"confidence_score": 0.9 - 0.1 * k, "ptm": 0.8, "iptm": 0.7,
                  "complex_plddt": 0.85, "complex_iplddt": 0.8}
        (folder / f"confidence_{config_id}_model_{k}.json").write_text(json.dumps(scores))
        if pae_for is None or k in pae_for:
            np.savez_compressed(folder / f"pae_{config_id}_model_{k}.npz", pae=np.full((6, 6), float(k)))


def _boltz_postprocess(root, ids, *flags, pae_args=True, provenance=None):
    for folder in ("_configuration", "structures", "tables", "sequences", "msas", "_extras", "pae"):
        (root / folder).mkdir(exist_ok=True)
    sequence_ids = root / "_configuration" / "sequence_ids.csv"
    pd.DataFrame({"id": ids}).to_csv(sequence_ids, index=False)
    if provenance is not None:
        (root / "_configuration" / "combinatorics_config.json").write_text(json.dumps(provenance))
    args = [root / "execution", root, sequence_ids,
            "--structures-map", root / "structures" / "structures_map.csv",
            "--structures-folder", root / "structures",
            "--msas-folder", root / "msas",
            "--confidence-csv", root / "tables" / "confidence.csv",
            "--affinity-csv", root / "tables" / "affinity.csv",
            "--sequences-csv", root / "sequences" / "sequences.csv",
            "--msas-csv", root / "msas" / "msas_map.csv",
            "--scores-info", root / "_extras" / "scores_info.txt"]
    if pae_args:
        args += ["--pae-map", root / "pae" / "pae_map.csv", "--pae-folder", root / "pae"]
    return _run(BOLTZ_POST, *args, *flags)


def _matrix(path):
    with np.load(path) as data:
        assert list(data.files) == ["pae"], "the array key is the one boltz writes"
        return data["pae"]


class TestBoltzPostprocessing:
    def test_top_only_keeps_the_top_models_matrix_as_the_bare_id(self, tmp_path):
        _boltz_prediction(tmp_path, "p1", samples=2)
        _boltz_postprocess(tmp_path, ["p1"])
        assert sorted(os.listdir(tmp_path / "pae")) == ["p1.npz", "pae_map.csv"]
        assert float(_matrix(tmp_path / "pae" / "p1.npz")[0, 0]) == 0.0

    def test_the_ids_mirror_the_structures_stream(self, tmp_path):
        _boltz_prediction(tmp_path, "c1", samples=1)
        _boltz_prediction(tmp_path, "c2", samples=1)
        _boltz_postprocess(tmp_path, ["c1", "c2"])
        structures = pd.read_csv(tmp_path / "structures" / "structures_map.csv")
        pae = pd.read_csv(tmp_path / "pae" / "pae_map.csv")
        assert sorted(pae["id"]) == sorted(structures["id"]) == ["c1", "c2"]
        assert all(os.path.exists(f) for f in pae["file"])

    def test_every_sample_has_its_own_matrix_when_all_samples_are_kept(self, tmp_path):
        _boltz_prediction(tmp_path, "p1", samples=3)
        _boltz_postprocess(tmp_path, ["p1"], "--all-samples")
        structures = pd.read_csv(tmp_path / "structures" / "structures_map.csv")
        pae = pd.read_csv(tmp_path / "pae" / "pae_map.csv")
        assert list(pae["id"]) == list(structures["id"]) == ["p1_1", "p1_2", "p1_3"]
        for k, sample_id in enumerate(pae["id"]):
            assert float(_matrix(tmp_path / "pae" / f"{sample_id}.npz")[0, 0]) == float(k)

    def test_a_single_kept_sample_still_gets_the_numbered_id(self, tmp_path):
        """top_only=False with one diffusion sample declares `<id>_<1..1>`, so the matrix must be `<id>_1`."""
        _boltz_prediction(tmp_path, "p1", samples=1)
        _boltz_postprocess(tmp_path, ["p1"], "--all-samples")
        assert list(pd.read_csv(tmp_path / "pae" / "pae_map.csv")["id"]) == ["p1_1"]

    def test_a_sample_boltz_wrote_no_matrix_for_is_left_out_of_the_map(self, tmp_path):
        _boltz_prediction(tmp_path, "p1", samples=3, pae_for={0, 2})
        out = _boltz_postprocess(tmp_path, ["p1"], "--all-samples")
        assert list(pd.read_csv(tmp_path / "pae" / "pae_map.csv")["id"]) == ["p1_1", "p1_3"]
        assert not (tmp_path / "pae" / "p1_2.npz").exists()
        assert "No PAE file found for p1_2" in out
        assert list(pd.read_csv(tmp_path / "structures" / "structures_map.csv")["id"]) == ["p1_1", "p1_2", "p1_3"]

    def test_no_matrix_at_all_writes_no_map_and_changes_nothing_else(self, tmp_path):
        _boltz_prediction(tmp_path, "p1", samples=1, pae_for=set())
        _boltz_postprocess(tmp_path, ["p1"])
        assert not (tmp_path / "pae" / "pae_map.csv").exists()
        assert list(pd.read_csv(tmp_path / "tables" / "confidence.csv")["id"]) == ["p1"]

    def test_the_existing_outputs_keep_their_shape(self, tmp_path):
        _boltz_prediction(tmp_path, "p1", samples=2)
        _boltz_postprocess(tmp_path, ["p1"])
        structures = pd.read_csv(tmp_path / "structures" / "structures_map.csv")
        assert list(structures.columns) == ["id", "file", "value"]
        assert (tmp_path / "structures" / "p1.pdb").exists()
        confidence = pd.read_csv(tmp_path / "tables" / "confidence.csv")
        assert list(confidence.columns) == BOLTZ_CONFIDENCE_COLUMNS
        assert list(pd.read_csv(tmp_path / "tables" / "affinity.csv").columns) == ["id", "input_file"]

    def test_the_pae_arguments_are_optional(self, tmp_path):
        """A caller that predates them keeps working; the matrices go to <output>/pae."""
        _boltz_prediction(tmp_path, "p1", samples=1)
        _boltz_postprocess(tmp_path, ["p1"], pae_args=False)
        assert (tmp_path / "pae" / "p1.npz").exists()
        assert list(pd.read_csv(tmp_path / "tables" / "confidence.csv")["id"]) == ["p1"]

    def test_the_map_carries_the_same_provenance_columns_as_structures(self, tmp_path):
        _boltz_prediction(tmp_path, "p1+l1", samples=2)
        provenance = {"predicted_ids": ["p1+l1"], "provenance": {"proteins": ["p1"], "ligands": ["l1"]}}
        _boltz_postprocess(tmp_path, ["p1+l1"], "--all-samples", provenance=provenance)
        structures = pd.read_csv(tmp_path / "structures" / "structures_map.csv")
        pae = pd.read_csv(tmp_path / "pae" / "pae_map.csv")
        for column in ("proteins.id", "ligands.id"):
            assert list(pae[column]) == list(structures[column])
        assert list(pae["proteins.id"]) == ["p1", "p1"]
        assert list(pae["ligands.id"]) == ["l1", "l1"]
        assert list(pae["id"]) == list(structures["id"])


# ── AlphaFold confidence step on a fake ColabFold folder ─────────────────────

MODEL = "alphafold2_ptm_model_2_seed_000"


def _scores(tag, with_pae=True):
    """A ColabFold scores JSON: `pae` (and `max_pae`) only exist for the ptm and multimer models."""
    scores = {"plddt": [90.0, 80.0, 70.0], "ptm": 0.81, "tag": tag}
    if with_pae:
        scores.update({"max_pae": 9.5, "pae": [[0.5, 1.0, 2.0], [1.0, 0.5, 3.0], [2.0, 3.0, 0.5]]})
    return scores


def _colabfold(folder, prefix, relaxed=False, rank=1, with_pae=True, write_scores=True):
    kind = "relaxed" if relaxed else "unrelaxed"
    folder.mkdir(parents=True, exist_ok=True)
    (folder / f"{prefix}_{kind}_rank_{rank:03d}_{MODEL}.pdb").write_text("ATOM\nEND\n")
    if write_scores:
        scores = _scores(f"rank{rank}", with_pae)
        (folder / f"{prefix}_scores_rank_{rank:03d}_{MODEL}.json").write_text(json.dumps(scores))


def _af_confidence(root, with_pae_args=True):
    (root / "structures").mkdir(exist_ok=True)
    (root / "tables").mkdir(exist_ok=True)
    args = [root / "execution", root / "structures", root / "tables" / "confidence.csv"]
    if with_pae_args:
        args += ["--pae-folder", root / "pae", "--pae-map", root / "pae" / "pae_map.csv"]
    return _run(AF_CONFIDENCE, *args)


class TestAlphaFoldScores:
    def test_the_best_models_scores_json_is_copied_per_id(self, tmp_path):
        _colabfold(tmp_path / "execution", "c1")
        _colabfold(tmp_path / "execution", "c1", rank=2)
        _af_confidence(tmp_path)
        copied = json.loads((tmp_path / "pae" / "c1.json").read_text())
        assert copied["tag"] == "rank1", "rank_002 must never be taken for the id"
        assert copied["pae"] == _scores("rank1")["pae"]
        pae = pd.read_csv(tmp_path / "pae" / "pae_map.csv")
        assert list(pae.columns) == ["id", "file"]
        assert list(pae["id"]) == ["c1"] and os.path.exists(pae["file"][0])

    def test_relaxed_in_the_folder_or_the_id_does_not_break_the_scores_lookup(self, tmp_path):
        root = tmp_path / "run_relaxed_1"
        _colabfold(root / "execution", "my_relaxed_design")
        _af_confidence(root)
        confidence = pd.read_csv(root / "tables" / "confidence.csv")
        assert list(confidence["id"]) == ["my_relaxed_design"]
        assert list(pd.read_csv(root / "pae" / "pae_map.csv")["id"]) == ["my_relaxed_design"]

    def test_two_variants_of_the_best_model_make_one_entry(self, tmp_path):
        _colabfold(tmp_path / "execution", "c1", relaxed=True)
        _colabfold(tmp_path / "execution", "c1", relaxed=False, write_scores=False)
        _af_confidence(tmp_path)
        assert list(pd.read_csv(tmp_path / "pae" / "pae_map.csv")["id"]) == ["c1"]
        assert sorted(os.listdir(tmp_path / "pae")) == ["c1.json", "pae_map.csv"]

    def test_the_ids_mirror_the_structures_that_get_a_confidence_row(self, tmp_path):
        for prefix in ("c1", "c2"):
            _colabfold(tmp_path / "execution", prefix)
        _af_confidence(tmp_path)
        confidence = pd.read_csv(tmp_path / "tables" / "confidence.csv")
        pae = pd.read_csv(tmp_path / "pae" / "pae_map.csv")
        assert sorted(pae["id"]) == sorted(confidence["id"]) == ["c1", "c2"]

    def test_a_structure_without_a_scores_json_gets_no_entry(self, tmp_path):
        _colabfold(tmp_path / "execution", "c1")
        _colabfold(tmp_path / "execution", "c2", write_scores=False)
        _af_confidence(tmp_path)
        assert list(pd.read_csv(tmp_path / "pae" / "pae_map.csv")["id"]) == ["c1"]
        assert not (tmp_path / "pae" / "c2.json").exists()

    def test_a_model_without_a_pae_matrix_still_hands_over_its_json(self, tmp_path):
        """`--model-type alphafold2` has no PAE head. The file is ColabFold's own scores JSON either way, so the consumer reports the missing key instead of the step failing."""
        _colabfold(tmp_path / "execution", "c1", with_pae=False)
        _af_confidence(tmp_path)
        assert "pae" not in json.loads((tmp_path / "pae" / "c1.json").read_text())
        assert list(pd.read_csv(tmp_path / "pae" / "pae_map.csv")["id"]) == ["c1"]

    def test_the_confidence_table_is_unchanged(self, tmp_path):
        _colabfold(tmp_path / "execution", "c1")
        _af_confidence(tmp_path)
        confidence = pd.read_csv(tmp_path / "tables" / "confidence.csv")
        assert list(confidence.columns) == ["id", "structure", "plddt", "max_pae", "ptm"]
        assert confidence.loc[0, "plddt"] == 80.0 and confidence.loc[0, "max_pae"] == 9.5

    def test_without_the_pae_arguments_nothing_new_is_written(self, tmp_path):
        _colabfold(tmp_path / "execution", "c1")
        _af_confidence(tmp_path, with_pae_args=False)
        assert not (tmp_path / "pae").exists()
        assert list(pd.read_csv(tmp_path / "tables" / "confidence.csv")["id"]) == ["c1"]

    def test_a_bundled_complex_id_survives_colabfolds_sanitizing(self, tmp_path):
        """ColabFold cannot write a `+`, so `p1+p2` is folded as `p1_p2`; the unsanitize pass the wrapper runs next must restore the file, the map id and the map path."""
        _colabfold(tmp_path / "execution", "p1_p2")
        _af_confidence(tmp_path)
        queries = tmp_path / "queries.csv"
        pd.DataFrame({"id": ["p1+p2"], "sequence": ["AAA:CCC"]}).to_csv(queries, index=False)
        _run(AF_UNSANITIZE, "--predicted-ids", queries,
             "--rename-files", tmp_path / "pae", "json",
             "--fix-csv", tmp_path / "tables" / "confidence.csv",
             "--fix-csv", tmp_path / "pae" / "pae_map.csv")
        assert sorted(os.listdir(tmp_path / "pae")) == ["p1+p2.json", "pae_map.csv"]
        pae = pd.read_csv(tmp_path / "pae" / "pae_map.csv")
        assert list(pae["id"]) == ["p1+p2"]
        assert os.path.basename(pae["file"][0]) == "p1+p2.json" and os.path.exists(pae["file"][0])


# ── the wrappers ─────────────────────────────────────────────────────────────

def _protein_pair():
    from biopipelines.sequence import Sequence
    return Sequence(seq=["MKTAYIAKQRQISFVKSHFSRQLEERLGL", "MEEKVEEIKKGAELLQSLVHPGEKDE"],
                    type="protein", ids=["p1", "p2"])


def _step_script(script_path, tool_name):
    runtime = os.path.dirname(script_path)
    matches = sorted(f for f in os.listdir(runtime) if f.endswith(f"_{tool_name}.sh"))
    assert matches, f"no {tool_name} step script under {runtime}"
    with open(os.path.join(runtime, matches[0]), encoding="utf-8", errors="replace") as handle:
        return handle.read()


def _build(new_pipeline, job, make):
    """The tool's own instance and its generated step script, built the way a user's pipeline builds it."""
    pipeline = new_pipeline(job)
    with pipeline:
        built = make()
        script_path = pipeline.save()
    return getattr(built, "producer", built), script_path


class TestBoltz2Wrapper:
    def _boltz(self, new_pipeline, job="pae_boltz", **kwargs):
        from biopipelines.boltz2 import Boltz2
        return _build(new_pipeline, job, lambda: Boltz2(proteins=_protein_pair(), **kwargs))

    def test_the_stream_is_declared_with_the_structure_ids(self, local_config, isolated_cwd, new_pipeline):
        tool, _ = self._boltz(new_pipeline)
        outputs = tool.get_output_files()
        pae, structures = outputs["pae"], outputs["structures"]
        assert pae.name == "pae" and pae.format == "npz"
        assert pae.ids == structures.ids and pae.ids_expanded == ["p1", "p2"]
        assert pae.files == [os.path.join(tool.output_folder, "pae", "<id>.npz")]
        assert pae.map_table == os.path.join(tool.output_folder, "pae", "pae_map.csv")

    @pytest.mark.parametrize("samples", [None, 1, 3])
    def test_every_diffusion_sample_has_a_matrix_when_top_only_is_off(
            self, local_config, isolated_cwd, new_pipeline, samples):
        tool, _ = self._boltz(new_pipeline, top_only=False, diffusion_samples=samples)
        outputs = tool.get_output_files()
        expected = [f"{pid}_{k}" for pid in ("p1", "p2") for k in range(1, (samples or 1) + 1)]
        assert outputs["structures"].ids_expanded == expected
        assert outputs["pae"].ids_expanded == expected

    def test_the_flag_reaches_boltz_predict(self, local_config, isolated_cwd, new_pipeline):
        _, script_path = self._boltz(new_pipeline)
        body = _step_script(script_path, "Boltz2")
        (predict,) = [l for l in body.splitlines() if "boltz predict" in l and not l.lstrip().startswith("echo")]
        assert " --write_full_pae" in predict and "--output_format pdb" in predict

    def test_the_postprocessing_step_is_told_where_to_put_them(self, local_config, isolated_cwd, new_pipeline):
        tool, script_path = self._boltz(new_pipeline)
        body = _step_script(script_path, "Boltz2")
        assert f'--pae-folder "{tool.output_folder}{os.sep}pae"' in body
        assert f"--pae-map {tool.pae_map_csv}" in body
        assert "--all-samples" not in body

    def test_top_only_off_is_still_forwarded_to_the_same_step(self, local_config, isolated_cwd, new_pipeline):
        _, script_path = self._boltz(new_pipeline, top_only=False, diffusion_samples=2)
        body = _step_script(script_path, "Boltz2")
        assert "--pae-map" in body and "--all-samples" in body

    def test_a_forwarded_write_full_pae_does_not_break_the_command(self, local_config, isolated_cwd, new_pipeline):
        """The flag is on by default; a user who also forwards it must still get a runnable command."""
        _, script_path = self._boltz(new_pipeline, write_full_pae=True)
        body = _step_script(script_path, "Boltz2")
        (predict,) = [l for l in body.splitlines() if "boltz predict" in l and not l.lstrip().startswith("echo")]
        assert predict.count("--write_full_pae") == 2

    def test_the_stream_is_registered_so_no_naming_warning_fires(self):
        from biopipelines import contract_enforcement
        assert contract_enforcement.check_stream_name("pae") is None


class TestAlphaFoldWrapper:
    def _alphafold(self, new_pipeline, job="pae_af", bundle=False):
        from biopipelines.alphafold import AlphaFold
        from biopipelines.combinatorics import Bundle

        def make():
            proteins = _protein_pair()
            return AlphaFold(proteins=Bundle(proteins) if bundle else proteins)
        return _build(new_pipeline, job, make)

    def test_the_stream_is_declared_with_the_structure_ids(self, local_config, isolated_cwd, new_pipeline):
        tool, _ = self._alphafold(new_pipeline)
        outputs = tool.get_output_files()
        pae, structures = outputs["pae"], outputs["structures"]
        assert pae.name == "pae" and pae.format == "json"
        assert pae.ids == structures.ids and pae.ids_expanded == ["p1", "p2"]
        assert pae.files == [os.path.join(tool.output_folder, "pae", "<id>.json")]
        assert pae.map_table == os.path.join(tool.output_folder, "pae", "pae_map.csv")

    def test_a_bundled_complex_declares_its_joined_id(self, local_config, isolated_cwd, new_pipeline):
        tool, _ = self._alphafold(new_pipeline, bundle=True)
        outputs = tool.get_output_files()
        assert outputs["pae"].ids == outputs["structures"].ids == ["p1+p2"]

    def test_the_confidence_step_is_told_where_to_put_them(self, local_config, isolated_cwd, new_pipeline):
        tool, script_path = self._alphafold(new_pipeline)
        body = _step_script(script_path, "AlphaFold")
        pae_dir = f"{tool.output_folder}{os.sep}pae"
        assert f'--pae-folder "{pae_dir}" --pae-map "{tool.pae_map}"' in body

    def test_the_unsanitize_pass_renames_and_fixes_them_too(self, local_config, isolated_cwd, new_pipeline):
        tool, script_path = self._alphafold(new_pipeline, bundle=True)
        body = _step_script(script_path, "AlphaFold")
        assert f'--rename-files "{tool.output_folder}{os.sep}pae" json' in body
        assert f'--fix-csv "{tool.pae_map}"' in body


# ── the completion check: the missing-id handling the stream rides on ────────

class TestCompletionCheck:
    """Both wrappers declare the ids and a `<id>` file template, so the generic check applies unchanged."""

    def _expected(self, tool, missing_rows):
        outputs = tool.get_output_files()
        tables = os.path.join(tool.output_folder, "tables")
        os.makedirs(tables, exist_ok=True)
        path = os.path.join(tables, "missing.csv")
        pd.DataFrame(missing_rows, columns=["id", "removed_by", "kind", "cause"]).to_csv(path, index=False)
        return {"output_folder": tool.output_folder,
                "structures": outputs["structures"].to_dict(),
                "pae": outputs["pae"].to_dict(),
                "tables": {"missing": {"path": path}}}

    def _write(self, tool, stream, ids, ext):
        folder = os.path.join(tool.output_folder, stream)
        for sid in ids:
            with open(os.path.join(folder, f"{sid}.{ext}"), "w") as handle:
                handle.write("x")

    @pytest.mark.parametrize("which", ["boltz", "alphafold"])
    def test_an_id_filtered_upstream_owes_neither_a_structure_nor_a_matrix(
            self, local_config, isolated_cwd, new_pipeline, which):
        if which == "boltz":
            tool, _ = TestBoltz2Wrapper()._boltz(new_pipeline)
            structure_ext, pae_ext, name = "pdb", "npz", "Boltz2"
        else:
            tool, _ = TestAlphaFoldWrapper()._alphafold(new_pipeline)
            structure_ext, pae_ext, name = "pdb", "json", "AlphaFold"
        expected = self._expected(tool, [("p2", "001_Panda", "filter", "dropped upstream")])
        self._write(tool, "structures", ["p1"], structure_ext)
        self._write(tool, "pae", ["p1"], pae_ext)
        ok, missing = pipe_check_completion.check_expected_outputs(expected, name, tool.output_folder)
        assert ok, missing

    def test_a_surviving_structure_without_its_matrix_is_reported(self, local_config, isolated_cwd, new_pipeline):
        tool, _ = TestBoltz2Wrapper()._boltz(new_pipeline)
        expected = self._expected(tool, [("p2", "001_Panda", "filter", "dropped upstream")])
        self._write(tool, "structures", ["p1"], "pdb")
        ok, missing = pipe_check_completion.check_expected_outputs(expected, "Boltz2", tool.output_folder)
        assert not ok and list(missing) == ["pae"]
