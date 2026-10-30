# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Pair ids and the rows that name them, for the tools that take a protein and a ligand or a pair table.

The repo rule is that an output id and its provenance are produced only by `combinatorics.py`.
Here: the Gnina/Vina `missing` rows named the wrong step, USalign and TMalign let two rows share an
id, and PocketGen joined a pair id by hand and recovered its parts with a split on `+`.
PLACER, Gnina's conformer choice and HBDesigner have files of their own.

Nothing here runs a model or a cluster job; the pipe scripts are driven on synthetic inputs.
"""

import importlib.util
import os
import re
import subprocess
import sys

import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPTS = os.path.join(REPO_ROOT, "pipe_scripts")


def _load(name, filename):
    sys.path.insert(0, SCRIPTS)
    spec = importlib.util.spec_from_file_location(name, os.path.join(SCRIPTS, filename))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _script_text(name):
    with open(os.path.join(SCRIPTS, name)) as handle:
        return handle.read()


# ── Gnina / Vina: a failure names the step that ran ───────────────────────────

def _vina_docking_config(step_folder):
    tables = step_folder / "tables"
    tables.mkdir(parents=True)
    return {
        "output_folder": str(step_folder / "_execution"), "backend": "vina", "gnina_binary": "vina",
        "exhaustiveness": 8, "num_modes": 2, "num_runs": 1, "seed": 1, "cnn_scoring": "rescore",
        "cnn_score_threshold": 0.0, "scoring": "vina", "box": {"center": "0,0,0", "size": 20},
        "extra_args": [], "dock_timeout": 0, "rmsd_threshold": 2.0,
        "docking_results_csv": str(tables / "docking_results.csv"),
        "docking_summary_csv": str(tables / "docking_summary.csv"),
    }


def _receptor_that_fails(pipe_gnina, monkeypatch, step_folder, tmp_path):
    """Run the docking stage with a receptor whose PDBQT conversion raises, as Vina's preparation does."""
    def refuse(src, dst):
        raise RuntimeError("receptor has no atoms")

    monkeypatch.setattr(pipe_gnina.vina_backend, "receptor_pdbqt", refuse)
    pdb = tmp_path / "prot.pdb"
    pdb.write_text("END\n")
    conformers = [{"ligand_id": "lig", "conformer_id": 0, "sdf_file": str(tmp_path / "lig.sdf"), "energy": ""}]
    _poses, failures = pipe_gnina.run_docking(
        {"prot": {"protein": str(pdb), "crystal_ligand": None}}, conformers, _vina_docking_config(step_folder))
    return failures


class TestDockingFailuresNameTheirStep:

    def test_a_vina_receptor_failure_names_the_vina_step(self, monkeypatch, tmp_path):
        """`removed_by` was the literal "Gnina", which no Vina step is, so the completion check read the row as one an upstream tool had left behind."""
        step = tmp_path / "004_Vina"
        failures = _receptor_that_fails(_load("pipe_gnina_under_test", "pipe_gnina.py"), monkeypatch, step, tmp_path)
        assert [(row["id"], row["removed_by"], row["kind"]) for row in failures] == [
            ("prot", "004_Vina", "failure")]

    def test_the_completion_check_does_not_excuse_it(self, monkeypatch, tmp_path):
        """The failed pair's file is absent, and a `failure` of this step must leave that a FAILED step instead of an excused one."""
        import pandas as pd

        checker = _load("pipe_check_completion_under_test", "pipe_check_completion.py")
        step = tmp_path / "004_Vina"
        failures = _receptor_that_fails(_load("pipe_gnina_under_test", "pipe_gnina.py"), monkeypatch, step, tmp_path)
        missing = step / "tables" / "missing.csv"
        pd.DataFrame(failures, columns=["id", "removed_by", "kind", "cause"]).to_csv(missing, index=False)
        structures = step / "structures"
        structures.mkdir()
        expected = {
            "output_folder": str(step),
            "tables": {"missing": {"path": str(missing)}},
            "structures": {"name": "structures", "format": "pdb", "ids": ["prot+lig"],
                           "files": [str(structures / "<id>.pdb")],
                           "map_table": str(structures / "structures_map.csv")},
        }
        all_exist, _by_category = checker.check_expected_outputs(expected, "Vina", str(step))
        assert not all_exist

    def test_an_upstream_row_is_still_excused(self, tmp_path):
        """Only this step's own rows lose the excuse; a row an earlier step stamped keeps it."""
        checker = _load("pipe_check_completion_under_test", "pipe_check_completion.py")
        missing = tmp_path / "missing.csv"
        missing.write_text("id,removed_by,kind,cause\nprot,002_Boltz2,failure,upstream\nother,004_Vina,failure,own\n")
        excused = checker._load_expected_missing_ids(
            {"tables": {"missing": {"path": str(missing)}}}, "004_Vina", "Vina")
        assert excused == ["prot"]


# ── USalign and TMalign: a row id is its pair with a mode appended ────────────

class TestPairAlignmentRowIds:

    MODES = [{"label": "sequential", "args": []}, {"label": "non_sequential", "args": ["-mm", "5"]}]

    def _align(self, pairs):
        from biopipelines import pair_alignment as pa

        def never(cmd):
            raise AssertionError("an aligner ran before the row ids were checked")

        files = {sid: f"{sid}.pdb" for pair in pairs for sid in pair}
        return pa.align_pairs(pairs, files, set(files), self.MODES, "USalign", "004_USalign",
                              run=never, log=lambda message: None)

    def test_two_rows_that_would_share_an_id_are_refused_before_anything_runs(self):
        """A mode name holds the `_` that separates an id's suffix, so (a, b_non) in `sequential` and (a, b) in `non_sequential` were both `a+b_non_sequential`."""
        with pytest.raises(ValueError, match=r"both have the id 'a\+b_non_sequential'"):
            self._align([("a", "b_non"), ("a", "b")])

    def test_no_pairs_is_no_rows(self):
        rows, missing = self._align([])
        assert rows == [] and missing == []


# ── PocketGen and PLACER: the pair id comes from the shared rule, its parts from provenance ──

def _stream(tmp_path, name, ids, fmt, with_files):
    from biopipelines.datastream import DataStream, create_map_table

    table = tmp_path / f"{name}_map.csv"
    files = [str(tmp_path / "<id>.pdb")] if with_files else []
    extra = {"code": ["LIG"] * len(ids), "smiles": ["C"] * len(ids)} if name == "compounds" else None
    create_map_table(str(table), ids=list(ids), files=files or None, additional_columns=extra)
    return DataStream(name=name, ids=list(ids), files=files, map_table=str(table), format=fmt)


# Ids as upstream tools make them: a protein that is itself a pair (`a+b`), a ligand with `+` and one with `_`.
PROTEINS = ["a+b", "c"]
LIGANDS = ["l+1", "m_2"]
PAIRS = {
    "a+b+l+1": ("a+b", "l+1"), "a+b+m_2": ("a+b", "m_2"), "c+l+1": ("c", "l+1"), "c+m_2": ("c", "m_2"),
}


def _run(script, *args):
    return subprocess.run([sys.executable, os.path.join(SCRIPTS, script), *args], capture_output=True, text=True)


class TestPocketGenPairs:

    HETATM_SCAFFOLD = (
        "ATOM      1  N   GLY A   1       0.000   0.000   0.000  1.00  0.00           N\n"
        "ATOM      2  CA  GLY A   1       1.458   0.000   0.000  1.00  0.00           C\n"
        "HETATM    5  C1  LIG A   2       5.000   5.000   5.000  1.00  0.00           C\n"
        "END\n")

    def _stage(self, tmp_path, proteins, ligands):
        pytest.importorskip("rdkit")
        for protein in proteins:
            (tmp_path / f"{protein}.pdb").write_text(self.HETATM_SCAFFOLD)
        structures = _stream(tmp_path, "structures", proteins, "pdb", True)
        compounds = _stream(tmp_path, "compounds", ligands, "csv", False)
        structures.save_json(str(tmp_path / "structures.json"))
        compounds.save_json(str(tmp_path / "ligand.json"))
        staging = tmp_path / "staging"
        result = _run("pipe_pocketgen_stage.py", "--structures-json", str(tmp_path / "structures.json"),
                      "--ligand-json", str(tmp_path / "ligand.json"), "--staging-folder", str(staging),
                      "--extras-dir", str(tmp_path / "extras"))
        return staging, result, structures, compounds

    def _postprocess(self, tmp_path, staging, samples=2, without=()):
        for pair in sorted(d for d in os.listdir(staging) if (staging / d).is_dir()):
            if pair in without:
                continue
            for sample in range(samples):
                (staging / pair / f"{sample}_whole_relaxed.pdb").write_text(
                    "ATOM      1  N   GLY A   1       0.000   0.000   0.000  1.00  0.00           N\nEND\n")
        out = tmp_path / "step" / "structures"
        result = _run("pipe_pocketgen_postprocess.py", "--staging-folder", str(staging),
                      "--structures-folder", str(out), "--sequences-folder", str(tmp_path / "step" / "sequences"),
                      "--structures-map", str(out / "structures_map.csv"),
                      "--sequences-csv", str(tmp_path / "step" / "sequences" / "sequences.csv"),
                      "--missing-csv", str(tmp_path / "step" / "tables" / "missing.csv"))
        return result, out

    def test_each_row_names_the_scaffold_and_ligand_it_was_made_from(self, tmp_path):
        """The post-processing split a pair id on its first `+`, so a scaffold id that holds one (`a+b`) gave `a` and `b+l+1`."""
        import pandas as pd

        staging, staged, _s, _c = self._stage(tmp_path, PROTEINS, LIGANDS)
        assert staged.returncode == 0, staged.stderr
        post, out = self._postprocess(tmp_path, staging)
        assert post.returncode == 0, post.stderr
        table = pd.read_csv(out / "structures_map.csv", dtype=str, keep_default_na=False)
        got = {row["id"]: (row["structures.id"], row["compounds.id"]) for _, row in table.iterrows()}
        assert got == {f"{pair}_{k}": parts for pair, parts in PAIRS.items() for k in (0, 1)}
        sequences = pd.read_csv(tmp_path / "step" / "sequences" / "sequences.csv", dtype=str, keep_default_na=False)
        assert {row["id"]: (row["structures.id"], row["compounds.id"]) for _, row in sequences.iterrows()} == got

    def test_the_ids_written_are_covered_by_the_ids_declared(self, tmp_path, local_config, isolated_cwd):
        import pandas as pd
        from biopipelines import id_patterns
        from biopipelines.pocketgen import PocketGen

        staging, staged, structures, compounds = self._stage(tmp_path, PROTEINS, LIGANDS)
        assert staged.returncode == 0, staged.stderr
        tool = PocketGen(structures=structures, ligand=compounds)
        tool.configure_inputs({"PocketGen": str(tmp_path / "repo")})
        declared = list(tool.get_output_files()["structures"].ids)
        post, out = self._postprocess(tmp_path, staging)
        assert post.returncode == 0, post.stderr
        written = list(pd.read_csv(out / "structures_map.csv", dtype=str, keep_default_na=False)["id"])
        assert len(written) == len(set(written)) == 8
        assert id_patterns.select_ids(declared, written) == written

    def test_a_pair_without_samples_is_a_failure_row_under_its_own_id(self, tmp_path):
        import pandas as pd

        staging, staged, _s, _c = self._stage(tmp_path, PROTEINS, LIGANDS)
        assert staged.returncode == 0, staged.stderr
        post, _out = self._postprocess(tmp_path, staging, without=("c+m_2",))
        assert post.returncode == 0, post.stderr
        missing = pd.read_csv(tmp_path / "step" / "tables" / "missing.csv", dtype=str, keep_default_na=False)
        assert [(row["id"], row["kind"]) for _, row in missing.iterrows()] == [("c+m_2", "failure")]

    def test_two_pairs_that_join_to_one_id_are_refused_before_anything_is_staged(self, tmp_path):
        staging, staged, _s, _c = self._stage(tmp_path, ["a+b", "a"], ["c", "b+c"])
        assert staged.returncode != 0
        assert "a+b+c" in staged.stderr
        assert not staging.exists() or not any((staging / d).is_dir() for d in os.listdir(staging))

    def test_neither_script_joins_or_splits_an_id_itself(self):
        for name in ("pipe_pocketgen_stage.py", "pipe_pocketgen_postprocess.py"):
            text = _script_text(name)
            assert not re.search(r"""f['"]\{\w+\}\+\{""", text), name
            assert not re.search(r"""split\(['"]\+['"]""", text), name
            assert not re.search(r"""f['"]\{\w*pair_id\}_\{""", text), name
