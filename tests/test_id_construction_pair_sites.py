# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Pair ids and the rows that name them, for the tools that take a protein and a ligand.

The repo rule is that an output id and its provenance are produced only by `combinatorics.py`.
Here: the Gnina/Vina `missing` rows named the wrong step, PocketGen and PLACER joined a pair id by
hand and recovered its parts with a split on `+`, and HBDesigner built its rank ids by hand in the
collector while the wrapper declares them through the shared multiplier.

Nothing here runs a model or a cluster job; the pipe scripts are driven on synthetic inputs.
"""

import importlib.util
import json
import os
import re
import subprocess
import sys

import pytest

from biopipelines import combinatorics as cb

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
