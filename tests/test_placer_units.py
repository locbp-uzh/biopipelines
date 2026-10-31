# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""PLACER: a run's id comes from the shared rule and its parents from the manifest the driver writes.

The post-processing split a run id on its first `+` to recover the structure and the ligand, so a structure id that holds one gave the wrong parents in `structures.id` and `compounds.id`. The driver now writes `units.json` and the post-processing looks a run up by exact id.

Nothing here runs PLACER; the driver's composition is called directly and the post-processing is driven on synthetic runs.
"""

import json
import re

import pytest

from tests._helpers import LIGANDS, PAIRS, PROTEINS, load_pipe_script, pipe_script_text, run_pipe_script


class TestPlacerUnits:

    def test_the_manifest_gives_each_run_its_parents_and_its_id_is_the_shared_composition(self):
        driver = load_pipe_script("pipe_placer_driver_under_test", "pipe_placer_driver.py")
        manifest = driver.unit_manifest(PROTEINS, LIGANDS)
        assert manifest == {pair: {"structures.id": s, "compounds.id": l} for pair, (s, l) in PAIRS.items()}
        assert list(manifest) == list(PAIRS), "a run's id is the structure then the ligand, joined by +"

    def test_sidechain_mode_has_the_structure_axis_alone(self):
        driver = load_pipe_script("pipe_placer_driver_under_test", "pipe_placer_driver.py")
        assert driver.unit_manifest(["a+b", "c"], [None]) == {
            "a+b": {"structures.id": "a+b"}, "c": {"structures.id": "c"}}

    def test_two_runs_that_join_to_one_id_are_refused(self):
        driver = load_pipe_script("pipe_placer_driver_under_test", "pipe_placer_driver.py")
        with pytest.raises(ValueError, match=r"a\+b\+c"):
            driver.unit_manifest(["a+b", "a"], ["c", "b+c"])

    def _run_postprocess(self, tmp_path, mode, units, manifest):
        runs = tmp_path / "runs"
        runs.mkdir()
        for unit in units:
            (runs / f"{unit}_model.pdb").write_text(
                "".join(f"MODEL {n}\nATOM      1  N   GLY A   1       0.000   0.000   0.000  1.00  0.00           N\nENDMDL\n"
                        for n in (1, 2)))
            (runs / f"{unit}.csv").write_text("model_idx,plddt\n1,0.5\n2,0.6\n")
        (runs / "units.json").write_text(json.dumps(manifest))
        step = tmp_path / "step"
        args = ["--mode", mode, "--runs-folder", str(runs), "--structures-folder", str(step / "structures"),
                "--structures-map", str(step / "structures" / "structures_map.csv"),
                "--scores-csv", str(step / "tables" / "scores.csv"),
                "--missing-csv", str(step / "tables" / "missing.csv")]
        if mode == "ligand":
            (tmp_path / "ligand.json").write_text(json.dumps({"name": "compounds", "ids": ["l+1", "m_2"], "files": [],
                                                              "map_table": str(tmp_path / "compounds.csv"),
                                                              "format": "csv"}))
            (tmp_path / "compounds.csv").write_text("id,code\nl+1,LIG\nm_2,LIG\n")
            args += ["--ligand-json", str(tmp_path / "ligand.json"),
                     "--compounds-map", str(step / "compounds" / "compounds_map.csv")]
        return run_pipe_script("pipe_placer_postprocess.py", *args), step

    def test_postprocess_reads_the_parents_from_the_manifest_not_from_the_id(self, tmp_path):
        import pandas as pd

        driver = load_pipe_script("pipe_placer_driver_under_test", "pipe_placer_driver.py")
        manifest = driver.unit_manifest(PROTEINS, LIGANDS)
        result, step = self._run_postprocess(tmp_path, "ligand", list(manifest), manifest)
        assert result.returncode == 0, result.stderr
        table = pd.read_csv(step / "structures" / "structures_map.csv", dtype=str, keep_default_na=False)
        assert {row["id"]: (row["structures.id"], row["compounds.id"]) for _, row in table.iterrows()} == {
            f"{pair}_{n}": parts for pair, parts in PAIRS.items() for n in (1, 2)}

    def test_a_run_the_manifest_does_not_know_stops_the_step(self, tmp_path):
        result, _step = self._run_postprocess(tmp_path, "ligand", ["a+b+l+1"], {})
        assert result.returncode != 0
        assert "a+b+l+1" in result.stderr

    def test_neither_script_joins_or_splits_an_id_itself(self):
        for name in ("pipe_placer_driver.py", "pipe_placer_postprocess.py"):
            text = pipe_script_text(name)
            assert not re.search(r"""f['"]\{\w+\}\+\{""", text), name
            assert not re.search(r"""split\(['"]\+['"]""", text), name
            assert not re.search(r"""f['"]\{unit_id\}_\{""", text), name
