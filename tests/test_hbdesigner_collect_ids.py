# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""HBDesigner: the ids its collector writes are the ids the wrapper declares, taken from one suffix pattern.

The wrapper declares `<parent>_<1..top_k>` through `generate_multiplied_ids_pattern` and hands the same pattern to the step that writes the structures map, while the collector built each rank id, and the id of each rank a run did not produce, with its own `f"{struct_id}_{rank}"` and a separate `--top-k`. The collector now takes the wrapper's suffix pattern and reads its ids from `combinatorics.multiplied_id_parents`, so a produced id and an unproduced one are drawn from the set the wrapper declares.

Nothing here runs HBDesigner; the collector is driven on synthetic run folders.
"""

import json
import os
import re
import subprocess
import sys

import pandas as pd
import pytest

from biopipelines import combinatorics as cb
from biopipelines import id_patterns
from biopipelines.datastream import DataStream, create_map_table

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
COLLECT = os.path.join(REPO_ROOT, "pipe_scripts", "pipe_hbdesigner.py")

PDB = ("ATOM      1  N   GLY A   1       0.000   0.000   0.000  1.00  0.00           N\n"
       "ATOM      2  CA  GLY A   1       1.458   0.000   0.000  1.00  0.00           C\n"
       "END\n")

# A design count per input: all of top_k, a partial run, and a run that crashed.
INPUTS = {"p_1": 3, "a+b": 1, "z": 0}
TOP_K = 3
PATTERN = f"<1..{TOP_K}>"


def _collect(tmp_path, pattern=PATTERN):
    exec_root = tmp_path / "execution"
    for struct_id, produced in INPUTS.items():
        run = exec_root / struct_id
        run.mkdir(parents=True)
        names = [f"{struct_id}_HBDes_rank_{rank}.pdb" for rank in range(1, produced + 1)]
        for name in names:
            (run / name).write_text(PDB)
        pd.DataFrame({"Rank": list(range(1, produced + 1)), "Output_PDB": names}).to_csv(
            run / f"{struct_id}_HBDes_stats.csv")
        if produced == 0:
            (run / ".hbdes_failed").write_text("")
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    for struct_id in INPUTS:
        (inputs / f"{struct_id}.pdb").write_text(PDB)
    table = tmp_path / "structures_map.csv"
    files = [str(inputs / "<id>.pdb")]
    create_map_table(str(table), ids=list(INPUTS), files=files)
    DataStream(name="structures", ids=list(INPUTS), files=files, map_table=str(table), format="pdb").save_json(
        str(tmp_path / "structures.json"))
    step = tmp_path / "step"
    result = subprocess.run(
        [sys.executable, COLLECT, "--structures-json", str(tmp_path / "structures.json"),
         "--exec-root", str(exec_root), "--structures-dir", str(step / "structures"),
         "--summary-csv", str(step / "tables" / "summary.csv"),
         "--sequences-csv", str(step / "sequences" / "sequences.csv"),
         "--missing-csv", str(step / "tables" / "missing.csv"), "--suffix-pattern", pattern],
        capture_output=True, text=True)
    return result, step


def _read(path):
    return pd.read_csv(path, dtype=str, keep_default_na=False)


class TestCollector:

    def test_the_ids_it_writes_and_the_ranks_it_reports_missing_are_the_declared_ids(self, tmp_path):
        result, step = _collect(tmp_path)
        assert result.returncode == 0, result.stderr
        declared = {sid: id_patterns.expand_ids(cb.generate_multiplied_ids_pattern([sid], PATTERN))
                    for sid in INPUTS}
        written = _read(step / "sequences" / "sequences.csv")
        missing = _read(step / "tables" / "missing.csv")
        for sid, ids in declared.items():
            produced = INPUTS[sid]
            assert [i for i in written["id"] if i in ids] == ids[:produced], sid
            assert [i for i in missing["id"] if i in ids] == ids[produced:], sid
        every = [i for ids in declared.values() for i in ids]
        assert sorted(list(written["id"]) + list(missing["id"])) == sorted(every), "an id is both or neither"
        assert os.path.exists(step / "structures" / "a+b_1.pdb")

    def test_each_row_names_the_input_it_was_made_from(self, tmp_path):
        result, step = _collect(tmp_path)
        assert result.returncode == 0, result.stderr
        written = _read(step / "sequences" / "sequences.csv")
        assert {row["id"]: row["structures.id"] for _, row in written.iterrows()} == {
            "p_1_1": "p_1", "p_1_2": "p_1", "p_1_3": "p_1", "a+b_1": "a+b"}
        summary = _read(step / "tables" / "summary.csv")
        assert {row["id"]: row["structures.id"] for _, row in summary.iterrows()} == {
            "p_1_1": "p_1", "p_1_2": "p_1", "p_1_3": "p_1", "a+b_1": "a+b"}

    def test_a_missing_rank_is_a_filter_when_the_run_produced_some_and_a_failure_when_it_produced_none(
            self, tmp_path):
        result, step = _collect(tmp_path)
        assert result.returncode == 0, result.stderr
        kinds = {row["id"]: row["kind"] for _, row in _read(step / "tables" / "missing.csv").iterrows()}
        assert kinds == {"a+b_2": "filter", "a+b_3": "filter", "z_1": "failure", "z_2": "failure", "z_3": "failure"}

    def test_the_count_comes_from_the_pattern_alone(self, tmp_path):
        """A pattern for two ranks keeps two designs of the three a run made, and reports no missing rank for it."""
        result, step = _collect(tmp_path, pattern="<1..2>")
        assert result.returncode == 0, result.stderr
        assert list(_read(step / "sequences" / "sequences.csv")["id"]) == ["p_1_1", "p_1_2", "a+b_1"]

    def test_the_wrapper_hands_the_collector_the_pattern_it_declares_with(self, local_config, isolated_cwd,
                                                                            new_pipeline):
        import pathlib
        from biopipelines.hbdesigner import HBDesigner

        pipeline = new_pipeline("hbdesigner_pattern")
        with pipeline:
            stream = DataStream(name="structures", ids=["p_1"], files=["<id>.pdb"], map_table="", format="pdb")
            out = HBDesigner(structures=stream, top_k=4)
            script_path = pipeline.save()
        body = next(pathlib.Path(script_path).parent.glob("*_HBDesigner.sh")).read_text(encoding="utf-8")
        collector = body[body.index("pipe_hbdesigner.py"):]
        collector = collector[:collector.index("\n\n")]
        assert '--suffix-pattern "<1..4>"' in collector
        assert "--top-k" not in collector
        assert list(out.streams.structures.ids) == ["p_1_<1..4>"]

    def test_the_collector_builds_no_id_itself(self):
        with open(COLLECT) as handle:
            text = handle.read()
        assert not re.search(r"""f['"]\{\w+\}_\{""", text)
