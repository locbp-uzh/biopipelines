"""GEMS: a pair that gets no prediction keeps its row of `affinity`, with a NaN `pkd_pred`.

The `affinity` table dropped the row of a pair GEMS did not predict and recorded it only in `missing`, so a `Panda`
merge of `affinity` against the docking table that fed it silently came back shorter. A result table has one row for
every id that entered; the measurement that is unavailable is NaN. The step still exits non-zero when no pair at all
got a prediction, which `rows` being empty used to decide.
"""

import importlib.util
import math
import os
import subprocess
import sys

import pytest

pd = pytest.importorskip("pandas")

from biopipelines.datastream import DataStream

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


@pytest.fixture
def gems(monkeypatch):
    spec = importlib.util.spec_from_file_location("pipe_gems_failed_pair", os.path.join(REPO_ROOT, "pipe_scripts", "pipe_gems.py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    monkeypatch.setattr(module, "write_ligand_sdf", lambda src, dst, smiles: open(dst, "w").write("sdf\n"))
    return module


def _meta(gems, tmp_path):
    proteins = []
    for pid in ("p1", "p2", "p3"):
        path = tmp_path / f"{pid}.pdb"
        path.write_text("ATOM      1  N   ALA A   1       0.0   0.0   0.0\n")
        proteins.append((pid, str(path)))
    return gems.stage_pairs(proteins, [("l1", str(tmp_path / "l1.sdf"))], str(tmp_path / "pairs"), {})


def test_a_pair_without_a_prediction_keeps_a_nan_row_beside_its_failure_row(gems, tmp_path):
    meta = _meta(gems, tmp_path)
    predictions = {meta[0]["stage_name"]: 5.0}

    rows, missing = gems.affinity_rows(meta, predictions, step_id="004")

    assert [r["id"] for r in rows] == ["p1+l1", "p2+l1", "p3+l1"], "every pair that entered has a row"
    assert rows[0]["pkd_pred"] == 5.0
    assert all(math.isnan(r["pkd_pred"]) for r in rows[1:])
    assert [(r["structures.id"], r["ligands.id"]) for r in rows] == [("p1", "l1"), ("p2", "l1"), ("p3", "l1")]
    assert [(m["id"], m["kind"], m["removed_by"]) for m in missing] == [
        ("p2+l1", "failure", "004"), ("p3+l1", "failure", "004")]


def _run_step(gems, monkeypatch, tmp_path, predicted):
    """main() on three proteins and one ligand; `predicted` are the staged names GEMS answers for."""
    for pid in ("p1", "p2", "p3"):
        (tmp_path / f"{pid}.pdb").write_text("ATOM      1  N   ALA A   1       0.0   0.0   0.0\n")
    (tmp_path / "l1.sdf").write_text("sdf\n")
    DataStream(name="structures", ids=["p1", "p2", "p3"], files=[str(tmp_path / "<id>.pdb")], format="pdb",
               map_table=str(tmp_path / "structures_map.csv")).save_json(str(tmp_path / "structures.json"))
    DataStream(name="compounds", ids=["l1"], files=[str(tmp_path / "<id>.sdf")], format="sdf",
               map_table=str(tmp_path / "compounds_map.csv")).save_json(str(tmp_path / "compounds.json"))
    repo = tmp_path / "gems_repo"
    repo.mkdir()

    def run(cmd, **kwargs):
        if cmd[1] == "GEMS_dataprep_workflow.py":
            (repo / "pairs_dataset.pt").write_text("dataset")
        else:
            names = [f"pair{i}" for i in predicted]
            pd.DataFrame({"id": names, "y_pred": [6.0 + i for i in predicted]}).to_csv(
                repo / "pairs_dataset_predictions.csv", index=False)
        return subprocess.CompletedProcess(cmd, 0, "", "")

    monkeypatch.setattr(gems.subprocess, "run", run)
    affinity = tmp_path / "tables" / "004_GEMS" / "affinity.csv"
    missing = tmp_path / "tables" / "004_GEMS" / "missing.csv"
    monkeypatch.setattr(sys, "argv", [
        "pipe_gems.py", "--structures-json", str(tmp_path / "structures.json"),
        "--compounds-json", str(tmp_path / "compounds.json"), "--gems-repo", str(repo),
        "--scratch-dir", str(tmp_path / "scratch"), "--affinity-csv", str(affinity), "--missing-csv", str(missing)])
    code = 0
    try:
        gems.main()
    except SystemExit as stopped:
        code = stopped.code
    return code, pd.read_csv(affinity), pd.read_csv(missing, dtype=str, keep_default_na=False)


def test_a_step_with_one_failed_pair_writes_every_row_and_exits_zero(gems, monkeypatch, tmp_path):
    code, affinity, missing = _run_step(gems, monkeypatch, tmp_path, predicted=[0, 2])

    assert code == 0
    assert list(affinity["id"]) == ["p1+l1", "p2+l1", "p3+l1"]
    assert affinity["pkd_pred"].isna().tolist() == [False, True, False]
    assert list(missing["id"]) == ["p2+l1"]


def test_a_step_in_which_no_pair_got_a_prediction_exits_non_zero(gems, monkeypatch, tmp_path):
    code, affinity, missing = _run_step(gems, monkeypatch, tmp_path, predicted=[])

    assert code == 1
    assert list(affinity["id"]) == ["p1+l1", "p2+l1", "p3+l1"], "the NaN rows are written before the step fails"
    assert affinity["pkd_pred"].isna().all()
    assert len(missing) == 3
