# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Runtime tests for the Consensus resi-csv aggregator (pipe_consensus.py)."""

import os
import sys
import json
import subprocess

import pandas as pd
import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PIPE = os.path.join(REPO_ROOT, "pipe_scripts", "pipe_consensus.py")


def _write_pose(directory, pid, rows):
    path = os.path.join(directory, f"{pid}.csv")
    pd.DataFrame(rows, columns=["id", "chain", "resi", "distance"]).to_csv(path, index=False)
    return path


def _run(tmp_path, poses, group_ids, operations):
    files = {pid: _write_pose(tmp_path, pid, rows) for pid, rows in poses.items()}
    stream = {"name": "distances", "ids": list(poses.keys()),
              "files": [files[p] for p in poses], "map_table": "", "format": "resi-csv"}
    groups = {"name": "structures", "ids": list(group_ids),
              "files": [], "map_table": "", "format": "csv"}
    sj = os.path.join(tmp_path, "stream.json")
    gj = os.path.join(tmp_path, "groups.json")
    json.dump(stream, open(sj, "w"))
    json.dump(groups, open(gj, "w"))
    out_dir = os.path.join(tmp_path, "out")
    cfg = {"stream_json": sj, "groups_json": gj, "operations": operations,
           "consensus_dir": out_dir,
           "consensus_map_csv": os.path.join(out_dir, "consensus_map.csv")}
    cj = os.path.join(tmp_path, "config.json")
    json.dump(cfg, open(cj, "w"))
    r = subprocess.run([sys.executable, PIPE, cj], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    return out_dir


def _read(out_dir, gid):
    path = os.path.join(out_dir, f"{gid}.csv")
    assert os.path.exists(path), (
        f"consensus file for group '{gid}' not found at {path}; "
        f"wrote: {sorted(os.listdir(out_dir)) if os.path.isdir(out_dir) else '(no out dir)'}"
    )
    return pd.read_csv(path)


def test_fraction_per_group(tmp_path):
    tmp_path = str(tmp_path)
    poses = {
        "x1_1": [("x1_1", "A", 10, 3.0), ("x1_1", "A", 11, 9.0)],
        "x1_2": [("x1_2", "A", 10, 4.0), ("x1_2", "A", 11, 2.0)],
        "x1_3": [("x1_3", "A", 10, 8.0), ("x1_3", "A", 11, 9.0)],
        "x2_1": [("x2_1", "A", 10, 12.0), ("x2_1", "A", 11, 3.0)],
        "x2_2": [("x2_2", "A", 10, 11.0), ("x2_2", "A", 11, 4.0)],
        "x2_3": [("x2_3", "A", 10, 10.0), ("x2_3", "A", 11, 2.0)],
    }
    out = _run(tmp_path, poses, ["x1", "x2"],
               [{"type": "fraction", "params": {"predicate": "distance<=6", "name": "frequency"}}])

    # One consensus file per GROUP id (the Consensus contract), not per pose.
    x1 = _read(out, "x1").set_index("resi")
    assert x1.loc[10, "frequency"] == pytest.approx(2 / 3)
    assert x1.loc[11, "frequency"] == pytest.approx(1 / 3)
    assert (x1["n_group"] == 3).all()

    x2 = _read(out, "x2").set_index("resi")
    assert x2.loc[10, "frequency"] == pytest.approx(0.0)
    assert x2.loc[11, "frequency"] == pytest.approx(1.0)


def test_group_aggregates_all_members(tmp_path):
    """Both poses of a group fold into the single group consensus file."""
    tmp_path = str(tmp_path)
    poses = {
        "x1_1": [("x1_1", "A", 10, 3.0)],
        "x1_2": [("x1_2", "A", 10, 4.0)],
    }
    out = _run(tmp_path, poses, ["x1"],
               [{"type": "fraction", "params": {"predicate": "distance<=6", "name": "frequency"}}])
    x1 = _read(out, "x1").set_index("resi")
    assert x1.loc[10, "frequency"] == pytest.approx(1.0)  # both <= 6
    assert (x1["n_group"] == 2).all()


def test_min_and_count_ops(tmp_path):
    tmp_path = str(tmp_path)
    poses = {
        "x1_1": [("x1_1", "A", 10, 3.0)],
        "x1_2": [("x1_2", "A", 10, 7.0)],
    }
    out = _run(tmp_path, poses, ["x1"],
               [{"type": "min", "params": {"column": "distance", "name": "distance_min"}},
                {"type": "count", "params": {"name": "count"}}])
    row = _read(out, "x1").set_index("resi").loc[10]
    assert row["distance_min"] == pytest.approx(3.0)
    assert row["count"] == 2


def test_one_global_group(tmp_path):
    """A single group id collapses all poses into one consensus."""
    tmp_path = str(tmp_path)
    poses = {
        "prot_1": [("prot_1", "A", 10, 3.0)],
        "prot_2": [("prot_2", "A", 10, 9.0)],
    }
    out = _run(tmp_path, poses, ["prot"],
               [{"type": "fraction", "params": {"predicate": "distance<=6", "name": "frequency"}}])
    assert _read(out, "prot").set_index("resi").loc[10, "frequency"] == pytest.approx(0.5)
