"""stream_ids defaults to the ids that actually have files.

A filtered Pool/Panda declares every original id but materializes only the survivors,
so looping over the declared list hands a tool paths that do not exist. LigandMPNN did
exactly that and died in prody's parsePDB after one design. Presence is therefore the
default, and a stream with no files at all falls back to its declared ids rather than
being reported empty.
"""

import json
import os
import subprocess
import sys

import pandas as pd
import pytest

RESOLVER = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                        "pipe_scripts", "resolve_stream_ids.py")


def _run(ds_json, valid_set=True):
    cmd = [sys.executable, RESOLVER, ds_json] + (["--valid-set"] if valid_set else [])
    proc = subprocess.run(cmd, capture_output=True, text=True)
    return proc.stdout.split(), proc.returncode


@pytest.fixture
def filtered_stream(tmp_path):
    """Ten declared ids, two survivors on disk — a Panda-filtered structures stream."""
    ids = [f"pose_{i}" for i in range(1, 11)]
    for sid in ("pose_3", "pose_7"):
        (tmp_path / f"{sid}.pdb").write_text(
            "ATOM      1  CA  ALA A   1       0.0   0.0   0.0  1.00  0.00           C\n")
    mp = tmp_path / "map.csv"
    pd.DataFrame({"id": ids,
                  "file_path": [str(tmp_path / f"{i}.pdb") for i in ids]}).to_csv(mp, index=False)
    p = tmp_path / "ds.json"
    p.write_text(json.dumps({"name": "structures", "ids": ids,
                             "files": [str(tmp_path / "<id>.pdb")],
                             "map_table": str(mp), "format": "pdb"}))
    return str(p)


@pytest.fixture
def value_stream(tmp_path):
    """A compounds stream: ids and a code column, no files anywhere."""
    mp = tmp_path / "compounds.csv"
    pd.DataFrame({"id": ["JF646closed"], "code": ["UNL"]}).to_csv(mp, index=False)
    p = tmp_path / "compounds_ds.json"
    p.write_text(json.dumps({"name": "compounds", "ids": ["JF646closed"], "files": [],
                             "map_table": str(mp), "format": "smiles"}))
    return str(p)


def test_filtered_stream_yields_only_survivors(filtered_stream):
    ids, rc = _run(filtered_stream)
    assert rc == 0
    assert ids == ["pose_3", "pose_7"], "must not name ids whose file was filtered away"


def test_declared_ids_still_reachable(filtered_stream):
    ids, rc = _run(filtered_stream, valid_set=False)
    assert rc == 0 and len(ids) == 10


def test_value_stream_falls_back_to_declared_ids(value_stream):
    """No files to check, so presence says nothing — iterate_files would raise here."""
    ids, rc = _run(value_stream)
    assert rc == 0, "a file-less stream must not error under the presence check"
    assert ids == ["JF646closed"]
