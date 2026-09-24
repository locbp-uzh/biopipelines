"""Selection accepts residue selection strings as literals, positionally or inside add/subtract/intersect."""
from __future__ import annotations

import glob
import json
import os
import subprocess
import sys

import pandas as pd
import pytest


def _mock(Mock):
    return Mock(
        ids=["s1", "s2"],
        streams={"structures": {"format": "pdb", "file": "<id>.pdb"}},
        tables={"sel": {"columns": ["within"], "fill": {"within": "A1-5"}}},
    )


def _saved_operations(root) -> list:
    matches = glob.glob(os.path.join(str(root), "**", "*_Selection", "_configuration", "config.json"), recursive=True)
    assert len(matches) == 1, matches
    with open(matches[0]) as f:
        return json.load(f)["operations"]


def test_implicit_literal_is_serialized(local_config, isolated_cwd, new_pipeline):
    from biopipelines.mock import Mock
    from biopipelines.selection import Selection

    pipeline = new_pipeline("sel_literal")
    with pipeline:
        m = _mock(Mock)
        Selection("A10-20+A45", Selection.expand(1), structures=m)
        pipeline.save()

    ops = _saved_operations(isolated_cwd)
    assert ops[0] == {"op": "add", "literals": ["A10-20+A45"]}


def test_literal_mixes_with_table_reference(local_config, isolated_cwd, new_pipeline):
    from biopipelines.mock import Mock
    from biopipelines.selection import Selection

    pipeline = new_pipeline("sel_literal_mixed")
    with pipeline:
        m = _mock(Mock)
        Selection(Selection.add(m.tables.sel.within, "A30"), Selection.subtract("A2"), structures=m)
        pipeline.save()

    ops = _saved_operations(isolated_cwd)
    assert ops[0]["literals"] == ["A30"] and ops[0]["refs"][0]["column"] == "within"
    assert ops[1] == {"op": "subtract", "literals": ["A2"]}


def test_literal_only_without_structures_is_rejected(local_config, isolated_cwd, new_pipeline):
    from biopipelines.selection import Selection

    pipeline = new_pipeline("sel_literal_no_ids")
    with pytest.raises(ValueError, match="structures="):
        with pipeline:
            Selection("A10-20")


@pytest.mark.parametrize("bad", ["A10-", "chainA", "10,20", ""])
def test_malformed_literal_is_rejected(bad):
    from biopipelines.selection import Selection

    with pytest.raises(ValueError, match="not a residue selection"):
        Selection.add(bad)


def test_non_reference_object_is_rejected():
    from biopipelines.selection import Selection

    with pytest.raises(TypeError, match="Selection.subtract"):
        Selection.subtract(42)


def test_literals_apply_to_every_structure_id(tmp_path):
    lines = [
        f"ATOM  {i:5d}  CA  ALA A{i:4d}      0.000   0.000   {i:.3f}  1.00  0.00           C"
        for i in range(1, 11)
    ] + ["END"]
    files = []
    for sid in ("p", "q"):
        pdb = tmp_path / f"{sid}.pdb"
        pdb.write_text("\n".join(lines) + "\n")
        files.append(str(pdb))
    ds_json = tmp_path / "structures.json"
    ds_json.write_text(json.dumps({"name": "structures", "ids": ["p", "q"], "files": files,
                                   "map_table": "", "format": "pdb"}))

    out_csv = tmp_path / "out.csv"
    cfg = tmp_path / "config.json"
    cfg.write_text(json.dumps({
        "operations": [
            {"op": "add", "literals": ["A2-6"]},
            {"op": "subtract", "literals": ["A4"]},
            {"op": "intersect", "literals": ["A1-5"]},
        ],
        "output_csv": str(out_csv),
        "structures_json": str(ds_json),
    }))

    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    script = os.path.join(repo_root, "pipe_scripts", "pipe_selection.py")
    r = subprocess.run([sys.executable, script, str(cfg)], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr

    rows = pd.read_csv(out_csv)
    assert rows["id"].tolist() == ["p", "q"]
    assert rows["selection"].tolist() == ["A2-3+A5", "A2-3+A5"]
    assert rows["n_residues"].tolist() == [3, 3]
