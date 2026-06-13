"""Config-time wiring tests for the Pocket Validation example and the two
framework changes it relies on:

  1. Selection now emits an ``n_residues`` column (count of residues in the
     resulting selection string), so a Selection result can be used directly
     as a size metric.
  2. DistanceSelector / SASA / Contacts / PoseChange take their ligand
     reference as a ``Ligand`` / compounds-stream only (per the Ligand
     Contract), resolving the residue code from the compounds map_table at
     runtime. A bare code string is rejected.

These exercise configuration-time wiring + emitted bash only (no GPU tool is
run), matching the suite's scope.
"""
from __future__ import annotations

import glob
import os

import pandas as pd


def _tool_script(pipeline_sh: str, tool_name: str) -> str:
    """Read the per-tool ``NNN_<Tool>.sh`` emitted next to pipeline.sh."""
    runtime_dir = os.path.dirname(pipeline_sh)
    matches = sorted(glob.glob(os.path.join(runtime_dir, f"*_{tool_name}.sh")))
    assert matches, f"no {tool_name} script in {runtime_dir}"
    # Scripts may carry non-ASCII (e.g. the "Å" in the distance echo) written
    # in the platform default encoding; decode leniently for marker checks.
    with open(matches[0], "rb") as f:
        return f.read().decode("utf-8", errors="replace")


# ── Selection n_residues column ───────────────────────────────────────────────

def test_selection_declares_n_residues_column(
    local_config, isolated_cwd, new_pipeline,
):
    """Selection's `selections` table declares id | selection | n_residues."""
    from biopipelines.mock import Mock
    from biopipelines.selection import Selection

    pipeline = new_pipeline("sel_n_residues")
    with pipeline:
        m = Mock(
            ids=["s1"],
            streams={"structures": {"format": "pdb", "file": "<id>.pdb"}},
            tables={"sel": {"columns": ["within"], "fill": {"within": "A1-5"}}},
        )
        sel = Selection(Selection.add(m.tables.sel.within), structures=m)
        pipeline.save()

    assert sel.tables.selections.info.columns == ["id", "selection", "n_residues"]


def test_selection_n_residues_counts_residues(tmp_path):
    """pipe_selection writes n_residues = number of residues in the selection.

    Drives pipe_selection.py directly: add "A2-5" then expand(1) -> "A1-6",
    which is 6 residues.
    """
    import json
    import os
    import subprocess
    import sys

    # Minimal CA-only PDB, chain A residues 1..10.
    lines = [
        f"ATOM  {i:5d}  CA  ALA A{i:4d}      0.000   0.000   {i:.3f}  1.00  0.00           C"
        for i in range(1, 11)
    ] + ["END"]
    pdb = tmp_path / "p.pdb"
    pdb.write_text("\n".join(lines) + "\n")

    ds = {"name": "structures", "ids": ["p"], "files": [str(pdb)],
          "map_table": "", "format": "pdb"}
    ds_json = tmp_path / "structures.json"
    ds_json.write_text(json.dumps(ds))

    tbl = tmp_path / "sel.csv"
    pd.DataFrame([{"id": "p", "within": "A2-5"}]).to_csv(tbl, index=False)

    out_csv = tmp_path / "out.csv"
    config = {
        "operations": [
            {"op": "add", "refs": [{"table": str(tbl), "column": "within"}]},
            {"op": "expand", "value": 1, "direction": "nc"},
        ],
        "output_csv": str(out_csv),
        "structures_json": str(ds_json),
    }
    cfg = tmp_path / "config.json"
    cfg.write_text(json.dumps(config))

    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    script = os.path.join(repo_root, "pipe_scripts", "pipe_selection.py")
    r = subprocess.run([sys.executable, script, str(cfg)],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr

    rows = pd.read_csv(out_csv)
    assert list(rows.columns) == ["id", "selection", "n_residues"]
    assert rows.iloc[0]["selection"] == "A1-6"
    assert int(rows.iloc[0]["n_residues"]) == 6


def test_selection_intersect_is_pairwise_not_union(tmp_path):
    """intersect(A, B) keeps current ∩ A ∩ B, not current ∩ (A ∪ B).

    current = A1-6; intersect against A2-5 and A4-8 must leave A4-5 (the
    common residues). A union of the operands would wrongly leave A2-6.
    """
    import json
    import os
    import subprocess
    import sys

    lines = [
        f"ATOM  {i:5d}  CA  ALA A{i:4d}      0.000   0.000   {i:.3f}  1.00  0.00           C"
        for i in range(1, 11)
    ] + ["END"]
    pdb = tmp_path / "p.pdb"
    pdb.write_text("\n".join(lines) + "\n")

    ds = {"name": "structures", "ids": ["p"], "files": [str(pdb)],
          "map_table": "", "format": "pdb"}
    ds_json = tmp_path / "structures.json"
    ds_json.write_text(json.dumps(ds))

    tbl = tmp_path / "sel.csv"
    pd.DataFrame([{"id": "p", "base": "A1-6", "a": "A2-5", "b": "A4-8"}]).to_csv(tbl, index=False)

    out_csv = tmp_path / "out.csv"
    config = {
        "operations": [
            {"op": "add", "refs": [{"table": str(tbl), "column": "base"}]},
            {"op": "intersect", "refs": [
                {"table": str(tbl), "column": "a"},
                {"table": str(tbl), "column": "b"},
            ]},
        ],
        "output_csv": str(out_csv),
        "structures_json": str(ds_json),
    }
    cfg = tmp_path / "config.json"
    cfg.write_text(json.dumps(config))

    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    script = os.path.join(repo_root, "pipe_scripts", "pipe_selection.py")
    r = subprocess.run([sys.executable, script, str(cfg)],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr

    rows = pd.read_csv(out_csv)
    assert rows.iloc[0]["selection"] == "A4-5"
    assert int(rows.iloc[0]["n_residues"]) == 2


def test_selection_empty_intersect_is_rejected(local_config, isolated_cwd, new_pipeline):
    """Selection.intersect() with no operands fails validation, not silently
    clears the running selection."""
    import pytest
    from biopipelines.mock import Mock
    from biopipelines.selection import Selection

    pipeline = new_pipeline("sel_empty_intersect")
    with pytest.raises(ValueError, match="intersect"):
        with pipeline:
            m = Mock(
                ids=["s1"],
                streams={"structures": {"format": "pdb", "file": "<id>.pdb"}},
                tables={"sel": {"columns": ["within"], "fill": {"within": "A1-5"}}},
            )
            Selection(Selection.add(m.tables.sel.within), Selection.intersect(), structures=m)


# ── DistanceSelector ligand reference: string vs compounds stream ─────────────

def test_distance_selector_rejects_string_ligand(
    local_config, isolated_cwd, new_pipeline,
):
    """A bare code string is no longer accepted — ligand must be a stream."""
    import pytest
    from biopipelines.mock import Mock
    from biopipelines.distance_selector import DistanceSelector

    pipeline = new_pipeline("dsel_string")
    with pytest.raises(ValueError, match="ligand"):
        with pipeline:
            m = Mock(ids=["c"], streams={"structures": {"format": "pdb", "file": "<id>.pdb"}})
            DistanceSelector(structures=m, ligand="ATP", distance=5.0)


def test_distance_selector_compounds_stream_resolves_code_at_runtime(
    local_config, isolated_cwd, new_pipeline,
):
    """A Ligand / compounds stream defers the residue code to runtime: the
    script resolves it from the compounds map_table `code` column and builds
    `ligand:$LIG_CODE`."""
    from biopipelines.entities import Ligand
    from biopipelines.mock import Mock
    from biopipelines.distance_selector import DistanceSelector

    pipeline = new_pipeline("dsel_stream")
    with pipeline:
        m = Mock(ids=["c"], streams={"structures": {"format": "pdb", "file": "<id>.pdb"}})
        lig = Ligand(code="LIG")
        ds = DistanceSelector(structures=m, ligand=lig, distance=5.0)
        script_path = pipeline.save()

    content = _tool_script(script_path, "DistanceSelector")
    assert 'LIG_CODE=$(resolve_stream_item' in content
    assert '"code"' in content
    assert '--reference "$LIG_CODE"' in content


def test_distance_selector_requires_a_ligand_reference(
    local_config, isolated_cwd, new_pipeline,
):
    """Neither a string nor a stream -> clear ValueError at construction."""
    import pytest
    from biopipelines.mock import Mock
    from biopipelines.distance_selector import DistanceSelector

    pipeline = new_pipeline("dsel_missing")
    with pytest.raises(ValueError, match="ligand"):
        with pipeline:
            m = Mock(ids=["c"], streams={"structures": {"format": "pdb", "file": "<id>.pdb"}})
            DistanceSelector(structures=m, distance=5.0)


# ── SASA ligand reference: string vs compounds stream ────────────────────────

def test_sasa_rejects_string_ligand(
    local_config, isolated_cwd, new_pipeline,
):
    """A bare code string is no longer accepted — ligand must be a stream."""
    import pytest
    from biopipelines.mock import Mock
    from biopipelines.sasa import SASA

    pipeline = new_pipeline("sasa_string")
    with pytest.raises(ValueError, match="ligand"):
        with pipeline:
            m = Mock(ids=["c"], streams={"structures": {"format": "pdb", "file": "<id>.pdb"}})
            SASA(structures=m, ligand="LIG")


def test_sasa_compounds_stream_resolves_code_at_runtime(
    local_config, isolated_cwd, new_pipeline,
):
    """A Ligand / compounds stream resolves the residue code from the compounds
    map_table at runtime and strips colons for the PyMOL selection."""
    from biopipelines.entities import Ligand
    from biopipelines.mock import Mock
    from biopipelines.sasa import SASA

    pipeline = new_pipeline("sasa_stream")
    with pipeline:
        m = Mock(ids=["c"], streams={"structures": {"format": "pdb", "file": "<id>.pdb"}})
        lig = Ligand(code="LIG")
        SASA(structures=m, ligand=lig)
        script_path = pipeline.save()

    content = _tool_script(script_path, "SASA")
    assert 'LIG_CODE_RAW=$(resolve_stream_item' in content
    assert '"code"' in content
    assert '--ligand "$LIGAND_RESN"' in content
