"""A moved run must rebase the paths INSIDE a map_table, not just the pointer to it.

`.expected_outputs.json` is rebased when a step folder moves, but a map_table is a
CSV on disk and its `file` column routinely points at a sibling step of the same
run. Rebasing only the JSON leaves those rows holding the old absolute paths, so a
run archived to another machine resolves every id to a file that is not there —
which is how a staged pose set failed with "Upstream file not found" against paths
on a cluster it was no longer running on.
"""

import json
import os

import pytest

pytest.importorskip("pandas")


def _make_run(root, run_name="run_001"):
    """A two-step run: 010_Pool holds the files, 020_Panda selects from them."""
    run = root / run_name
    pool = run / "010_Pool" / "structures"
    pool.mkdir(parents=True)
    panda = run / "020_Panda"
    (panda / "structures").mkdir(parents=True)
    (panda / "tables").mkdir(parents=True)

    ids = ["a", "b"]
    for i in ids:
        (pool / f"pose_{i}.pdb").write_text("ATOM\nEND\n")

    # the map_table points back at the POOL step, not at Panda's own folder
    mt = panda / "tables" / "structures_map.csv"
    mt.write_text(
        "id,file\n" + "".join(f"{i},{pool / f'pose_{i}.pdb'}\n" for i in ids)
    )

    result = {
        "tool_name": "Panda",
        "tool_class": "Panda",
        "output_structure": {
            "output_folder": str(panda),
            "structures": {
                "name": "structures",
                "ids": ids,
                "files": str(panda / "structures" / "<id>.pdb"),
                "map_table": str(mt),
            },
            "tables": {},
        },
    }
    (panda / ".expected_outputs.json").write_text(json.dumps(result))
    return run, panda, mt


def test_map_table_paths_follow_the_moved_run(tmp_path):
    """After moving the whole run, map_table rows must point into the new location."""
    from biopipelines.load import Load

    src = tmp_path / "src"
    src.mkdir()
    run, panda, mt = _make_run(src)

    # move the run somewhere else, as archiving to another filesystem would
    dst = tmp_path / "dst"
    dst.mkdir()
    moved_run = dst / "run_001"
    os.rename(run, moved_run)
    moved_panda = moved_run / "020_Panda"

    loader = Load(path=str(moved_panda), validate_files=False)
    ds = loader.loaded_result["output_structure"]["structures"]
    paths = loader._map_table_file_paths(ds)

    assert paths, "map_table produced no paths"
    for p in paths:
        assert str(dst) in p, f"path was not rebased into the new location: {p}"
        assert os.path.exists(p), f"rebased path does not exist: {p}"


def test_unmoved_run_is_left_alone(tmp_path):
    """No move, no rewriting — the stored paths are already correct."""
    from biopipelines.load import Load

    src = tmp_path / "src"
    src.mkdir()
    run, panda, mt = _make_run(src)

    loader = Load(path=str(panda), validate_files=False)
    ds = loader.loaded_result["output_structure"]["structures"]
    paths = loader._map_table_file_paths(ds)

    assert paths
    for p in paths:
        assert os.path.exists(p)
        assert str(src) in p


def test_paths_outside_the_run_are_not_rewritten(tmp_path):
    """Only the moved run's own prefix is rebased.

    A map_table may legitimately reference a file elsewhere on the filesystem;
    rewriting that would break a path that was never broken.
    """
    from biopipelines.load import Load

    src = tmp_path / "src"
    src.mkdir()
    run, panda, mt = _make_run(src)

    outside = tmp_path / "elsewhere" / "shared.pdb"
    outside.parent.mkdir()
    outside.write_text("ATOM\nEND\n")
    mt.write_text(mt.read_text() + f"c,{outside}\n")

    dst = tmp_path / "dst"
    dst.mkdir()
    os.rename(run, dst / "run_001")

    loader = Load(path=str(dst / "run_001" / "020_Panda"), validate_files=False)
    ds = loader.loaded_result["output_structure"]["structures"]
    paths = loader._map_table_file_paths(ds)

    assert str(outside) in paths, "a path outside the run must be left untouched"


def test_downstream_reads_a_rebased_csv(tmp_path):
    """The fix must reach tools that read the map_table CSV themselves.

    Reading it as data would hand back a stream with no files at all while the
    recovery report still said "recovered N of N via map_table", so it is named
    and refused instead.
    """
    import pandas as pd
    from biopipelines.load import Load

    src = tmp_path / "src"
    src.mkdir()
    run, panda, mt = _make_run(src)
    # A pre-1.4.0 run folder: the stream file lives under the retired `file_path`.
    rows = pd.read_csv(mt).rename(columns={"file": "file_path"})
    rows.to_csv(mt, index=False)

    dst = tmp_path / "dst"
    dst.mkdir()
    os.rename(run, dst / "run_001")
    moved_panda = dst / "run_001" / "020_Panda"

    # validate_files=True (the default): validate_files=False is an explicit
    # "do not resolve files", and then the map's file column is never read.
    with pytest.raises(ValueError) as excinfo:
        Load(path=str(moved_panda))

    message = str(excinfo.value)
    assert "file_path" in message
    assert "'file'" in message
    assert str(mt.name) in message or "map_table" in message.lower()


def test_tool_specific_path_column_is_rebased(tmp_path):
    """A tool's own path column must be rebased, not just file/file_path.

    Gnina uses `best_pose_file`, PyMOL `session_file`, and so on. These name a
    file the tool produced, not the stream's own file, so they keep their names
    and still have to be rebased. Rebasing only the canonical `file` left such a
    map pointing at /shares paths after the run was copied to another cluster,
    so the consumer found nothing and failed downstream of the real cause.
    """
    import pandas as pd
    from biopipelines.load import Load

    src = tmp_path / "src"
    src.mkdir()
    run, panda, mt = _make_run(src)

    # Gnina records a pose path beside the stream file; both must be rebased.
    rows = pd.read_csv(mt)
    rows["best_pose_file"] = rows["file"]
    rows.to_csv(mt, index=False)

    dst = tmp_path / "dst"
    dst.mkdir()
    os.rename(run, dst / "run_001")
    moved = dst / "run_001" / "020_Panda"

    loader = Load(path=str(moved), validate_files=False)
    ds = loader.loaded_result["output_structure"]["structures"]

    on_disk = pd.read_csv(ds["map_table"])
    for v in on_disk["best_pose_file"]:
        assert str(dst) in v, f"best_pose_file column was not rebased: {v}"
        assert os.path.exists(v), f"rebased best_pose_file does not exist: {v}"
