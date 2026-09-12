"""A map_table must survive Load unaltered except for the paths.

Companion to test_load_rebases_map_table.py, which covers rebasing's happy path.
Three defects hid from it because all of its cases use plain-string ids, one
column layout, and validate_files=False:

  * pd.read_csv without dtype=str retyped the map. Since 24d551d made the map
    authoritative for stream ids, "0001" loaded as 1 while its file stayed
    0001.pdb -- silently renaming the id and breaking every downstream join. The
    rebasing path wrote that damage back to disk.
  * keep_default_na left off turned real strings into blanks: a "NA" or "nan"
    cell (a chain id, a flag, a ligand code) came back empty.
  * The prefix test was a bare startswith, so a row pointing at a sibling run
    Job_001_retry was rewritten when Job_001 moved -- corrupting a path that was
    never broken.
"""

import json

import pytest

pytest.importorskip("pandas")

import pandas as pd  # noqa: E402

from biopipelines.load import _rebase_path, _under_prefix  # noqa: E402


# --- the prefix test -------------------------------------------------------

@pytest.mark.parametrize("path,inside", [
    ("/jobs/Job_001/010_Pool/a.pdb", True),
    ("/jobs/Job_001", True),
    ("/jobs/Job_001_retry/010_Pool/a.pdb", False),
    ("/jobs/Job_0010/a.pdb", False),
    ("/jobs/Job_00/a.pdb", False),
    ("/elsewhere/a.pdb", False),
    ("", False),
    (None, False),
    (3, False),
])
def test_prefix_test_matches_whole_components(path, inside):
    assert _under_prefix(path, "/jobs/Job_001") is inside


def test_sibling_run_is_not_rewritten():
    """The corruption case: a path that was never broken must be left alone."""
    sibling = "/jobs/Job_001_retry/010_Pool/pose.pdb"
    assert _rebase_path(sibling, "/jobs/Job_001", "/jobs/Job_001_v2") == sibling


def test_rebase_preserves_posix_separators():
    """Paths recorded on a cluster must not come back Windows-flavoured."""
    out = _rebase_path("/jobs/Job_001/010_Pool/a.pdb", "/jobs/Job_001", "/jobs/Job_001_v2")
    assert out == "/jobs/Job_001_v2/010_Pool/a.pdb"
    assert "\\" not in out


def test_rebase_is_idempotent_when_the_new_prefix_contains_the_old():
    """Job_001 -> Job_001_v2 must not become Job_001_v2_v2 on a second pass."""
    once = _rebase_path("/jobs/Job_001/a.pdb", "/jobs/Job_001", "/jobs/Job_001_v2")
    assert once == "/jobs/Job_001_v2/a.pdb"
    assert _rebase_path(once, "/jobs/Job_001", "/jobs/Job_001_v2") == once


# --- id and value fidelity through a real Load -----------------------------

def _make_run(root, ids, extra_cols=""):
    """A one-step run whose map_table carries the given ids."""
    step = root / "run_001" / "010_Mock"
    (step / "structures").mkdir(parents=True)
    (step / "tables").mkdir(parents=True)
    for i in ids:
        (step / "structures" / f"{i}.pdb").write_text("ATOM\nEND\n")

    header = "id,file" + (f",{extra_cols}" if extra_cols else "")
    rows = []
    for i in ids:
        row = f"{i},{step / 'structures' / f'{i}.pdb'}"
        if extra_cols:
            row += ",NA"
        rows.append(row)
    mt = step / "tables" / "structures_map.csv"
    mt.write_text(header + "\n" + "\n".join(rows) + "\n")

    (step / ".expected_outputs.json").write_text(json.dumps({
        "tool_name": "Mock",
        "tool_class": "Mock",
        "output_structure": {
            "output_folder": str(step),
            "structures": {
                "name": "structures",
                "ids": list(ids),
                "files": [str(step / "structures" / "<id>.pdb")],
                "map_table": str(mt),
            },
            "tables": {},
        },
    }))
    return step, mt


def test_zero_padded_ids_survive_the_map_table_read(tmp_path):
    """The headline defect: "0001" must not load as "1"."""
    ids = ["0001", "0002", "0010"]
    step, mt = _make_run(tmp_path, ids)

    from biopipelines.load import Load
    loader = Load.__new__(Load)
    rows = loader._read_map_table({"map_table": str(mt)})

    assert rows is not None
    assert list(rows["id"]) == ids, (
        f"ids were retyped: {list(rows['id'])} != {ids}. A numeric-looking id "
        f"loaded as an int no longer matches the file it names, nor any "
        f"downstream provenance column.")


def test_na_strings_are_not_blanked(tmp_path):
    """"NA" and "nan" are real values in a chain id or flag column."""
    step = tmp_path / "run" / "010_Mock" / "tables"
    step.mkdir(parents=True)
    mt = step / "structures_map.csv"
    mt.write_text("id,file,chain\na,/x/a.pdb,NA\nb,/x/b.pdb,nan\nc,/x/c.pdb,B\n")

    from biopipelines.load import Load
    loader = Load.__new__(Load)
    rows = loader._read_map_table({"map_table": str(mt)})

    assert list(rows["chain"]) == ["NA", "nan", "B"], (
        f"real strings were eaten by NA inference: {list(rows['chain'])}")


def test_rebased_copy_keeps_ids_and_values_verbatim(tmp_path):
    """Rebasing writes to disk, so any retyping it does is permanent."""
    old_root = tmp_path / "old"
    old_root.mkdir()
    mt = old_root / "structures_map.csv"
    mt.write_text(
        "id,file,replicate,flag\n"
        "0001,/old/run/010_Mock/structures/0001.pdb,1,NA\n"
        "0007,/old/run/010_Mock/structures/0007.pdb,2,nan\n"
    )

    from biopipelines.load import Load
    loader = Load.__new__(Load)
    loader._rebased_map_tables = set()
    structure = {"structures": {"name": "structures", "map_table": str(mt)}}
    loader._rebase_map_tables(structure, "/old/run", "/new/run")

    rebased = structure["structures"]["map_table"]
    assert rebased != str(mt), "no rebased copy was written"
    out = pd.read_csv(rebased, dtype=str, keep_default_na=False)

    assert list(out["id"]) == ["0001", "0007"], f"ids retyped on disk: {list(out['id'])}"
    assert list(out["replicate"]) == ["1", "2"], f"ints reformatted: {list(out['replicate'])}"
    assert list(out["flag"]) == ["NA", "nan"], f"strings blanked: {list(out['flag'])}"
    assert list(out["file"]) == ["/new/run/010_Mock/structures/0001.pdb",
                                "/new/run/010_Mock/structures/0007.pdb"]


def test_reading_a_rebased_copy_does_not_rebase_again(tmp_path):
    """The on-disk copy is already correct; a second pass would double it."""
    old_root = tmp_path / "old"
    old_root.mkdir()
    mt = old_root / "structures_map.csv"
    mt.write_text("id,file\na,/jobs/Job_001/010_Mock/a.pdb\n")

    from biopipelines.load import Load
    loader = Load.__new__(Load)
    loader._rebased_map_tables = set()
    structure = {"structures": {"name": "structures", "map_table": str(mt)}}
    # new prefix contains the old one -- the pathological rename
    loader._path_rebase = ("/jobs/Job_001", "/jobs/Job_001_v2")
    loader._rebase_map_tables(structure, "/jobs/Job_001", "/jobs/Job_001_v2")

    rows = loader._read_map_table({"map_table": structure["structures"]["map_table"]})
    assert list(rows["file"]) == ["/jobs/Job_001_v2/010_Mock/a.pdb"], (
        f"path was rebased twice: {list(rows['file'])}")


def test_unwritable_destination_raises_instead_of_keeping_stale_paths(tmp_path,
                                                                     monkeypatch):
    """Silently keeping the old map hands downstream tools dead paths."""
    old_root = tmp_path / "old"
    old_root.mkdir()
    mt = old_root / "structures_map.csv"
    mt.write_text("id,file\na,/old/run/010_Mock/a.pdb\n")

    def _boom(*args, **kwargs):
        raise OSError("read-only file system")

    monkeypatch.setattr(pd.DataFrame, "to_csv", _boom)

    from biopipelines.load import Load
    loader = Load.__new__(Load)
    loader._rebased_map_tables = set()
    structure = {"structures": {"name": "structures", "map_table": str(mt)}}

    with pytest.raises(RuntimeError, match="could not write the rebased map_table"):
        loader._rebase_map_tables(structure, "/old/run", "/new/run")
