"""Load resolves `<id>` file templates through the map_table.

A `files` entry containing `<id>` is a pattern, not a path — statting it always
fails, which used to report every file-based stream as missing regardless of what
was on disk. These tests pin the four distinct outcomes: a template resolves
against the map_table, ids narrow to the rows the producer actually wrote, a
genuinely absent file is still reported, and an unresolvable template is
reported as such rather than as a missing file.
"""

import csv
import json
import os


def _write_loadable(tool_folder, declared_ids, produced_ids, *,
                    with_map=True, delete_files=()):
    """Build a tool folder with a templated structures stream and its map_table."""
    stream_folder = os.path.join(tool_folder, "structures")
    os.makedirs(stream_folder, exist_ok=True)

    for sid in produced_ids:
        if sid not in delete_files:
            with open(os.path.join(stream_folder, f"{sid}.pdb"), "w") as f:
                f.write("ATOM\n")

    ds = {
        "name": "structures",
        "ids": [declared_ids],
        "files": [os.path.join(stream_folder, "<id>.pdb")],
        "map_table": "",
        "format": "pdb",
        "metadata": {},
    }

    if with_map:
        map_table = os.path.join(stream_folder, "structures_map.csv")
        with open(map_table, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["id", "file", "value"])
            for sid in produced_ids:
                writer.writerow([sid, os.path.join(stream_folder, f"{sid}.pdb"), ""])
        ds["map_table"] = map_table

    with open(os.path.join(tool_folder, ".expected_outputs.json"), "w") as f:
        json.dump({
            "tool_name": "Mock",
            "tool_class": "Mock",
            "output_structure": {
                "structures": ds,
                "tables": {},
                "output_folder": tool_folder,
            },
        }, f)

    return tool_folder


def test_template_is_resolved_not_statted(tmp_path, record_case):
    """A `<id>.pdb` template with all files present reports nothing missing."""
    from biopipelines.load import Load

    folder = _write_loadable(str(tmp_path / "001_Mock"),
                             "design_<1..5>",
                             [f"design_{i}" for i in range(1, 6)])
    loaded = Load(folder)

    record_case(input="files=['<id>.pdb'], 5/5 present",
                expected=([], []),
                actual=(loaded.missing_files, loaded.unresolved_streams))
    assert loaded.missing_files == []
    assert loaded.unresolved_streams == []


def test_ids_narrow_to_map_table_rows(tmp_path, record_case):
    """Declared ids are a prediction; the map_table rows are what exists."""
    from biopipelines.load import Load

    produced = [f"design_{i}" for i in range(1, 291)]
    folder = _write_loadable(str(tmp_path / "001_Mock"), "design_<1..300>", produced)

    structures = Load(folder).get_output_files()["structures"]

    record_case(input="declared design_<1..300>, map_table has 290 rows",
                expected=(290, 290, True),
                actual=(len(structures.ids), len(structures.files),
                        all("<id>" not in f for f in structures.files)))
    assert structures.ids == produced
    assert len(structures.files) == 290
    assert all("<id>" not in f for f in structures.files)


def test_absent_file_still_reported(tmp_path, record_case):
    """A row the map_table lists but disk lacks is a genuine missing file."""
    from biopipelines.load import Load

    produced = [f"design_{i}" for i in range(1, 11)]
    folder = _write_loadable(str(tmp_path / "001_Mock"), "design_<1..10>", produced,
                             delete_files={"design_1", "design_2", "design_3"})
    loaded = Load(folder)

    record_case(input="map_table lists 10, 3 pdbs deleted",
                expected=(3, 0),
                actual=(len(loaded.missing_files), len(loaded.unresolved_streams)))
    assert len(loaded.missing_files) == 3
    assert loaded.unresolved_streams == []


def test_unresolved_template_is_distinct_from_absent(tmp_path, record_case):
    """No map_table behind a template means unresolved, not missing."""
    from biopipelines.load import Load

    folder = _write_loadable(str(tmp_path / "001_Mock"), "design_<1..10>",
                             [f"design_{i}" for i in range(1, 11)], with_map=False)
    loaded = Load(folder)

    record_case(input="files=['<id>.pdb'] with no map_table",
                expected=(0, 1, "structures"),
                actual=(len(loaded.missing_files), len(loaded.unresolved_streams),
                        loaded.unresolved_streams[0][0]))
    assert loaded.missing_files == []
    assert len(loaded.unresolved_streams) == 1
    assert loaded.unresolved_streams[0][0] == "structures"


def test_config_display_counts_expanded_items(tmp_path, record_case):
    """A compact pattern is one list entry standing for many items."""
    from biopipelines.load import Load

    produced = [f"design_{i}" for i in range(1, 291)]
    folder = _write_loadable(str(tmp_path / "001_Mock"), "design_<1..300>", produced)

    lines = Load(folder).get_config_display()
    structures_line = next(l for l in lines if l.startswith("Loaded structures:"))

    record_case(input="declared design_<1..300>, 290 map_table rows",
                expected="Loaded structures: 290",
                actual=structures_line)
    assert structures_line == "Loaded structures: 290"


def test_validate_files_false_propagates_declared_ids(tmp_path, record_case):
    """With validation off, ids are not reconciled — the disclosed behaviour."""
    from biopipelines.load import Load

    folder = _write_loadable(str(tmp_path / "001_Mock"), "design_<1..10>",
                             [f"design_{i}" for i in range(1, 6)])

    structures = Load(folder, validate_files=False).get_output_files()["structures"]

    record_case(input="validate_files=False, declared design_<1..10>",
                expected=["design_<1..10>"],
                actual=structures.ids)
    assert structures.ids == ["design_<1..10>"]
