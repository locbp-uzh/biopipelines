"""Load must recover whatever a partially-completed step produced.

Every generated tool script exits 0 by design, so a framework-only failure
cannot kill a pipeline that would otherwise run. The price of that decision is
paid here: ``Load`` is the component that has to get results back out of a
half-written output folder. These tests pin that guarantee across the shapes a
dying step actually leaves behind — some ids present and some absent, a
truncated map_table, an absent map_table, a ``_FAILED`` marker, a populated
``missing.csv`` — and in each one ``Load`` must report precisely what it
recovered and what it could not, never raise, and never hand back an empty
stream that reads as success.

The fixture layout mirrors a real ``Mock`` run byte for byte (verified against
``pipe_mock.py`` output): a ``<id>``-templated ``files`` entry, an ``id/file/value``
map_table under the stream folder, and a ``tables`` dict. One test runs the real
``pipe_mock.py`` so the layout stays anchored to the producer.
"""

import csv
import io
import json
import os
import subprocess
import sys

import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
IDS = [f"design_{i}" for i in range(1, 11)]


# ── fixture builders ──────────────────────────────────────────────────────────

def _build(tool_folder, *, declared_ids=None, produced_ids=None,
           map_rows="all", tables=None, suffix=""):
    """A Mock-shaped output folder.

    map_rows: "all" (one row per produced id) | "none" (no map_table file) |
    a list of ids to list | "empty" (header only) | "zero" (zero-byte file).

    suffix: text the producer puts between the id and the extension, so the
    declared template reads ``<id>_best.pdb`` rather than ``<id>.pdb``. Gnina,
    Aggrescan3D and CABSflex all ship such a template.
    """
    declared_ids = list(IDS if declared_ids is None else declared_ids)
    produced_ids = list(declared_ids if produced_ids is None else produced_ids)
    stream_folder = os.path.join(tool_folder, "structures")
    tables_folder = os.path.join(tool_folder, "tables")
    os.makedirs(stream_folder, exist_ok=True)
    os.makedirs(tables_folder, exist_ok=True)

    for sid in produced_ids:
        with open(os.path.join(stream_folder, f"{sid}{suffix}.pdb"), "w") as f:
            f.write("ATOM\n")

    map_table = os.path.join(stream_folder, "structures_map.csv")
    if map_rows == "none":
        map_table_field = map_table  # declared, but never written
    elif map_rows == "zero":
        open(map_table, "w").close()
        map_table_field = map_table
    else:
        listed = [] if map_rows == "empty" else (
            produced_ids if map_rows == "all" else list(map_rows))
        with open(map_table, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["id", "file", "value"])
            for sid in listed:
                w.writerow([sid, os.path.join(stream_folder, f"{sid}{suffix}.pdb"), ""])
        map_table_field = map_table

    table_decls = {}
    for name, rows in (tables or {"scores": [[sid, "1.0"] for sid in produced_ids]}).items():
        path = os.path.join(tables_folder, f"{name}.csv")
        if rows is not None:
            columns = (["id", "removed_by", "kind", "cause"] if name == "missing"
                       else ["id", "score"])
            with open(path, "w", newline="") as f:
                w = csv.writer(f)
                w.writerow(columns)
                for row in rows:
                    w.writerow(row)
        else:
            columns = ["id", "score"]  # declared but never written
        table_decls[name] = {"name": name, "path": path,
                             "columns": columns, "description": ""}

    with open(os.path.join(tool_folder, ".expected_outputs.json"), "w") as f:
        json.dump({
            "tool_name": "Mock",
            "tool_class": "Mock",
            "output_structure": {
                "structures": {
                    "name": "structures",
                    "ids": declared_ids,
                    "files": [os.path.join(stream_folder, f"<id>{suffix}.pdb")],
                    "map_table": map_table_field,
                    "format": "pdb",
                    "metadata": {},
                },
                "tables": table_decls,
                "output_folder": tool_folder,
            },
        }, f)
    return tool_folder


def _delete_files(tool_folder, ids, suffix=""):
    for sid in ids:
        os.remove(os.path.join(tool_folder, "structures", f"{sid}{suffix}.pdb"))


def _truncate_map_mid_row(tool_folder, at_id):
    """Cut the map_table off inside a row, as a killed writer does."""
    path = os.path.join(tool_folder, "structures", "structures_map.csv")
    raw = open(path).read()
    open(path, "w").write(raw[: raw.index(at_id)] + f"{at_id},partial")


def _mark(tool_folder, status):
    parent = os.path.dirname(tool_folder)
    name = os.path.basename(tool_folder)
    open(os.path.join(parent, f"{name.split('_')[0]}_Mock_{status}"), "w").write(
        f"Status: {status}\n")


@pytest.fixture
def loadable_env(local_config):
    """Pin the ConfigManager to the local fixture config.

    ``Load`` is constructed without a Pipeline, so nothing else sets the config
    variant for it and auto-detection raises on any machine whose username no
    committed config claims.
    """
    from biopipelines.config_manager import ConfigManager

    ConfigManager("local")
    yield
    ConfigManager._instance = None
    ConfigManager._config = None
    ConfigManager._variant = None


def _load(folder, **kw):
    from biopipelines.load import Load
    loaded = Load(folder, **kw)
    return loaded, loaded.get_output_files()


# ── the shapes ────────────────────────────────────────────────────────────────

def test_partial_files_keeps_map_ids_and_names_the_absent_paths(loadable_env, tmp_path, record_case):
    """The map is authoritative for ids; absent files are reported, not dropped.

    The producer recorded 10 rows and then lost 4 files. Load keeps the 10 ids
    the map claims (a downstream tool needs them to excuse the gap) while
    stating that 4 of the recovered paths are gone.
    """
    folder = _build(str(tmp_path / "001_Mock"))
    _delete_files(folder, IDS[6:])

    loaded, out = _load(folder)
    rec = loaded.recovery["structures"]

    record_case(input="map lists 10, 6 pdbs on disk",
                expected=(10, 10, 4),
                actual=(len(out["structures"].ids), rec["recovered"],
                        len(rec["absent_files"])))
    assert len(out["structures"].ids) == 10
    assert rec == {"declared": 10, "recovered": 10, "source": "map_table",
                   "unrecovered_ids": [], "absent_files": rec["absent_files"]}
    assert len(rec["absent_files"]) == 4
    assert len(loaded.missing_files) == 4


def test_truncated_map_table_recovers_the_rows_that_were_written(loadable_env, tmp_path, record_case):
    """A writer killed mid-row leaves a short map; Load takes what is there."""
    folder = _build(str(tmp_path / "001_Mock"))
    _truncate_map_mid_row(folder, "design_7")
    _delete_files(folder, IDS[6:])

    loaded, out = _load(folder)
    ids = list(out["structures"].ids)

    record_case(input="map truncated inside the design_7 row",
                expected=(7, IDS[:6]),
                actual=(len(ids), ids[:6]))
    assert ids[:6] == IDS[:6]
    assert len(ids) == 7  # the half-written design_7 row still parses
    assert loaded.recovery["structures"]["declared"] == 10
    assert loaded.recovery["structures"]["recovered"] == 7
    # Its file path was cut off, so it is reported absent rather than trusted.
    assert len(loaded.recovery["structures"]["absent_files"]) == 1


def test_absent_map_table_recovers_by_glob_without_inventing_duplicate_ids(loadable_env, 
        tmp_path, record_case):
    """No map_table: Load globs, and every recovered id is distinct and real.

    Matching in declared order let the substring rules pair ``design_10`` with
    ``design_1.pdb``: the recovery came back listing ``design_1`` twice and
    ``design_10`` not at all, which reads as 7 recovered designs when only 6
    exist. Exact basenames are claimed first and one file now satisfies at most
    one id.
    """
    folder = _build(str(tmp_path / "001_Mock"), map_rows="none")
    _delete_files(folder, IDS[6:])

    loaded, out = _load(folder)
    ids = list(out["structures"].ids)
    rec = loaded.recovery["structures"]

    record_case(input="no map_table, design_1..6 on disk, design_7..10 gone",
                expected=(IDS[:6], ["design_7", "design_8", "design_9", "design_10"]),
                actual=(ids, rec["unrecovered_ids"]))
    assert ids == IDS[:6]
    assert len(ids) == len(set(ids))
    assert [os.path.basename(f) for f in out["structures"].files] == \
        [f"{sid}.pdb" for sid in IDS[:6]]
    assert rec["source"] == "glob"
    assert rec["unrecovered_ids"] == ["design_7", "design_8", "design_9", "design_10"]


def test_absent_map_table_inverts_a_suffixed_template_instead_of_guessing(
        loadable_env, tmp_path, record_case):
    """A template with a suffix must be inverted, not substring-matched.

    Claiming exact basenames first only covers ``<id>.pdb``, where the basename
    IS the id. Gnina declares ``<id>_best.pdb``, Aggrescan3D ``<id>_A3D.csv``:
    no basename ever equals an id, every id fell through to the substring rules,
    and ``design_1`` was handed ``design_10_best.pdb`` — another id's file —
    while ``design_1_best.pdb`` went unused and ``design_10`` was reported
    unrecovered. The declared template already states how the name is built, so
    the id comes back out of it.
    """
    folder = _build(str(tmp_path / "001_Mock"), map_rows="none", suffix="_best")
    _delete_files(folder, IDS[1:9], suffix="_best")  # design_1 and design_10 survive

    loaded, out = _load(folder)
    ds = out["structures"]
    pairs = list(zip(ds.ids, [os.path.basename(f) for f in ds.files]))

    record_case(input="no map_table, template <id>_best.pdb, "
                      "design_1_best.pdb and design_10_best.pdb on disk",
                expected=[("design_1", "design_1_best.pdb"),
                          ("design_10", "design_10_best.pdb")],
                actual=pairs)
    assert pairs == [("design_1", "design_1_best.pdb"),
                     ("design_10", "design_10_best.pdb")]
    assert loaded.recovery["structures"]["unrecovered_ids"] == IDS[1:9]


def test_absent_map_table_keeps_the_declared_id_not_the_basename(
        loadable_env, tmp_path, record_case):
    """Recovery must not rename ids to whatever the filename happens to say.

    Reading the id off the basename turned every ``design_N`` into
    ``design_N_A3D``, which joins against no table and no upstream map_table.
    """
    folder = _build(str(tmp_path / "001_Mock"), map_rows="none", suffix="_A3D")

    loaded, out = _load(folder)
    ids = list(out["structures"].ids)

    record_case(input="no map_table, template <id>_A3D.pdb, all 10 on disk",
                expected=IDS, actual=ids)
    assert ids == IDS


def test_absent_map_table_leaves_another_ids_file_alone(
        loadable_env, tmp_path, record_case):
    """A file the template names for an undeclared id is still that id's.

    ``design_20_best.pdb`` inverts to ``design_20``; the substring rules would
    otherwise let ``design_2`` claim it once the template tier found nothing.
    """
    folder = _build(str(tmp_path / "001_Mock"), declared_ids=["design_2"],
                    map_rows="none", suffix="_best")
    _delete_files(folder, ["design_2"], suffix="_best")
    open(os.path.join(folder, "structures", "design_20_best.pdb"), "w").write("ATOM\n")

    loaded, out = _load(folder)
    ids = list(out["structures"].ids)

    record_case(input="declared design_2, only design_20_best.pdb on disk",
                expected=[], actual=ids)
    assert ids == []
    assert loaded.recovery["structures"]["unrecovered_ids"] == ["design_2"]


def test_empty_map_table_says_it_recovered_nothing(loadable_env, tmp_path, record_case):
    """A step that ran and produced nothing must not read as success.

    An empty stream is indistinguishable from a successful load of zero items,
    so the shortfall is stated outright and ``recovered_nothing`` is True.
    """
    folder = _build(str(tmp_path / "001_Mock"), map_rows="empty")
    _delete_files(folder, IDS)

    loaded, out = _load(folder)
    summary = loaded.recovery_summary()

    record_case(input="map_table header-only, no files produced",
                expected=(0, True, True),
                actual=(len(out["structures"].ids), loaded.recovered_nothing,
                        any("RECOVERED NOTHING" in line for line in summary)))
    assert list(out["structures"].ids) == []
    assert loaded.recovered_nothing is True
    assert any("RECOVERED NOTHING" in line for line in summary)
    assert any("0 of 10" in line for line in summary)
    # A readable map with no rows is an answer, not an unresolvable template.
    assert loaded.unresolved_streams == []


def test_absent_map_table_is_distinct_from_an_empty_one(loadable_env, tmp_path, record_case):
    """"Nothing was produced" and "I cannot tell" must not report the same."""
    empty = _build(str(tmp_path / "e" / "001_Mock"), map_rows="empty")
    _delete_files(empty, IDS)
    absent = _build(str(tmp_path / "a" / "001_Mock"), map_rows="none")
    _delete_files(absent, IDS)

    e_loaded, _ = _load(empty)
    a_loaded, _ = _load(absent)

    record_case(input="header-only map vs absent map, no files either way",
                expected=("map_table", "glob"),
                actual=(e_loaded.recovery["structures"]["source"],
                        a_loaded.recovery["structures"]["source"]))
    assert e_loaded.recovery["structures"]["source"] == "map_table"
    assert e_loaded.unresolved_streams == []
    assert a_loaded.recovery["structures"]["source"] == "glob"
    assert [n for n, _p in a_loaded.unresolved_streams] == ["structures"]


def test_failed_marker_is_reported_and_the_partial_output_still_loads(loadable_env, 
        tmp_path, record_case):
    """A _FAILED step's output is still worth loading — but say where it came from."""
    folder = _build(str(tmp_path / "001_Mock"), map_rows=IDS[:6])
    _delete_files(folder, IDS[6:])
    _mark(folder, "FAILED")

    loaded, out = _load(folder)

    record_case(input="_FAILED marker, 6 of 10 produced",
                expected=("FAILED", 6),
                actual=(loaded.completion_status, len(out["structures"].ids)))
    assert loaded.completion_status == "FAILED"
    assert list(out["structures"].ids) == IDS[:6]
    assert any("FAILED" in line for line in loaded.recovery_summary())
    assert loaded.get_loaded_metadata()["completion_status"] == "FAILED"


def test_completed_marker_is_read_too(loadable_env, tmp_path, record_case):
    folder = _build(str(tmp_path / "001_Mock"))
    _mark(folder, "COMPLETED")
    loaded, _ = _load(folder)
    record_case(input="_COMPLETED marker", expected="COMPLETED",
                actual=loaded.completion_status)
    assert loaded.completion_status == "COMPLETED"


def test_no_marker_reads_as_none_not_as_failure(loadable_env, tmp_path, record_case):
    folder = _build(str(tmp_path / "001_Mock"))
    loaded, _ = _load(folder)
    record_case(input="no completion marker", expected=None,
                actual=loaded.completion_status)
    assert loaded.completion_status is None
    assert any("none found" in line for line in loaded.recovery_summary())


def test_populated_missing_csv_excuses_its_ids(loadable_env, tmp_path, record_case):
    """Ids the producer recorded as lost are dropped and counted as excused."""
    folder = _build(
        str(tmp_path / "001_Mock"), produced_ids=IDS[:8], map_rows="all",
        tables={"scores": [[sid, "1.0"] for sid in IDS[:8]],
                "missing": [[sid, "001_Mock", "failure", "mock"] for sid in IDS[8:]]},
    )

    loaded, out = _load(folder)

    record_case(input="missing.csv lists design_9/10, map lists 8",
                expected=(IDS[:8], ["design_10", "design_9"]),
                actual=(list(out["structures"].ids), loaded.excused_ids))
    assert list(out["structures"].ids) == IDS[:8]
    assert loaded.excused_ids == ["design_10", "design_9"]
    assert any("excused by tables/missing.csv" in line
               for line in loaded.recovery_summary())


def test_declared_missing_csv_that_was_never_written_does_not_raise(loadable_env, 
        tmp_path, record_case):
    """The exact shape a step killed before writing missing.csv leaves behind.

    Raising here denied the user every id that step *did* produce, which is the
    one thing Load exists to prevent.
    """
    folder = _build(
        str(tmp_path / "001_Mock"),
        tables={"scores": [[sid, "1.0"] for sid in IDS], "missing": None},
    )

    loaded, out = _load(folder)

    record_case(input="tables.missing declared, file never written",
                expected=(10, ["missing"], []),
                actual=(len(out["structures"].ids), loaded.absent_tables,
                        loaded.excused_ids))
    assert len(out["structures"].ids) == 10
    assert loaded.absent_tables == ["missing"]
    assert loaded.excused_ids == []
    assert any("missing" in line for line in loaded.recovery_summary())


def test_missing_csv_without_an_id_column_excuses_nothing(loadable_env, tmp_path, record_case):
    folder = _build(str(tmp_path / "001_Mock"))
    path = os.path.join(folder, "tables", "missing.csv")
    with open(path, "w", newline="") as f:
        csv.writer(f).writerow(["removed_by", "kind"])
    doc_path = os.path.join(folder, ".expected_outputs.json")
    doc = json.load(open(doc_path))
    doc["output_structure"]["tables"]["missing"] = {
        "name": "missing", "path": path, "columns": ["removed_by", "kind"],
        "description": "",
    }
    json.dump(doc, open(doc_path, "w"))

    loaded, out = _load(folder)

    record_case(input="missing.csv with no 'id' column",
                expected=(10, []),
                actual=(len(out["structures"].ids), loaded.excused_ids))
    assert len(out["structures"].ids) == 10
    assert loaded.excused_ids == []


def test_declared_table_absent_from_disk_is_named(loadable_env, tmp_path, record_case):
    folder = _build(str(tmp_path / "001_Mock"),
                    tables={"scores": None})

    loaded, out = _load(folder)

    record_case(input="tables/scores.csv declared but absent",
                expected=(["scores"], True),
                actual=(loaded.absent_tables, "scores" in out["tables"]))
    assert loaded.absent_tables == ["scores"]
    # Still declared, so a caller can see the path it was meant to be at.
    assert "scores" in out["tables"]
    assert any("scores" in line for line in loaded.recovery_summary())


# ── the guarantee, across every shape ─────────────────────────────────────────

def _every_shape(root):
    """(label, folder) for each partially-completed shape."""
    shapes = []

    f = _build(str(root / "b" / "001_Mock"))
    _delete_files(f, IDS[6:])
    shapes.append(("some ids absent", f))

    f = _build(str(root / "c" / "001_Mock"))
    _truncate_map_mid_row(f, "design_7")
    _delete_files(f, IDS[6:])
    shapes.append(("truncated map_table", f))

    f = _build(str(root / "d" / "001_Mock"), map_rows="empty")
    _delete_files(f, IDS)
    shapes.append(("header-only map_table", f))

    f = _build(str(root / "e" / "001_Mock"), map_rows="none")
    _delete_files(f, IDS[6:])
    shapes.append(("absent map_table", f))

    f = _build(str(root / "f" / "001_Mock"), map_rows="zero")
    shapes.append(("zero-byte map_table", f))

    f = _build(str(root / "g" / "001_Mock"), map_rows=IDS[:6])
    _delete_files(f, IDS[6:])
    _mark(f, "FAILED")
    shapes.append(("_FAILED marker", f))

    f = _build(str(root / "h" / "001_Mock"), produced_ids=IDS[:8],
               tables={"scores": [[s, "1.0"] for s in IDS[:8]],
                       "missing": [[s, "001_Mock", "failure", "m"] for s in IDS[8:]]})
    shapes.append(("populated missing.csv", f))

    f = _build(str(root / "i" / "001_Mock"),
               tables={"scores": [[s, "1.0"] for s in IDS], "missing": None})
    shapes.append(("declared missing.csv absent", f))

    f = _build(str(root / "j" / "001_Mock"), tables={"scores": None})
    shapes.append(("declared table absent", f))

    return shapes


def test_every_partial_shape_loads_and_reports(loadable_env, tmp_path, record_case):
    """Across every shape: no raise, a non-empty report, and honest counts.

    The report is the deliverable — an empty stream returned quietly is exactly
    the failure mode exit-0 makes possible.
    """
    import contextlib

    outcomes = {}
    for label, folder in _every_shape(tmp_path):
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            loaded, out = _load(folder)
        summary = loaded.recovery_summary()
        printed = buf.getvalue()
        outcomes[label] = {
            "ids": len(out["structures"].ids),
            "declared": loaded.recovery["structures"]["declared"],
            "summary_lines": len(summary),
            "summary_printed": "recovery summary" in printed,
        }

    record_case(input=f"{len(outcomes)} partial-output shapes",
                expected="every shape loads, reports, and states a count",
                actual=outcomes)

    for label, o in outcomes.items():
        assert o["declared"] == 10, label
        assert o["summary_lines"] >= 3, label
        assert o["summary_printed"], f"{label}: recovery summary was not printed"

    # A recovery of zero is stated as such, never left to look like success.
    assert outcomes["header-only map_table"]["ids"] == 0


def test_load_never_raises_on_a_partial_folder_even_without_validation(loadable_env, tmp_path):
    """validate_files=False is the documented escape hatch; it must not raise either."""
    for _label, folder in _every_shape(tmp_path):
        loaded, out = _load(folder, validate_files=False)
        assert out["structures"] is not None


# ── anchored to the real producer ─────────────────────────────────────────────

def test_real_mock_output_recovers_after_files_and_map_rows_are_lost(loadable_env, 
        tmp_path, record_case):
    """Same guarantee against a folder written by ``pipe_mock.py`` itself.

    The hand-built fixtures above copy the Mock layout; this one is the layout,
    so the fixtures cannot drift away from what a producer actually writes.
    """
    out_dir = tmp_path / "001_Mock"
    stream_dir = out_dir / "structures"
    stream_dir.mkdir(parents=True)
    map_table = stream_dir / "structures_map.csv"
    cfg = {
        "output_folder": str(out_dir),
        "parent_ids": IDS,
        "output_ids": IDS,
        "provenance": {},
        "axis_names": [],
        "children": None,
        "produce": None,
        "streams": {"structures": {"format": "pdb", "file": "<id>.pdb",
                                   "values": None, "map_table": str(map_table),
                                   "stream_folder": str(stream_dir)}},
        "tables": {},
        "map_table_strategy": "runtime",
        "missing": [],
        "source_streams": [],
    }
    cfg_path = tmp_path / "mock_config.json"
    json.dump(cfg, open(cfg_path, "w"))
    subprocess.run([sys.executable,
                    os.path.join(REPO_ROOT, "pipe_scripts", "pipe_mock.py"),
                    str(cfg_path)], check=True, capture_output=True)

    assert map_table.exists(), "pipe_mock.py did not write the map_table"

    with open(os.path.join(str(out_dir), ".expected_outputs.json"), "w") as f:
        json.dump({
            "tool_name": "Mock", "tool_class": "Mock",
            "output_structure": {
                "structures": {"name": "structures", "ids": IDS,
                               "files": [str(stream_dir / "<id>.pdb")],
                               "map_table": str(map_table), "format": "pdb",
                               "metadata": {}},
                "tables": {}, "output_folder": str(out_dir),
            },
        }, f)

    # Now break it the way a killed step does: drop the last 4 map rows and files.
    rows = map_table.read_text().splitlines()
    map_table.write_text("\n".join(rows[:7]) + "\n")
    for sid in IDS[6:]:
        os.remove(str(stream_dir / f"{sid}.pdb"))

    loaded, out = _load(str(out_dir))

    record_case(input="real pipe_mock.py output, last 4 rows and files removed",
                expected=(IDS[:6], 10, 6),
                actual=(list(out["structures"].ids),
                        loaded.recovery["structures"]["declared"],
                        loaded.recovery["structures"]["recovered"]))
    assert list(out["structures"].ids) == IDS[:6]
    assert loaded.recovery["structures"] == {
        "declared": 10, "recovered": 6, "source": "map_table",
        "unrecovered_ids": [], "absent_files": [],
    }
