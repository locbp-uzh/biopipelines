"""``table_utils`` against every container a caller actually holds.

These five functions are the only supported way to read a finished run's tables
back into Python, so each one has to accept a tool's ``StandardizedOutput``, the
``TableContainer`` behind its ``.tables``, and the plain ``{name: TableInfo}``
dict ``Load(...).get_output_files()["tables"]`` returns — and reject anything
else with a message naming what was passed.

Three defects are pinned here. ``get_table_path`` and ``get_indexed_table`` used to branch on ``hasattr(table, '_entries')``, which ``TableInfo.__getattr__`` made true for every attribute name, so ``get_table_path`` raised ``TypeError: object of type 'TableInfo' has no len()`` for every ordinary table. ``table_exists`` caught the ``TypeError`` that a non-dict source raised and returned ``False`` — a wrong answer ("this run has no such table") in place of an error. And the ``__getattr__`` behind the first of those answered to any name at all; it now refuses private ones, which is what the last three tests in this file cover.
"""

import pytest

from biopipelines import table_utils as tu
from biopipelines.base_config import (
    IndexedTableContainer,
    StandardizedOutput,
    TableContainer,
    TableInfo,
    ToolOutput,
)


ENTRY_POINTS = ["get_table", "get_table_path", "get_indexed_table",
                "list_tables", "table_exists"]


def _sources():
    """One ``scores`` table, reached through each accepted container type."""
    info = TableInfo(name="scores", path="/run/tables/scores.csv",
                     columns=["id", "score"])
    return {
        "dict": {"scores": info},
        "TableContainer": TableContainer({"scores": info}),
        "StandardizedOutput": StandardizedOutput({
            "tables": {"scores": {"path": "/run/tables/scores.csv",
                                  "columns": ["id", "score"]}},
            "output_folder": "/run",
        }),
    }


def _call(name, source):
    """Invoke one entry point against ``source`` for the ``scores`` table."""
    if name == "list_tables":
        return tu.list_tables(source)
    if name == "get_indexed_table":
        return tu.get_indexed_table(source, "scores", "any")
    return getattr(tu, name)(source, "scores")


# ── every entry point x every accepted container ──────────────────────────────

@pytest.mark.parametrize("container", sorted(_sources()))
def test_get_table_returns_the_table_info(container, record_case):
    source = _sources()[container]
    table = tu.get_table(source, "scores")
    record_case(input=f"get_table({container})",
                expected="/run/tables/scores.csv", actual=table.info.path)
    assert isinstance(table, TableInfo)
    assert table.info.path == "/run/tables/scores.csv"


@pytest.mark.parametrize("container", sorted(_sources()))
def test_get_table_path_returns_the_path_for_an_ordinary_table(container, record_case):
    """The headline defect: this raised for every ordinary table, in every container."""
    source = _sources()[container]
    path = tu.get_table_path(source, "scores")
    record_case(input=f"get_table_path({container})",
                expected="/run/tables/scores.csv", actual=path)
    assert path == "/run/tables/scores.csv"


@pytest.mark.parametrize("container", sorted(_sources()))
def test_list_tables_names_the_tables(container, record_case):
    source = _sources()[container]
    names = tu.list_tables(source)
    record_case(input=f"list_tables({container})", expected=["scores"], actual=names)
    assert names == ["scores"]


@pytest.mark.parametrize("container", sorted(_sources()))
def test_table_exists_answers_both_ways(container, record_case):
    source = _sources()[container]
    result = (tu.table_exists(source, "scores"), tu.table_exists(source, "absent"))
    record_case(input=f"table_exists({container}, present/absent)",
                expected=(True, False), actual=result)
    assert result == (True, False)


@pytest.mark.parametrize("container", sorted(_sources()))
def test_get_indexed_table_rejects_an_ordinary_table_by_naming_its_type(
        container, record_case):
    """Same latent flaw: the ``_entries`` probe let a TableInfo through, and it
    died on ``table[entry_id]`` with 'not subscriptable' instead of saying why."""
    source = _sources()[container]
    with pytest.raises(TypeError) as excinfo:
        tu.get_indexed_table(source, "scores", "any")
    message = str(excinfo.value)
    record_case(input=f"get_indexed_table({container}) on an ordinary table",
                expected="TypeError naming TableInfo",
                actual=message)
    assert "TableInfo" in message
    assert "IndexedTableContainer" in message


@pytest.mark.parametrize("container", sorted(_sources()))
def test_missing_name_raises_key_error_listing_what_is_available(container):
    source = _sources()[container]
    for fn in ("get_table", "get_table_path"):
        with pytest.raises(KeyError) as excinfo:
            getattr(tu, fn)(source, "nope")
        assert "scores" in str(excinfo.value)


# ── ToolOutput reaches its tables through the same path ───────────────────────

def test_tool_output_is_accepted(record_case):
    """``ToolOutput.tables`` is a TableContainer, and every docstring names the type."""
    class _FakeConfig:
        TOOL_NAME = "Fake"
        environments = []
        output_folder = "/run"
        execution_order = 1
        tables = TableContainer({"scores": TableInfo(
            name="scores", path="/run/tables/scores.csv", columns=["id"])})

    out = ToolOutput(_FakeConfig())
    result = (tu.list_tables(out), tu.get_table_path(out, "scores"),
              tu.table_exists(out, "scores"))
    record_case(input="ToolOutput", expected=(["scores"], "/run/tables/scores.csv", True),
                actual=result)
    assert result == (["scores"], "/run/tables/scores.csv", True)


# ── indexed tables ────────────────────────────────────────────────────────────

def _indexed_source():
    container = IndexedTableContainer(name="rmsf", columns=["id", "rmsf"])
    container.add("1a2j", path="/run/tables/1a2j_RMSF.csv")
    container.add("1gfl", path="/run/tables/1gfl_RMSF.csv")
    return {"rmsf": container}


def test_get_indexed_table_returns_one_entry(record_case):
    entry = tu.get_indexed_table(_indexed_source(), "rmsf", "1gfl")
    record_case(input="get_indexed_table(dict, 'rmsf', '1gfl')",
                expected="/run/tables/1gfl_RMSF.csv", actual=entry.info.path)
    assert entry.info.path == "/run/tables/1gfl_RMSF.csv"


def test_get_indexed_table_unknown_id_raises_key_error():
    with pytest.raises(KeyError):
        tu.get_indexed_table(_indexed_source(), "rmsf", "nope")


def test_get_table_path_on_an_indexed_table_points_at_get_indexed_table(record_case):
    """An indexed table has one path per entry, so there is no single answer —
    and the error has to say which call does have one."""
    with pytest.raises(TypeError) as excinfo:
        tu.get_table_path(_indexed_source(), "rmsf")
    message = str(excinfo.value)
    record_case(input="get_table_path on an IndexedTableContainer",
                expected="TypeError naming get_indexed_table", actual=message)
    assert "IndexedTableContainer" in message
    assert "get_indexed_table" in message


def test_indexed_container_survives_list_and_exists(record_case):
    source = _indexed_source()
    result = (tu.list_tables(source), tu.table_exists(source, "rmsf"))
    record_case(input="list_tables/table_exists over an indexed table",
                expected=(["rmsf"], True), actual=result)
    assert result == (["rmsf"], True)


# ── rejection: an error, never a wrong answer ────────────────────────────────

@pytest.mark.parametrize("bad", [
    pytest.param(["scores"], id="list"),
    pytest.param("scores", id="str"),
    pytest.param(None, id="None"),
    pytest.param(42, id="int"),
])
@pytest.mark.parametrize("entry_point", ENTRY_POINTS)
def test_unsupported_source_is_rejected_by_every_entry_point(entry_point, bad):
    with pytest.raises(TypeError):
        _call(entry_point, bad)


def test_table_exists_does_not_convert_the_rejection_into_false(record_case):
    """``False`` here reads as "the run has no such table", which is a lie.

    ``table_exists`` used to swallow the ``TypeError`` that every other entry
    point raised, so passing a ``StandardizedOutput`` — the container a user is
    most likely to be holding — silently answered that the table was absent.
    """
    with pytest.raises(TypeError):
        tu.table_exists(object(), "scores")

    # The container that used to trigger the wrong answer now answers correctly.
    out = _sources()["StandardizedOutput"]
    record_case(input="table_exists(StandardizedOutput, 'scores')",
                expected=True, actual=tu.table_exists(out, "scores"))
    assert tu.table_exists(out, "scores") is True


@pytest.mark.parametrize("entry_point", ENTRY_POINTS)
def test_rejection_message_names_what_was_passed_and_what_to_pass(entry_point):
    class Widget:
        pass

    with pytest.raises(TypeError) as excinfo:
        _call(entry_point, Widget())
    message = str(excinfo.value)
    assert "Widget" in message
    for expected in ("StandardizedOutput", "TableContainer", "dict"):
        assert expected in message, f"{entry_point}: message omits {expected}"


# ── the root cause ────────────────────────────────────────────────────────────

def test_table_info_does_not_claim_attributes_it_has_no_business_having():
    """A private name is never a column, so ``hasattr`` must answer honestly.

    ``TableInfo.__getattr__`` used to hand back a TableReference for any name at all, which made every duck-type probe against a TableInfo true -- the ``hasattr(table, '_entries')`` branch in ``table_utils`` being the one that broke ``get_table_path`` for every ordinary table.
    """
    info = TableInfo(name="scores", path="/x.csv", columns=["id"])
    assert not hasattr(info, "_entries")
    assert not hasattr(info, "_tables")
    assert not hasattr(info, "__setstate__")


def test_table_info_still_references_columns_by_dot_notation(record_case):
    """The permissive part is the documented feature, and has to survive.

    Declared columns are set as real attributes by ``__init__``, while an undeclared public name still resolves through ``__getattr__`` -- a tool may know a table's path long before its schema (``columns=[]`` at construction, as plot.py and the Install sentinel do).
    """
    from biopipelines.biopipelines_io import TableReference

    declared = TableInfo(name="scores", path="/x.csv", columns=["id", "score"])
    undeclared_schema = TableInfo(name="plot_data", path="/p.csv", columns=[])

    result = (str(declared.score), str(undeclared_schema.whatever_column))
    record_case(input="declared column / column on a schema-less table",
                expected=("TABLE_REFERENCE:/x.csv:score",
                          "TABLE_REFERENCE:/p.csv:whatever_column"),
                actual=result)
    assert isinstance(declared.score, TableReference)
    assert isinstance(undeclared_schema.whatever_column, TableReference)
    assert result == ("TABLE_REFERENCE:/x.csv:score",
                      "TABLE_REFERENCE:/p.csv:whatever_column")


def test_table_info_survives_copy_instead_of_recursing():
    """The old ``__getattr__`` answered the copy protocol's own dunder lookups.

    ``copy.copy`` asks for ``__copy__``/``__reduce_ex__``, and the permissive version answered by reaching for ``self._info`` on an instance whose ``__init__`` had not run, which came straight back through ``__getattr__``: measured as ``RecursionError`` for ``copy.copy`` and ``TypeError: 'TableReference' object is not callable`` for ``copy.deepcopy``.
    """
    import copy

    for clone in (copy.copy(TableInfo(name="scores", path="/x.csv", columns=["id"])),
                  copy.deepcopy(TableInfo(name="scores", path="/x.csv", columns=["id"]))):
        assert clone.info.path == "/x.csv"
        assert clone.info.columns == ["id"]
