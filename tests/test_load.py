"""Load's narrowing of a value-based stream, and what it does with ids it cannot expand.

``_reconcile_value_stream`` expands the declared ids so it can intersect them with the rows
its map_table lists. A malformed pattern in the declaration made that expansion raise, and the
method used to swallow every exception — a programming error included — and carry on with the
unexpanded ids, with nothing said. The project forbids a guessed fallback; the fallback that
remains is narrowed to the expansion's own errors and announces itself.
"""

import pandas as pd
import pytest


def _stub_load():
    """A Load with only the state ``_reconcile_value_stream`` touches."""
    from biopipelines.load import Load

    load = Load.__new__(Load)
    load.recovery = {}
    load.absent_files_by_stream = {}
    return load


@pytest.fixture
def map_table(tmp_path):
    path = tmp_path / "map.csv"
    pd.DataFrame({"id": ["a", "b"], "sequence": ["AA", "BB"]}).to_csv(path, index=False)
    return str(path)


def _ds(ids, map_table):
    return {"name": "sequences", "ids": list(ids), "files": [],
            "map_table": map_table, "format": "csv"}


def test_expandable_ids_are_narrowed_to_the_map(map_table, capsys):
    load = _stub_load()
    out = load._reconcile_value_stream("sequences", _ds(["a", "b", "c"], map_table))

    assert out is not None
    assert out["ids"] == ["a", "b"]
    assert "cannot expand" not in capsys.readouterr().out


def test_an_unexpandable_id_is_reported_and_the_rest_still_narrow(map_table, capsys):
    """`<1..>` raises ValueError inside the expansion; the fallback must say so, not hide it."""
    load = _stub_load()
    out = load._reconcile_value_stream("sequences", _ds(["a", "b", "<1..>"], map_table))

    printed = capsys.readouterr().out
    assert "cannot expand the declared ids of 'sequences'" in printed
    assert "using them unexpanded" in printed
    assert out is not None
    assert out["ids"] == ["a", "b"]


def test_the_warning_is_printed_once_per_stream(map_table, capsys):
    load = _stub_load()
    for _ in range(3):
        load._reconcile_value_stream("sequences", _ds(["a", "b", "<1..>"], map_table))

    assert capsys.readouterr().out.count("cannot expand the declared ids") == 1


def test_a_programming_error_in_the_expansion_is_no_longer_swallowed(map_table, monkeypatch):
    """The bare `except Exception` hid any bug behind the unexpanded ids; only expansion errors are caught."""
    from biopipelines import datastream

    def boom(_cls, _data):
        raise AttributeError("a typo, not a malformed pattern")

    monkeypatch.setattr(datastream.DataStream, "from_dict", classmethod(boom))
    with pytest.raises(AttributeError):
        _stub_load()._reconcile_value_stream("sequences", _ds(["a", "b"], map_table))
