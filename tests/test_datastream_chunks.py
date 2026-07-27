"""Tests for ``DataStream.chunks()`` — splitting a stream across N workers.

``__iter__`` yields one stream per id, which is the wrong granularity for
fanning work over a handful of parallel tasks. ``chunks()`` groups instead.
"""
from __future__ import annotations

import pytest

from biopipelines.datastream import DataStream


def _stream(ids, files=None, name="s", fmt="pdb"):
    return DataStream(name=name, ids=ids, files=files or [], format=fmt)


# ── count / size ──────────────────────────────────────────────────────────────

def test_count_spreads_remainder_over_leading_chunks():
    """10 into 4 is 3/3/2/2, not 3/3/3/1."""
    chunks = _stream(["p_<0..9>"]).chunks(4)
    assert [len(c) for c in chunks] == [3, 3, 2, 2]
    assert chunks[0].ids_expanded == ["p_0", "p_1", "p_2"]
    assert chunks[-1].ids_expanded == ["p_8", "p_9"]


def test_count_divides_evenly():
    assert [len(c) for c in _stream(["p_<0..7>"]).chunks(4)] == [2, 2, 2, 2]


def test_size_fixes_items_per_chunk():
    chunks = _stream(["p_<0..9>"]).chunks(size=3)
    assert [len(c) for c in chunks] == [3, 3, 3, 1]


def test_more_chunks_than_items_drops_empties():
    chunks = _stream(["a", "b"]).chunks(5)
    assert [c.ids_expanded for c in chunks] == [["a"], ["b"]]


def test_empty_stream_yields_no_chunks():
    assert _stream([]).chunks(4) == []


def test_chunks_partition_the_stream_exactly():
    ids = _stream(["p_<0..9>"]).ids_expanded
    for spec in ({"count": 3}, {"count": 7}, {"size": 4}):
        flat = [i for c in _stream(["p_<0..9>"]).chunks(**spec) for i in c.ids_expanded]
        assert flat == ids, spec


# ── carried stream properties ─────────────────────────────────────────────────

def test_files_follow_their_ids():
    chunks = _stream(["p_<0..3>"], ["p_<0..3>.pdb"]).chunks(2)
    assert chunks[0].files_expanded == ["p_0.pdb", "p_1.pdb"]
    assert chunks[1].files_expanded == ["p_2.pdb", "p_3.pdb"]


def test_chunks_keep_name_and_format():
    for chunk in _stream(["p_<0..3>"], name="designs", fmt="cif").chunks(2):
        assert chunk.name == "designs"
        assert chunk.format == "cif"


def test_shared_file_stream_keeps_the_shared_path():
    """Every chunk of a multi-record artifact points at the same file."""
    ds = DataStream(name="s", ids=["p_<0..3>"], files="all.fasta", format="fasta")
    chunks = ds.chunks(2)
    assert all(c.is_shared_file and c.files == "all.fasta" for c in chunks)


# ── lazy ids ──────────────────────────────────────────────────────────────────

def test_lazy_ids_split_on_their_deterministic_axis():
    """A lazy id has a known outer axis and an unknown bracket fan-out; the
    split happens on the former and the latter is preserved for runtime."""
    ds = _stream(["prot_<0..9>[_<N><S A L K>]"])
    chunks = ds.chunks(4)

    assert [len(c) for c in chunks] == [3, 3, 2, 2]
    assert chunks[0].ids_expanded == [
        "prot_0[_<N><S A L K>]",
        "prot_1[_<N><S A L K>]",
        "prot_2[_<N><S A L K>]",
    ]
    assert all(c.is_lazy for c in chunks)


def test_top_level_lazy_id_raises():
    """No deterministic prefix means no config-time axis to split on."""
    with pytest.raises(ValueError, match="lazy at the top level"):
        _stream(["[_<N><S A L K>]"]).chunks(4)


# ── argument validation ───────────────────────────────────────────────────────

def test_count_and_size_are_mutually_exclusive():
    for kwargs in ({}, {"count": 2, "size": 2}):
        with pytest.raises(ValueError, match="exactly one"):
            _stream(["a", "b"]).chunks(**kwargs)


def test_rejects_non_positive_ints():
    for bad in (0, -1, 2.5, True, "4"):
        with pytest.raises(ValueError, match="positive integer"):
            _stream(["a", "b"]).chunks(bad)
        with pytest.raises(ValueError, match="positive integer"):
            _stream(["a", "b"]).chunks(size=bad)
