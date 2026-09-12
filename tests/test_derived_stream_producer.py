"""A derived stream keeps the tool that produced it.

`_producer` is the back-pointer the dataflow recovery walks at `save()` time
(`base_config._walk_producer_refs`) to work out which step feeds which.
`StandardizedOutput` preserves it across indexing, iteration and chunking, but
`DataStream`'s own derivations built a fresh object and dropped it.

The visible consequence was a false statement rather than a missing one: with
no edge recovered, the run page falls into the `not in_edges` branch and
renders a chip reading "no dataflow recovered -- this is the preceding step,
not a dependency". The steps are wired; the page says they are not.

Nothing schedules off `_producer` -- dependency directives come from batch
order -- so this is a reporting defect. It is still the page people read to
find out what a run did.
"""
from __future__ import annotations

import pytest

from biopipelines.datastream import DataStream

SENTINEL = object()


def _stream(ids=("a", "b", "c", "d"), producer=SENTINEL):
    ds = DataStream(
        name="structures",
        ids=list(ids),
        files=[f"/tmp/{i}.pdb" for i in ids],
        format="pdb",
    )
    if producer is not SENTINEL:
        ds._producer = producer
    return ds


def test_filter_by_ids_keeps_the_producer():
    """The shape that hits this most: a consumer taking a narrowed stream."""
    ds = _stream(producer="STEP1")

    assert ds.filter_by_ids(["a", "c"])._producer == "STEP1"


def test_slicing_keeps_the_producer():
    assert _stream(producer="STEP1")[0:2]._producer == "STEP1"


def test_indexing_keeps_the_producer():
    assert _stream(producer="STEP1")[0]._producer == "STEP1"


def test_chunking_keeps_the_producer():
    """`for chunk in stream.chunks(n)` is how a fan-out is written."""
    chunks = _stream(producer="STEP1").chunks(2)

    assert len(chunks) == 2
    for chunk in chunks:
        assert chunk._producer == "STEP1"


def test_iteration_keeps_the_producer():
    for item in _stream(producer="STEP1"):
        assert item._producer == "STEP1"


def test_iteration_of_a_shared_file_stream_keeps_the_producer():
    """The shared-file branch builds its sub-streams separately."""
    ds = DataStream(name="fasta", ids=["a", "b"], files="/tmp/all.fasta", format="fasta")
    ds._producer = "STEP2"

    for item in ds:
        assert item._producer == "STEP2"


@pytest.mark.parametrize("derive", [
    lambda ds: ds.filter_by_ids(["a"]),
    lambda ds: ds[0:1],
    lambda ds: ds[0],
    lambda ds: ds.chunks(2)[0],
    lambda ds: next(iter(ds)),
])
def test_a_stream_with_no_producer_stays_without_one(derive):
    """Inheriting must not invent a producer where the parent had none.

    A stream built directly in a pipe script has no producing tool, and
    stamping one would fabricate a dataflow edge rather than recover it.
    """
    derived = derive(_stream(producer=SENTINEL))

    assert getattr(derived, "_producer", None) is None
