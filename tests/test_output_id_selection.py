# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Selecting one id out of a StandardizedOutput must return that id's file or refuse.

`__getitem__` resolved the file by position: `idx = ds.ids_expanded.index(key)` indexing into `ds.files`. When a stream declares explicit files against unexpanded pattern ids the two lists live in different index spaces, and the mismatch is silent in both directions -- most ids fall off the end and come back with no file at all, and the ones that do not come back holding a *different* id's file. The second is the dangerous one: nothing downstream can tell it from a correct answer.
"""

import pytest

from biopipelines.datastream import DataStream
from biopipelines.outputs import StandardizedOutput


def _output(**stream_kwargs):
    return StandardizedOutput({"structures": DataStream(name="structures", **stream_kwargs)})


# ── the index-space mismatch ──────────────────────────────────────────────────

MISMATCHED = dict(
    ids=["a_<0..2>", "b_<0..2>"],
    files=["/x/a_<0..2>.pdb", "/x/b_<0..2>.pdb"],
    map_table="/x/m.csv",
    format="pdb",
)


@pytest.mark.parametrize("key", ["a_0", "a_1", "a_2", "b_0", "b_1", "b_2"])
def test_an_unresolvable_file_refuses_instead_of_guessing(key):
    output = _output(**MISMATCHED)
    with pytest.raises(KeyError, match="Cannot resolve a file"):
        output[key]


def test_the_refusal_names_both_index_spaces():
    """The message has to say 2 files against 6 expanded ids, or the caller cannot tell which list to fix."""
    output = _output(**MISMATCHED)
    with pytest.raises(KeyError) as excinfo:
        output["a_0"]
    message = str(excinfo.value)
    assert "2 explicit file(s)" in message
    assert "6 expanded id(s)" in message
    assert "2 before expansion" in message


def test_no_id_can_ever_receive_another_ids_file():
    """`so['a_1']` returned `/x/b_<0..2>.pdb` -- b's file for an a id, with no error. This is the whole reason the mismatch is a refusal and not a best effort."""
    output = _output(**MISMATCHED)
    for key in DataStream(name="structures", **MISMATCHED).ids_expanded:
        try:
            files = output[key].streams.structures.files
        except KeyError:
            continue
        for path in files:
            assert key in path, f"id {key!r} was handed {path!r}"


# ── the shapes that do resolve ────────────────────────────────────────────────

def test_a_file_template_covers_every_id():
    output = _output(ids=["a_<0..2>"], files=["/x/<id>.pdb"], map_table="/x/m.csv", format="pdb")
    assert output["a_1"].streams.structures.files == ["/x/<id>.pdb"]
    assert output["a_1"].streams.structures.ids == ["a_1"]


def test_one_explicit_file_per_expanded_id_resolves_by_id():
    output = _output(ids=["a_0", "a_1"], files=["/x/a_0.pdb", "/x/a_1.pdb"],
                     map_table="/x/m.csv", format="pdb")
    assert output["a_0"].streams.structures.files == ["/x/a_0.pdb"]
    assert output["a_1"].streams.structures.files == ["/x/a_1.pdb"]


def test_an_absent_id_still_reports_the_ids_that_exist():
    output = _output(ids=["a_<0..2>"], files=["/x/<id>.pdb"], map_table="/x/m.csv", format="pdb")
    with pytest.raises(KeyError, match="not found in stream"):
        output["a_9"]
