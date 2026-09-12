# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""IdSet owns the id rules that used to be remembered at each call site.

An id was a `str` and an id set a `List[str]`, so nothing carried the rules: a slot value never crosses the `_` suffix delimiter, `[...]` is a zero-or-more repetition, `+` composes axes while `_` separates parent from child. Those rules lived as free functions and were re-derived, or forgotten, wherever ids were handled. These tests pin the type's behavior against the free functions it wraps, so introducing it cannot quietly change what an id means.
"""

import pytest

from biopipelines import id_patterns
from biopipelines.idset import AXIS_SEPARATOR, IdSet, bundle, product


# ── enumeration: all the ids, or the ids that are still patterns ──────────────

def test_an_expandable_set_enumerates_to_concrete_ids():
    assert IdSet(["a_<0..2>"]).enumerated().ids == ["a_0", "a_1", "a_2"]
    assert IdSet(["a_<0..2>"]).can_expand is True


def test_a_lazy_set_enumerates_to_ids_that_keep_their_unexpandable_parts():
    """A lazy bracket's values exist only once the upstream tool has run, so enumerating cannot invent them; it expands the deterministic slots and keeps the bracket."""
    lazy = IdSet(["a_<1 2>[_<?>]"])
    assert lazy.enumerated().ids == ["a_1[_<?>]", "a_2[_<?>]"]
    assert lazy.is_lazy is True
    assert lazy.can_expand is False


def test_enumeration_agrees_with_the_free_function_it_replaces():
    ids = ["a_<0..2>", "b_<1 2>[_<?>]", "literal"]
    assert IdSet(ids).enumerated().ids == id_patterns.partial_expand_ids(ids)


# ── counting: len is declared, count is expanded ──────────────────────────────

def test_len_counts_declared_ids_and_count_counts_expanded_ones():
    """These differ on purpose. At configuration time a pattern's expansion is often not knowable, and reporting a count you cannot compute is worse than reporting the one you can."""
    one_pattern = IdSet(["a_<0..2>"])
    assert len(one_pattern) == 1
    assert one_pattern.count() == 3


def test_count_of_a_lazy_set_counts_the_deterministic_prefix_only():
    assert IdSet(["a_<0..4>[_<?>]"]).count() == 5


# ── bundling ──────────────────────────────────────────────────────────────────

def test_bundling_collapses_a_set_into_one_prefix():
    assert IdSet(["l1", "l2"]).bundled().ids == ["l1+l2"]
    assert len(IdSet(["l1", "l2"]).bundled()) == 1


def test_bundling_drops_duplicates_but_keeps_order():
    assert IdSet(["l2", "l1", "l2"]).bundled().ids == ["l2+l1"]


def test_bundling_several_sets_yields_one_entity():
    assert bundle(IdSet(["l1"]), IdSet(["l2"])).ids == ["l1+l2"]


def test_a_bundled_axis_does_not_multiply_a_product():
    """This is the distinction a hand-rolled cartesian product loses: a bundled axis is one entity however many ids it holds."""
    bundled = IdSet(["LIG1", "LIG2"]).bundled()
    assert bundled.product(IdSet(["seedA_1", "seedA_2"])).ids == [
        "LIG1+LIG2+seedA_1", "LIG1+LIG2+seedA_2"]


# ── cartesian product ─────────────────────────────────────────────────────────

def test_product_is_row_major_left_to_right():
    assert IdSet(["p1", "p2"]).product(IdSet(["l1", "l2"])).ids == [
        "p1+l1", "p1+l2", "p2+l1", "p2+l2"]


def test_product_of_three_axes():
    got = product(IdSet(["a"]), IdSet(["b1", "b2"]), IdSet(["c"]))
    assert got.ids == ["a+b1+c", "a+b2+c"]


def test_product_with_nothing_is_the_set_itself():
    assert IdSet(["p1", "p2"]).product().ids == ["p1", "p2"]


def test_product_refuses_a_bare_list():
    """The point of the type is that a list which skipped the rule cannot be passed where an IdSet belongs."""
    with pytest.raises(TypeError):
        IdSet(["a"]).product(["b"])


# ── suffix multiplication ─────────────────────────────────────────────────────

def test_a_suffix_turns_every_parent_into_a_child():
    assert IdSet(["p1", "p2"]).multiplied_by_suffix("1").ids == ["p1_1", "p2_1"]


def test_a_pattern_suffix_keeps_the_result_a_pattern():
    assert IdSet(["5HG6_<0..4>"]).multiplied_by_suffix("<1..3>").ids == ["5HG6_<0..4>_<1..3>"]


def test_suffix_multiplication_agrees_with_the_free_function():
    ids = ["5HG6_<0..4>", "other"]
    assert IdSet(ids).multiplied_by_suffix("<1..3>").ids == id_patterns.append_suffix(ids, "<1..3>")


def test_a_suffix_may_not_carry_the_axis_separator():
    """`+` composes axes and `_` separates parent from child; an id mixing them cannot be decomposed back."""
    with pytest.raises(ValueError, match="compose"):
        IdSet(["a"]).multiplied_by_suffix(f"x{AXIS_SEPARATOR}y")


# ── renaming ──────────────────────────────────────────────────────────────────

def test_renaming_leaves_unnamed_ids_alone():
    assert IdSet(["a", "b"]).renamed({"a": "x"}).ids == ["x", "b"]


def test_renaming_accepts_a_callable():
    assert IdSet(["a", "b"]).renamed(lambda i: i.upper()).ids == ["A", "B"]


def test_renaming_preserves_cardinality():
    assert len(IdSet(["a", "b", "c"]).renamed(lambda i: i + "z")) == 3


def test_a_rename_that_would_collapse_two_ids_raises():
    """A collapse silently loses a row, and provenance columns join on the id, so it must fail loudly."""
    with pytest.raises(ValueError, match="collapse"):
        IdSet(["a", "b"]).renamed({"a": "b"})


# ── selection ─────────────────────────────────────────────────────────────────

def test_selection_is_exact_outside_the_unresolved_slots():
    rows = [f"design_{i}_1A" for i in range(1, 13)]
    assert IdSet(["design_<1..2>[_<?><A>]"]).select(rows).ids == ["design_1_1A", "design_2_1A"]


def test_selection_agrees_with_the_free_function():
    rows = ["a_1", "a_2", "b_1"]
    patterns = ["a_[<?>]"]
    assert IdSet(patterns).select(rows).ids == id_patterns.select_ids(patterns, rows)


# ── value semantics ───────────────────────────────────────────────────────────

def test_equality_is_over_the_ids():
    assert IdSet(["a", "b"]) == IdSet(["a", "b"])
    assert IdSet(["a", "b"]) != IdSet(["b", "a"])


def test_an_id_set_is_hashable_so_it_can_be_deduplicated():
    assert len({IdSet(["a"]), IdSet(["a"]), IdSet(["b"])}) == 2


def test_a_single_string_is_accepted_as_one_id():
    assert IdSet("a_0").ids == ["a_0"]


def test_a_non_string_id_is_refused():
    with pytest.raises(TypeError):
        IdSet([1])


def test_slicing_yields_an_id_set_and_indexing_yields_an_id():
    ids = IdSet(["a", "b", "c"])
    assert ids[1] == "b"
    assert ids[1:] == IdSet(["b", "c"])
