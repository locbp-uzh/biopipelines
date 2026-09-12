"""Regression tests for :func:`biopipelines.id_patterns.select_ids` — exact selection of map_table rows.

``select_ids`` used to build a glob (``[...]`` → ``*``) and match with ``fnmatchcase``. The wildcard was unbounded on the numeric side, so ``design_<1..2>[_<?>]`` also selected ``design_10_*``, ``design_11_*`` and ``design_12_*``: every count and statistic downstream was then computed on the wrong row set. Selection must fail in both directions — no fabricated ids (rows the patterns cover but that are absent stay absent) and no extra ids (rows the patterns do not cover).
"""

import pytest

from biopipelines import id_patterns as idp


DESIGNS_1_TO_12 = [f"design_{i}_1A" for i in range(1, 13)]
PADDED_1_TO_12 = [f"design_{i:02d}_1A" for i in range(1, 13)]


# ── the 1-vs-10 over-selection (the bug) ──────────────────────────────────────

@pytest.mark.parametrize("patterns, rows, expected", [
    # Asking for designs 1-2 out of 12 must give exactly two rows, not five.
    (["design_<1..2>[_<?><A>]"], DESIGNS_1_TO_12, ["design_1_1A", "design_2_1A"]),
    (["design_<1..2>[_<?>]"],    DESIGNS_1_TO_12, ["design_1_1A", "design_2_1A"]),
    # A single deterministic id must not reach its longer-numbered neighbors.
    (["design_1[_<?>]"],         DESIGNS_1_TO_12, ["design_1_1A"]),
    # Nor when the bracket is followed by more literal text.
    (["prot_<0..1>[_<?><A I L V>]+9DP"],
     ["prot_0_1A+9DP", "prot_1_2V+9DP", "prot_10_1A+9DP", "prot_11_3L+9DP"],
     ["prot_0_1A+9DP", "prot_1_2V+9DP"]),
    # A shared literal prefix is a prefix, not a match: 'ab' is not 'abc'.
    (["ab[_<?>]"], ["ab_1", "abc_1", "ab_2", "abcd_1"], ["ab_1", "ab_2"]),
])
def test_no_over_selection_across_digit_boundary(record_case, patterns, rows, expected):
    actual = idp.select_ids(patterns, rows)
    record_case(input=(patterns, rows), expected=expected, actual=actual)
    assert actual == expected


def test_selection_count_matches_pattern_count(record_case):
    """A deterministic prefix of N ids selects at most N rows per lazy child."""
    patterns = ["design_<1..2>[_<?>]"]
    actual = idp.select_ids(patterns, DESIGNS_1_TO_12)
    record_case(input=(patterns, DESIGNS_1_TO_12), expected=2, actual=len(actual))
    assert len(actual) == idp.count_ids(patterns) == 2


# ── zero-padded ids ───────────────────────────────────────────────────────────

@pytest.mark.parametrize("patterns, rows, expected", [
    # '<1..2>' expands to unpadded '1','2', which must not match padded rows at all.
    (["design_<1..2>[_<?>]"], PADDED_1_TO_12, []),
    # The padding has to be written into the pattern to select padded rows.
    (["design_0<1..2>[_<?>]"], PADDED_1_TO_12, ["design_01_1A", "design_02_1A"]),
    # '01' must not match '1' either — the deterministic part is exact both ways.
    (["design_<01 02>"], ["design_1", "design_2", "design_01", "design_02"],
     ["design_01", "design_02"]),
])
def test_zero_padded_ids(record_case, patterns, rows, expected):
    actual = idp.select_ids(patterns, rows)
    record_case(input=(patterns, rows), expected=expected, actual=actual)
    assert actual == expected


# ── multi-axis patterns ───────────────────────────────────────────────────────

@pytest.mark.parametrize("patterns, rows, expected", [
    # Two deterministic axes: the cartesian product selects exactly, and 'p_1_A' must not reach 'p_10_A' or 'p_1_AB'.
    (["p_<1..2>_<A B>"], ["p_1_A", "p_1_B", "p_2_A", "p_2_B", "p_10_A", "p_1_AB"],
     ["p_1_A", "p_1_B", "p_2_A", "p_2_B"]),
    # Combinatorial '+' ids with a lazy multiplier suffix on the second axis.
    (["prot_<1..2>+lig_<1..2>[_<?>]"],
     ["prot_1+lig_1_1", "prot_1+lig_2_1", "prot_2+lig_1_1",
      "prot_10+lig_1_1", "prot_1+lig_10_1", "prot_2+lig_20_3"],
     ["prot_1+lig_1_1", "prot_1+lig_2_1", "prot_2+lig_1_1"]),
    # A union of patterns is the union of their exact selections, in row order.
    (["design_<1..2>[_<?>]", "design_5[_<?>]"], DESIGNS_1_TO_12,
     ["design_1_1A", "design_2_1A", "design_5_1A"]),
])
def test_multi_axis_patterns(record_case, patterns, rows, expected):
    actual = idp.select_ids(patterns, rows)
    record_case(input=(patterns, rows), expected=expected, actual=actual)
    assert actual == expected


# ── enumeration slots ─────────────────────────────────────────────────────────

@pytest.mark.parametrize("patterns, rows, expected", [
    (["pos_<42A 42V 42W>"], ["pos_42A", "pos_42V", "pos_42W", "pos_42WX", "pos_142A"],
     ["pos_42A", "pos_42V", "pos_42W"]),
    (["design_<1 2>"], [f"design_{i}" for i in range(1, 13)], ["design_1", "design_2"]),
])
def test_enumeration_slot_matches_exactly(record_case, patterns, rows, expected):
    actual = idp.select_ids(patterns, rows)
    record_case(input=(patterns, rows), expected=expected, actual=actual)
    assert actual == expected


# ── lazy suffixes ─────────────────────────────────────────────────────────────

@pytest.mark.parametrize("patterns, rows, expected", [
    # The documented case: the bracket covers whatever the multiplier produced.
    (["a_<1 2>[_<?>]"], ["a_1_x", "a_2_y", "b_1"], ["a_1_x", "a_2_y"]),
    (["a_1[_<?>]"], ["a_1_p", "a_1_q", "a_2"], ["a_1_p", "a_1_q"]),
    # The bracket is a wildcard over its slots, so suffix length is free.
    (["prot_<0..1>[_<?><S A L K>]"],
     ["prot_0_1S", "prot_0_12A", "prot_1_1L", "prot_10_1K", "prot_11_2S"],
     ["prot_0_1S", "prot_0_12A", "prot_1_1L"]),
    # A bracket is optional, so the bare parent row is covered too — the group may resolve to nothing.
    (["a_1[_<?>]"], ["a_1", "a_1_x"], ["a_1", "a_1_x"]),
    # And repeatable: each occurrence carries its own separator.
    (["a_1[_<?>]"], ["a_1", "a_1_x", "a_1_x_y", "a_10_x"], ["a_1", "a_1_x", "a_1_x_y"]),
])
def test_lazy_suffix_selection(record_case, patterns, rows, expected):
    actual = idp.select_ids(patterns, rows)
    record_case(input=(patterns, rows), expected=expected, actual=actual)
    assert actual == expected


# ── zero matches, and never fabricating ids ───────────────────────────────────

@pytest.mark.parametrize("patterns, rows", [
    (["zzz_<1..2>[_<?>]"], DESIGNS_1_TO_12),       # no row shares the prefix
    (["design_<1..2>[_<?>]"], []),                 # no rows at all
    ([], DESIGNS_1_TO_12),                         # no patterns
    (["design_<20..21>[_<?>]"], DESIGNS_1_TO_12),  # rows exist, these ids do not
])
def test_pattern_matching_zero_rows(record_case, patterns, rows):
    actual = idp.select_ids(patterns, rows)
    record_case(input=(patterns, rows), expected=[], actual=actual)
    assert actual == []


def test_never_fabricates_absent_ids(record_case):
    """The rows are the source of truth: a pattern covering absent ids yields fewer rows, never invented ones."""
    patterns = ["a_<0..4>"]
    rows = ["a_0", "a_2", "a_4"]
    actual = idp.select_ids(patterns, rows)
    record_case(input=(patterns, rows), expected=rows, actual=actual)
    assert actual == rows


def test_selection_is_a_subsequence_of_the_rows(record_case):
    """Every result is a row, in row order, without duplicates from overlapping patterns."""
    patterns = ["design_<1..3>[_<?>]", "design_2[_<?>]"]
    rows = ["design_3_1A", "design_1_1A", "design_2_1A", "design_10_1A"]
    actual = idp.select_ids(patterns, rows)
    expected = ["design_3_1A", "design_1_1A", "design_2_1A"]
    record_case(input=(patterns, rows), expected=expected, actual=actual)
    assert actual == expected
    assert len(actual) == len(set(actual))
    assert all(a in rows for a in actual)


# ── the two bracket spellings, and what they cost to match ────────────────────

@pytest.mark.parametrize("pattern, expected", [
    # The delimiter inside the bracket: genuinely optional and repeatable.
    ("parent[_<?>]", ["parent", "parent_4E", "parent_4E_6U", "parent_1", "parent_4E_6U_9Z"]),
    # The delimiter outside it: the value cannot cross a delimiter, so at most one suffix. Zero
    # repetitions would give the bare `parent_`, which is never a real id, so this means exactly one.
    ("parent_[<?>]", ["parent_4E", "parent_", "parent_1"]),
])
def test_where_the_delimiter_sits_decides_optional_and_repeatable(record_case, pattern, expected):
    rows = ["parent", "parent_4E", "parent_4E_6U", "parent_", "parent_1", "parentX",
            "parent_4E_6U_9Z"]
    actual = idp.select_ids([pattern], rows)
    record_case(input=f"select_ids([{pattern!r}], {len(rows)} rows)", expected=expected, actual=actual)
    assert actual == expected


@pytest.mark.parametrize("pattern", ["parent_[<?>]", "parent[_<?>]", "parent_[<?><?>]"])
def test_a_non_matching_id_is_rejected_in_linear_time(pattern):
    """`parent_[<?>]` compiled to `(?:[^_]+)*`, whose every split of a delimiter-free run is a candidate, so rejecting an id cost 2**n: 0.36 s at 24 characters and double for each one after. The equivalent `[^_]*` cannot backtrack, and Python 3.10 offers no atomic group to express the same intent directly.

    The run is kept to 26 characters on purpose. Longer and the unfixed regex does not return at all, so this test would hang instead of failing — at 26 it takes about 1.4 s unfixed against a few microseconds fixed, which the threshold below separates cleanly.
    """
    import time

    compiled = idp._regex_from_lazy(pattern)
    # A trailing delimiter with nothing after it cannot match either spelling, so the regex must try
    # every split of the run before giving up. `..._x` would match `parent[_<?>]` as two repetitions.
    subject = "parent_" + "A" * 26 + "_"

    start = time.perf_counter()
    assert compiled.fullmatch(subject) is None
    assert time.perf_counter() - start < 0.05, "matching backtracks exponentially"
