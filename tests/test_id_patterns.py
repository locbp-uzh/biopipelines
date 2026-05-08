"""Unit tests for biopipelines.id_patterns — pure pattern/ID logic.

Parametrized so the XLSX report shows the input pattern, the expected
expansion, and the actual expansion for every case.
"""

import pytest

from biopipelines import id_patterns as idp


# ── expand_pattern ────────────────────────────────────────────────────────────

@pytest.mark.parametrize("pattern, expected", [
    ("literal",              ["literal"]),
    ("base_<0..2>",          ["base_0", "base_1", "base_2"]),
    ("pos_<A B C>",          ["pos_A", "pos_B", "pos_C"]),
    ("<0..1>_<A B>",         ["0_A", "0_B", "1_A", "1_B"]),
    ("5HG6_<0..1>_<A B>",    ["5HG6_0_A", "5HG6_0_B", "5HG6_1_A", "5HG6_1_B"]),
])
def test_expand_pattern(record_case, pattern, expected):
    actual = idp.expand_pattern(pattern)
    record_case(input=pattern, expected=expected, actual=actual)
    assert actual == expected


def test_expand_pattern_lazy_raises(record_case):
    pattern = "prot_<0..2>[_<N><A V>]"
    record_case(input=pattern, expected="LazyPatternError", actual="LazyPatternError")
    with pytest.raises(idp.LazyPatternError):
        idp.expand_pattern(pattern)


# ── expand_ids ────────────────────────────────────────────────────────────────

@pytest.mark.parametrize("ids, expected", [
    (["a_<0..1>", "literal"],      ["a_0", "a_1", "literal"]),
    (["<0..1>_<A B>"],             ["0_A", "0_B", "1_A", "1_B"]),
    (["x", "y", "z"],              ["x", "y", "z"]),
])
def test_expand_ids(record_case, ids, expected):
    actual = idp.expand_ids(ids)
    record_case(input=ids, expected=expected, actual=actual)
    assert actual == expected


def test_expand_ids_lazy_raises(record_case):
    ids = ["prot[_<N><A V>]"]
    record_case(input=ids, expected="LazyPatternError", actual="LazyPatternError")
    with pytest.raises(idp.LazyPatternError):
        idp.expand_ids(ids)


# ── try_expand ────────────────────────────────────────────────────────────────

@pytest.mark.parametrize("pattern, expected_ids, expected_complete", [
    ("base_<0..2>",              ["base_0", "base_1", "base_2"], True),
    ("prot_<0..1>[_<N><A V>]",   ["prot_0", "prot_1"],            False),
    ("base[_<N><A V>]",          ["base"],                        False),
])
def test_try_expand(record_case, pattern, expected_ids, expected_complete):
    ids, complete = idp.try_expand(pattern)
    actual = (ids, complete)
    expected = (expected_ids, expected_complete)
    record_case(input=pattern, expected=expected, actual=actual)
    assert actual == expected


def test_try_expand_ids_mixed(record_case):
    ids = ["a_<0..1>", "b[_<N><X>]"]
    expected = (["a_0", "a_1", "b"], False)
    actual = idp.try_expand_ids(ids)
    record_case(input=ids, expected=expected, actual=actual)
    assert actual == expected


# ── count_pattern / count_ids ─────────────────────────────────────────────────

@pytest.mark.parametrize("pattern, expected", [
    ("base_<0..49>",              50),
    ("<0..2>_<A B>",              6),
    ("literal",                    1),
    ("prot_<0..4>[_<N><A V>]",    5),  # lazy: counts prefix only
])
def test_count_pattern(record_case, pattern, expected):
    actual = idp.count_pattern(pattern)
    record_case(input=pattern, expected=expected, actual=actual)
    assert actual == expected


def test_count_ids_sums(record_case):
    ids = ["a_<0..2>", "lit", "b_<0..1>"]
    expected = 3 + 1 + 2
    actual = idp.count_ids(ids)
    record_case(input=ids, expected=expected, actual=actual)
    assert actual == expected


# ── bracket / lazy handling ───────────────────────────────────────────────────

@pytest.mark.parametrize("value, expected", [
    ("x[_<N><A V>]",  True),
    ("x_<0..2>",      False),
    ("plain",         False),
])
def test_is_lazy(record_case, value, expected):
    actual = idp.is_lazy(value)
    record_case(input=value, expected=expected, actual=actual)
    assert actual is expected


@pytest.mark.parametrize("value, expected", [
    ("prot_<0..4>[_<N><S A L K>]",  "prot_<0..4>"),
    ("base[_<N><A V>]",             "base"),
    ("literal",                     "literal"),
])
def test_strip_brackets(record_case, value, expected):
    actual = idp.strip_brackets(value)
    record_case(input=value, expected=expected, actual=actual)
    assert actual == expected


@pytest.mark.parametrize("value, expected", [
    ("prot_<0..2>[_<N><A V>]+9DP",  "prot_<0..2>*+9DP"),
    ("prot_<0..2>[_<N><A I L V>]",  "prot_<0..2>*"),
    ("literal",                     "literal"),
])
def test_glob_from_lazy(record_case, value, expected):
    actual = idp.glob_from_lazy(value)
    record_case(input=value, expected=expected, actual=actual)
    assert actual == expected


@pytest.mark.parametrize("ids, expected", [
    (["prot_<0..1>[_<N><A V>]+X"],  ["prot_0*+X", "prot_1*+X"]),
    (["prot_<0..1>"],               ["prot_0", "prot_1"]),
])
def test_glob_from_lazy_ids(record_case, ids, expected):
    actual = idp.glob_from_lazy_ids(ids)
    record_case(input=ids, expected=expected, actual=actual)
    assert actual == expected


# ── expand_at (regression for Reviewer A3 bug) ────────────────────────────────

@pytest.mark.parametrize("pattern, index, expected", [
    # The exact A3 regression case the reviewer cited: expand_at("base_<0..2>", 1)
    # used to IndexError because the buggy code wrote result[m.end()] (single
    # char) instead of result[m.end():] (slice).
    ("base_<0..2>",          0, "base_0"),
    ("base_<0..2>",          1, "base_1"),
    ("base_<0..2>",          2, "base_2"),
    ("base_<0..2>_suffix",   1, "base_1_suffix"),
    ("literal",              0, "literal"),
    ("<0..1>_<A B>",         0, "0_A"),
    ("<0..1>_<A B>",         1, "0_B"),
    ("<0..1>_<A B>",         2, "1_A"),
    ("<0..1>_<A B>",         3, "1_B"),
])
def test_expand_at_regression(record_case, pattern, index, expected):
    actual = idp.expand_at(pattern, index)
    record_case(
        input=f"expand_at({pattern!r}, {index})",
        expected=expected,
        actual=actual,
    )
    assert actual == expected


@pytest.mark.parametrize("pattern, index", [
    ("base_<0..2>",  3),
    ("literal",      1),
    ("base_<0..2>", -1),
])
def test_expand_at_out_of_range(record_case, pattern, index):
    record_case(
        input=f"expand_at({pattern!r}, {index})",
        expected="IndexError",
        actual="IndexError",
    )
    with pytest.raises(IndexError):
        idp.expand_at(pattern, index)


# ── dedup_parent_children ─────────────────────────────────────────────────────

@pytest.mark.parametrize("ids, expected", [
    (["prot_<0..2>", "prot_0", "prot_1"],   ["prot_<0..2>"]),
    (["prot_0", "prot_1", "other"],         ["prot_0", "prot_1", "other"]),
    (["prot_<0..2>", "other"],              ["prot_<0..2>", "other"]),
])
def test_dedup_parent_children(record_case, ids, expected):
    actual = idp.dedup_parent_children(ids)
    record_case(input=ids, expected=expected, actual=actual)
    assert actual == expected


# ── composition helpers ───────────────────────────────────────────────────────

def test_expand_file_pattern(record_case):
    inp = ("<id>.pdb", "5HG6_0")
    expected = "5HG6_0.pdb"
    actual = idp.expand_file_pattern(*inp)
    record_case(input=inp, expected=expected, actual=actual)
    assert actual == expected


@pytest.mark.parametrize("template, expected", [
    ("foo*.pdb", True),
    ("foo.pdb",  False),
])
def test_file_has_glob(record_case, template, expected):
    actual = idp.file_has_glob(template)
    record_case(input=template, expected=expected, actual=actual)
    assert actual is expected


def test_make_range(record_case):
    inp = ("design", 0, 49)
    expected = "design_<0..49>"
    actual = idp.make_range(*inp)
    record_case(input=inp, expected=expected, actual=actual)
    assert actual == expected


def test_make_set(record_case):
    inp = ("pos", ["42A", "42V", "42W"])
    expected = "pos_<42A 42V 42W>"
    actual = idp.make_set(*inp)
    record_case(input=inp, expected=expected, actual=actual)
    assert actual == expected


def test_append_suffix(record_case):
    inp = (["5HG6_<0..4>"], "<1..3>")
    expected = ["5HG6_<0..4>_<1..3>"]
    actual = idp.append_suffix(*inp)
    record_case(input=inp, expected=expected, actual=actual)
    assert actual == expected


@pytest.mark.parametrize("value, expected", [
    ("a_<0..1>",  True),
    ("literal",   False),
])
def test_contains_pattern(record_case, value, expected):
    actual = idp.contains_pattern(value)
    record_case(input=value, expected=expected, actual=actual)
    assert actual is expected


@pytest.mark.parametrize("value, expected", [
    ("literal",    True),
    ("a_<0..1>",   False),
])
def test_is_literal(record_case, value, expected):
    actual = idp.is_literal(value)
    record_case(input=value, expected=expected, actual=actual)
    assert actual is expected


# ── Mock-driven: lazy pattern end-to-end ──────────────────────────────────────

def test_id_patterns_lazy_children_stays_lazy_via_mock(
    local_config, isolated_cwd, new_pipeline, record_case,
):
    """Mock with lazy children keeps ids bracketed at config time (no expansion)."""
    from biopipelines.mock import Mock

    pipeline = new_pipeline("idp_lazy_mock")
    with pipeline:
        m = Mock(
            ids=["p0", "p1"],
            children="[_<N><A V>]",
            produce=["_1A", "_2V"],
            streams={"s": {"format": "pdb", "file": "<id>.pdb"}},
        )
        pipeline.save()

    ids = list(m.streams.s.ids)
    all_lazy = all(idp.is_lazy(i) for i in ids)
    record_case(input="Mock lazy children → is_lazy(all)",
                expected=True, actual=all_lazy)
    assert all_lazy
    assert idp.strip_brackets(ids[0]).startswith("p0")


def test_id_patterns_deterministic_children_append_suffix_via_mock(
    local_config, isolated_cwd, new_pipeline, record_case,
):
    """Deterministic children is equivalent to append_suffix on parents."""
    from biopipelines.mock import Mock

    pipeline = new_pipeline("idp_det_mock")
    with pipeline:
        m = Mock(
            ids=["p0", "p1"],
            children="<1..2>",
            streams={"s": {"format": "pdb", "file": "<id>.pdb"}},
        )
        pipeline.save()

    ids = list(m.streams.s.ids)
    expected = idp.append_suffix(["p0", "p1"], "<1..2>")
    record_case(input="Mock children=<1..2>",
                expected=list(expected), actual=ids)
    assert ids == list(expected)
