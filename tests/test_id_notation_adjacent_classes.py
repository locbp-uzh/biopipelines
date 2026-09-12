"""Adjacent match classes inside one bracket: same language, linear time.

``[<#><?>]`` compiled to ``(?:\\d+[^_]+)*`` -- two unbounded quantifiers
competing for the same characters, since ``<?>`` is ``[^_]+`` and that
includes digits. Rejecting an id then costs one attempt per ordered partition
of its run: measured 395 ms at 22 characters, doubling with each one.

Only ``[<?><?>]`` was fast, because it was the one shape the old collapse
caught -- and that collapse was wrong, rewriting it to ``[^_]*`` and so
dropping the two-character minimum that ``[<#><#>]`` correctly kept.

Committing to the first class (an atomic group) is linear but changes what
matches. The rewrite here does not: these tests pin the equivalence
exhaustively rather than on examples, because "same language" is the whole
claim.
"""
from __future__ import annotations

import itertools
import re
import time

import pytest

from biopipelines.id_patterns import (
    DIGIT_REGEX,
    SEGMENT_REGEX,
    _regex_from_lazy,
    _translate_bracket,
)

CLASSES = ["<#>", "<?>"]
RAW = {"<#>": DIGIT_REGEX, "<?>": SEGMENT_REGEX}

# Small alphabet, both digit and non-digit, plus the suffix delimiter, which
# neither class matches and which is what forces the backtracking to exhaust.
ALPHABET = "01A_"


def _all_strings(max_len):
    out = [""]
    for n in range(1, max_len + 1):
        out += ["".join(t) for t in itertools.product(ALPHABET, repeat=n)]
    return out


SUBJECTS = _all_strings(6)


def _bodies(n):
    return ["".join(c) for c in itertools.product(CLASSES, repeat=n)]


@pytest.mark.parametrize("body", _bodies(2) + _bodies(3))
def test_collapsed_form_matches_exactly_what_the_naive_form_matches(body):
    """The rewrite is a rewrite, not a narrowing: identical over every string."""
    naive = re.compile("^(?:" + "".join(RAW[c] for c in re.findall(r"<.>", body)) + ")*$")
    actual = re.compile("^" + _translate_bracket(body) + "$")

    mismatched = [s for s in SUBJECTS if bool(naive.match(s)) != bool(actual.match(s))]
    assert not mismatched, f"[{body}] differs on {mismatched[:8]}"


@pytest.mark.parametrize("body", _bodies(2))
def test_adjacent_classes_reject_in_linear_time(body):
    """A rejecting id must not cost exponential time.

    The subject is all digits followed by the delimiter, which neither class
    matches, so every one of these is a genuine rejection -- the case the
    engine has to prove by exhausting the partitions.
    """
    rx = re.compile("^p_" + _translate_bracket(body) + "$")

    start = time.perf_counter()
    assert not rx.match("p_" + "1" * 200 + "_")
    elapsed = time.perf_counter() - start

    # The old form needed ~0.4 s for 22 characters; this is 200.
    assert elapsed < 0.05, f"[{body}] took {elapsed * 1000:.1f} ms to reject 200 characters"


def test_two_digit_classes_still_mean_two_or_more_digits():
    """``[<#><#>]`` keeps its minimum; it is not the same as ``[<#>]``."""
    two = re.compile("^" + _regex_from_lazy("p_[<#><#>]").pattern + "$")
    one = re.compile("^" + _regex_from_lazy("p_[<#>]").pattern + "$")

    assert not two.match("p_7")
    assert two.match("p_77")
    assert two.match("p_")
    # The single class accepts the one-digit case the pair refuses.
    assert one.match("p_7")


def test_two_segment_classes_keep_their_minimum_too():
    """``[<?><?>]`` used to collapse to ``[^_]*``, dropping the minimum.

    That made it match a one-character suffix, which is exactly what the pair
    spelling says it should not, and what the digit pair already refused.
    """
    rx = re.compile("^" + _regex_from_lazy("p_[<?><?>]").pattern + "$")

    assert not rx.match("p_A")
    assert rx.match("p_AB")
    assert rx.match("p_")


@pytest.mark.parametrize(
    "pattern,expected",
    [
        ("p_[<#><#>]", r"p_(?:\d\d+)?"),
        ("p_[<?><?>]", "p_(?:[^_][^_]+)?"),
        ("p_[<#><?>]", r"p_(?:\d[^_]+)?"),
        ("p_[<?><#>]", r"p_(?:[^_]+\d)?"),
        # Both segment classes keep their quantifier; only the digit between
        # them has a wider neighbour to give its characters to.
        ("p_[<?><#><?>]", r"p_(?:[^_]+\d[^_]+)?"),
        # Unchanged: a single class, and anything a literal separates.
        ("p_[<?>]", "p_[^_]*"),
        ("p_[<#>]", r"p_\d*"),
        ("p[_<#>]", r"p(?:_\d+)*"),
        ("p_[<#>_<#>]", r"p_(?:\d+_\d+)*"),
    ],
)
def test_compiled_shape_is_pinned(pattern, expected):
    """Pin the emitted regex: the timing test alone would pass on a fast box."""
    assert _regex_from_lazy(pattern).pattern == expected


def test_delimiter_separated_slots_are_untouched():
    """A literal between the slots already bounds each repetition.

    ``[<#>_<#>]`` genuinely repeats -- two iterations carry two delimiters and
    are not subsumed by one -- so it must keep its ``*``.
    """
    rx = re.compile("^" + _regex_from_lazy("p_[<#>_<#>]").pattern + "$")

    assert rx.match("p_1_2")
    assert rx.match("p_1_23_4")
    assert not rx.match("p_1")
