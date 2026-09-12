# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Angle brackets hold literal values, and a match class is punctuation.

`id_patterns` uses `<...>` to compact a *predictable* id, so a slot is a finite set of literal values: `<0..2>` is three ids and `<N>` is the one-element set containing the letter N. That reading now holds inside a `[...]` bracket too, where the set cannot be enumerated but each value is still matched literally, so `prot_[<N>]` repeats the letter.

A lettered class cannot survive in this domain: `<N><S E>` reads as serine or glutamate at position N, while `<15 16><N>` reads as asparagine at positions 15 and 16. So the classes are punctuation and nothing else -- `<#>` for digits, `<?>` for one segment. The letter spellings are gone rather than deprecated, because a slot that means a literal everywhere cannot also mean a class anywhere.

`[...]` marks an *unpredictable* id: its content says what shape the runtime value will take, which is what stops `design_1[_<#>]` from covering `design_1_ZZZ`.
"""

import time

import pytest

from biopipelines import id_map_utils as imu
from biopipelines import id_patterns as ip


# ── one vocabulary: a slot is literals, a class is punctuation ────────────────

def test_an_id_patterns_slot_is_a_set_of_literals_not_a_class():
    """`<N>` compacts a predictable id, so it is the letter N. This is the behavior the notation is for and it does not change."""
    assert ip.expand_pattern("prot_<N>") == ["prot_N"]
    assert ip.expand_pattern("prot_<0..2>") == ["prot_0", "prot_1", "prot_2"]
    assert ip.expand_pattern("prot_<A V>") == ["prot_A", "prot_V"]


def test_an_id_map_class_is_punctuation_so_it_cannot_read_as_a_residue():
    assert imu.parse_id_map_pattern({"*": "*_<#>"}).pattern == r"^(.+)_\d+$"
    assert imu.parse_id_map_pattern({"*": "*_<?>"}).pattern == r"^(.+)_[^_]+$"


@pytest.mark.parametrize("retired", ["*_<N>", "*_<S>", "*-<N>", "*_seq_<N>"])
def test_a_retired_letter_spelling_refuses_instead_of_matching_nothing(retired):
    """A leftover slot would be escaped into the literal text `<N>`, which no real id contains, so the pattern would strip nothing and the caller would see unmapped ids with no error. An id_map can arrive from a user's JSON, so it has to say what is wrong."""
    with pytest.raises(ValueError) as excinfo:
        imu.parse_id_map_pattern({"*": retired})
    assert "<#>" in str(excinfo.value) and "<?>" in str(excinfo.value)


def test_a_canonical_spelling_still_strips():
    assert imu.map_table_ids_to_ids("protein_1", {"*": "*_<#>"}) == ["protein_1", "protein"]
    assert imu.map_table_ids_to_ids("protein_19A", {"*": "*_<?>"}) == ["protein_19A", "protein"]


def test_the_class_tuples_hold_punctuation_only():
    assert imu.DIGIT_CLASSES == ("<#>",)
    assert imu.SEGMENT_CLASSES == ("<?>",)


def test_the_default_id_map_uses_the_canonical_spelling():
    assert imu.DEFAULT_ID_MAP == {"*": "*_<?>"}


# ── the delimiter and the segment rule have one definition ────────────────────

def test_both_modules_take_the_segment_rule_from_id_patterns():
    """The rule that a segment cannot contain the delimiter is one line now, not two: changing `_SUFFIX_DELIMITER` moves the bracket translation and the id_map classes together."""
    assert ip.SEGMENT_REGEX == f"[^{ip.SUFFIX_DELIMITER}]+"
    assert ip.SEGMENT_REGEX in imu.parse_id_map_pattern({"*": "*_<?>"}).pattern
    assert ip.SEGMENT_REGEX in ip._regex_from_lazy("parent[_<?>]").pattern


# ── a bracket says what shape the runtime value takes ─────────────────────────

@pytest.mark.parametrize("pattern, covers, refuses", [
    ("prot_[<N>]", ["prot_", "prot_N", "prot_NN"], ["prot_1", "prot_S", "prot_NS"]),
    ("prot[_<N>]", ["prot", "prot_N", "prot_N_N"], ["prot_1", "prot_NN"]),
])
def test_a_lettered_slot_in_a_bracket_is_that_letter(pattern, covers, refuses):
    """A bare word used to fall through to a segment wildcard, so `prot_[<N>]` silently matched every row. A slot is a literal value set everywhere now, so it matches the letters it spells and nothing else."""
    assert ip.select_ids([pattern], covers + refuses) == covers


@pytest.mark.parametrize("pattern, covers, refuses", [
    ("design_1[_<#>]", ["design_1_99"], ["design_1_1A", "design_1_ZZZ"]),
    ("design_1[_<#><A V>]", ["design_1_1A", "design_1_2V"], ["design_1_99", "design_1_ZZZ"]),
    ("design_1[_<0..9>]", ["design_1_7"], ["design_1_77", "design_1_A"]),
    ("design_1[_<?>]", ["design_1_99", "design_1_ZZZ"], []),
])
def test_an_explicit_shape_in_a_bracket_constrains(pattern, covers, refuses):
    """A bracket used to discard its content entirely, so `design_1[_<#>]` covered `design_1_ZZZ` as readily as `design_1_99` -- it accepted a constraint and honoured none of it. A slot cannot be *enumerated* inside a bracket, since the values arrive at runtime, but it can still say what shape they will have."""
    rows = covers + refuses
    assert ip.select_ids([pattern], rows) == covers


@pytest.mark.parametrize("pattern", ["design_1_[<#>]", "design_1_[<?>]", "design_1[_<#>]"])
def test_a_single_slot_bracket_rejects_in_linear_time(pattern):
    """`(?:X+)*` re-splits every unbroken run, so rejecting an n-character subject costs 2**n. `prot_[<#>]` compiled to `prot_(?:\\d+)*` and took seconds on a 24-digit id; both single-slot shapes must collapse to `X*`.

    The delimiter-outside digit shape is the one that regressed, and it is the one six tools emit.
    """
    compiled = ip._regex_from_lazy(pattern)
    subject = "design_1_" + "1" * 26 + "x"
    start = time.perf_counter()
    compiled.fullmatch(subject)
    assert time.perf_counter() - start < 0.05


@pytest.mark.parametrize("pattern, expected", [
    ("design_1_[<#>]", r"design_1_\d*"),
    ("design_1_[<?>]", "design_1_[^_]*"),
])
def test_a_single_slot_bracket_compiles_without_a_nested_quantifier(pattern, expected):
    """The timing test above proves the symptom is gone; this pins the shape, so a future edit cannot reintroduce the nested quantifier and pass on a fast machine."""
    assert ip._regex_from_lazy(pattern).pattern == expected
