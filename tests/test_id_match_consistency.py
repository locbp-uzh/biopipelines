# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Tests for the id-match tier tags and the lookup-consistency score.

``get_mapped_ids`` answers each id from the first tier of its ladder that fits, and the last three tiers are string operations on the id, so a lookup can return a DIFFERENT row's value. The permissiveness is deliberate — an upstream tool that renames ids depends on it — so instead of tightening the match, the tier that answered is kept and a whole lookup is scored by how CONSISTENTLY its ids resolved.

Covers, in order:

* every tier of the ladder produces its own tag, including the ``closest_siblings_only`` branch and the no-match fall-off;
* ``get_mapped_ids`` returns exactly what it returned before the tiers existed — the tags are additive;
* the score's defining properties: uniform is 1.0 at any depth, a mix never reaches 1.0, and a RARER minority scores WORSE than a bigger one;
* unmatched ids are counted but kept out of the score, since they return ``None`` rather than a wrong value;
* the reported line, its id-naming cap, and that reporting is log-only;
* the five real cases the id audit found, each with the score it produces.
"""

import os

import pandas as pd
import pytest

from biopipelines import contract_enforcement as ce
from biopipelines import id_map_utils as im
from biopipelines.id_map_utils import (
    TIER_CHILD, TIER_EXACT, TIER_GROUP, TIER_NONE, TIER_PARENT,
    TIER_PROVENANCE, TIER_SIBLING, get_mapped_ids, get_mapped_ids_with_tiers,
    score_id_match_tiers,
)


@pytest.fixture(autouse=True)
def _clean_state(monkeypatch):
    """The provenance/table caches and the contract dedup cache are process-global."""
    from biopipelines import biopipelines_io

    for var in list(os.environ):
        if var.startswith("BIOPIPELINES_ENFORCE"):
            monkeypatch.delenv(var, raising=False)
    im._provenance_lookup_cache.clear()
    biopipelines_io._table_cache.clear()
    biopipelines_io._provenance_index_cache.clear()
    ce.reset_reported()
    yield
    im._provenance_lookup_cache.clear()
    biopipelines_io._table_cache.clear()
    biopipelines_io._provenance_index_cache.clear()
    ce.reset_reported()


def tiers_of(sources, targets, **kwargs):
    return {sid: tier for sid, (_v, tier) in
            get_mapped_ids_with_tiers(sources, targets, **kwargs).items()}


def write_map(tmp_path, name, ids, provenance):
    path = tmp_path / name
    pd.DataFrame({"id": ids, "structures.id": provenance}).to_csv(path, index=False)
    return str(path)


# -- the tiers themselves ----------------------------------------------------

def test_exact_match_is_tagged_exact(record_case):
    actual = get_mapped_ids_with_tiers(["design_1"], ["design_1"])
    record_case(input="design_1 -> [design_1]", expected=("design_1", TIER_EXACT),
                actual=actual["design_1"])
    assert actual["design_1"] == ("design_1", TIER_EXACT)


def test_provenance_match_is_tagged_provenance(tmp_path, record_case):
    path = write_map(tmp_path, "map.csv", ["Panda_1"], ["LID_001_1"])
    actual = get_mapped_ids_with_tiers(["Panda_1"], ["LID_001_1"], map_table_paths=[path])
    record_case(input="Panda_1 -> [LID_001_1] with provenance",
                expected=("LID_001_1", TIER_PROVENANCE), actual=actual["Panda_1"])
    assert actual["Panda_1"] == ("LID_001_1", TIER_PROVENANCE)


def test_child_match_is_tagged_child():
    actual = get_mapped_ids_with_tiers(["design"], ["design_1", "design_2"])
    assert actual["design"] == ("design_1", TIER_CHILD)


def test_parent_match_is_tagged_parent():
    actual = get_mapped_ids_with_tiers(["design_1_1"], ["design_1"])
    assert actual["design_1_1"] == ("design_1", TIER_PARENT)


def test_sibling_match_is_tagged_sibling():
    """The tier that can hand back another row's value, and the reason this file exists."""
    actual = get_mapped_ids_with_tiers(["design_1_1"], ["design_1_2"])
    assert actual["design_1_1"] == ("design_1_2", TIER_SIBLING)


def test_no_match_is_tagged_none():
    actual = get_mapped_ids_with_tiers(["alpha"], ["beta"])
    assert actual["alpha"] == (None, TIER_NONE)


def test_closest_siblings_branch_has_its_own_tier():
    """`closest_siblings_only` replaces the ladder, so it is not one of its tiers."""
    actual = get_mapped_ids_with_tiers(
        ["Panda_29_2"], ["Panda_29_1", "Panda_2_1"], unique=False,
        closest_siblings_only=True)
    assert actual["Panda_29_2"] == (["Panda_29_1"], TIER_GROUP)
    empty = get_mapped_ids_with_tiers(
        ["Panda_29_2"], ["Panda_2_1"], unique=False, closest_siblings_only=True)
    assert empty["Panda_29_2"] == ([], TIER_NONE)


def test_non_unique_lookup_is_tagged_too():
    actual = get_mapped_ids_with_tiers(["design"], ["design_1", "design_2"], unique=False)
    assert actual["design"] == (["design_1", "design_2"], TIER_CHILD)


# -- the tags are additive: matching itself is unchanged ---------------------

@pytest.mark.parametrize("sources,targets,kwargs", [
    (["design_1"], ["design_1"], {}),
    (["design"], ["design_1", "design_2"], {}),
    (["design_1_1"], ["design_1"], {}),
    (["design_1_1"], ["design_1_2"], {}),
    (["alpha"], ["beta"], {}),
    (["design"], ["design_1", "design_2"], {"unique": False}),
    (["prot1+lig3"], ["prot1+lig1", "prot1+lig2"], {}),
    (["Panda_29_2"], ["Panda_29_1", "Panda_2_1"],
     {"unique": False, "closest_siblings_only": True}),
])
def test_get_mapped_ids_returns_the_values_the_tier_core_matched(sources, targets, kwargs):
    """The public entry point must be the tier core with the tier projected away."""
    with_tiers = get_mapped_ids_with_tiers(sources, targets, **kwargs)
    assert get_mapped_ids(sources, targets, **kwargs) == {
        sid: value for sid, (value, _tier) in with_tiers.items()
    }


# -- the score ---------------------------------------------------------------

def test_uniform_exact_lookup_scores_perfectly(record_case):
    ids = [f"design_{i}" for i in range(1, 501)]
    summary = score_id_match_tiers(tiers_of(ids, ids))
    record_case(input="500 ids, all exact", expected=1.0, actual=summary.score)
    assert summary.score == 1.0
    assert summary.tier_counts == {TIER_EXACT: 500}


@pytest.mark.parametrize("tier", [TIER_PROVENANCE, TIER_CHILD, TIER_PARENT, TIER_SIBLING])
def test_uniform_degradation_scores_perfectly_at_any_depth(tier):
    """Depth is not the signal: an upstream rename downgrades every row identically, and that is the case the permissive ladder exists for."""
    summary = score_id_match_tiers({f"design_{i}": tier for i in range(500)})
    assert summary.score == 1.0
    assert summary.majority_tier == tier


def test_a_rare_minority_scores_worse_than_a_big_one(record_case):
    """The counter-intuitive ordering the maintainer asked for: fewer rows affected is WORSE, because a tier answering a handful of rows is an accident while one answering half the lookup is a second rule."""
    three = score_id_match_tiers(
        {f"d_{i}": (TIER_SIBLING if i < 3 else TIER_EXACT) for i in range(500)}).score
    fifty = score_id_match_tiers(
        {f"d_{i}": (TIER_SIBLING if i < 50 else TIER_EXACT) for i in range(500)}).score
    half = score_id_match_tiers(
        {f"d_{i}": (TIER_SIBLING if i < 250 else TIER_EXACT) for i in range(500)}).score
    record_case(input="3/500 vs 50/500 vs 250/500 sibling",
                expected="3 < 50 < 250", actual=f"{three:.3f} < {fifty:.3f} < {half:.3f}")
    assert three < fifty < half
    assert three < 0.05


def test_no_mixed_lookup_reaches_a_perfect_score():
    """A mixed lookup is a different object from a uniform one, so the top score is reserved."""
    half = score_id_match_tiers(
        {f"d_{i}": (TIER_SIBLING if i < 250 else TIER_EXACT) for i in range(500)})
    assert half.score == pytest.approx(im.MIXED_LOOKUP_CEILING)
    assert half.score < 1.0


def test_minority_ids_are_the_ids_that_left_the_majority():
    tiers = {"a": TIER_EXACT, "b": TIER_EXACT, "c": TIER_SIBLING, "d": TIER_PARENT}
    summary = score_id_match_tiers(tiers)
    assert summary.majority_tier == TIER_EXACT
    assert set(summary.minority_ids) == {"c", "d"}


def test_tier_counts_render_along_the_ladder():
    """Ladder order, not insertion order, so the line reads as the story of the lookup."""
    summary = score_id_match_tiers({"a": TIER_SIBLING, "b": TIER_EXACT, "c": TIER_PARENT})
    assert list(summary.tier_counts) == [TIER_EXACT, TIER_PARENT, TIER_SIBLING]


def test_a_tie_between_tiers_gives_the_majority_to_the_shallower_one():
    summary = score_id_match_tiers({"a": TIER_SIBLING, "b": TIER_EXACT})
    assert summary.majority_tier == TIER_EXACT
    assert summary.minority_ids == ("a",)


def test_unmatched_ids_are_counted_but_kept_out_of_the_score(record_case):
    """An unmatched id returns None — a KeyError or a filtered row, never a wrong value — so it is not evidence of an inconsistent match."""
    tiers = {f"d_{i}": TIER_EXACT for i in range(10)}
    tiers.update({"x": TIER_NONE, "y": TIER_NONE})
    summary = score_id_match_tiers(tiers)
    record_case(input="10 exact + 2 unmatched", expected=(1.0, 2, 10),
                actual=(summary.score, len(summary.unmatched), summary.resolved))
    assert summary.score == 1.0
    assert set(summary.unmatched) == {"x", "y"}
    assert summary.resolved == 10
    assert summary.total == 12


def test_a_lookup_that_matched_nothing_has_no_consistency_to_judge():
    summary = score_id_match_tiers({"x": TIER_NONE, "y": TIER_NONE})
    assert summary.score == 1.0
    assert summary.majority_tier is None
    assert summary.resolved == 0


# -- the reported line -------------------------------------------------------

def _check(tiers, **kwargs):
    summary = score_id_match_tiers(tiers)
    return ce.check_id_match_consistency(
        tier_counts=summary.tier_counts, minority_ids=summary.minority_ids,
        score=summary.score, unmatched=len(summary.unmatched), **kwargs)


def test_a_uniform_lookup_reports_nothing():
    assert _check({f"d_{i}": TIER_SIBLING for i in range(500)}) is None


def test_an_empty_lookup_reports_nothing():
    assert _check({}) is None


def test_the_line_names_the_ids_that_degraded_differently(record_case):
    tiers = {f"design_{i}": TIER_EXACT for i in range(1, 501)}
    for odd in ("design_9", "design_11", "design_12"):
        tiers[odd] = TIER_SIBLING
    violation = _check(tiers)
    assert violation is not None
    record_case(input="497 exact + 3 sibling", expected="497 exact, 3 via sibling (ids named)",
                actual=violation.message)
    assert "500 rows: 497 exact, 3 via sibling" in violation.message
    assert "(design_9, design_11, design_12)" in violation.message
    assert violation.check == "id_match_consistency"


def test_the_line_reports_unmatched_ids_alongside_the_tiers():
    tiers = {f"d_{i}": TIER_EXACT for i in range(10)}
    tiers["odd"] = TIER_SIBLING
    tiers["gone"] = TIER_NONE
    violation = _check(tiers)
    assert "12 rows: 10 exact, 1 via sibling, 1 unmatched" in violation.message


def test_naming_is_capped_so_a_pathological_lookup_stays_one_line():
    tiers = {f"d_{i}": TIER_SIBLING for i in range(2000)}
    tiers.update({f"e_{i}": TIER_EXACT for i in range(3000)})
    violation = _check(tiers)
    assert violation.message.count("d_") == ce.MAX_NAMED_IDS
    assert f"+{2000 - ce.MAX_NAMED_IDS} more" in violation.message


def test_the_line_names_the_calling_context():
    violation = _check({"a": TIER_EXACT, "b": TIER_SIBLING}, where="DistanceSelector")
    assert "(in DistanceSelector)" in violation.message


def test_the_line_picks_up_an_ambient_tool_context():
    with ce.tool_context("Boltz2"):
        violation = _check({"a": TIER_EXACT, "b": TIER_SIBLING})
    assert "(in Boltz2)" in violation.message


def test_reporting_is_log_only(capsys):
    """No gate, no refusal: an inconsistent lookup still returns its matches."""
    sources = [f"design_{i}" for i in range(1, 21)]
    targets = [s for s in sources if s != "design_9"]
    result = get_mapped_ids(sources, targets)
    assert result["design_9"] == "design_1"
    assert result["design_10"] == "design_10"
    assert "id_match_consistency" in capsys.readouterr().err


def test_the_check_obeys_the_severity_ladder(monkeypatch):
    """It is a normal contract row, so a maintainer can promote or silence it."""
    assert ce.severity_for("id_match_consistency") == ce.SEVERITY["id_match_consistency"] == ce.WARN
    monkeypatch.setenv("BIOPIPELINES_ENFORCE_ID_MATCH_CONSISTENCY", "raise")
    with pytest.raises(ce.ContractViolation):
        ce.report(_check({"a": TIER_EXACT, "b": TIER_SIBLING}))


def test_report_helper_returns_the_score_it_reported(capsys):
    summary = im.report_id_match_consistency({"a": TIER_EXACT, "b": TIER_SIBLING})
    assert summary.score == pytest.approx(im.MIXED_LOOKUP_CEILING)
    assert "1 exact, 1 via sibling" in capsys.readouterr().err


# -- the five cases the audit found ------------------------------------------

def test_case_clean_exact_lookup_is_silent():
    ids = [f"design_{i}" for i in range(1, 501)]
    assert score_id_match_tiers(tiers_of(ids, ids)).score == 1.0


def test_case_upstream_rename_degrading_uniformly_scores_well(record_case):
    """Every row loses a suffix, so every row answers at the parent tier: exactly the case the permissiveness exists for."""
    sources = [f"design_{i}_1" for i in range(1, 501)]
    targets = [f"design_{i}" for i in range(1, 501)]
    summary = score_id_match_tiers(tiers_of(sources, targets))
    record_case(input="500 renamed ids, all parent", expected=1.0, actual=summary.score)
    assert summary.score == 1.0
    assert summary.tier_counts == {TIER_PARENT: 500}


def test_case_design_9_taking_design_1s_row_scores_badly(record_case):
    """The audit's headline case: design_9 is absent from the table, so the sibling tier hands it design_1's value while its neighbors match exactly."""
    sources = [f"design_{i}" for i in range(1, 501)]
    targets = [s for s in sources if s not in {"design_9", "design_11", "design_12"}]
    resolved = get_mapped_ids_with_tiers(sources, targets)
    assert resolved["design_9"] == ("design_1", TIER_SIBLING)
    summary = score_id_match_tiers({s: t for s, (_v, t) in resolved.items()})
    record_case(input="497 exact + 3 sibling", expected="< 0.05", actual=summary.score)
    assert summary.score < 0.05
    assert set(summary.minority_ids) == {"design_9", "design_11", "design_12"}


def test_case_a_multi_axis_id_degrading_on_one_axis():
    """`prot1+lig3` has no row, so the shared `prot1` axis pulls in `prot1+lig1`'s box."""
    sources = [f"prot1+lig{i}" for i in range(1, 501)]
    targets = [s for s in sources if s != "prot1+lig3"]
    resolved = get_mapped_ids_with_tiers(sources, targets)
    assert resolved["prot1+lig3"] == ("prot1+lig1", TIER_SIBLING)
    summary = score_id_match_tiers({s: t for s, (_v, t) in resolved.items()})
    assert summary.score < 0.05
    assert summary.minority_ids == ("prot1+lig3",)


def test_case_numeric_ids_do_reach_the_provenance_tier(tmp_path, record_case):
    """A map_table whose provenance cells are ints still resolves: the index coerces to str, so int64 `1` indexes as `'1'`."""
    sources = [f"Panda_{i}" for i in range(1, 11)]
    targets = [str(i) for i in range(1, 11)]
    path = write_map(tmp_path, "ints.csv", sources, list(range(1, 11)))
    summary = score_id_match_tiers(tiers_of(sources, targets, map_table_paths=[path]))
    record_case(input="int provenance column", expected={TIER_PROVENANCE: 10},
                actual=summary.tier_counts)
    assert summary.tier_counts == {TIER_PROVENANCE: 10}
    assert summary.score == 1.0


def test_case_one_blank_provenance_cell_scatters_the_degradation(tmp_path, record_case):
    """The real numeric hazard: a single missing cell drops one row out of the provenance tier and into the sibling tier, which is the scattered pattern the score is built to catch."""
    sources = [f"{i}_2" for i in range(1, 11)]
    targets = [f"{i}_1" for i in range(1, 11)]
    provenance = list(targets)
    provenance[9] = None
    path = write_map(tmp_path, "partial.csv", sources, provenance)
    summary = score_id_match_tiers(tiers_of(sources, targets, map_table_paths=[path]))
    record_case(input="9 provenance + 1 blank cell",
                expected={TIER_PROVENANCE: 9, TIER_SIBLING: 1}, actual=summary.tier_counts)
    assert summary.tier_counts == {TIER_PROVENANCE: 9, TIER_SIBLING: 1}
    assert summary.minority_ids == ("10_2",)
    assert summary.score < 0.25


def test_case_nan_coerced_numeric_provenance_is_unreachable(tmp_path):
    """With a blank cell pandas reads the whole column as float64, so every cell indexes as `'1.0'` and matches nothing — the tier is unreachable for ALL rows, not just the blank one. Uniformly unmatched, so the consistency score has nothing to say; the KeyError each id raises is the report."""
    sources = [f"Panda_{i}" for i in range(1, 11)]
    targets = [str(i) for i in range(1, 11)]
    path = write_map(tmp_path, "floats.csv", sources, list(range(1, 10)) + [None])
    summary = score_id_match_tiers(tiers_of(sources, targets, map_table_paths=[path]))
    assert summary.resolved == 0
    assert len(summary.unmatched) == 10


# -- through a real table lookup ---------------------------------------------

def test_iterate_table_values_scores_the_whole_walk(capsys):
    """One id's tier says nothing; the set's consistency does, so the line is emitted once the walk finishes."""
    from biopipelines.biopipelines_io import iterate_table_values

    ids = [f"design_{i}" for i in range(1, 21)]
    rows = [i for i in ids if i != "design_9"]
    table = pd.DataFrame({"id": rows, "within": [f"val_{i}" for i in range(len(rows))]})
    values = dict(iterate_table_values(table, ids, "within"))
    assert values["design_9"] == values["design_1"]
    err = capsys.readouterr().err
    assert "20 rows: 19 exact, 1 via sibling (design_9)" in err
    assert "table lookup of column 'within'" in err


def test_a_single_lookup_cannot_be_inconsistent(capsys):
    """`lookup_table_value` resolves one id, so its own tier mix is uniform by construction and it stays quiet."""
    from biopipelines.biopipelines_io import lookup_table_value

    table = pd.DataFrame({"id": ["design_1"], "within": ["A"]})
    assert lookup_table_value(table, "design_9", "within") == "A"
    assert "id_match_consistency" not in capsys.readouterr().err
