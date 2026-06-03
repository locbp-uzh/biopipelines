"""Unit tests for biopipelines.id_map_utils.

These cover get_mapped_ids and friends directly. Until this file existed
the matcher was only exercised indirectly through test_panda.py and
test_remap.py, so refactors of the inner loops (e.g. inverted-index
optimization) had no targeted regression coverage.

The test cases mirror the priority order documented in get_mapped_ids
docstring: exact > provenance > child/parent > sibling > none.
"""
from __future__ import annotations

import os

import pandas as pd
import pytest


@pytest.fixture(autouse=True)
def _reset_caches():
    """Clear the per-process provenance + table caches between tests so
    one test's CSV doesn't leak into another's lookup."""
    from biopipelines import id_map_utils
    from biopipelines import biopipelines_io

    id_map_utils._provenance_lookup_cache.clear()
    biopipelines_io._table_cache.clear()
    biopipelines_io._provenance_index_cache.clear()
    yield
    id_map_utils._provenance_lookup_cache.clear()
    biopipelines_io._table_cache.clear()
    biopipelines_io._provenance_index_cache.clear()


# ── map_table_ids_to_ids: suffix-stripping helper ─────────────────────────────

def test_map_table_ids_to_ids_basic_numeric_suffix():
    from biopipelines.id_map_utils import map_table_ids_to_ids

    assert map_table_ids_to_ids("rifampicin_1_2", {"*": "*_<N>"}) == [
        "rifampicin_1_2", "rifampicin_1", "rifampicin",
    ]


def test_map_table_ids_to_ids_segment_suffix():
    from biopipelines.id_map_utils import map_table_ids_to_ids

    # <S> matches any segment, including alphanumeric suffixes.
    assert map_table_ids_to_ids("protein_1_19A", {"*": "*_<S>"}) == [
        "protein_1_19A", "protein_1", "protein",
    ]


def test_map_table_ids_to_ids_identity_pattern():
    from biopipelines.id_map_utils import map_table_ids_to_ids

    assert map_table_ids_to_ids("no_change", {"*": "*"}) == ["no_change"]


def test_map_table_ids_to_ids_multi_axis_plus_separator():
    from biopipelines.id_map_utils import map_table_ids_to_ids

    # Multi-axis IDs (joined with '+') generate sub-sequences and
    # suffix-stripped variants.
    out = map_table_ids_to_ids("prot1+lig1_2", {"*": "*_<S>"})
    # Order is most-specific-first; we don't pin the full enumeration
    # but the most-specific must come first and key prefixes must be
    # present.
    assert out[0] == "prot1+lig1_2"
    assert "prot1+lig1" in out
    assert "prot1" in out
    assert "lig1_2" in out
    assert "lig1" in out


# ── get_mapped_ids: priority 1 — exact match ──────────────────────────────────

def test_get_mapped_ids_exact_match_unique():
    from biopipelines.id_map_utils import get_mapped_ids

    out = get_mapped_ids(["protein_1"], ["protein_1"])
    assert out == {"protein_1": "protein_1"}


def test_get_mapped_ids_exact_match_includes_other_descendants_when_not_unique():
    """unique=False at exact match collects exact + targets whose base list
    contains the source id (i.e. children of the exact match)."""
    from biopipelines.id_map_utils import get_mapped_ids

    out = get_mapped_ids(
        ["protein_1"],
        ["protein_1", "protein_1_a", "protein_1_b"],
        unique=False,
    )
    matches = out["protein_1"]
    assert matches[0] == "protein_1"  # exact wins first slot
    assert set(matches) == {"protein_1", "protein_1_a", "protein_1_b"}


# ── get_mapped_ids: priority 3a — child match (target = source + suffix) ──────

def test_get_mapped_ids_child_match_unique_picks_closest():
    from biopipelines.id_map_utils import get_mapped_ids

    # protein_1 has two descendants; closer (shorter distance) wins.
    out = get_mapped_ids(
        ["protein_1"],
        ["protein_1_1", "protein_1_1_a"],
    )
    # protein_1 sits at distance 1 in protein_1_1's base list and at
    # distance 2 in protein_1_1_a's base list — protein_1_1 wins.
    assert out == {"protein_1": "protein_1_1"}


def test_get_mapped_ids_child_match_not_unique_returns_all():
    from biopipelines.id_map_utils import get_mapped_ids

    out = get_mapped_ids(
        ["protein_1"],
        ["protein_1_1", "protein_1_2"],
        unique=False,
    )
    assert set(out["protein_1"]) == {"protein_1_1", "protein_1_2"}


# ── get_mapped_ids: priority 3b — parent match (source = target + suffix) ─────

def test_get_mapped_ids_parent_match():
    from biopipelines.id_map_utils import get_mapped_ids

    out = get_mapped_ids(["protein_1_1"], ["protein_1"])
    assert out == {"protein_1_1": "protein_1"}


def test_get_mapped_ids_parent_match_picks_closest_when_two_ancestors():
    from biopipelines.id_map_utils import get_mapped_ids

    out = get_mapped_ids(["protein_1_1"], ["protein_1", "protein"])
    # protein_1 is at distance 1 from protein_1_1, protein at distance 2
    assert out == {"protein_1_1": "protein_1"}


# ── get_mapped_ids: priority 4 — sibling match (common ancestor) ──────────────

def test_get_mapped_ids_sibling_match():
    """Source and target share an ancestor but neither contains the other."""
    from biopipelines.id_map_utils import get_mapped_ids

    out = get_mapped_ids(["protein_1_1"], ["protein_1_2"])
    assert out == {"protein_1_1": "protein_1_2"}


def test_get_mapped_ids_sibling_prefers_closer_common_ancestor():
    """Among multiple siblings, the one with the smallest combined distance wins."""
    from biopipelines.id_map_utils import get_mapped_ids

    # Source: protein_1_1 -> bases [protein_1_1, protein_1, protein]
    # Targets:
    #   protein_1_2  -> [protein_1_2, protein_1, protein]   (1+1=2 via protein_1)
    #   protein_2_1  -> [protein_2_1, protein_2, protein]   (2+2=4 via protein)
    out = get_mapped_ids(["protein_1_1"], ["protein_1_2", "protein_2_1"])
    assert out == {"protein_1_1": "protein_1_2"}


# ── get_mapped_ids: closest_siblings_only (design-group matching) ─────────────

def test_closest_siblings_only_returns_same_parent_group():
    """A source returns ALL targets sharing its immediate parent (design group),
    including itself, and nothing from distant lineages."""
    from biopipelines.id_map_utils import get_mapped_ids

    subs = ["protein_29_1", "protein_29_2", "protein_2_1", "protein_2_2"]
    out = get_mapped_ids(["protein_29_1"], subs, unique=False,
                         closest_siblings_only=True)
    assert sorted(out["protein_29_1"]) == ["protein_29_1", "protein_29_2"]


def test_closest_siblings_only_no_cross_design_overmatch():
    """A source whose exact id is absent matches only its same-parent siblings,
    not the multi-segment-collapsed lineages (29_2 must not pull 2_*)."""
    from biopipelines.id_map_utils import get_mapped_ids

    # protein_29 made only one substitution sample (29_1); 29_2 has no exact match.
    subs = ["protein_29_1", "protein_2_1", "protein_2_2", "protein_49_1"]
    out = get_mapped_ids(["protein_29_2"], subs, unique=False,
                         closest_siblings_only=True)
    assert out["protein_29_2"] == ["protein_29_1"]


def test_closest_siblings_only_empty_when_no_same_parent():
    """No target shares the source's immediate parent -> empty, never a distant id."""
    from biopipelines.id_map_utils import get_mapped_ids

    subs = ["protein_2_1", "protein_2_2"]
    out = get_mapped_ids(["protein_99_1"], subs, unique=False,
                         closest_siblings_only=True)
    assert out["protein_99_1"] == []


def test_closest_siblings_only_default_off_keeps_legacy_overmatch():
    """Without the flag, the sibling tier still returns the broad match set."""
    from biopipelines.id_map_utils import get_mapped_ids

    subs = ["protein_29_1", "protein_2_1", "protein_2_2"]
    out = get_mapped_ids(["protein_29_2"], subs, unique=False)
    # legacy behavior: sibling tier returns more than just the same-parent match
    assert len(out["protein_29_2"]) > 1
    assert "protein_29_1" in out["protein_29_2"]


# ── get_mapped_ids: no match ──────────────────────────────────────────────────

def test_get_mapped_ids_no_match_returns_none():
    from biopipelines.id_map_utils import get_mapped_ids

    out = get_mapped_ids(["alpha_1"], ["beta_1"])
    assert out == {"alpha_1": None}


def test_get_mapped_ids_no_match_returns_empty_list_when_not_unique():
    from biopipelines.id_map_utils import get_mapped_ids

    out = get_mapped_ids(["alpha_1"], ["beta_1"], unique=False)
    assert out == {"alpha_1": []}


# ── get_mapped_ids: multi-source preserves per-source independence ────────────

def test_get_mapped_ids_multiple_sources_independent_results():
    from biopipelines.id_map_utils import get_mapped_ids

    out = get_mapped_ids(
        ["protein_1", "protein_2"],
        ["protein_1_a", "protein_2_b"],
    )
    assert out == {"protein_1": "protein_1_a", "protein_2": "protein_2_b"}


# ── _get_provenance_ids + provenance-based matching ──────────────────────────

def test_provenance_ids_reads_dotid_columns(tmp_path):
    from biopipelines.id_map_utils import _get_provenance_ids

    csv = tmp_path / "map.csv"
    pd.DataFrame({
        "id": ["Panda_1", "Panda_2"],
        "structures.id": ["LID_001", "LID_002"],
        "sequences.id": ["seq_001", "seq_002"],
    }).to_csv(csv, index=False)

    assert sorted(_get_provenance_ids("Panda_1", [str(csv)])) == [
        "LID_001", "seq_001",
    ]
    assert sorted(_get_provenance_ids("Panda_2", [str(csv)])) == [
        "LID_002", "seq_002",
    ]
    # Unknown id -> empty list
    assert _get_provenance_ids("Panda_3", [str(csv)]) == []


def test_provenance_ids_skips_self_and_nan(tmp_path):
    from biopipelines.id_map_utils import _get_provenance_ids

    csv = tmp_path / "map.csv"
    pd.DataFrame({
        "id": ["a"],
        "structures.id": ["a"],   # self → skipped
        "sequences.id": [None],   # NaN → skipped
    }).to_csv(csv, index=False)

    assert _get_provenance_ids("a", [str(csv)]) == []


def test_provenance_ids_cached_across_calls(tmp_path):
    """The lookup index is built once per CSV — re-reading the file (with
    different contents) after the first call should NOT change the answer
    because the cached index is authoritative."""
    from biopipelines.id_map_utils import _get_provenance_ids

    csv = tmp_path / "map.csv"
    pd.DataFrame({
        "id": ["x"],
        "structures.id": ["original_x"],
    }).to_csv(csv, index=False)

    first = _get_provenance_ids("x", [str(csv)])
    assert first == ["original_x"]

    # Mutate the file: a fresh read would yield a different value, but the
    # cached lookup must keep returning the original.
    pd.DataFrame({
        "id": ["x"],
        "structures.id": ["MUTATED"],
    }).to_csv(csv, index=False)

    second = _get_provenance_ids("x", [str(csv)])
    assert second == ["original_x"]


def test_get_mapped_ids_provenance_match_via_dotid(tmp_path):
    """Panda auto-renames source to Panda_1; provenance column ties it back to
    the underlying LID_001 family which has children in target_ids."""
    from biopipelines.id_map_utils import get_mapped_ids

    csv = tmp_path / "panda_map.csv"
    pd.DataFrame({
        "id": ["Panda_1"],
        "structures.id": ["LID_001"],
    }).to_csv(csv, index=False)

    out = get_mapped_ids(
        source_ids=["Panda_1"],
        target_ids=["LID_001_1"],
        map_table_paths=[str(csv)],
    )
    assert out == {"Panda_1": "LID_001_1"}


# ── post-refactor invariants ──────────────────────────────────────────────────

def test_get_mapped_ids_inverted_index_handles_repeated_bases():
    """Multiple targets may share a base id. The inverted index must list
    them all, and tie-break by distance (most specific first)."""
    from biopipelines.id_map_utils import get_mapped_ids

    # Both protein_1_a and protein_1_b have protein_1 as their direct
    # parent (distance 1), so unique=True picks the first inserted (which
    # is protein_1_a since it appears first in target_ids).
    out_unique = get_mapped_ids(
        ["protein_1"],
        ["protein_1_a", "protein_1_b"],
    )
    assert out_unique == {"protein_1": "protein_1_a"}

    out_all = get_mapped_ids(
        ["protein_1"],
        ["protein_1_a", "protein_1_b"],
        unique=False,
    )
    assert set(out_all["protein_1"]) == {"protein_1_a", "protein_1_b"}


def test_get_mapped_ids_distance_tiebreak_is_deterministic():
    """When two targets are equally close, the first listed in target_ids wins
    under unique=True. (Stable insertion order preserved through the
    inverted index.)"""
    from biopipelines.id_map_utils import get_mapped_ids

    out_first = get_mapped_ids(["root"], ["root_x", "root_y"])
    assert out_first == {"root": "root_x"}

    out_swapped = get_mapped_ids(["root"], ["root_y", "root_x"])
    assert out_swapped == {"root": "root_y"}


def test_get_mapped_ids_handles_empty_inputs():
    from biopipelines.id_map_utils import get_mapped_ids

    assert get_mapped_ids([], ["a"]) == {}
    # Empty target list -> no match for the source.
    assert get_mapped_ids(["a"], []) == {"a": None}


# ── prune_redundant_provenance_columns ────────────────────────────────────────
#
# Redundancy is defined through the matcher itself: a provenance cell p in
# column <axis>.id is redundant for row id=x iff get_mapped_ids([x],[p],
# map_table_paths=None) resolves x to p — i.e. id-semantics alone recover it
# with no provenance lookup. A whole column is dropped only when every
# non-empty cell is redundant.


def test_prune_drops_multiplier_provenance():
    """structure_1 -> structure_1_3: structures.id is recoverable, dropped."""
    from biopipelines.id_map_utils import prune_redundant_provenance_columns

    df = pd.DataFrame({
        "id": ["structure_1_1", "structure_1_2", "structure_1_3"],
        "structures.id": ["structure_1", "structure_1", "structure_1"],
        "score": [1, 2, 3],
    })
    out = prune_redundant_provenance_columns(df)
    assert list(out.columns) == ["id", "score"]


def test_prune_keeps_unrecoverable_rename():
    """structure_1 -> Panda_1: rename is not id-recoverable, column kept."""
    from biopipelines.id_map_utils import prune_redundant_provenance_columns

    df = pd.DataFrame({
        "id": ["Panda_1", "Panda_2"],
        "structures.id": ["a", "b"],
    })
    out = prune_redundant_provenance_columns(df)
    assert "structures.id" in out.columns


def test_prune_drops_multi_axis_when_matcher_resolves():
    """a+b with a.id,b.id: both recoverable from the joined id, both dropped.

    This is the case the naive 'differs from id' rule mishandles: neither
    'a' nor 'b' equals 'a+b' as a string, yet the matcher recovers both via
    +-component expansion, so both columns are genuine no-ops.
    """
    from biopipelines.id_map_utils import prune_redundant_provenance_columns

    df = pd.DataFrame({
        "id": ["a+b", "a+c"],
        "proteins.id": ["a", "a"],
        "ligands.id": ["b", "c"],
    })
    out = prune_redundant_provenance_columns(df)
    assert list(out.columns) == ["id"]


def test_prune_keeps_mixed_column():
    """One row recoverable, one renamed -> whole column kept (per-column rule)."""
    from biopipelines.id_map_utils import prune_redundant_provenance_columns

    df = pd.DataFrame({
        "id": ["structure_1_1", "Panda_9"],
        "structures.id": ["structure_1", "zzz"],
    })
    out = prune_redundant_provenance_columns(df)
    assert "structures.id" in out.columns


def test_prune_ignores_empty_cells():
    """An empty provenance cell does not block pruning an otherwise-redundant
    column (every NON-empty value is recoverable)."""
    from biopipelines.id_map_utils import prune_redundant_provenance_columns

    df = pd.DataFrame({
        "id": ["structure_1_1", "structure_1_2"],
        "structures.id": ["structure_1", ""],
    })
    out = prune_redundant_provenance_columns(df)
    assert "structures.id" not in out.columns


def test_prune_excludes_generation_pool_and_original_columns():
    """Generation (.-N.id), Bundle (.0.id), pool.id/pool.path and original.id
    are never pruned, even when id-recoverable — they encode lineage the
    matcher cannot reconstruct from the current id."""
    from biopipelines.id_map_utils import prune_redundant_provenance_columns

    df = pd.DataFrame({
        "id": ["a_1"],
        "structures.-1.id": ["a"],
        "structures.0.id": ["a"],
        "pool.id": ["a"],
        "pool.path": ["/some/origin"],
        "original.id": ["a"],
        "structures.id": ["a"],   # plain + recoverable -> the only drop
    })
    out = prune_redundant_provenance_columns(df)
    assert list(out.columns) == [
        "id", "structures.-1.id", "structures.0.id",
        "pool.id", "pool.path", "original.id",
    ]


def test_prune_is_pure_no_mutation():
    """The helper returns a new frame and never mutates its input."""
    from biopipelines.id_map_utils import prune_redundant_provenance_columns

    df = pd.DataFrame({
        "id": ["structure_1_1"],
        "structures.id": ["structure_1"],
    })
    before = list(df.columns)
    prune_redundant_provenance_columns(df)
    assert list(df.columns) == before
