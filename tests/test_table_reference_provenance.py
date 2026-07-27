"""Regression tests for table-reference id matching and its failure modes.

Two related concerns:

1. ``lookup_table_value`` must resolve an id whose string relationship to the
   table's ids was destroyed by an upstream rename (Panda_1 vs design_3). The
   mapping that recovers it lives in the *consumer's* stream map_table, not in
   the referenced table, so provenance only works when those paths are supplied.
2. ``pipe_distance_selector.parse_restrict_spec`` must raise on an unresolvable
   restriction rather than returning a permissive or empty restriction — the
   silent forms turned a matching failure into structurally wrong selections.
"""

import os
import sys

import pandas as pd
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "pipe_scripts"))

from biopipelines.biopipelines_io import load_table, lookup_table_value, clear_table_cache
from biopipelines.id_map_utils import get_mapped_ids


@pytest.fixture(autouse=True)
def _clear_cache():
    clear_table_cache()
    yield
    clear_table_cache()


def _renamed_setup(tmp_path):
    """A selections table on original ids + a stream map that renamed them."""
    table = tmp_path / "selections.csv"
    pd.DataFrame({
        "id": ["design_1", "design_2", "design_3"],
        "within": ["A10+A11", "A20+A21", "A30+A31"],
    }).to_csv(table, index=False)

    stream_map = tmp_path / "panda_map.csv"
    pd.DataFrame({
        "id": ["Panda_1", "Panda_2", "Panda_3"],
        "file": ["/x/1.pdb", "/x/2.pdb", "/x/3.pdb"],
        "structures.id": ["design_3", "design_1", "design_2"],
    }).to_csv(stream_map, index=False)
    return str(table), str(stream_map)


# ── id matching across a rename ──────────────────────────────────────────────

def test_renamed_ids_unmatchable_without_provenance(tmp_path, record_case):
    """Suffix matching alone cannot relate Panda_1 to design_3."""
    matched = get_mapped_ids(["Panda_1"], ["design_1", "design_2", "design_3"])
    record_case(input="get_mapped_ids(['Panda_1'], ['design_1..3'])",
                expected={"Panda_1": None}, actual=matched)
    assert matched == {"Panda_1": None}


def test_stream_map_provenance_recovers_renamed_ids(tmp_path, record_case):
    """The consumer's stream map routes each renamed id to its original."""
    _, stream_map = _renamed_setup(tmp_path)
    matched = get_mapped_ids(
        ["Panda_1", "Panda_2", "Panda_3"],
        ["design_1", "design_2", "design_3"],
        map_table_paths=[stream_map],
    )
    expected = {"Panda_1": "design_3", "Panda_2": "design_1", "Panda_3": "design_2"}
    record_case(input="get_mapped_ids(renamed, originals, map_table_paths=[stream_map])",
                expected=expected, actual=matched)
    assert matched == expected


def test_referenced_table_provenance_is_the_wrong_direction(tmp_path, record_case):
    """Provenance ON the referenced table does not resolve a downstream rename.

    The table's own ``<stream>.id`` column records where its rows came from
    (design_1 -> parent_a); it says nothing about ids minted after it was
    written, so it cannot map Panda_1 back into the table's id space.
    """
    table = tmp_path / "selections_prov.csv"
    pd.DataFrame({
        "id": ["design_1", "design_2", "design_3"],
        "structures.id": ["parent_a", "parent_b", "parent_c"],
        "within": ["A10+A11", "A20+A21", "A30+A31"],
    }).to_csv(table, index=False)

    matched = get_mapped_ids(["Panda_1"], ["design_1", "design_2", "design_3"],
                             map_table_paths=[str(table)])
    record_case(input="get_mapped_ids(['Panda_1'], originals, map_table_paths=[the table itself])",
                expected={"Panda_1": None}, actual=matched)
    assert matched == {"Panda_1": None}


def test_lookup_table_value_raises_on_unmatched_id(tmp_path, record_case):
    """A renamed id that cannot be matched raises, naming id and available ids."""
    table_path, _ = _renamed_setup(tmp_path)
    table, column = load_table(f"TABLE_REFERENCE:{table_path}:within")

    with pytest.raises(KeyError) as exc:
        lookup_table_value(table, "Panda_1", column)
    msg = str(exc.value)
    record_case(input="lookup_table_value(table, 'Panda_1', 'within')",
                expected="KeyError naming id + available ids", actual=msg)
    assert "Panda_1" in msg and "design_1" in msg
    # The message must say provenance was never consulted, so the fix is findable.
    assert "map_table_paths" in msg


def test_lookup_table_value_resolves_renamed_id_with_map_tables(tmp_path, record_case):
    """The fix: forwarding the stream's map_table makes the renamed id resolve.

    This is the regression for the original defect — lookup_table_value did not
    accept map_table_paths at all, so the provenance tier of get_mapped_ids was
    unreachable and every renamed id raised.
    """
    table_path, stream_map = _renamed_setup(tmp_path)
    table, column = load_table(f"TABLE_REFERENCE:{table_path}:within")

    got = {sid: lookup_table_value(table, sid, column, map_table_paths=[stream_map])
           for sid in ("Panda_1", "Panda_2", "Panda_3")}
    expected = {"Panda_1": "A30+A31", "Panda_2": "A10+A11", "Panda_3": "A20+A21"}
    record_case(input="lookup_table_value(..., map_table_paths=[stream_map])",
                expected=expected, actual=got)
    assert got == expected


def test_restrict_spec_resolves_renamed_id_via_stream_map(tmp_path, record_case):
    """DistanceSelector's restriction resolves a renamed id through provenance."""
    from pipe_distance_selector import parse_restrict_spec

    table = tmp_path / "restrict.csv"
    pd.DataFrame({"id": ["design_1", "design_2", "design_3"],
                  "within": ["A1", "A2", "A1+A2"]}).to_csv(table, index=False)
    stream_map = tmp_path / "structures_map.csv"
    pd.DataFrame({"id": ["Panda_1"], "file": ["/x/Panda_1.pdb"],
                  "structures.id": ["design_3"]}).to_csv(stream_map, index=False)

    ref = f"TABLE_REFERENCE:{table}:within"
    chained, chainless, chains = parse_restrict_spec(
        ref, "Panda_1", {"A"}, [str(stream_map)])
    record_case(input="parse_restrict_spec('Panda_1', map_table_paths=[stream_map])",
                expected={("A", 1), ("A", 2)}, actual=chained)
    assert chained == {("A", 1), ("A", 2)}

    # Without the map the same call must raise rather than silently mis-restrict.
    with pytest.raises(KeyError):
        parse_restrict_spec(ref, "Panda_1", {"A"})


# ── the restriction-spec error branches ──────────────────────────────────────

def test_restrict_spec_raises_on_missing_table(tmp_path, record_case):
    """A missing restriction table must raise, not mean 'no restriction'."""
    from pipe_distance_selector import parse_restrict_spec

    ref = f"TABLE_REFERENCE:{tmp_path / 'nope.csv'}:within"
    with pytest.raises(FileNotFoundError) as exc:
        parse_restrict_spec(ref, "design_1", {"A"})
    record_case(input="parse_restrict_spec(<missing table>)",
                expected="FileNotFoundError", actual=type(exc.value).__name__)


def test_restrict_spec_raises_on_missing_id(tmp_path, record_case):
    """A missing id must raise, not mean 'nothing passes'.

    Returning empty sets here is what destroyed the catalytic residues: every
    residue failed the restriction, so the designed set was structurally wrong
    while the run reported success.
    """
    from pipe_distance_selector import parse_restrict_spec

    table = tmp_path / "restrict.csv"
    pd.DataFrame({"id": ["somethingelse"], "within": ["A1+A2"]}).to_csv(table, index=False)

    with pytest.raises(KeyError) as exc:
        parse_restrict_spec(f"TABLE_REFERENCE:{table}:within", "design_1", {"A"})
    record_case(input="parse_restrict_spec(<id absent from table>)",
                expected="KeyError", actual=type(exc.value).__name__)
    assert "design_1" in str(exc.value)


def test_restrict_spec_resolves_matching_id(tmp_path, record_case):
    """The success path still returns the parsed restriction."""
    from pipe_distance_selector import parse_restrict_spec

    table = tmp_path / "restrict.csv"
    pd.DataFrame({"id": ["design_1"], "within": ["A1+A2"]}).to_csv(table, index=False)

    chained, chainless, chains = parse_restrict_spec(
        f"TABLE_REFERENCE:{table}:within", "design_1", {"A"})
    record_case(input="parse_restrict_spec(<matching id>)",
                expected=({("A", 1), ("A", 2)}, None, None),
                actual=(chained, chainless, chains))
    assert chained == {("A", 1), ("A", 2)}
    assert chainless is None and chains is None


def test_empty_restrict_spec_means_no_restriction(tmp_path, record_case):
    """An unset restriction is still legitimately 'no restriction'."""
    from pipe_distance_selector import parse_restrict_spec

    result = parse_restrict_spec("", "design_1", {"A"})
    record_case(input="parse_restrict_spec('')", expected=(None, None, None), actual=result)
    assert result == (None, None, None)
