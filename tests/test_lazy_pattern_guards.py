"""Regression tests for the "is a pattern" vs "can be expanded now" split.

`contains_pattern` matches a `<..>` slot wherever it sits, brackets included, so it answers *"is this a pattern?"*. Three pipe scripts used it to guard `expand_ids`, which asks the other question — *"can I turn this into concrete ids right now?"* — and a lazy id passed the guard only to raise `LazyPatternError` inside it. `can_expand` / `can_expand_ids` answer that second question; `is_pattern` answers the first including a bracket with no slot in it.
"""

import importlib.util
import os

import pandas as pd
import pytest

from biopipelines import contract_enforcement
from biopipelines import id_patterns as idp


def _load_pipe_script(name):
    repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    spec = importlib.util.spec_from_file_location(
        name, os.path.join(repo, "pipe_scripts", f"{name}.py")
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ── the predicates ───────────────────────────────────────────────────────────

@pytest.mark.parametrize("value, contains_pattern, is_pattern, can_expand, is_literal", [
    ("literal",         False, False, False, True),
    ("a_<0..1>",        True,  True,  True,  False),
    ("prot[_<?>]",      True,  True,  False, False),
    ("prot_[<?>]",      True,  True,  False, False),
    ("a_<1 2>[_<?>]",   True,  True,  False, False),
    ("p[_x]",           False, True,  False, False),
])
def test_pattern_predicates(record_case, value, contains_pattern, is_pattern,
                            can_expand, is_literal):
    actual = (idp.contains_pattern(value), idp.is_pattern(value),
              idp.can_expand(value), idp.is_literal(value))
    expected = (contains_pattern, is_pattern, can_expand, is_literal)
    record_case(input=value, expected=expected, actual=actual)
    assert actual == expected


@pytest.mark.parametrize("ids, expected", [
    (["a_<0..1>"],                  True),
    (["a_<0..1>", "literal"],       True),
    (["literal"],                   False),
    (["prot[_<?>]"],                False),
    (["a_<0..1>", "prot[_<?>]"],    False),   # mixed: expanding would still raise
])
def test_can_expand_ids_is_the_guard_for_expand_ids(record_case, ids, expected):
    actual = idp.can_expand_ids(ids)
    record_case(input=ids, expected=expected, actual=actual)
    assert actual is expected
    if actual:
        idp.expand_ids(ids)  # the guard's promise: this does not raise


def test_contains_pattern_alone_is_not_a_safe_expand_guard(record_case):
    """The exact confusion the three guard sites had."""
    lazy = "prot_<1 2>[_<?>]"
    record_case(input=lazy, expected="contains_pattern=True, can_expand=False",
                actual=f"contains_pattern={idp.contains_pattern(lazy)}, "
                       f"can_expand={idp.can_expand(lazy)}")
    assert idp.contains_pattern(lazy) is True
    assert idp.can_expand(lazy) is False
    with pytest.raises(idp.LazyPatternError):
        idp.expand_ids([lazy])


# ── guard site 1 + 2: pipe_fuse_queries ──────────────────────────────────────

def _fuse_slot(tmp_path):
    table = tmp_path / "seqs.csv"
    pd.DataFrame([
        {"id": "prot_1_5N", "sequence": "AAAA"},
        {"id": "prot_2_7N", "sequence": "CCCC"},
        {"id": "other_1",   "sequence": "GGGG"},
    ]).to_csv(table, index=False)
    return {"ids": ["prot_<1 2>[_<?>]"], "map_table": str(table), "files": []}


def test_fuse_queries_lazy_slot_selects_map_rows(tmp_path, record_case):
    mod = _load_pipe_script("pipe_fuse_queries")
    actual = mod.load_sequences_from_slot(_fuse_slot(tmp_path))
    expected = {"prot_1_5N": "AAAA", "prot_2_7N": "CCCC"}
    record_case(input="lazy slot ids", expected=expected, actual=actual)
    assert actual == expected


def test_fuse_queries_lazy_slot_ids_reach_the_product(tmp_path, record_case):
    mod = _load_pipe_script("pipe_fuse_queries")
    slot = _fuse_slot(tmp_path)
    results = mod.generate_fusion_sequences([slot, slot], "GS", ["2"], "fus")
    actual = sorted(r["id"] for r in results)
    expected = ["prot_1_5N+2+prot_1_5N", "prot_1_5N+2+prot_2_7N",
                "prot_2_7N+2+prot_1_5N", "prot_2_7N+2+prot_2_7N"]
    record_case(input="lazy slots x lazy slots", expected=expected, actual=actual)
    assert actual == expected


# ── guard site 3: pipe_pdb ───────────────────────────────────────────────────

def test_pipe_pdb_lazy_ids_resolve_from_the_upstream_map(tmp_path, record_case):
    mod = _load_pipe_script("pipe_pdb")
    src = tmp_path / "upstream"
    src.mkdir()
    pdb = src / "prot_1_5N.pdb"
    pdb.write_text(
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\n"
        "END\n"
    )
    upstream_map = tmp_path / "upstream_map.csv"
    pd.DataFrame([{"id": "prot_1_5N", "file": str(pdb)}]).to_csv(upstream_map, index=False)
    out = tmp_path / "out"

    failures = mod.fetch_structures({
        "pdb_ids": ["prot_<1>[_<?>]"],
        "custom_ids": ["prot_<1>[_<?>]"],
        "from_upstream": True,
        "upstream_map_table": str(upstream_map),
        "repo_pdbs_folder": str(tmp_path / "pdbs"),
        "output_folder": str(out),
        "structures_table": str(tmp_path / "structures.csv"),
        "sequences_table": str(tmp_path / "sequences.csv"),
        "failed_table": str(tmp_path / "failed.csv"),
        "fetch_compounds": False,
    })
    actual = pd.read_csv(tmp_path / "structures.csv")["id"].astype(str).tolist()
    record_case(input="lazy pdb_ids + upstream map", expected=["prot_1_5N"], actual=actual)
    assert failures == 0
    assert actual == ["prot_1_5N"]


# ── guard site 4: pipe_remap ─────────────────────────────────────────────────

def test_pipe_remap_lazy_files_keep_their_brackets(record_case):
    mod = _load_pipe_script("pipe_remap")
    stream = {"ids": ["a_1"], "files": ["d/x[_<?>].pdb", "d/y[_<?>].pdb"]}
    actual = mod.expand_stream_ids(dict(stream))["files"]
    record_case(input=stream["files"], expected=stream["files"], actual=actual)
    assert actual == stream["files"]


def test_pipe_remap_lazy_ids_select_only_covered_map_rows(tmp_path, record_case):
    """resolve_pattern_ids wired in: the pattern selects, it does not take every row."""
    mod = _load_pipe_script("pipe_remap")
    table = tmp_path / "map.csv"
    pd.DataFrame([
        {"id": "prot_1_5N"}, {"id": "prot_2_7N"}, {"id": "other_1"},
    ]).to_csv(table, index=False)
    actual = mod.expand_stream_ids(
        {"ids": ["prot_1[_<?>]"], "files": [], "map_table": str(table)}
    )["ids"]
    record_case(input="lazy ids vs 3 map rows", expected=["prot_1_5N"], actual=actual)
    assert actual == ["prot_1_5N"]


# ── the zero-match report ────────────────────────────────────────────────────

@pytest.mark.parametrize("patterns, row_ids, expect_violation", [
    (["design_1"],  ["design_01", "design_02"], True),   # zero padding: exact match, zero hits
    (["design_1"],  ["design_1", "design_2"],   False),
    ([],            ["design_1"],               False),  # nothing asked for
    (["design_1"],  [],                         False),  # nothing to select from
])
def test_check_pattern_selection(record_case, patterns, row_ids, expect_violation):
    selected = [r for r in row_ids if r in patterns]
    violation = contract_enforcement.check_pattern_selection(patterns, row_ids, selected)
    actual = violation is not None
    record_case(input=(patterns, row_ids), expected=expect_violation, actual=actual)
    assert actual is expect_violation
    if violation is not None:
        assert violation.check == "pattern_selection"
        assert contract_enforcement.severity_for("pattern_selection") == contract_enforcement.WARN


def test_select_ids_reports_a_zero_match_selection(capsys, record_case):
    contract_enforcement.reset_reported()
    actual = idp.select_ids(["design_1"], ["design_01", "design_02"], where="unit test")
    err = capsys.readouterr().err
    record_case(input="design_1 vs design_01/design_02", expected="[] + warning",
                actual=f"{actual} + {'warning' if 'pattern_selection' in err else 'silence'}")
    assert actual == []
    assert "pattern_selection" in err
    assert "unit test" in err


def test_select_ids_stays_quiet_when_it_matches(capsys, record_case):
    contract_enforcement.reset_reported()
    actual = idp.select_ids(["design_<1 2>[_<?>]"], ["design_1_5", "design_2_7", "design_10_5"])
    err = capsys.readouterr().err
    expected = ["design_1_5", "design_2_7"]
    record_case(input="design_<1 2>[_<?>]", expected=expected, actual=actual)
    assert actual == expected
    assert "pattern_selection" not in err
