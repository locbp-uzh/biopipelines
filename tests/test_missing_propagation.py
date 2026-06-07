"""Unit tests for the missing-manifest propagation cascade.

Covers the framework-level pieces that excuse upstream-filtered ids downstream:

  * ``ligand_utils._read_codes`` reads codes from the rows the map_table
    actually carries, tolerating a compounds stream that statically declares
    more ids than survived an upstream filter (the XTB crash).
  * ``BaseConfig._collect_upstream_missing_paths`` gathers every input axis's
    ``missing`` table (not just the first), so a filter on more than one axis
    propagates fully.
  * ``pipe_propagate_missing`` merges multiple upstream manifests, de-duping
    by id.
"""

import json
import os
import subprocess
import sys

import pandas as pd
import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


# ── ligand_utils._read_codes tolerates an over-declared compounds stream ──────

def _write_compounds_stream(tmp_path, declared_ids, map_rows):
    map_csv = tmp_path / "compounds_map.csv"
    pd.DataFrame(map_rows).to_csv(map_csv, index=False)
    js = tmp_path / "compounds.json"
    json.dump(
        {
            "name": "compounds",
            "ids": list(declared_ids),
            "files": [],
            "map_table": str(map_csv),
            "format": "csv",
        },
        open(js, "w"),
    )
    return str(js)


def test_resolve_ligand_code_ignores_filtered_ids(tmp_path, record_case):
    """A compounds stream may declare ids an upstream filter dropped; the code
    resolver reads only the rows present in the map_table (no KeyError)."""
    from biopipelines.ligand_utils import resolve_ligand_code

    js = _write_compounds_stream(
        tmp_path,
        declared_ids=["dasatinib", "imatinib", "nilotinib", "bosutinib"],
        map_rows=[
            {"id": "dasatinib", "code": "LIG", "smiles": "CC"},
            {"id": "imatinib", "code": "LIG", "smiles": "CN"},
        ],
    )
    code = resolve_ligand_code(js)
    record_case(input="resolve_ligand_code(over-declared stream)",
                expected="LIG", actual=code)
    assert code == "LIG"


def test_resolve_ligand_codes_multiple(tmp_path):
    from biopipelines.ligand_utils import resolve_ligand_codes

    js = _write_compounds_stream(
        tmp_path,
        declared_ids=["a", "b", "c"],
        map_rows=[
            {"id": "a", "code": "LIG", "smiles": ""},
            {"id": "b", "code": "STI", "smiles": ""},
        ],
    )
    assert resolve_ligand_codes(js) == ["LIG", "STI"]


def test_resolve_ligand_code_rejects_multiple(tmp_path):
    from biopipelines.ligand_utils import resolve_ligand_code

    js = _write_compounds_stream(
        tmp_path,
        declared_ids=["a", "b"],
        map_rows=[
            {"id": "a", "code": "LIG", "smiles": ""},
            {"id": "b", "code": "STI", "smiles": ""},
        ],
    )
    with pytest.raises(ValueError):
        resolve_ligand_code(js)


def test_read_codes_empty_raises(tmp_path):
    from biopipelines.ligand_utils import resolve_ligand_code

    js = _write_compounds_stream(
        tmp_path, declared_ids=["a"], map_rows=[{"id": "a", "code": "", "smiles": ""}]
    )
    with pytest.raises(ValueError):
        resolve_ligand_code(js)


@pytest.mark.parametrize("code", ["DROX", "A7ZK"])
def test_resolve_ligand_code_accepts_extended_ccd(tmp_path, code, record_case):
    """4-5 char (extended CCD) codes resolve at the consumer boundary, matching
    Ligand's 1-5 acceptance — so an mmCIF-extracted or code=... extended ligand
    flows into PLACER/PLIP/PocketGen/PoseBusters/ProLIF/XTB/LigandMPNN."""
    from biopipelines.ligand_utils import resolve_ligand_code

    js = _write_compounds_stream(
        tmp_path, declared_ids=["x"], map_rows=[{"id": "x", "code": code, "smiles": ""}]
    )
    out = resolve_ligand_code(js)
    record_case(input=f"resolve_ligand_code(code={code!r})",
                expected=code, actual=out)
    assert out == code


def test_resolve_ligand_code_rejects_over_five_chars(tmp_path):
    """A code longer than the extended-CCD limit is still rejected clearly."""
    from biopipelines.ligand_utils import resolve_ligand_code

    js = _write_compounds_stream(
        tmp_path, declared_ids=["x"], map_rows=[{"id": "x", "code": "TOOLONG", "smiles": ""}]
    )
    with pytest.raises(ValueError, match="1-5 alphanumeric"):
        resolve_ligand_code(js)


# ── BaseConfig collects every input axis's missing table ──────────────────────

class _FakeTableInfo:
    def __init__(self, path):
        self.path = path


class _FakeTables:
    def __init__(self, missing_path):
        self._tables = {"missing": type("M", (), {"info": _FakeTableInfo(missing_path)})()}


class _FakeSource:
    def __init__(self, missing_path):
        self.tables = _FakeTables(missing_path)


def _bare_baseconfig():
    """A concrete BaseConfig instance (the collect/path helpers don't touch
    pipeline state, so the abstract methods can be no-ops)."""
    from biopipelines.base_config import BaseConfig

    class _Bare(BaseConfig):
        TOOL_NAME = "Bare"

        def validate_params(self):
            pass

        def configure_inputs(self, pipeline_folders):
            pass

        def generate_script(self, script_path):
            return ""

        def get_output_files(self):
            return {}

    return object.__new__(_Bare)


def test_collect_upstream_missing_paths_dedupes_and_orders():
    src_a = _FakeSource("/a/tables/missing.csv")
    src_b = _FakeSource("/b/tables/missing.csv")
    src_dup = _FakeSource("/a/tables/missing.csv")
    plain_stream = object()  # no .tables -> contributes nothing

    paths = _bare_baseconfig()._collect_upstream_missing_paths(
        src_a, plain_stream, src_b, None, src_dup
    )
    assert paths == ["/a/tables/missing.csv", "/b/tables/missing.csv"]


def test_get_upstream_missing_table_path_returns_first():
    src_a = _FakeSource("/a/tables/missing.csv")
    src_b = _FakeSource("/b/tables/missing.csv")
    self_ = _bare_baseconfig()
    assert (self_._get_upstream_missing_table_path(None, src_a, src_b)
            == "/a/tables/missing.csv")
    assert self_._get_upstream_missing_table_path(None, object()) is None


def test_upstream_missing_flag_emits_all_axes():
    """The pipe-script flag carries EVERY axis's path, not just the first —
    this is what fixes "multiple missing in, only one propagated"."""
    self_ = _bare_baseconfig()
    src_a = _FakeSource("/a/tables/missing.csv")
    src_b = _FakeSource("/b/tables/missing.csv")
    flag = self_.upstream_missing_flag(src_a, object(), src_b, None)
    assert flag == ' --upstream-missing "/a/tables/missing.csv" "/b/tables/missing.csv"'
    assert self_.upstream_missing_flag(object(), None) == ""
    # custom flag name (mmseqs-style underscore variant)
    assert self_.upstream_missing_flag(src_a, flag="--um").startswith(' --um ')


# ── biopipelines_io.read_upstream_missing merges + de-dupes ───────────────────

def test_read_upstream_missing_merges_dedupes(tmp_path):
    from biopipelines.biopipelines_io import read_upstream_missing

    a = tmp_path / "a.csv"
    b = tmp_path / "b.csv"
    pd.DataFrame([
        {"id": "x", "removed_by": "P", "kind": "filter", "cause": "1"},
        {"id": "y", "removed_by": "P", "kind": "filter", "cause": "2"},
    ]).to_csv(a, index=False)
    pd.DataFrame([
        {"id": "y", "removed_by": "Q", "kind": "filter", "cause": "dup"},
        {"id": "z", "removed_by": "Q", "kind": "filter", "cause": "3"},
    ]).to_csv(b, index=False)

    assert read_upstream_missing(None) == []
    assert read_upstream_missing([]) == []
    assert read_upstream_missing("/no/such.csv") == []
    assert [r["id"] for r in read_upstream_missing(str(a))] == ["x", "y"]
    # str and list both accepted; list merges + de-dupes by id (first position kept)
    merged = read_upstream_missing([str(a), str(b)])
    assert [r["id"] for r in merged] == ["x", "y", "z"]


# ── pipe_propagate_missing merges multiple manifests ──────────────────────────

def test_propagate_missing_merges_and_dedupes(tmp_path):
    up1 = tmp_path / "up1" / "tables"
    up2 = tmp_path / "up2" / "tables"
    up1.mkdir(parents=True)
    up2.mkdir(parents=True)
    pd.DataFrame([
        {"id": "nilotinib", "removed_by": "Panda", "kind": "filter", "cause": "pLDDT"},
        {"id": "bosutinib", "removed_by": "Panda", "kind": "filter", "cause": "pLDDT"},
    ]).to_csv(up1 / "missing.csv", index=False)
    pd.DataFrame([
        {"id": "bosutinib", "removed_by": "Panda", "kind": "filter", "cause": "dup"},
        {"id": "greasy_decoy", "removed_by": "Panda", "kind": "filter", "cause": "x"},
    ]).to_csv(up2 / "missing.csv", index=False)

    dest = tmp_path / "out" / "tables" / "missing.csv"
    r = subprocess.run(
        [sys.executable, os.path.join(REPO_ROOT, "pipe_scripts", "pipe_propagate_missing.py"),
         "--upstream-folders", str(up1), str(up2),
         "--output-folder", str(tmp_path / "out"),
         "--missing-csv", str(dest)],
        capture_output=True, text=True,
    )
    assert r.returncode == 0, r.stderr
    df = pd.read_csv(dest)
    assert sorted(df["id"]) == ["bosutinib", "greasy_decoy", "nilotinib"]


def test_propagate_missing_no_upstream_writes_empty(tmp_path):
    dest = tmp_path / "out" / "tables" / "missing.csv"
    r = subprocess.run(
        [sys.executable, os.path.join(REPO_ROOT, "pipe_scripts", "pipe_propagate_missing.py"),
         "--upstream-folders", str(tmp_path / "nonexistent"),
         "--output-folder", str(tmp_path / "out"),
         "--missing-csv", str(dest)],
        capture_output=True, text=True,
    )
    assert r.returncode == 0, r.stderr
    df = pd.read_csv(dest)
    assert list(df.columns) == ["id", "removed_by", "kind", "cause"]
    assert df.empty


# ── completion-check remaps input-axis ids into output id space ───────────────

def _load_pcc():
    import importlib.util
    path = os.path.join(REPO_ROOT, "pipe_scripts", "pipe_check_completion.py")
    spec = importlib.util.spec_from_file_location("pcc", path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def test_completion_remap_product_ids():
    """A raw ligand-axis missing id excuses the protein x ligand product files."""
    pcc = _load_pcc()
    expected = {"structures": {"ids": ["prot+lig1", "prot+lig2"], "files": ["<id>.cif"]},
                "tables": {}}
    assert "prot+lig2" in pcc._remap_missing_to_output_ids(["lig2"], expected)


def test_completion_remap_diffusion_fanout():
    pcc = _load_pcc()
    expected = {"structures": {"ids": ["prot+lig1_<1..2>", "prot+lig2_<1..2>"],
                               "files": ["<id>.cif"]}, "tables": {}}
    out = pcc._remap_missing_to_output_ids(["lig2"], expected)
    assert set(out) >= {"prot+lig2_1", "prot+lig2_2"}


def test_completion_remap_multiplier_ids():
    pcc = _load_pcc()
    expected = {"sequences": {"ids": ["p1_1", "p1_2", "p2_1", "p2_2"], "files": ["<id>.fa"]},
                "tables": {}}
    assert set(pcc._remap_missing_to_output_ids(["p2"], expected)) >= {"p2_1", "p2_2"}


def test_completion_remap_one_to_one_kept():
    """Exact-id (1:1 / single-axis) tools keep the id unchanged."""
    pcc = _load_pcc()
    expected = {"structures": {"ids": ["a", "b", "c"], "files": ["<id>.pdb"]}, "tables": {}}
    assert pcc._remap_missing_to_output_ids(["c"], expected) == ["c"]


def test_completion_remap_unmatched_axis_kept_harmless():
    """A missing id on an axis that is not an output dimension is kept as-is
    (excuses nothing, fails nothing)."""
    pcc = _load_pcc()
    expected = {"distances": {"ids": ["s1", "s2"], "files": ["<id>.csv"]}, "tables": {}}
    assert pcc._remap_missing_to_output_ids(["ligX"], expected) == ["ligX"]


# ── step id derivation + same-tool-chain excusal ──────────────────────────────

def test_step_id_from_table_path():
    from biopipelines.biopipelines_io import step_id_from_table_path

    assert step_id_from_table_path("/job/005_XTB/tables/missing.csv") == "005_XTB"
    assert step_id_from_table_path("") == ""


def _write_missing_table(tmp_path, step_folder, rows):
    """Build an expected_outputs dict pointing at a written missing.csv under
    <step_folder>/tables/missing.csv."""
    tdir = tmp_path / step_folder / "tables"
    tdir.mkdir(parents=True)
    path = tdir / "missing.csv"
    pd.DataFrame(rows, columns=["id", "removed_by", "kind", "cause"]).to_csv(path, index=False)
    return {"tables": {"missing": {"path": str(path)}}}, str(tmp_path / step_folder)


def test_same_tool_chain_excuses_upstream_failure(tmp_path):
    """Two steps of the same tool: an upstream-step failure row (002_XTB) is
    excused at the later step (005_XTB) instead of being read as a local
    failure of the bare tool name."""
    pcc = _load_pcc()
    expected, folder = _write_missing_table(
        tmp_path, "005_XTB",
        [{"id": "x", "removed_by": "002_XTB", "kind": "failure", "cause": "upstream"}],
    )
    step_id = pcc._step_id(folder, "XTB")
    excused = pcc._load_expected_missing_ids(expected, step_id, "XTB")
    assert excused == ["x"]


def test_local_failure_not_excused(tmp_path):
    """A failure row stamped with THIS step's id is a real local failure and
    must NOT be excused (would otherwise mask a genuine FAILED)."""
    pcc = _load_pcc()
    expected, folder = _write_missing_table(
        tmp_path, "005_XTB",
        [{"id": "x", "removed_by": "005_XTB", "kind": "failure", "cause": "local"}],
    )
    step_id = pcc._step_id(folder, "XTB")
    assert pcc._load_expected_missing_ids(expected, step_id, "XTB") == []


def test_bare_name_local_failure_not_excused(tmp_path):
    """Back-compat: a manifest written before step-prefixing (removed_by is the
    bare tool name) still counts as local and is not excused."""
    pcc = _load_pcc()
    expected, folder = _write_missing_table(
        tmp_path, "005_XTB",
        [{"id": "x", "removed_by": "XTB", "kind": "failure", "cause": "local"}],
    )
    step_id = pcc._step_id(folder, "XTB")
    assert pcc._load_expected_missing_ids(expected, step_id, "XTB") == []


def test_filter_kind_always_excused(tmp_path):
    """A kind=filter row is excused even when removed_by is this very step."""
    pcc = _load_pcc()
    expected, folder = _write_missing_table(
        tmp_path, "005_Panda",
        [{"id": "x", "removed_by": "005_Panda", "kind": "filter", "cause": "pLDDT"}],
    )
    step_id = pcc._step_id(folder, "Panda")
    assert pcc._load_expected_missing_ids(expected, step_id, "Panda") == ["x"]


def test_legacy_no_kind_branch_removed(tmp_path):
    """A manifest missing the kind column excuses nothing (the legacy
    excuse-everything branch is gone)."""
    pcc = _load_pcc()
    tdir = tmp_path / "005_XTB" / "tables"
    tdir.mkdir(parents=True)
    path = tdir / "missing.csv"
    pd.DataFrame([{"id": "x", "removed_by": "XTB", "cause": "c"}]).to_csv(path, index=False)
    expected = {"tables": {"missing": {"path": str(path)}}}
    step_id = pcc._step_id(str(tmp_path / "005_XTB"), "XTB")
    assert pcc._load_expected_missing_ids(expected, step_id, "XTB") == []


# ── _filter_expected_missing matches on the path's owner id ───────────────────

def test_filter_expected_missing_owner_id_not_basename():
    """Excusal matches the path's OWNER id, not its filename stem. A template
    that puts <id> anywhere (rank1_<id>.pdb, <id>/model.pdb) is excused by the
    same id, even though the basename is not the id."""
    pcc = _load_pcc()
    missing_pairs = [
        ("protA", "/o/rank1_protA.pdb"),   # template rank1_<id>.pdb
        ("protB", "/o/protB/model.pdb"),   # template <id>/model.pdb
        ("protC", "/o/rank1_protC.pdb"),   # not excused
    ]
    unexpected, expected = pcc._filter_expected_missing(missing_pairs, ["protA", "protB"])
    assert set(expected) == {"/o/rank1_protA.pdb", "/o/protB/model.pdb"}
    assert unexpected == ["/o/rank1_protC.pdb"]


def test_filter_expected_missing_descendant_lazy_fanout():
    """A lazy fan-out collapses to its prefix in the excused set; the concrete
    runtime child is excused because it descends from that prefix. A sibling
    design that shares only a top-level base is NOT excused."""
    pcc = _load_pcc()
    missing_pairs = [
        ("Panda_5_rank001", "/o/Panda_5_rank001.pdb"),   # child of excused Panda_5
        ("Panda_6_rank001", "/o/Panda_6_rank001.pdb"),   # different design
    ]
    unexpected, expected = pcc._filter_expected_missing(missing_pairs, ["Panda_5"])
    assert expected == ["/o/Panda_5_rank001.pdb"]
    assert unexpected == ["/o/Panda_6_rank001.pdb"]


def test_filter_expected_missing_no_owner_never_excused():
    """A path with no owner id (legacy/shared-file entry) is never excused."""
    pcc = _load_pcc()
    missing_pairs = [(None, "/o/shared.fasta")]
    unexpected, expected = pcc._filter_expected_missing(missing_pairs, ["shared"])
    assert unexpected == ["/o/shared.fasta"]
    assert expected == []


def test_extract_id_file_pairs_carries_owner_id():
    """extract_id_file_pairs threads the id through templates with <id> in a
    non-stem position, and yields owner_id None for bare/shared entries."""
    pcc = _load_pcc()
    # <id> as a non-stem prefix
    pairs = pcc.extract_id_file_pairs(
        {"ids": ["a", "b"], "files": ["rank1_<id>.pdb"]})
    assert pairs == [("a", "rank1_a.pdb"), ("b", "rank1_b.pdb")]
    # shared-file form -> no owner id
    assert pcc.extract_id_file_pairs({"ids": ["a"], "files": "shared.fa"}) == [(None, "shared.fa")]
    # legacy bare list -> no owner id
    assert pcc.extract_id_file_pairs(["x.pdb"]) == [(None, "x.pdb")]


# ── EnsembleAnalysis requires groups ──────────────────────────────────────────

def test_ensemble_analysis_requires_groups():
    from biopipelines.datastream import DataStream
    from biopipelines.ensemble_analysis import EnsembleAnalysis

    structures = DataStream(name="structures", ids=["c_1", "c_2"],
                            files=["<id>.pdb"], map_table="m.csv", format="pdb")
    tool = EnsembleAnalysis.__new__(EnsembleAnalysis)
    tool.structures_stream = structures
    tool.groups_stream = DataStream(name="sequences", ids=[], files=[],
                                    map_table="g.csv", format="csv")
    tool.selection = "CA"
    tool.reference = "mean"
    with pytest.raises(ValueError, match="groups is required"):
        tool.validate_params()
