"""Regression tests for the id-filter / lazy-id refactor.

A filter that restricts a stream's ids while leaving its map_table intact must
not let a bulk consumer reprocess the filtered-out rows. Streams carry id
*patterns*; runtime consumers select the map_table rows the patterns cover
(id_patterns.select_ids), so the rows are the source of truth and filtering is
honored without a config-time flag. Covers the generic materializer, both
combinatorics consumers (Boltz + AlphaFold query builders), and stitch.
"""

import json
import os

import pandas as pd
import pytest

from biopipelines.datastream import DataStream
from biopipelines.biopipelines_io import write_filtered_map_table


def _seq_stream(tmp_path, rows):
    table = tmp_path / "seqs.csv"
    pd.DataFrame(rows).to_csv(table, index=False)
    ds = DataStream(
        name="sequences",
        ids=[r["id"] for r in rows],
        files=[],
        map_table=str(table),
        format="csv",
    )
    return ds, str(table)


# ── write_filtered_map_table ───────────────────────────────────────────────

def test_filtered_subset_preserves_order(tmp_path, record_case):
    ds, _ = _seq_stream(tmp_path, [
        {"id": "a", "sequence": "AAA"},
        {"id": "b", "sequence": "BBB"},
        {"id": "c", "sequence": "CCC"},
    ])
    filt = ds.filter_by_ids(["c", "a"])
    out = tmp_path / "out.csv"
    write_filtered_map_table(filt, str(out))
    actual = pd.read_csv(out)["id"].tolist()
    # filter_by_ids preserves original expanded order (a before c)
    expected = ["a", "c"]
    record_case(input="filter [c,a]", expected=expected, actual=actual)
    assert actual == expected


def test_filtered_numeric_like_ids_match_as_strings(tmp_path, record_case):
    table = tmp_path / "seqs.csv"
    table.write_text("id,sequence\n001,AAA\n002,BBB\n")
    ds = DataStream(
        name="sequences",
        ids=["001", "002"],
        files=[],
        map_table=str(table),
        format="csv",
    )
    filt = ds.filter_by_ids(["002"])
    out = tmp_path / "out.csv"

    write_filtered_map_table(filt, str(out))

    actual = pd.read_csv(out, dtype={"id": str})["id"].tolist()
    expected = ["002"]
    record_case(input="numeric-looking id filter", expected=expected, actual=actual)
    assert actual == expected


def test_unfiltered_fast_path_byte_identical(tmp_path, record_case):
    ds, table = _seq_stream(tmp_path, [
        {"id": "a", "sequence": "AAA"},
        {"id": "b", "sequence": "BBB"},
    ])
    out = tmp_path / "out.csv"
    write_filtered_map_table(ds, str(out))
    actual = open(out, "rb").read()
    expected = open(table, "rb").read()
    record_case(input="unfiltered copy", expected="byte-identical", actual="byte-identical" if actual == expected else "differs")
    assert actual == expected


def test_column_subset(tmp_path, record_case):
    ds, _ = _seq_stream(tmp_path, [
        {"id": "a", "sequence": "AAA", "extra": "x"},
        {"id": "b", "sequence": "BBB", "extra": "y"},
    ])
    out = tmp_path / "out.csv"
    write_filtered_map_table(ds, str(out), columns=["sequence"])
    cols = list(pd.read_csv(out).columns)
    record_case(input="columns=[sequence]", expected=["id", "sequence"], actual=cols)
    assert cols == ["id", "sequence"]


def test_missing_id_raises(tmp_path, record_case):
    ds, _ = _seq_stream(tmp_path, [{"id": "a", "sequence": "AAA"}])
    ghost = DataStream(name="sequences", ids=["a", "ghost"], files=[],
                       map_table=ds.map_table, format="csv")
    record_case(input="id not in table", expected="KeyError", actual="KeyError")
    with pytest.raises(KeyError):
        write_filtered_map_table(ghost, str(tmp_path / "out.csv"))


def test_required_column_raises(tmp_path, record_case):
    ds, _ = _seq_stream(tmp_path, [{"id": "a", "sequence": "AAA"}])
    record_case(input="require missing col", expected="KeyError", actual="KeyError")
    with pytest.raises(KeyError):
        write_filtered_map_table(ds, str(tmp_path / "out.csv"), required_columns=["msa_file"])


# ── combinatorics carries filtered ids ──────────────────────────────────────

def test_combinatorics_source_carries_filtered_ids(tmp_path, record_case):
    from biopipelines.combinatorics import (
        generate_combinatorics_config,
        load_ids_from_sources,
    )

    ds, table = _seq_stream(tmp_path, [
        {"id": "a", "sequence": "AAA"},
        {"id": "b", "sequence": "BBB"},
        {"id": "c", "sequence": "CCC"},
    ])
    filt = ds.filter_by_ids(["a", "c"])
    cfg = tmp_path / "combo.json"
    generate_combinatorics_config(str(cfg), proteins=(filt, "sequences"))

    data = json.loads(cfg.read_text())
    src = data["axes"]["proteins"]["sources"][0]
    # Raw table path preserved; the stream's id patterns are carried and the
    # runtime selects matching rows (materialization is runtime, not config-time).
    same_path = os.path.abspath(src["path"]) == os.path.abspath(table)
    ids = load_ids_from_sources([src])
    record_case(input="filtered combinatorics", expected=("raw path", ["a", "c"]),
                actual=("raw path" if same_path else "rewritten", ids))
    assert same_path
    assert src["ids"] == ["a", "c"]
    assert ids == ["a", "c"]


def test_combinatorics_bare_filtered_datastream_source(tmp_path, record_case):
    """A bare filtered DataStream passed directly as a combinatorics axis (e.g.
    proteins=output.streams.sequences.filter_by_ids(...)) must yield its ids at
    config time — it has no .streams attribute, so the id collector must read
    .ids directly instead of raising 'No proteins found'."""
    from biopipelines.combinatorics import (
        generate_combinatorics_config,
        load_ids_from_sources,
    )

    ds, _ = _seq_stream(tmp_path, [
        {"id": "p_1", "sequence": "AAA"},
        {"id": "p_2", "sequence": "BBB"},
        {"id": "p_3", "sequence": "CCC"},
    ])
    kept = ds.filter_by_ids(["p_1", "p_3"])  # bare DataStream, not StandardizedOutput
    cfg = tmp_path / "combo.json"
    generate_combinatorics_config(str(cfg), proteins=(kept, "sequences"))  # must not raise
    src = json.loads(cfg.read_text())["axes"]["proteins"]["sources"][0]
    ids = load_ids_from_sources([src])
    record_case(input="bare filtered DS axis", expected=["p_1", "p_3"], actual=ids)
    assert src["ids"] == ["p_1", "p_3"]
    assert ids == ["p_1", "p_3"]


def test_combinatorics_lazy_stream_carries_pattern(tmp_path, record_case):
    """A lazy stream carries its id pattern; the runtime selects matching rows
    against the real map_table, honoring any upstream filter."""
    from biopipelines.combinatorics import (
        generate_combinatorics_config,
        load_ids_from_sources,
    )

    table = tmp_path / "seqs.csv"
    pd.DataFrame([
        {"id": "prot_1A", "sequence": "AAA"},
        {"id": "prot_2V", "sequence": "VVV"},
        {"id": "other_1", "sequence": "CCC"},
    ]).to_csv(table, index=False)
    lazy = DataStream(
        name="sequences",
        ids=["prot[_<N><A V>]"],
        files=[],
        map_table=str(table),
        format="csv",
    )
    cfg = tmp_path / "combo.json"
    generate_combinatorics_config(str(cfg), proteins=(lazy, "sequences"))
    src = json.loads(cfg.read_text())["axes"]["proteins"]["sources"][0]
    ids = load_ids_from_sources([src])
    record_case(input="lazy stream pattern select", expected=["prot_1A", "prot_2V"], actual=ids)
    assert src["ids"] == ["prot[_<N><A V>]"]
    assert ids == ["prot_1A", "prot_2V"]  # 'other_1' excluded by the pattern


def _load_boltz_helper():
    import importlib.util

    repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    spec = importlib.util.spec_from_file_location(
        "pipe_boltz_config_unified",
        os.path.join(repo, "pipe_scripts", "pipe_boltz_config_unified.py"),
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_boltz_unified_loader_selects_by_pattern(tmp_path, record_case):
    """Patterns select; an id the patterns do not cover is simply absent (no raise)."""
    mod = _load_boltz_helper()
    table = tmp_path / "compounds.csv"
    pd.DataFrame([
        {"id": "lig1", "smiles": "CCO"},
        {"id": "lig2", "smiles": "CCN"},
    ]).to_csv(table, index=False)
    # 'ghost' matches no row → contributes nothing, no error.
    axis = {"sources": [{"path": str(table), "iterate": True, "ids": ["lig1", "ghost"]}]}
    iterated, _static, _ = mod.load_axis_data(axis)
    ids = [r["id"] for r in iterated]
    record_case(input="boltz pattern select", expected=["lig1"], actual=ids)
    assert ids == ["lig1"]


def _load_stitch_helper():
    import importlib.util

    repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    spec = importlib.util.spec_from_file_location(
        "pipe_stitch_sequences",
        os.path.join(repo, "pipe_scripts", "pipe_stitch_sequences.py"),
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_stitch_selects_by_pattern(tmp_path, record_case):
    mod = _load_stitch_helper()
    csv = tmp_path / "seqs.csv"
    pd.DataFrame([
        {"id": "a", "sequence": "AAA"},
        {"id": "b", "sequence": "BBB"},
        {"id": "c", "sequence": "CCC"},
    ]).to_csv(csv, index=False)

    info = {"type": "tool_output", "sequences_file": str(csv), "id_patterns": ["a", "c"]}
    kept = mod.load_sequences_from_info(info)
    record_case(input="stitch concrete select", expected=["a", "c"], actual=list(kept))
    assert list(kept) == ["a", "c"]


def test_stitch_lazy_pattern_selects_matching(tmp_path, record_case):
    mod = _load_stitch_helper()
    csv = tmp_path / "seqs.csv"
    pd.DataFrame([
        {"id": "prot_1A", "sequence": "AAA"},
        {"id": "prot_2V", "sequence": "VVV"},
        {"id": "other_1", "sequence": "CCC"},
    ]).to_csv(csv, index=False)

    info = {"type": "tool_output", "sequences_file": str(csv), "id_patterns": ["prot[_<N><A V>]"]}
    kept = mod.load_sequences_from_info(info)
    record_case(input="stitch lazy pattern", expected=["prot_1A", "prot_2V"], actual=list(kept))
    assert list(kept) == ["prot_1A", "prot_2V"]


def test_combinatorics_unfiltered_no_ids_key_filter(tmp_path, record_case):
    from biopipelines.combinatorics import (
        generate_combinatorics_config,
        load_ids_from_sources,
    )

    ds, _ = _seq_stream(tmp_path, [
        {"id": "a", "sequence": "AAA"},
        {"id": "b", "sequence": "BBB"},
    ])
    cfg = tmp_path / "combo.json"
    generate_combinatorics_config(str(cfg), proteins=(ds, "sequences"))
    src = json.loads(cfg.read_text())["axes"]["proteins"]["sources"][0]
    ids = load_ids_from_sources([src])
    record_case(input="unfiltered combinatorics", expected=["a", "b"], actual=ids)
    assert ids == ["a", "b"]


# ── runtime materializer pipe script ────────────────────────────────────────

def test_materialize_pipe_script_runtime(tmp_path, record_case):
    import subprocess
    import sys

    ds, _ = _seq_stream(tmp_path, [
        {"id": "a", "sequence": "AAA"},
        {"id": "b", "sequence": "BBB"},
        {"id": "c", "sequence": "CCC"},
    ])
    filt = ds.filter_by_ids(["a", "c"])
    js = tmp_path / "seq.json"
    filt.save_json(str(js))
    out = tmp_path / "filtered.csv"

    repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    script = os.path.join(repo, "pipe_scripts", "materialize_filtered_map_table.py")
    r = subprocess.run(
        [sys.executable, script, str(js), str(out), "--require", "id,sequence"],
        capture_output=True, text=True,
    )
    actual = pd.read_csv(out)["id"].tolist() if r.returncode == 0 else r.stderr
    record_case(input="runtime materialize filtered", expected=["a", "c"], actual=actual)
    assert r.returncode == 0, r.stderr
    assert actual == ["a", "c"]


# ── multi-stream consumer honors ids (boltz unified config builder) ─────────

def test_boltz_unified_loader_honors_ids(tmp_path, record_case):
    import importlib.util

    repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    spec = importlib.util.spec_from_file_location(
        "pipe_boltz_config_unified",
        os.path.join(repo, "pipe_scripts", "pipe_boltz_config_unified.py"),
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)

    table = tmp_path / "compounds.csv"
    pd.DataFrame([
        {"id": "lig1", "smiles": "CCO"},
        {"id": "lig2", "smiles": "CCN"},
        {"id": "lig3", "smiles": "CCC"},
    ]).to_csv(table, index=False)

    axis = {"sources": [{"path": str(table), "iterate": True, "ids": ["lig1", "lig3"]}]}
    iterated, static, _ = mod.load_axis_data(axis)
    ids = [r["id"] for r in iterated]
    record_case(input="boltz unified filtered ids", expected=["lig1", "lig3"], actual=ids)
    assert ids == ["lig1", "lig3"]


def test_boltz_unified_loader_preserves_numeric_like_ids(tmp_path, record_case):
    mod = _load_boltz_helper()
    table = tmp_path / "compounds.csv"
    table.write_text("id,smiles\n001,CCO\n002,CCN\n")

    axis = {"sources": [{"path": str(table), "iterate": True, "ids": ["002"]}]}
    iterated, _static, _ = mod.load_axis_data(axis)

    ids = [r["id"] for r in iterated]
    record_case(input="boltz numeric-looking ids", expected=["002"], actual=ids)
    assert ids == ["002"]


# ── AlphaFold query builder is a combinatorics consumer too ──────────────────

def _load_af_query_helper():
    import importlib.util

    repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    spec = importlib.util.spec_from_file_location(
        "pipe_alphafold_queries",
        os.path.join(repo, "pipe_scripts", "pipe_alphafold_queries.py"),
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _af_seq_csv(tmp_path):
    p = tmp_path / "af_seqs.csv"
    pd.DataFrame([
        {"id": "prot_1A", "sequence": "MAAA"},
        {"id": "prot_2V", "sequence": "MVVV"},
        {"id": "other_1", "sequence": "MCCC"},
    ]).to_csv(p, index=False)
    return str(p)


def test_af_queries_lazy_pattern_selects(tmp_path, record_case):
    """A lazy source pattern selects matching rows instead of raising KeyError."""
    mod = _load_af_query_helper()
    got = mod._load_id_to_seq(_af_seq_csv(tmp_path), patterns=["prot[_<N><A V>]"])
    record_case(input="af lazy pattern", expected=["prot_1A", "prot_2V"], actual=list(got))
    assert list(got) == ["prot_1A", "prot_2V"]  # 'other_1' excluded


def test_af_queries_concrete_filter_and_nonmatch(tmp_path, record_case):
    mod = _load_af_query_helper()
    csv = _af_seq_csv(tmp_path)
    concrete = mod._load_id_to_seq(csv, patterns=["prot_1A"])
    nonmatch = mod._load_id_to_seq(csv, patterns=["ghost_<0..2>"])
    unfiltered = mod._load_id_to_seq(csv)
    record_case(input="af concrete/nonmatch/all",
                expected=(["prot_1A"], [], 3),
                actual=(list(concrete), list(nonmatch), len(unfiltered)))
    assert list(concrete) == ["prot_1A"]
    assert nonmatch == {}                       # non-matching pattern: no rows, no raise
    assert len(unfiltered) == 3                  # no patterns: every row


def test_af_queries_preserves_numeric_like_ids(tmp_path, record_case):
    mod = _load_af_query_helper()
    csv = tmp_path / "af_numeric.csv"
    csv.write_text("id,sequence\n001,MAAA\n002,MVVV\n")

    got = mod._load_id_to_seq(str(csv), patterns=["002"])

    record_case(input="af numeric-looking ids", expected=["002"], actual=list(got))
    assert list(got) == ["002"]


# ── materializer reads concrete rows in runtime mode for a lazy stream ───────

def test_materializer_lazy_runtime_uses_map_rows(tmp_path, record_case):
    """A lazy stream in runtime mode resolves ids_expanded from the map_table, so
    write_filtered_map_table emits exactly the rows the upstream wrote."""
    table = tmp_path / "lazy_map.csv"
    pd.DataFrame([
        {"id": "prot_1A", "sequence": "MAAA"},
        {"id": "prot_2V", "sequence": "MVVV"},
    ]).to_csv(table, index=False)
    ds = DataStream(name="sequences", ids=["prot[_<N><A V>]"], files=[],
                    map_table=str(table), format="csv", _runtime_mode=True)
    out = tmp_path / "out.csv"
    write_filtered_map_table(ds, str(out))
    actual = pd.read_csv(out)["id"].tolist()
    record_case(input="lazy runtime materialize", expected=["prot_1A", "prot_2V"], actual=actual)
    assert actual == ["prot_1A", "prot_2V"]
