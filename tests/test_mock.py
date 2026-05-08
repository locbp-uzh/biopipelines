"""Tests for the Mock tool — stub outputs for pipeline-wiring tests.

Covers:
- explicit-IDs mode: file/value streams, tables with fill / explicit rows
- fan-out via `children` (deterministic and lazy + `produce`)
- combinatorics via Bundle/Each passed as source
- map_table_strategy variants (runtime / config / both)
- `missing` post-expansion filtering
- runtime script actually emits files, map_tables, and CSV tables
"""

import csv
import json
import os
import subprocess
import sys

import pytest


# ── config-time: explicit IDs ─────────────────────────────────────────────────

def test_mock_explicit_ids_file_stream(
    local_config, isolated_cwd, new_pipeline, assert_valid_script, record_case,
):
    from biopipelines.mock import Mock

    pipeline = new_pipeline("mock_ids_files")
    with pipeline:
        m = Mock(
            ids=["a", "b", "c"],
            streams={"structures": {"format": "pdb", "file": "<id>.pdb"}},
        )
        script_path = pipeline.save()

    ids = list(m.streams.structures.ids)
    record_case(input="Mock(ids=[a,b,c])",
                expected=["a", "b", "c"], actual=ids)
    assert_valid_script(script_path, "Mock")
    assert ids == ["a", "b", "c"]


def test_mock_value_stream_and_table_fill(
    local_config, isolated_cwd, new_pipeline, assert_valid_script, record_case,
):
    from biopipelines.mock import Mock

    pipeline = new_pipeline("mock_values_table")
    with pipeline:
        m = Mock(
            ids=["x", "y"],
            streams={"seqs": {"format": "csv", "values": ["MKT", "AET"]}},
            tables={"scores": {"columns": ["score"], "fill": {"score": 0.0}}},
        )
        script_path = pipeline.save()

    cols = list(m.tables.scores.info.columns)
    record_case(input="Mock value-stream + table fill",
                expected=(["x", "y"], ["id", "score"]),
                actual=(list(m.streams.seqs.ids), cols))
    assert_valid_script(script_path, "Mock")
    assert list(m.streams.seqs.ids) == ["x", "y"]
    assert list(m.streams.seqs.files) == []
    assert cols == ["id", "score"]


def test_mock_rejects_ids_and_source_together(record_case):
    from biopipelines.mock import Mock

    record_case(input="Mock(ids=..., source=...)",
                expected="ValueError", actual="ValueError")
    with pytest.raises(ValueError):
        Mock(ids=["a"], source="anything")


def test_mock_rejects_lazy_children_without_produce(record_case):
    from biopipelines.mock import Mock

    record_case(input="Mock(children='[_<N><A V>]') without produce",
                expected="ValueError", actual="ValueError")
    with pytest.raises(ValueError):
        Mock(ids=["a"], children="[_<N><A V>]")


# ── config-time: deterministic children ───────────────────────────────────────

def test_mock_deterministic_children_keeps_pattern(
    local_config, isolated_cwd, new_pipeline, record_case,
):
    """Deterministic children keep the output IDs compact (no expansion)."""
    from biopipelines.mock import Mock

    pipeline = new_pipeline("mock_children_det")
    with pipeline:
        m = Mock(
            ids=["prot_0", "prot_1"],
            children="<1..3>",
            streams={"designs": {"format": "pdb", "file": "<id>.pdb"}},
        )
        pipeline.save()

    ids = list(m.streams.designs.ids)
    record_case(input="ids=[prot_0,prot_1], children=<1..3>",
                expected=["prot_0_<1..3>", "prot_1_<1..3>"],
                actual=ids)
    assert ids == ["prot_0_<1..3>", "prot_1_<1..3>"]


# ── config-time: lazy children ────────────────────────────────────────────────

def test_mock_lazy_children_keeps_lazy_pattern(
    local_config, isolated_cwd, new_pipeline, record_case,
):
    """Lazy children keep the output IDs lazy at config time."""
    from biopipelines.mock import Mock
    from biopipelines import id_patterns

    pipeline = new_pipeline("mock_children_lazy")
    with pipeline:
        m = Mock(
            ids=["prot_0", "prot_1"],
            children="[_<N><A V>]",
            produce=["_1A", "_1V"],
            streams={"mut": {"format": "pdb", "file": "<id>.pdb"}},
        )
        pipeline.save()

    ids = list(m.streams.mut.ids)
    record_case(input="lazy children=[_<N><A V>]",
                expected=("all lazy", True),
                actual=("all lazy", all(id_patterns.is_lazy(i) for i in ids)))
    assert all(id_patterns.is_lazy(i) for i in ids)


# ── runtime: deterministic fan-out materializes files ────────────────────────

def _run_pipe_mock(config_path):
    import biopipelines.mock as _m
    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(_m.__file__)))
    pipe_mock = os.path.join(repo_root, "pipe_scripts", "pipe_mock.py")
    subprocess.run([sys.executable, pipe_mock, str(config_path)], check=True)


def test_mock_runtime_deterministic_children_files(tmp_path, record_case):
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    cfg = {
        "output_folder": str(out_dir),
        "parent_ids": ["p0", "p1"],
        "output_ids": ["p0_<1..2>", "p1_<1..2>"],
        "provenance": {},
        "axis_names": [],
        "children": "<1..2>",
        "produce": None,
        "streams": {
            "designs": {"format": "pdb", "file": "<id>.pdb", "values": None,
                        "map_table": str(out_dir / "designs_map.csv")},
        },
        "tables": {},
        "map_table_strategy": "runtime",
        "missing": [],
        "source_streams": [],
    }
    cfg_path = tmp_path / "cfg.json"
    json.dump(cfg, open(cfg_path, "w"))
    _run_pipe_mock(cfg_path)

    produced = sorted(f for f in os.listdir(out_dir) if f.endswith(".pdb"))
    with open(out_dir / "designs_map.csv") as f:
        rows = list(csv.DictReader(f))
    ids = [r["id"] for r in rows]

    record_case(input="runtime children=<1..2>, parents=[p0,p1]",
                expected=(["p0_1.pdb", "p0_2.pdb", "p1_1.pdb", "p1_2.pdb"],
                          ["p0_1", "p0_2", "p1_1", "p1_2"]),
                actual=(produced, ids))
    assert produced == ["p0_1.pdb", "p0_2.pdb", "p1_1.pdb", "p1_2.pdb"]
    assert ids == ["p0_1", "p0_2", "p1_1", "p1_2"]


def test_mock_runtime_lazy_children_with_produce(tmp_path, record_case):
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    cfg = {
        "output_folder": str(out_dir),
        "parent_ids": ["prot_0", "prot_1"],
        "output_ids": ["prot_0[_<N><A V>]", "prot_1[_<N><A V>]"],
        "provenance": {},
        "axis_names": [],
        "children": "[_<N><A V>]",
        "produce": ["_1A", "_1V"],
        "streams": {
            "mut": {"format": "pdb", "file": "<id>.pdb", "values": None,
                    "map_table": str(out_dir / "mut_map.csv")},
        },
        "tables": {},
        "map_table_strategy": "runtime",
        "missing": [],
        "source_streams": [],
    }
    cfg_path = tmp_path / "cfg.json"
    json.dump(cfg, open(cfg_path, "w"))
    _run_pipe_mock(cfg_path)

    with open(out_dir / "mut_map.csv") as f:
        rows = list(csv.DictReader(f))
    ids = [r["id"] for r in rows]
    parents = [r["mut.parent"] for r in rows]

    record_case(input="lazy children + produce=[_1A,_1V]",
                expected=(["prot_0_1A", "prot_0_1V", "prot_1_1A", "prot_1_1V"],
                          ["prot_0", "prot_0", "prot_1", "prot_1"]),
                actual=(ids, parents))
    assert ids == ["prot_0_1A", "prot_0_1V", "prot_1_1A", "prot_1_1V"]
    assert parents == ["prot_0", "prot_0", "prot_1", "prot_1"]


def test_mock_runtime_missing_drops_ids(tmp_path, record_case):
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    cfg = {
        "output_folder": str(out_dir),
        "parent_ids": ["p0", "p1"],
        "output_ids": ["p0_<1..2>", "p1_<1..2>"],
        "provenance": {},
        "axis_names": [],
        "children": "<1..2>",
        "produce": None,
        "streams": {
            "designs": {"format": "pdb", "file": "<id>.pdb", "values": None,
                        "map_table": str(out_dir / "designs_map.csv")},
        },
        "tables": {},
        "map_table_strategy": "runtime",
        "missing": ["p0_1", "p1"],  # one expanded + one parent
        "source_streams": [],
    }
    cfg_path = tmp_path / "cfg.json"
    json.dump(cfg, open(cfg_path, "w"))
    _run_pipe_mock(cfg_path)

    with open(out_dir / "designs_map.csv") as f:
        rows = list(csv.DictReader(f))
    ids = [r["id"] for r in rows]
    record_case(input="missing=['p0_1','p1'] (drops 1 child + 1 parent)",
                expected=["p0_2"], actual=ids)
    assert ids == ["p0_2"]


# ── config-time map_table strategy ────────────────────────────────────────────

def test_mock_map_table_config_strategy(
    local_config, isolated_cwd, new_pipeline, record_case,
):
    """map_table_strategy='config' writes the map_table in get_output_files."""
    from biopipelines.mock import Mock

    pipeline = new_pipeline("mock_map_config")
    with pipeline:
        m = Mock(
            ids=["a", "b"],
            streams={"s": {"format": "pdb", "file": "<id>.pdb"}},
            map_table_strategy="config",
        )
        pipeline.save()

    map_path = m.streams.s.map_table
    exists = os.path.exists(map_path)
    record_case(input="map_table_strategy='config'",
                expected=("map_table exists at config time", True),
                actual=("map_table exists at config time", exists))
    assert exists


# ── combinatorics via Bundle/Each sources ─────────────────────────────────────

def test_mock_source_each_cartesian(
    local_config, isolated_cwd, new_pipeline, record_case,
):
    """Mock(source=[Each(a), Each(b)]) → cartesian parent IDs with provenance."""
    from biopipelines.datastream import DataStream, create_map_table
    from biopipelines.combinatorics import Each
    from biopipelines.mock import Mock

    a_map = isolated_cwd / "a_map.csv"
    b_map = isolated_cwd / "b_map.csv"
    create_map_table(str(a_map), ids=["a1", "a2"])
    create_map_table(str(b_map), ids=["b1", "b2"])

    a = DataStream(name="sequences", ids=["a1", "a2"], map_table=str(a_map))
    b = DataStream(name="compounds", ids=["b1", "b2"], map_table=str(b_map))

    pipeline = new_pipeline("mock_each_cartesian")
    with pipeline:
        m = Mock(
            source=[Each(a), Each(b)],
            streams={"out": {"format": "pdb", "file": "<id>.pdb"}},
        )
        pipeline.save()

    ids = list(m.streams.out.ids)
    record_case(input="source=[Each(a[2]), Each(b[2])]",
                expected=["a1+b1", "a1+b2", "a2+b1", "a2+b2"],
                actual=ids)
    assert ids == ["a1+b1", "a1+b2", "a2+b1", "a2+b2"]


# ── integration: Mock + downstream pipeline.save ─────────────────────────────

def test_mock_integration_script_emitted(
    local_config, isolated_cwd, new_pipeline, assert_valid_script, record_case,
):
    """Mock inside a Pipeline produces a runnable bash script."""
    from biopipelines.mock import Mock

    pipeline = new_pipeline("mock_integration")
    with pipeline:
        Mock(
            ids=["a", "b"],
            streams={"s": {"format": "pdb", "file": "<id>.pdb"}},
            tables={"t": {"columns": ["score"], "fill": 1.0}},
        )
        script_path = pipeline.save()

    record_case(input="Mock + pipeline.save()",
                expected="Mock in script",
                actual=os.path.basename(script_path))
    assert_valid_script(script_path, "Mock", "mock_integration")


# ── developer_manual.md coverage: additional patterns ────────────────────────

# -- combinatorics: Bundle / single-Each / bare DataStream --------------------

def test_mock_source_bundle_groups_ids(
    local_config, isolated_cwd, new_pipeline, record_case,
):
    """Bundle(a, b) keeps elements together as one parent (no cartesian)."""
    from biopipelines.datastream import DataStream, create_map_table
    from biopipelines.combinatorics import Bundle
    from biopipelines.mock import Mock

    a_map = isolated_cwd / "a_map.csv"
    b_map = isolated_cwd / "b_map.csv"
    create_map_table(str(a_map), ids=["a1", "a2"])
    create_map_table(str(b_map), ids=["b1", "b2"])
    a = DataStream(name="sequences", ids=["a1", "a2"], map_table=str(a_map))
    b = DataStream(name="sequences", ids=["b1", "b2"], map_table=str(b_map))

    pipeline = new_pipeline("mock_bundle")
    with pipeline:
        m = Mock(
            source=Bundle(a, b),
            streams={"out": {"format": "pdb", "file": "<id>.pdb"}},
        )
        pipeline.save()

    ids = list(m.streams.out.ids)
    record_case(input="source=Bundle(a[2], b[2])",
                expected=(1, True),
                actual=(len(ids), all("+" in i for i in ids)))
    # Bundle produces a single '+'-joined parent ID (not cartesian)
    assert len(ids) == 1
    assert "+" in ids[0]


def test_mock_source_single_each_axis_provenance(
    isolated_cwd, record_case,
):
    """Each(a) alone yields flat parent IDs and single-axis provenance.

    Uses a direct Mock._resolve_from_source call (no Pipeline) because inside a
    Pipeline context Mock(...) returns a StandardizedOutput which does not
    expose the raw axis_names/provenance attributes.
    """
    from biopipelines.datastream import DataStream, create_map_table
    from biopipelines.combinatorics import Each
    from biopipelines.mock import Mock

    a_map = isolated_cwd / "a_map.csv"
    create_map_table(str(a_map), ids=["a1", "a2", "a3"])
    a = DataStream(name="sequences", ids=["a1", "a2", "a3"], map_table=str(a_map))

    parents, prov, axes = Mock._resolve_from_source(Each(a))
    record_case(input="source=Each(a[3])",
                expected=(["a1", "a2", "a3"], ["sequences"],
                          {"sequences": ["a1", "a2", "a3"]}),
                actual=(parents, axes, prov))
    assert parents == ["a1", "a2", "a3"]
    assert axes == ["sequences"]
    assert prov == {"sequences": ["a1", "a2", "a3"]}


def test_mock_source_bare_datastream(
    isolated_cwd, record_case,
):
    """A bare DataStream (no Each/Bundle wrapper) is accepted as a single axis."""
    from biopipelines.datastream import DataStream, create_map_table
    from biopipelines.mock import Mock

    a_map = isolated_cwd / "a_map.csv"
    create_map_table(str(a_map), ids=["p1", "p2"])
    a = DataStream(name="structures", ids=["p1", "p2"], map_table=str(a_map))

    parents, prov, axes = Mock._resolve_from_source(a)
    record_case(input="source=<bare DataStream>",
                expected=(["p1", "p2"], ["structures"]),
                actual=(parents, axes))
    assert parents == ["p1", "p2"]
    assert axes == ["structures"]


# -- map_table_strategy = "both" ----------------------------------------------

def test_mock_map_table_both_strategy(
    local_config, isolated_cwd, new_pipeline, record_case,
):
    """'both' writes the map_table at config time AND the pipe script rewrites it."""
    from biopipelines.mock import Mock

    pipeline = new_pipeline("mock_map_both")
    with pipeline:
        m = Mock(
            ids=["a", "b"],
            streams={"s": {"format": "pdb", "file": "<id>.pdb"}},
            map_table_strategy="both",
        )
        pipeline.save()

    map_path = m.streams.s.map_table
    config_time_exists = os.path.exists(map_path)

    # Now run the runtime pipe script — it should overwrite/re-create the map
    cfg_path = os.path.join(m.output_folder, "_configuration", "mock_config.json")
    _run_pipe_mock(cfg_path)

    with open(map_path) as f:
        rows = list(csv.DictReader(f))

    record_case(input="map_table_strategy='both'",
                expected=(True, ["a", "b"]),
                actual=(config_time_exists, [r["id"] for r in rows]))
    assert config_time_exists
    assert [r["id"] for r in rows] == ["a", "b"]


# -- table variants: explicit rows, dict fill, multi-column fill --------------

def test_mock_table_explicit_rows(tmp_path, record_case):
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    cfg = {
        "output_folder": str(out_dir),
        "parent_ids": ["a", "b"],
        "output_ids": ["a", "b"],
        "provenance": {}, "axis_names": [],
        "children": None, "produce": None,
        "streams": {},
        "tables": {
            "scores": {
                "path": str(out_dir / "scores.csv"),
                "columns": ["id", "score"],
                "fill": None,
                "rows": [["a", 0.9], ["b", 0.7]],
            }
        },
        "map_table_strategy": "runtime", "missing": [], "source_streams": [],
    }
    cfg_path = tmp_path / "cfg.json"
    json.dump(cfg, open(cfg_path, "w"))
    _run_pipe_mock(cfg_path)

    with open(out_dir / "scores.csv") as f:
        rows = list(csv.reader(f))

    record_case(input="table rows=[[a,0.9],[b,0.7]]",
                expected=[["id", "score"], ["a", "0.9"], ["b", "0.7"]],
                actual=rows)
    assert rows == [["id", "score"], ["a", "0.9"], ["b", "0.7"]]


def test_mock_table_dict_fill_per_column(tmp_path, record_case):
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    cfg = {
        "output_folder": str(out_dir),
        "parent_ids": ["x", "y"],
        "output_ids": ["x", "y"],
        "provenance": {}, "axis_names": [],
        "children": None, "produce": None,
        "streams": {},
        "tables": {
            "scores": {
                "path": str(out_dir / "scores.csv"),
                "columns": ["id", "plddt", "rmsd"],
                "fill": {"plddt": 0.85, "rmsd": 1.2},
                "rows": None,
            }
        },
        "map_table_strategy": "runtime", "missing": [], "source_streams": [],
    }
    cfg_path = tmp_path / "cfg.json"
    json.dump(cfg, open(cfg_path, "w"))
    _run_pipe_mock(cfg_path)

    with open(out_dir / "scores.csv") as f:
        rows = list(csv.DictReader(f))
    record_case(input="fill={plddt:0.85, rmsd:1.2}",
                expected=[("x", "0.85", "1.2"), ("y", "0.85", "1.2")],
                actual=[(r["id"], r["plddt"], r["rmsd"]) for r in rows])
    assert [(r["id"], r["plddt"], r["rmsd"]) for r in rows] == [
        ("x", "0.85", "1.2"), ("y", "0.85", "1.2")
    ]


def test_mock_table_scalar_fill_applied_to_all_columns(tmp_path, record_case):
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    cfg = {
        "output_folder": str(out_dir),
        "parent_ids": ["i"],
        "output_ids": ["i"],
        "provenance": {}, "axis_names": [],
        "children": None, "produce": None,
        "streams": {},
        "tables": {
            "t": {
                "path": str(out_dir / "t.csv"),
                "columns": ["id", "a", "b"],
                "fill": 42,
                "rows": None,
            }
        },
        "map_table_strategy": "runtime", "missing": [], "source_streams": [],
    }
    cfg_path = tmp_path / "cfg.json"
    json.dump(cfg, open(cfg_path, "w"))
    _run_pipe_mock(cfg_path)

    with open(out_dir / "t.csv") as f:
        rows = list(csv.reader(f))
    record_case(input="scalar fill=42",
                expected=[["id", "a", "b"], ["i", "42", "42"]],
                actual=rows)
    assert rows == [["id", "a", "b"], ["i", "42", "42"]]


# -- compact ID pattern input (config time) -----------------------------------

def test_mock_compact_ids_preserved_at_config_time(
    local_config, isolated_cwd, new_pipeline, record_case,
):
    """Compact patterns in `ids` must be preserved, not expanded at config time
    (matches the developer manual rule: don't call ids_expanded at config time)."""
    from biopipelines.mock import Mock

    pipeline = new_pipeline("mock_compact_ids")
    with pipeline:
        m = Mock(
            ids=["prot_<0..2>"],
            streams={"s": {"format": "pdb", "file": "<id>.pdb"}},
        )
        pipeline.save()

    ids = list(m.streams.s.ids)
    record_case(input="ids=['prot_<0..2>']",
                expected=["prot_<0..2>"], actual=ids)
    assert ids == ["prot_<0..2>"]


# -- runtime provenance columns from multi-axis source ------------------------

def test_mock_runtime_emits_axis_provenance_columns(tmp_path, record_case):
    """Cartesian source should materialize {axis}.id provenance columns in the map_table."""
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    cfg = {
        "output_folder": str(out_dir),
        "parent_ids": ["a1+b1", "a1+b2", "a2+b1", "a2+b2"],
        "output_ids": ["a1+b1", "a1+b2", "a2+b1", "a2+b2"],
        "provenance": {
            "sequences": ["a1", "a1", "a2", "a2"],
            "compounds": ["b1", "b2", "b1", "b2"],
        },
        "axis_names": ["sequences", "compounds"],
        "children": None, "produce": None,
        "streams": {
            "out": {"format": "pdb", "file": "<id>.pdb", "values": None,
                    "map_table": str(out_dir / "out_map.csv")},
        },
        "tables": {},
        "map_table_strategy": "runtime", "missing": [], "source_streams": [],
    }
    cfg_path = tmp_path / "cfg.json"
    json.dump(cfg, open(cfg_path, "w"))
    _run_pipe_mock(cfg_path)

    with open(out_dir / "out_map.csv") as f:
        rows = list(csv.DictReader(f))
    cols = list(rows[0].keys())
    seq_prov = [r["sequences.id"] for r in rows]
    cmp_prov = [r["compounds.id"] for r in rows]

    record_case(input="multi-axis provenance @ runtime",
                expected=(["sequences.id", "compounds.id"],
                          ["a1", "a1", "a2", "a2"],
                          ["b1", "b2", "b1", "b2"]),
                actual=([c for c in cols if c.endswith(".id")],
                        seq_prov, cmp_prov))
    assert "sequences.id" in cols and "compounds.id" in cols
    assert seq_prov == ["a1", "a1", "a2", "a2"]
    assert cmp_prov == ["b1", "b2", "b1", "b2"]


# -- validation rejections ----------------------------------------------------

def test_mock_rejects_produce_without_children(record_case):
    from biopipelines.mock import Mock
    record_case(input="Mock(ids=['a'], produce=['_1A'])",
                expected="ValueError", actual="ValueError")
    with pytest.raises(ValueError):
        Mock(ids=["a"], produce=["_1A"])


def test_mock_rejects_invalid_map_table_strategy(record_case):
    from biopipelines.mock import Mock
    record_case(input="map_table_strategy='nope'",
                expected="ValueError", actual="ValueError")
    with pytest.raises(ValueError):
        Mock(ids=["a"], map_table_strategy="nope")


def test_mock_rejects_stream_values_length_mismatch(record_case):
    from biopipelines.mock import Mock
    record_case(input="values length != parent_ids length",
                expected="ValueError", actual="ValueError")
    with pytest.raises(ValueError):
        Mock(
            ids=["a", "b", "c"],
            streams={"s": {"format": "csv", "values": ["only_one"]}},
        )


def test_mock_rejects_table_rows_length_mismatch(record_case):
    from biopipelines.mock import Mock
    record_case(input="table rows length != parent_ids length",
                expected="ValueError", actual="ValueError")
    with pytest.raises(ValueError):
        Mock(
            ids=["a", "b"],
            tables={"t": {"columns": ["score"], "rows": [["a", 0.1]]}},
        )


def test_mock_rejects_neither_ids_nor_source(record_case):
    from biopipelines.mock import Mock
    record_case(input="Mock() — no ids/no source",
                expected="ValueError", actual="ValueError")
    with pytest.raises(ValueError):
        Mock(streams={"s": {"format": "pdb", "file": "<id>.pdb"}})


# -- StandardizedOutput dot-notation (developer_manual §Tool Development) -----

def test_mock_standardized_output_dot_notation(
    local_config, isolated_cwd, new_pipeline, record_case,
):
    """Per developer_manual.md line 481 and user_manual.md §DataStream vs Tables:
    a tool assignment inside a Pipeline context yields a StandardizedOutput
    exposing dot-notation access to streams, tables, and output_folder."""
    from biopipelines.mock import Mock
    from biopipelines.base_config import StandardizedOutput

    pipeline = new_pipeline("mock_stdout_dot")
    with pipeline:
        m = Mock(
            ids=["a", "b"],
            streams={"structures": {"format": "pdb", "file": "<id>.pdb"}},
            tables={"scores": {"columns": ["score"], "fill": 0.0}},
        )
        pipeline.save()

    record_case(
        input="StandardizedOutput dot-notation",
        expected=("StandardizedOutput", ["a", "b"], "scores", True),
        actual=(type(m).__name__,
                list(m.streams.structures.ids),
                m.tables.scores.info.name,
                bool(m.output_folder)),
    )
    assert isinstance(m, StandardizedOutput)
    assert list(m.streams.structures.ids) == ["a", "b"]
    assert m.tables.scores.info.name == "scores"
    assert m.output_folder  # populated


# -- runtime file-template <id> substitution ----------------------------------

def test_mock_runtime_file_template_substitution(tmp_path, record_case):
    """Each expanded ID gets its own file at <output_folder>/<id>.pdb."""
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    cfg = {
        "output_folder": str(out_dir),
        "parent_ids": ["prot1", "prot2", "prot3"],
        "output_ids": ["prot1", "prot2", "prot3"],
        "provenance": {}, "axis_names": [],
        "children": None, "produce": None,
        "streams": {
            "s": {"format": "pdb", "file": "<id>.pdb", "values": None,
                  "map_table": str(out_dir / "s_map.csv")},
        },
        "tables": {},
        "map_table_strategy": "runtime", "missing": [], "source_streams": [],
    }
    cfg_path = tmp_path / "cfg.json"
    json.dump(cfg, open(cfg_path, "w"))
    _run_pipe_mock(cfg_path)

    produced = sorted(f for f in os.listdir(out_dir) if f.endswith(".pdb"))
    # Each map_table 'file' column points to the corresponding <id>.pdb
    with open(out_dir / "s_map.csv") as f:
        rows = list(csv.DictReader(f))
    file_cells = [os.path.basename(r["file"]) for r in rows]

    record_case(
        input="<id>.pdb template × 3 IDs",
        expected=(["prot1.pdb", "prot2.pdb", "prot3.pdb"],
                  ["prot1.pdb", "prot2.pdb", "prot3.pdb"]),
        actual=(produced, file_cells),
    )
    assert produced == ["prot1.pdb", "prot2.pdb", "prot3.pdb"]
    assert file_cells == ["prot1.pdb", "prot2.pdb", "prot3.pdb"]


# -- value-stream runtime with missing filter --------------------------------

def test_mock_runtime_value_stream_with_missing(tmp_path, record_case):
    """Value-based streams (no file template) still get a map_table, and
    `missing` drops rows."""
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    cfg = {
        "output_folder": str(out_dir),
        "parent_ids": ["s1", "s2", "s3"],
        "output_ids": ["s1", "s2", "s3"],
        "provenance": {}, "axis_names": [],
        "children": None, "produce": None,
        "streams": {
            "seqs": {"format": "csv", "file": None,
                     "values": ["MKT", "AET", "GFT"],
                     "map_table": str(out_dir / "seqs_map.csv")},
        },
        "tables": {},
        "map_table_strategy": "runtime",
        "missing": ["s2"], "source_streams": [],
    }
    cfg_path = tmp_path / "cfg.json"
    json.dump(cfg, open(cfg_path, "w"))
    _run_pipe_mock(cfg_path)

    with open(out_dir / "seqs_map.csv") as f:
        rows = list(csv.DictReader(f))
    record_case(
        input="value stream, missing=['s2']",
        expected=(["s1", "s3"], ["MKT", "GFT"]),
        actual=([r["id"] for r in rows], [r["value"] for r in rows]),
    )
    assert [r["id"] for r in rows] == ["s1", "s3"]
    assert [r["value"] for r in rows] == ["MKT", "GFT"]


# -- multi-stream + multi-table composition ----------------------------------

def test_mock_multiple_streams_and_tables(
    local_config, isolated_cwd, new_pipeline, record_case,
):
    from biopipelines.mock import Mock

    pipeline = new_pipeline("mock_multi")
    with pipeline:
        m = Mock(
            ids=["a", "b"],
            streams={
                "structures": {"format": "pdb", "file": "<id>.pdb"},
                "sequences":  {"format": "csv", "values": ["MKT", "AET"]},
            },
            tables={
                "scores":    {"columns": ["plddt"], "fill": 0.9},
                "distances": {"columns": ["rmsd"],  "fill": 1.1},
            },
        )
        pipeline.save()

    record_case(
        input="2 streams + 2 tables",
        expected=({"structures", "sequences"}, {"scores", "distances"}),
        actual=(set(m.streams.keys()),
                set(m.tables._tables.keys())),
    )
    assert set(m.streams.keys()) == {"structures", "sequences"}
    assert set(m.tables._tables.keys()) == {"scores", "distances"}


# ── user_manual.md coverage: DataStream iteration + column references ───────

def test_mock_stream_iteration_yields_single_item_datastreams(
    local_config, isolated_cwd, new_pipeline, record_case,
):
    """Per user_manual.md: iterating a stream yields single-ID DataStream items
    that can be passed directly to downstream tools."""
    from biopipelines.datastream import DataStream
    from biopipelines.mock import Mock

    pipeline = new_pipeline("mock_stream_iter")
    with pipeline:
        m = Mock(
            ids=["a", "b", "c"],
            streams={"structures": {"format": "pdb", "file": "<id>.pdb"}},
        )
        pipeline.save()

    items = list(m.streams.structures)
    per_item_ids = [list(item.ids) for item in items]
    all_datastreams = all(isinstance(item, DataStream) for item in items)

    record_case(
        input="for x in m.streams.structures",
        expected=(3, [["a"], ["b"], ["c"]], True),
        actual=(len(items), per_item_ids, all_datastreams),
    )
    assert len(items) == 3
    assert per_item_ids == [["a"], ["b"], ["c"]]
    assert all_datastreams


def test_mock_tool_iteration_yields_single_item_standardized_outputs(
    local_config, isolated_cwd, new_pipeline, record_case,
):
    """Per user_manual.md: `for item in tool` yields single-item
    StandardizedOutputs when all streams share the same IDs."""
    from biopipelines.base_config import StandardizedOutput
    from biopipelines.mock import Mock

    pipeline = new_pipeline("mock_tool_iter")
    with pipeline:
        m = Mock(
            ids=["a", "b"],
            streams={
                "structures": {"format": "pdb", "file": "<id>.pdb"},
                "sequences":  {"format": "csv", "values": ["MKT", "AET"]},
            },
        )
        pipeline.save()

    items = list(m)
    types_ok = all(isinstance(it, StandardizedOutput) for it in items)
    per_item_struct_ids = [list(it.streams.structures.ids) for it in items]
    per_item_seq_ids = [list(it.streams.sequences.ids) for it in items]

    record_case(
        input="for item in m (shared-ID streams)",
        expected=(2, True, [["a"], ["b"]], [["a"], ["b"]]),
        actual=(len(items), types_ok, per_item_struct_ids, per_item_seq_ids),
    )
    assert len(items) == 2
    assert types_ok
    assert per_item_struct_ids == [["a"], ["b"]]
    assert per_item_seq_ids == [["a"], ["b"]]


def test_mock_table_column_reference_string(
    local_config, isolated_cwd, new_pipeline, record_case,
):
    """Per user_manual.md §Table Column References: `tool.tables.x.col` yields
    a TABLE_REFERENCE string that downstream tools resolve at runtime via
    biopipelines_io.load_table."""
    from biopipelines.mock import Mock
    from biopipelines.biopipelines_io import load_table

    pipeline = new_pipeline("mock_table_ref")
    with pipeline:
        m = Mock(
            ids=["a", "b"],
            streams={"structures": {"format": "pdb", "file": "<id>.pdb"}},
            tables={"scores": {"columns": ["plddt"], "fill": 0.5}},
        )
        pipeline.save()

    ref = m.tables.scores.plddt
    ref_str = str(ref)
    # Also exercise the full runtime round-trip: emit the table, then parse.
    _run_pipe_mock(os.path.join(m.output_folder, "_configuration", "mock_config.json"))
    table, column = load_table(ref_str)

    record_case(
        input="m.tables.scores.plddt → TABLE_REFERENCE + load_table",
        expected=("TABLE_REFERENCE:", ":plddt", "plddt", 2),
        actual=(ref_str[:16], ref_str[-6:], column, len(table)),
    )
    assert ref_str.startswith("TABLE_REFERENCE:")
    assert ref_str.endswith(":plddt")
    assert ref.column == "plddt"
    assert column == "plddt"
    assert len(table) == 2
