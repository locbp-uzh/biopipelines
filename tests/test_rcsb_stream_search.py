"""RCSB compounds-stream runtime regressions."""

import importlib.util
import json
import sys

import pandas as pd

from biopipelines.datastream import DataStream


def _load_helper():
    spec = importlib.util.spec_from_file_location(
        "pipe_rcsb_search", "pipe_scripts/pipe_rcsb_search.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _stream_json(tmp_path, rows, ids):
    table = tmp_path / "compounds.csv"
    pd.DataFrame(rows).to_csv(table, index=False)
    stream = DataStream(
        name="compounds",
        ids=ids,
        files=[],
        map_table=str(table),
        format="smiles",
    )
    source = tmp_path / "compounds.json"
    stream.save_json(str(source))
    return source


def test_rcsb_load_smiles_resolves_lazy_ids_from_runtime_map(tmp_path):
    module = _load_helper()
    source = _stream_json(
        tmp_path,
        [
            {"id": "6_Panda_1", "smiles": "CCO"},
            {"id": "6_Panda_2", "smiles": "CCN"},
        ],
        ["6_Panda_[<N>]"],
    )

    assert module.load_smiles(str(source)) == [
        ("6_Panda_1", "CCO"),
        ("6_Panda_2", "CCN"),
    ]


def test_rcsb_load_smiles_honours_concrete_stream_filter(tmp_path):
    module = _load_helper()
    source = _stream_json(
        tmp_path,
        [
            {"id": "dye_1", "smiles": "CCO"},
            {"id": "dye_2", "smiles": "CCN"},
        ],
        ["dye_2"],
    )

    assert module.load_smiles(str(source)) == [("dye_2", "CCN")]


def test_rcsb_zero_hits_writes_concrete_empty_fetch_config(tmp_path, monkeypatch):
    module = _load_helper()
    source = _stream_json(
        tmp_path,
        [{"id": "dye_1", "smiles": "CCO"}],
        ["dye_[<N>]"],
    )
    fetch_config = tmp_path / "fetch.json"
    fetch_config.write_text(json.dumps({
        "pdb_ids": [],
        "custom_ids": ["[<hits>]"],
    }))
    hits_table = tmp_path / "hits.csv"
    config = tmp_path / "search.json"
    config.write_text(json.dumps({
        "query_nodes": [{
            "type": "terminal",
            "service": "chemical",
            "parameters": {"value": ""},
        }],
        "stream_query_index": 0,
        "logical_operator": "and",
        "return_type": "entry",
        "max_results": 10,
        "total_max_results": None,
        "sort_field": "score",
        "compounds_stream": str(source),
        "fetch_config": str(fetch_config),
        "hits_table": str(hits_table),
    }))
    monkeypatch.setattr(module, "run_search", lambda request: [])
    monkeypatch.setattr(sys, "argv", ["pipe_rcsb_search.py", "--config", str(config)])

    module.main()

    assert pd.read_csv(hits_table).empty
    resolved = json.loads(fetch_config.read_text())
    assert resolved["pdb_ids"] == []
    assert resolved["custom_ids"] == []
