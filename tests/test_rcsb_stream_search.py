"""RCSB compounds-stream runtime regressions."""

import importlib.util
import json
import os
import subprocess
import sys

import pandas as pd

from biopipelines.datastream import DataStream

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HELPER = os.path.join(REPO_ROOT, "pipe_scripts", "pipe_rcsb_search.py")


def _load_helper():
    spec = importlib.util.spec_from_file_location("pipe_rcsb_search", HELPER)
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


def test_helper_imports_when_run_outside_the_repository(tmp_path):
    """The helper runs on a node whose cwd is not the checkout and whose env does
    not have biopipelines installed, so it must bootstrap sys.path itself."""
    env = {k: v for k, v in os.environ.items() if k != "PYTHONPATH"}
    result = subprocess.run(
        [sys.executable, HELPER, "--help"],
        cwd=str(tmp_path), env=env, capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr


def test_runtime_search_uses_the_shared_request_builder(tmp_path, monkeypatch):
    """A non-score sort must emit sort_by=<attribute path>, not the rejected
    sort_by="attribute" form, and must honour per-field direction."""
    module = _load_helper()
    source = _stream_json(tmp_path, [{"id": "dye_1", "smiles": "CCO"}], ["dye_1"])
    fetch_config = tmp_path / "fetch.json"
    fetch_config.write_text(json.dumps({"pdb_ids": [], "custom_ids": []}))
    config = tmp_path / "search.json"
    config.write_text(json.dumps({
        "query_nodes": [{"type": "terminal", "service": "chemical",
                         "parameters": {"value": ""}}],
        "stream_query_index": 0,
        "logical_operator": "and",
        "return_type": "entry",
        "max_results": 10,
        "total_max_results": None,
        "sort_field": "release_date",
        "compounds_stream": str(source),
        "fetch_config": str(fetch_config),
        "hits_table": str(tmp_path / "hits.csv"),
    }))

    seen = []
    monkeypatch.setattr(module, "run_search", lambda request: seen.append(request) or [])
    monkeypatch.setattr(sys, "argv", ["pipe_rcsb_search.py", "--config", str(config)])
    module.main()

    sort = seen[0]["request_options"]["sort"][0]
    assert sort == {"sort_by": "rcsb_accession_info.initial_release_date",
                    "direction": "desc"}
    assert "attribute" not in sort


def test_entry_matched_by_several_compounds_keeps_every_provenance(tmp_path, monkeypatch):
    """One entry matched by two compounds writes one canonical provenance row per
    relationship while the structure fetch remains deduplicated."""
    module = _load_helper()
    source = _stream_json(
        tmp_path,
        [{"id": "dye_1", "smiles": "CCO"}, {"id": "dye_2", "smiles": "CCN"}],
        ["dye_1", "dye_2"],
    )
    fetch_config = tmp_path / "fetch.json"
    fetch_config.write_text(json.dumps({"pdb_ids": [], "custom_ids": []}))
    hits_table = tmp_path / "hits.csv"
    config = tmp_path / "search.json"
    config.write_text(json.dumps({
        "query_nodes": [{"type": "terminal", "service": "chemical",
                         "parameters": {"value": ""}}],
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

    # Both compounds hit 4HHB; the second scores higher.
    calls = []

    def fake_search(request):
        calls.append(request)
        score = 0.5 if len(calls) == 1 else 0.9
        return [{"identifier": "4HHB", "score": score}]

    monkeypatch.setattr(module, "run_search", fake_search)
    monkeypatch.setattr(sys, "argv", ["pipe_rcsb_search.py", "--config", str(config)])
    module.main()

    hits = pd.read_csv(hits_table)
    assert len(hits) == 2
    assert set(hits["compounds.id"]) == {"dye_1", "dye_2"}
    assert set(hits["id"]) == {"4hhb"}
    scores = hits.set_index("compounds.id")["score"].to_dict()
    assert scores == {"dye_1": 0.5, "dye_2": 0.9}
    # The structure fetch stays unique even though the relationship table is long.
    assert json.loads(fetch_config.read_text())["custom_ids"] == ["4hhb"]
