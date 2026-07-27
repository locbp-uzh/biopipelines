"""Ligand PubChem vendor-enrichment regressions."""

import importlib.util
import json
import os
import sys

import pandas as pd

from biopipelines.datastream import DataStream


REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HELPER = os.path.join(REPO_ROOT, "pipe_scripts", "pipe_ligand.py")


def _load_helper():
    spec = importlib.util.spec_from_file_location("pipe_ligand_vendor", HELPER)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


class _Response:
    def __init__(self, payload):
        self.payload = payload

    def raise_for_status(self):
        return None

    def json(self):
        return self.payload


def test_vendor_enrichment_preserves_rows_and_appends_evidence(tmp_path, monkeypatch):
    source = tmp_path / "compounds.csv"
    pd.DataFrame([
        {"id": "dye_a", "format": "csv", "smiles": "CCO", "cid": "",
         "brightness": 12.5},
    ]).to_csv(source, index=False)
    ds = DataStream(
        name="compounds", ids=["dye_a"], files=[],
        map_table=str(source), format="csv",
    )
    ds_json = tmp_path / "compounds.json"
    ds.save_json(str(ds_json))

    module = _load_helper()

    def fake_get(url, **kwargs):
        if "/cids/JSON" in url:
            return _Response({"IdentifierList": {"CID": [702]}})
        return _Response({
            "SourceCategories": {
                "Categories": [{
                    "Category": "Chemical Vendors",
                    "Sources": [
                        {"SourceName": "Vendor A",
                         "SourceURL": "https://vendor-a.example"},
                        {"SourceName": "Vendor B",
                         "SourceURL": "https://vendor-b.example"},
                    ],
                }],
            },
        })

    monkeypatch.setattr("requests.get", fake_get)
    monkeypatch.setattr(module.time, "sleep", lambda _seconds: None)
    output = tmp_path / "enriched.csv"
    failed = tmp_path / "failed.csv"
    module.enrich_vendor_evidence({
        "enrich_compounds_json": str(ds_json),
        "compounds_table": str(output),
        "failed_table": str(failed),
    })

    row = pd.read_csv(output).iloc[0]
    assert row["id"] == "dye_a"
    assert row["brightness"] == 12.5
    assert row["cid"] == 702
    assert row["vendor_status"] == "available"
    assert row["vendor_count"] == 2
    assert json.loads(row["vendors"]) == ["Vendor A", "Vendor B"]
    assert json.loads(row["vendor_urls"]) == [
        "https://vendor-a.example",
        "https://vendor-b.example",
    ]


def test_vendor_lookup_wrapper_emits_runtime_enrichment_config(
    local_config, isolated_cwd, new_pipeline,
):
    from biopipelines.ligand import Ligand

    source = isolated_cwd / "library.csv"
    pd.DataFrame([{"id": "dye_a", "smiles": "CCO"}]).to_csv(source, index=False)
    compounds = DataStream(
        name="compounds", ids=["dye_a"], files=[],
        map_table=str(source), format="csv",
    )
    pipeline = new_pipeline("ligand_vendor")
    with pipeline:
        enriched = Ligand(compounds=compounds, vendor_lookup=True)
        pipeline.save()

    config = json.loads(open(pipeline.tools[-1].config_file, encoding="utf-8").read())
    assert config["enrich_compounds_json"]
    assert enriched.tables.compounds.info.columns[-6:] == [
        "vendor_status", "vendor_count", "vendors", "vendor_urls",
        "vendor_checked_at", "vendor_error",
    ]
