# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Execution-time RCSB search for queries whose descriptor comes from a compounds stream.

Resolves the stream's SMILES, runs one search per compound, unions the hits, then
writes the resolved PDB ids into the fetch config that pipe_pdb.py consumes.
"""

import argparse
import csv
import json
import os
import sys

from biopipelines.biopipelines_io import iterate_values, load_datastream

SEARCH_API_URL = "https://search.rcsb.org/rcsbsearch/v2/query"


def load_smiles(stream_source):
    """Read canonical SMILES through the runtime DataStream resolver."""
    stream = load_datastream(stream_source)
    pairs = []
    for compound_id, values in iterate_values(stream, columns=["smiles"]):
        value = values["smiles"]
        if value is None or value != value:
            continue
        smiles = str(value).strip()
        if not smiles:
            continue
        pairs.append((compound_id, smiles))
    return pairs


def build_request(query_nodes, logical_operator, return_type, max_results, sort_field):
    nodes = []
    for index, node in enumerate(query_nodes):
        node = dict(node)
        node["node_id"] = index
        nodes.append(node)

    query_node = nodes[0] if len(nodes) == 1 else {
        "type": "group",
        "logical_operator": logical_operator,
        "nodes": nodes,
    }

    request = {
        "query": query_node,
        "return_type": return_type,
        "request_options": {"paginate": {"start": 0, "rows": max_results}},
    }
    if sort_field != "score":
        request["request_options"]["sort"] = [
            {"sort_by": "attribute", "attribute": sort_field, "direction": "asc"}
        ]
    return request


def run_search(request):
    """Return the result_set for one request; [] when RCSB reports no matches."""
    import requests

    response = requests.post(
        SEARCH_API_URL, json=request,
        headers={"Content-Type": "application/json"}, timeout=30,
    )
    if response.status_code == 204:
        return []
    if response.status_code != 200:
        raise SystemExit(
            f"RCSB search failed with status {response.status_code}: {response.text[:200]}"
        )
    return response.json().get("result_set", [])


def extract_pdb_id(identifier, return_type):
    if return_type == "assembly":
        return identifier.split("-")[0]
    if return_type == "polymer_entity":
        return identifier.split("_")[0]
    if return_type == "polymer_instance":
        return identifier.split(".")[0]
    return identifier


def main():
    parser = argparse.ArgumentParser(description="RCSB search with per-compound fan-out")
    parser.add_argument("--config", required=True, help="Search configuration JSON")
    args = parser.parse_args()

    with open(args.config) as handle:
        config = json.load(handle)

    pairs = load_smiles(config["compounds_stream"])

    slot = config["stream_query_index"]
    return_type = config["return_type"]
    print(f"RCSB search: {len(pairs)} compounds, one search each", flush=True)

    hits = {}
    per_compound = []
    for compound_id, smiles in pairs:
        nodes = [dict(node) for node in config["query_nodes"]]
        parameters = dict(nodes[slot]["parameters"])
        parameters["value"] = smiles
        nodes[slot] = dict(nodes[slot], parameters=parameters)

        request = build_request(
            nodes, config["logical_operator"], return_type,
            config["max_results"], config["sort_field"],
        )
        results = run_search(request)

        found = []
        for result in results:
            identifier = result.get("identifier", "")
            pdb_id = extract_pdb_id(identifier, return_type)
            found.append(pdb_id)
            if pdb_id not in hits:
                hits[pdb_id] = {
                    "pdb_id": pdb_id,
                    "result_id": identifier,
                    "score": result.get("score", 0.0),
                    "compound_id": compound_id,
                    "smiles": smiles,
                }
        per_compound.append((compound_id, smiles, found))
        print(f"  {compound_id}: {len(found)} hits", flush=True)

    total_cap = config.get("total_max_results")
    ordered = sorted(hits.values(), key=lambda h: h["score"], reverse=True)
    if total_cap is not None and len(ordered) > total_cap:
        print(f"  capping {len(ordered)} unique entries to total_max_results={total_cap}")
        ordered = ordered[:total_cap]

    pdb_ids = [hit["pdb_id"] for hit in ordered]
    output_ids = [pdb_id.lower() for pdb_id in pdb_ids]
    print(f"RCSB search: {len(pdb_ids)} unique entries across {len(pairs)} compounds", flush=True)

    with open(config["hits_table"], "w", newline='', encoding='utf-8') as handle:
        writer = csv.writer(handle)
        writer.writerow(["id", "pdb_id", "result_id", "score", "compounds.id", "query_smiles"])
        for output_id, hit in zip(output_ids, ordered):
            writer.writerow([
                output_id, hit["pdb_id"], hit["result_id"],
                hit["score"], hit["compound_id"], hit["smiles"],
            ])

    with open(config["fetch_config"]) as handle:
        fetch_config = json.load(handle)
    fetch_config["pdb_ids"] = pdb_ids
    fetch_config["custom_ids"] = output_ids
    with open(config["fetch_config"], "w") as handle:
        json.dump(fetch_config, handle, indent=2)


if __name__ == "__main__":
    main()
