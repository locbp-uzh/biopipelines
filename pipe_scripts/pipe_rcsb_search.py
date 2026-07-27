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

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from biopipelines.biopipelines_io import iterate_values, load_datastream
from biopipelines.rcsb import SEARCH_API_URL, build_search_request, extract_pdb_id


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
    for compound_id, smiles in pairs:
        nodes = [dict(node) for node in config["query_nodes"]]
        parameters = dict(nodes[slot]["parameters"])
        parameters["value"] = smiles
        nodes[slot] = dict(nodes[slot], parameters=parameters)

        request = build_search_request(
            nodes, config["logical_operator"], return_type,
            config["max_results"], config["sort_field"],
        )
        results = run_search(request)

        found = []
        for result in results:
            identifier = result.get("identifier", "")
            pdb_id = extract_pdb_id(identifier, return_type)
            found.append(pdb_id)
            score = result.get("score", 0.0)
            hit = hits.get(pdb_id)
            if hit is None:
                hits[pdb_id] = {
                    "pdb_id": pdb_id,
                    "result_id": identifier,
                    "score": score,
                    # Several compounds can match one entry. Keep one relationship
                    # record per match; the search-results table writes these in
                    # long form so compounds.id remains canonical and joinable.
                    "matches": [{
                        "compound_id": compound_id,
                        "smiles": smiles,
                        "result_id": identifier,
                        "score": score,
                    }],
                }
            else:
                # Rank the entry by its best match across compounds.
                if score > hit["score"]:
                    hit["score"] = score
                    hit["result_id"] = identifier
                hit["matches"].append({
                    "compound_id": compound_id,
                    "smiles": smiles,
                    "result_id": identifier,
                    "score": score,
                })
        print(f"  {compound_id}: {len(found)} hits", flush=True)

    total_cap = config.get("total_max_results")
    ordered = sorted(hits.values(), key=lambda h: h["score"], reverse=True)
    if total_cap is not None and len(ordered) > total_cap:
        print(f"  capping {len(ordered)} unique entries to total_max_results={total_cap}")
        ordered = ordered[:total_cap]

    pdb_ids = [hit["pdb_id"] for hit in ordered]
    output_ids = [pdb_id.lower() for pdb_id in pdb_ids]
    print(f"RCSB search: {len(pdb_ids)} unique entries across {len(pairs)} compounds", flush=True)

    # One row per entry/compound relationship. The fetch config below remains
    # deduplicated by entry, while this table keeps canonical single-valued
    # compounds.id provenance for filtering and equality joins.
    with open(config["hits_table"], "w", newline='', encoding='utf-8') as handle:
        writer = csv.writer(handle)
        writer.writerow([
            "id", "pdb_id", "result_id", "score", "compounds.id", "query_smiles",
        ])
        for output_id, hit in zip(output_ids, ordered):
            for match in hit["matches"]:
                writer.writerow([
                    output_id, hit["pdb_id"], match["result_id"], match["score"],
                    match["compound_id"], match["smiles"],
                ])

    with open(config["fetch_config"]) as handle:
        fetch_config = json.load(handle)
    fetch_config["pdb_ids"] = pdb_ids
    fetch_config["custom_ids"] = output_ids
    with open(config["fetch_config"], "w") as handle:
        json.dump(fetch_config, handle, indent=2)


if __name__ == "__main__":
    main()
