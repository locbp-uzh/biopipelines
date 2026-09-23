#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Collect OpenFold3's predictions into the framework's streams and tables.

`run_openfold` writes `<out>/<query>/seed_<N>/<query>_seed_<N>_sample_<K>_model.cif` plus the
confidence JSON beside it. This turns that into one `<id>` per complex (or `<id>_1..K` when the
user asked for every sample), a structures map_table, and a confidence table.

Three artifacts, kept separate on purpose: the confidence table has a row for every id that
ENTERED, with NaN where a number is unavailable; the stream map has only the ids whose file was
really written; the local missing table says why an id has no structure. A query that failed is
in all three.
"""

import argparse
import csv
import glob
import json
import os
import re
import sys
from typing import Dict, List, Optional, Tuple

_biopipelines_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'biopipelines')
sys.path.insert(0, _biopipelines_dir)
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from combinatorics import CombinatoricsConfig

CONFIDENCE_COLUMNS = ["id", "plddt", "ptm", "iptm", "gpde", "has_clash",
                      "disorder", "ranking_score", "seed", "sample"]

# The aggregated JSON's own spellings. Read off a real run rather than the docs: `avg_plddt`,
# `disorder` and `sample_ranking_score` are not what a reader would guess, and guessing left
# three columns silently empty in the first cluster run.
CONFIDENCE_KEYS = {
    "plddt": "avg_plddt",
    "ptm": "ptm",
    "iptm": "iptm",
    "gpde": "gpde",
    "has_clash": "has_clash",
    "disorder": "disorder",
    "ranking_score": "sample_ranking_score",
}

# `<query>_seed_<N>_sample_<K>_model.<ext>`
SAMPLE_RE = re.compile(r"_seed_(?P<seed>\d+)_sample_(?P<sample>\d+)_model\.(?:cif|pdb)$")


def _rank(scores: Dict[str, Optional[float]]) -> float:
    """A sample's ranking score, with an unscored sample sorting last rather than crashing."""
    value = scores.get("ranking_score")
    return value if isinstance(value, (int, float)) else float("-inf")


def parse_bool(text: str) -> bool:
    return str(text).strip().lower() in ("1", "true", "yes", "on")


def discover(prediction_folder: str, query: str, extension: str) -> List[Tuple[int, int, str]]:
    """(seed, sample, path) for one query's structures, ordered by seed then sample."""
    pattern = os.path.join(prediction_folder, query, "seed_*", f"*_model.{extension}")
    found = []
    for path in glob.glob(pattern):
        match = SAMPLE_RE.search(os.path.basename(path))
        if match:
            found.append((int(match.group("seed")), int(match.group("sample")), path))
    return sorted(found)


def read_confidence(structure_path: str) -> Dict[str, Optional[float]]:
    """The aggregated scores written beside a structure, or empty when absent.

    Absent is not an error: the structure is the result, and a confidence file that failed to
    write should leave NaN in a column rather than discard a prediction that exists.
    """
    base = structure_path.rsplit("_model.", 1)[0]
    path = f"{base}_confidences_aggregated.json"
    if not os.path.exists(path):
        return {}
    try:
        with open(path) as handle:
            data = json.load(handle)
    except (OSError, json.JSONDecodeError):
        return {}
    return {column: data.get(key) for column, key in CONFIDENCE_KEYS.items()}


def main():
    parser = argparse.ArgumentParser(description="Collect OpenFold3 predictions")
    parser.add_argument("--prediction-folder", required=True)
    parser.add_argument("--combinatorics-config", required=True)
    parser.add_argument("--queries-json", default=None)
    parser.add_argument("--structures-dir", required=True)
    parser.add_argument("--structures-map-csv", required=True)
    parser.add_argument("--confidence-csv", required=True)
    parser.add_argument("--local-missing-csv", required=True)
    parser.add_argument("--output-format", default="cif")
    parser.add_argument("--top-only", default="True")
    args = parser.parse_args()

    top_only = parse_bool(args.top_only)
    extension = "pdb" if args.output_format == "pdb" else "cif"

    # The queries that entered are the ones queries.json named; a query run_openfold never made a
    # folder for (a featurization or MSA failure) must still get a row, so the folders are not the list.
    if args.queries_json and os.path.exists(args.queries_json):
        with open(args.queries_json) as handle:
            queries = sorted((json.load(handle).get("queries") or {}).keys())
    else:
        queries = sorted(name for name in os.listdir(args.prediction_folder)
                         if os.path.isdir(os.path.join(args.prediction_folder, name))
                         and not name.startswith("_") and name != "msas")
    if not queries:
        print(f"ERROR: no prediction folders under {args.prediction_folder} — "
              f"run_openfold produced nothing")
        sys.exit(1)

    for directory in (args.structures_dir, os.path.dirname(args.structures_map_csv),
                      os.path.dirname(args.confidence_csv),
                      os.path.dirname(args.local_missing_csv)):
        if directory:
            os.makedirs(directory, exist_ok=True)

    map_rows, confidence_rows, missing_rows = [], [], []
    produced = 0

    for query in queries:
        samples = discover(args.prediction_folder, query, extension)
        if not samples:
            confidence_rows.append({"id": query, **{c: "" for c in CONFIDENCE_COLUMNS[1:]}})
            missing_rows.append({
                "id": query, "removed_by": "OpenFold3", "kind": "failure",
                "cause": f"run_openfold wrote no {extension} for this query"})
            continue

        scored = [(seed, sample, path, read_confidence(path)) for seed, sample, path in samples]
        if top_only:
            # "Best", not "first". run_openfold writes every sample and ranks them by
            # sample_ranking_score; taking samples[0] took whichever had the lowest sample
            # number, which is an arbitrary structure presented as the chosen one.
            chosen = [max(scored, key=lambda row: _rank(row[3]))]
        else:
            chosen = scored

        for index, (seed, sample, path, scores) in enumerate(chosen, start=1):
            out_id = query if top_only else f"{query}_{index}"
            destination = os.path.join(args.structures_dir, f"{out_id}.{extension}")
            with open(path, "rb") as source, open(destination, "wb") as target:
                target.write(source.read())
            map_rows.append({"id": out_id, "file": destination})
            row = {"id": out_id, "seed": seed, "sample": sample}
            row.update({k: v for k, v in scores.items() if v is not None})
            confidence_rows.append(row)
            produced += 1

    with open(args.structures_map_csv, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["id", "file"])
        writer.writeheader()
        writer.writerows(map_rows)

    with open(args.confidence_csv, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=CONFIDENCE_COLUMNS)
        writer.writeheader()
        for row in confidence_rows:
            writer.writerow({c: row.get(c, "") for c in CONFIDENCE_COLUMNS})

    with open(args.local_missing_csv, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["id", "removed_by", "kind", "cause"])
        writer.writeheader()
        writer.writerows(missing_rows)

    print(f"Collected {produced} structure(s) for {len(queries)} quer(y/ies); "
          f"{len(missing_rows)} produced nothing")

    # Every query failed: the step did not work, and exiting 0 would rubber-stamp it.
    if produced == 0:
        print("ERROR: no OpenFold3 prediction produced a structure")
        sys.exit(1)


if __name__ == "__main__":
    main()
