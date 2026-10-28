#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Collect Protenix's predictions into the framework's streams and tables.

`protenix pred` writes `<out>/<query>/seed_<N>/predictions/<query>_sample_<K>.cif` and `<query>_summary_confidence_sample_<K>.json` beside it.
Within a seed the samples are numbered by ranking score, so sample 0 is the best of that seed.
With --need_atom_confidence it also writes `<query>_full_data_sample_<K>.json`, whose `token_pair_pae` is the PAE matrix.
A query that fails is recorded under `<out>/ERR/<query>.txt` and the command still exits 0.

This turns that into one `<id>` per complex (or `<id>_1..N` when the user asked for every structure), a structures map_table, a confidence table and, with --pae-dir, one `<id>.npz` per structure.

Three artifacts, kept separate on purpose.
The confidence table has a row for every id that ENTERED, with NaN where a number is unavailable.
The stream maps have only the ids whose file was really written.
The local missing table says why an id has no structure; a query that failed or was skipped is in all three.
"""

import argparse
import csv
import json
import os
import re
import shutil
import sys
from typing import Dict, List, Optional, Tuple

import numpy as np

CONFIDENCE_COLUMNS = ["id", "plddt", "ptm", "iptm", "gpde", "has_clash",
                      "disorder", "ranking_score", "seed", "sample"]

# The keys of upstream's summary_confidence JSON that the table carries; plddt is on a 0 to 100 scale there.
SCORE_COLUMNS = CONFIDENCE_COLUMNS[1:-2]

CAUSE_LIMIT = 300


def _rank(scores: Dict[str, Optional[float]]) -> float:
    """A sample's ranking score, with an unscored sample sorting last rather than crashing."""
    value = scores.get("ranking_score")
    return value if isinstance(value, (int, float)) and not isinstance(value, bool) else float("-inf")


def parse_bool(text: str) -> bool:
    return str(text).strip().lower() in ("1", "true", "yes", "on")


def discover(prediction_folder: str, query: str) -> List[Tuple[int, int, str, str]]:
    """(seed, sample, structure path, file prefix) for one query, ordered by seed then sample.

    Files are `<query>_sample_<K>.cif` under `seed_<N>/predictions/`, or `<query>_<seed>_sample_<K>.cif` in the layout upstream's docs describe.
    Walks the query's own folder, so a query whose name is a prefix of another's cannot pick up its files.
    """
    pattern = re.compile(rf"^(?P<prefix>{re.escape(query)}(?:_(?P<seed>\d+))?)_sample_(?P<sample>\d+)\.cif$")
    found = []
    for directory, _, files in os.walk(os.path.join(prediction_folder, query)):
        for name in files:
            match = pattern.match(name)
            if not match:
                continue
            seed = match.group("seed")
            if seed is None:
                parent = re.search(r"seed_(\d+)", directory)
                seed = parent.group(1) if parent else 0
            found.append((int(seed), int(match.group("sample")),
                          os.path.join(directory, name), match.group("prefix")))
    return sorted(found)


def sibling(structure_path: str, prefix: str, kind: str, sample: int) -> str:
    """The `<prefix>_<kind>_sample_<K>.json` written beside a structure."""
    return os.path.join(os.path.dirname(structure_path), f"{prefix}_{kind}_sample_{sample}.json")


def read_confidence(path: str) -> Dict[str, Optional[float]]:
    """The scores of one sample, or empty when the file is absent.

    Absent is not an error: the structure is the result, and a confidence file that failed to write should leave NaN in a column rather than discard a prediction that exists.
    """
    if not os.path.exists(path):
        return {}
    try:
        with open(path) as handle:
            data = json.load(handle)
    except (OSError, json.JSONDecodeError):
        return {}
    return {column: data.get(column) for column in SCORE_COLUMNS}


def read_pae(path: str) -> Tuple[Optional[np.ndarray], str]:
    """(token x token PAE matrix, reason it is missing) from a full-data JSON; the matrix is None when it cannot be used."""
    if not os.path.exists(path):
        return None, f"{os.path.basename(path)} was not written"
    try:
        with open(path) as handle:
            data = json.load(handle)
        matrix = np.asarray(data["token_pair_pae"], dtype=np.float32)
    except (OSError, ValueError, KeyError, TypeError) as exc:
        return None, f"{os.path.basename(path)} has no readable token_pair_pae ({exc})"
    while matrix.ndim > 2 and matrix.shape[0] == 1:
        matrix = matrix[0]
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        return None, f"token_pair_pae in {os.path.basename(path)} is {matrix.shape}, not a square matrix"
    return matrix, ""


def read_error(prediction_folder: str, query: str) -> Optional[str]:
    """First line of the error upstream recorded for a query, if it recorded one."""
    path = os.path.join(prediction_folder, "ERR", f"{query}.txt")
    if not os.path.exists(path):
        return None
    with open(path, errors="replace") as handle:
        for line in handle:
            if line.strip():
                return line.strip()[:CAUSE_LIMIT]
    return None


def remove_full_data(prediction_folder: str, query: str) -> None:
    """Delete a query's full-data JSONs; each holds three N_token x N_token matrices and the PAE is already extracted."""
    for directory, _, files in os.walk(os.path.join(prediction_folder, query)):
        for name in files:
            if "_full_data_sample_" in name and name.endswith(".json"):
                os.remove(os.path.join(directory, name))


def read_ids(path: str) -> List[str]:
    with open(path, newline="") as handle:
        return [row["id"] for row in csv.DictReader(handle)]


def read_skipped(path: Optional[str]) -> Dict[str, str]:
    if not path or not os.path.exists(path):
        return {}
    with open(path, newline="") as handle:
        return {row["id"]: row["cause"] for row in csv.DictReader(handle)}


def write_rows(path: str, fieldnames: List[str], rows: List[Dict]) -> None:
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({name: row.get(name, "") for name in fieldnames})


def main():
    parser = argparse.ArgumentParser(description="Collect Protenix predictions")
    parser.add_argument("--prediction-folder", required=True)
    parser.add_argument("--query-ids-csv", required=True)
    parser.add_argument("--skipped-csv", default=None)
    parser.add_argument("--structures-dir", required=True)
    parser.add_argument("--structures-map-csv", required=True)
    parser.add_argument("--confidence-csv", required=True)
    parser.add_argument("--pae-dir", default=None)
    parser.add_argument("--pae-map-csv", default=None)
    parser.add_argument("--local-missing-csv", required=True)
    parser.add_argument("--top-only", default="True")
    args = parser.parse_args()

    top_only = parse_bool(args.top_only)
    queries = read_ids(args.query_ids_csv)
    skipped = read_skipped(args.skipped_csv)
    if not queries:
        print(f"ERROR: {args.query_ids_csv} lists no queries")
        sys.exit(1)

    for directory in (args.structures_dir, os.path.dirname(args.structures_map_csv),
                      os.path.dirname(args.confidence_csv), os.path.dirname(args.local_missing_csv),
                      args.pae_dir, os.path.dirname(args.pae_map_csv or "")):
        if directory:
            os.makedirs(directory, exist_ok=True)

    map_rows, confidence_rows, missing_rows, pae_rows = [], [], [], []
    produced = 0

    def fail(query: str, cause: str) -> None:
        confidence_rows.append({"id": query})
        missing_rows.append({"id": query, "removed_by": "Protenix", "kind": "failure", "cause": cause})

    for query in queries:
        if query in skipped:
            fail(query, skipped[query])
            continue

        samples = discover(args.prediction_folder, query)
        if not samples:
            fail(query, read_error(args.prediction_folder, query)
                 or "protenix wrote no structure for this query; see the step log")
            continue

        scored = [(seed, sample, path, prefix,
                   read_confidence(sibling(path, prefix, "summary_confidence", sample)))
                  for seed, sample, path, prefix in samples]
        # "Best", not "first": upstream numbers samples per seed, so across seeds the first file is arbitrary.
        chosen = [max(scored, key=lambda row: _rank(row[4]))] if top_only else scored

        for index, (seed, sample, path, prefix, scores) in enumerate(chosen, start=1):
            out_id = query if top_only else f"{query}_{index}"
            destination = os.path.join(args.structures_dir, f"{out_id}.cif")
            shutil.copyfile(path, destination)
            map_rows.append({"id": out_id, "file": destination})

            if args.pae_dir:
                matrix, reason = read_pae(sibling(path, prefix, "full_data", sample))
                if matrix is None:
                    print(f"WARNING: no PAE for {out_id}, left out of the pae stream: {reason}")
                else:
                    pae_destination = os.path.join(args.pae_dir, f"{out_id}.npz")
                    np.savez_compressed(pae_destination, pae=matrix)
                    pae_rows.append({"id": out_id, "file": pae_destination})

            row = {"id": out_id, "seed": seed, "sample": sample}
            row.update({k: v for k, v in scores.items() if v is not None})
            confidence_rows.append(row)
            produced += 1

        if args.pae_dir:
            remove_full_data(args.prediction_folder, query)

    write_rows(args.structures_map_csv, ["id", "file"], map_rows)
    if args.pae_map_csv:
        write_rows(args.pae_map_csv, ["id", "file"], pae_rows)
    write_rows(args.confidence_csv, CONFIDENCE_COLUMNS, confidence_rows)
    write_rows(args.local_missing_csv, ["id", "removed_by", "kind", "cause"], missing_rows)

    print(f"Collected {produced} structure(s) for {len(queries)} quer(y/ies); "
          f"{len(missing_rows)} produced nothing")

    # Every query failed: the step did not work, and exiting 0 would rubber-stamp it.
    if produced == 0:
        print("ERROR: no Protenix prediction produced a structure")
        sys.exit(1)


if __name__ == "__main__":
    main()
