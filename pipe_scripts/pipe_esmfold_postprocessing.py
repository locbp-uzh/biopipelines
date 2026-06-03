#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Post-processing helper for ESMFold.

Reads the per-sequence PDB files and ESMFold_scores.json produced by
pipe_esmfold_inference.py, then builds the map table and merged metrics CSV.

This script runs under the 'biopipelines' conda environment.

Usage:
    python pipe_esmfold.py \\
        --sequences-json <path> \\
        --predictions-dir <path> \\
        --map-csv <path> \\
        --merged-csv <path>
"""

import argparse
import json
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from biopipelines.biopipelines_io import load_datastream

MERGED_COLUMNS = ["id", "file", "plddt", "ptm"]
MAP_COLUMNS    = ["id", "file", "sequences.id"]


def load_scores(predictions_dir: str) -> dict:
    """Load ESMFold_scores.json, return a dict keyed by sequence id."""
    scores_path = os.path.join(predictions_dir, "ESMFold_scores.json")
    if not os.path.exists(scores_path):
        print(f"Warning: ESMFold_scores.json not found at {scores_path}", file=sys.stderr)
        return {}
    with open(scores_path) as f:
        raw = json.load(f)
    # raw schema: {description: [...], plddt: [...], ptm: [...], perresidue_plddt: [...]}
    return {
        desc: {"plddt": plddt, "ptm": ptm}
        for desc, plddt, ptm in zip(raw["description"], raw["plddt"], raw["ptm"])
    }


def find_pdb(predictions_dir: str, seq_id: str) -> str:
    """Return path to <seq_id>.pdb in predictions_dir, or empty string if absent."""
    path = os.path.join(predictions_dir, f"{seq_id}.pdb")
    return path if os.path.exists(path) else ""


def build_outputs(sequences_ds, predictions_dir: str, map_csv: str, merged_csv: str) -> int:
    """Build map table and merged metrics CSV. Returns number of structures found."""
    scores = load_scores(predictions_dir)

    map_rows    = []
    merged_rows = []

    for seq_id in sequences_ds.ids_expanded:
        pdb_path = find_pdb(predictions_dir, seq_id)
        score    = scores.get(seq_id, {})

        map_rows.append({
            "id":           seq_id,
            "file":         pdb_path,
            "sequences.id": seq_id,
        })

        if pdb_path:
            plddt = score.get("plddt")
            ptm   = score.get("ptm")
            merged_rows.append({
                "id":    seq_id,
                "file":  pdb_path,
                "plddt": round(plddt, 2) if plddt is not None else None,
                "ptm":   round(ptm, 3) if ptm is not None else None,
            })
            plddt_str = f"{plddt:.1f}" if plddt is not None else "N/A"
            ptm_str   = f"{ptm:.3f}"   if ptm is not None else "N/A"
            print(f"  {seq_id}: pLDDT={plddt_str}  pTM={ptm_str}")
        else:
            print(f"Warning: PDB not found for {seq_id}", file=sys.stderr)

    pd.DataFrame(map_rows, columns=MAP_COLUMNS).to_csv(map_csv, index=False)
    print(f"Map table: {map_csv} ({len(map_rows)} entries)")

    if merged_rows:
        pd.DataFrame(merged_rows, columns=MERGED_COLUMNS).to_csv(merged_csv, index=False)
        print(f"Merged CSV: {merged_csv} ({len(merged_rows)} structures)")
    else:
        pd.DataFrame(columns=MERGED_COLUMNS).to_csv(merged_csv, index=False)
        print("Warning: No ESMFold structures found — wrote empty merged CSV", file=sys.stderr)

    return len(merged_rows)


def main():
    parser = argparse.ArgumentParser(description="ESMFold post-processor")
    parser.add_argument("--sequences-json", required=True)
    parser.add_argument("--predictions-dir", required=True)
    parser.add_argument("--map-csv",         required=True)
    parser.add_argument("--merged-csv",      required=True)
    args = parser.parse_args()

    sequences_ds = load_datastream(args.sequences_json)
    if not sequences_ds.ids_expanded:
        print("Error: No sequences found in input DataStream", file=sys.stderr)
        sys.exit(1)

    n = len(sequences_ds.ids_expanded)
    print(f"Post-processing ESMFold output for {n} sequence(s)")

    success = build_outputs(
        sequences_ds=sequences_ds,
        predictions_dir=args.predictions_dir,
        map_csv=args.map_csv,
        merged_csv=args.merged_csv,
    )

    if success == 0:
        print("Error: No ESMFold structures found", file=sys.stderr)
        sys.exit(1)

    if success < n:
        print(f"Warning: {n - success}/{n} sequences have no output structure", file=sys.stderr)


if __name__ == "__main__":
    main()
