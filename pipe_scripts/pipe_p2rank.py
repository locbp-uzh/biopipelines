#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""P2Rank runner. For each input structure, runs `prank predict -f <file>`,
then parses P2Rank's `<basename>_predictions.csv` (pockets) and
`<basename>_residues.csv` (per-residue scores), and writes per-pocket,
per-residue, and summary CSVs with biopipelines IDs.

P2Rank CSV format (comma-separated, values space-padded — headers and cells
are stripped on read):

  <name>_predictions.csv:
    name, rank, score, probability, sas_points, surf_atoms,
    center_x, center_y, center_z, residue_ids, surf_atom_ids
    (residue_ids is a space-separated list of "<chain>_<resnum>" tokens)

  <name>_residues.csv:
    chain, residue_label, residue_name, score, zscore, probability, pocket
    (pocket = rank of the pocket the residue belongs to; 0 = unassigned)
"""

import argparse
import os
import shutil
import subprocess
import sys

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from biopipelines.biopipelines_io import load_datastream, iterate_files, step_id_from_table_path, container_argv_prefix  # noqa: E402
from biopipelines.sele_utils import chain_aware_sele  # noqa: E402


POCKET_COLS = ["id", "pocket_idx", "rank", "score", "probability",
               "n_residues", "residues", "center_x", "center_y", "center_z"]
RESIDUE_COLS = ["id", "chain", "resi", "resn", "pocket_idx", "score", "probability"]
SUMMARY_COLS = ["id", "n_pockets", "top_score", "top_probability", "top_residues"]


def _read_p2rank_csv(path: str) -> pd.DataFrame:
    """Read a P2Rank CSV, stripping whitespace from headers and string cells."""
    df = pd.read_csv(path, skipinitialspace=True)
    df.columns = [c.strip() for c in df.columns]
    for col in df.columns:
        if df[col].dtype == object:
            df[col] = df[col].astype(str).str.strip()
    return df


def _residue_tokens_to_sele(residue_ids: str):
    """Convert P2Rank 'A_12 A_45 B_3' tokens to (chain, resnum) tuples."""
    residues = []
    for tok in str(residue_ids).split():
        if "_" not in tok:
            continue
        chain, _, resnum = tok.rpartition("_")
        try:
            residues.append((chain, int(resnum)))
        except ValueError:
            continue
    return residues


def run_prank(pdb_path: str, config: str, work_dir: str,
              threads: int = 1, visualizations: bool = False,
              container_prefix: str = "") -> str:
    """Run `prank predict` and return the output directory."""
    os.makedirs(work_dir, exist_ok=True)
    cmd = container_argv_prefix(container_prefix) + ["prank", "predict", "-f", pdb_path, "-c", config,
           "-o", work_dir,
           "-threads", str(threads),
           "-visualizations", "1" if visualizations else "0"]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"prank failed: {result.stderr.strip() or result.stdout.strip()}")
    return work_dir


def parse_pockets(pred_csv: str, sid: str):
    """Parse predictions CSV -> list of pocket row dicts."""
    rows = []
    if not os.path.exists(pred_csv):
        return rows
    df = _read_p2rank_csv(pred_csv)
    for _, r in df.iterrows():
        residues = _residue_tokens_to_sele(r.get("residue_ids", ""))
        rows.append({
            "id": sid,
            "pocket_idx": int(r["rank"]),
            "rank": int(r["rank"]),
            "score": float(r["score"]),
            "probability": float(r["probability"]),
            "n_residues": len(residues),
            "residues": chain_aware_sele(residues),
            "center_x": float(r["center_x"]),
            "center_y": float(r["center_y"]),
            "center_z": float(r["center_z"]),
        })
    return rows


def parse_residues(res_csv: str, sid: str):
    """Parse residues CSV -> list of residue row dicts."""
    rows = []
    if not os.path.exists(res_csv):
        return rows
    df = _read_p2rank_csv(res_csv)
    for _, r in df.iterrows():
        rows.append({
            "id": sid,
            "chain": r["chain"],
            "resi": int(r["residue_label"]),
            "resn": r["residue_name"],
            "pocket_idx": int(r["pocket"]),
            "score": float(r["score"]),
            "probability": float(r["probability"]),
        })
    return rows


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--structures-json", required=True)
    p.add_argument("--config", default="default")
    p.add_argument("--threads", type=int, default=1)
    p.add_argument("--visualizations", type=int, default=0)
    p.add_argument("--scratch-dir", required=True)
    p.add_argument("--pockets-csv", required=True)
    p.add_argument("--residues-csv", required=True)
    p.add_argument("--residues-dir", required=True)
    p.add_argument("--residues-map-csv", required=True)
    p.add_argument("--summary-csv", required=True)
    p.add_argument("--missing-csv", required=True)
    p.add_argument("--container-prefix", default="")
    p.add_argument("--upstream-missing", default=None)
    args = p.parse_args()

    os.makedirs(args.scratch_dir, exist_ok=True)
    os.makedirs(args.residues_dir, exist_ok=True)
    ds = load_datastream(args.structures_json)

    pocket_rows, residue_rows, summary_rows, missing_rows, residue_map_rows = [], [], [], [], []
    step_id = step_id_from_table_path(args.missing_csv)

    for sid, pdb_path in iterate_files(ds):
        try:
            work_dir = os.path.join(args.scratch_dir, sid)
            run_prank(pdb_path, args.config, work_dir,
                      threads=args.threads, visualizations=bool(args.visualizations),
                      container_prefix=args.container_prefix)
        except Exception as e:
            print(f"WARNING: {sid} P2Rank failed: {e}", file=sys.stderr)
            missing_rows.append({"id": sid, "removed_by": step_id, "kind": "failure", "cause": str(e)[:200]})
            continue

        base = os.path.splitext(os.path.basename(pdb_path))[0]
        pred_csv = os.path.join(work_dir, f"{base}.pdb_predictions.csv")
        res_csv = os.path.join(work_dir, f"{base}.pdb_residues.csv")
        # P2Rank names outputs after the full input filename. Fall back to a
        # glob match if the .pdb-suffixed name is not found (e.g. .cif inputs).
        if not os.path.exists(pred_csv):
            pred_csv = _find_one(work_dir, "_predictions.csv")
        if not os.path.exists(res_csv):
            res_csv = _find_one(work_dir, "_residues.csv")

        pockets = parse_pockets(pred_csv, sid) if pred_csv else []
        residues = parse_residues(res_csv, sid) if res_csv else []

        if not pockets and not residues:
            missing_rows.append({"id": sid, "removed_by": step_id, "kind": "failure", "cause": "P2Rank produced no parseable output"})
            continue

        pocket_rows.extend(pockets)
        residue_rows.extend(residues)

        # One resi-csv file per structure for the `residues` stream.
        res_path = os.path.join(args.residues_dir, f"{sid}.csv")
        pd.DataFrame(residues, columns=RESIDUE_COLS).to_csv(res_path, index=False)
        residue_map_rows.append({"id": sid, "file": res_path})

        top = pockets[0] if pockets else None
        summary_rows.append({
            "id": sid,
            "n_pockets": len(pockets),
            "top_score": top["score"] if top else 0.0,
            "top_probability": top["probability"] if top else 0.0,
            "top_residues": top["residues"] if top else "",
        })
        print(f"  {sid}: {len(pockets)} pockets, top score {top['score'] if top else 0:.3f}")

    upstream_rows = []
    if args.upstream_missing and os.path.exists(args.upstream_missing):
        try:
            up_df = pd.read_csv(args.upstream_missing)
            if not up_df.empty:
                upstream_rows = up_df.to_dict("records")
        except Exception as e:
            print(f"Warning: could not read upstream missing.csv: {e}", file=sys.stderr)

    all_missing = upstream_rows + missing_rows

    for d in (args.pockets_csv, args.residues_csv, args.residues_map_csv,
              args.summary_csv, args.missing_csv):
        os.makedirs(os.path.dirname(d), exist_ok=True)
    pd.DataFrame(pocket_rows, columns=POCKET_COLS).to_csv(args.pockets_csv, index=False)
    pd.DataFrame(residue_rows, columns=RESIDUE_COLS).to_csv(args.residues_csv, index=False)
    pd.DataFrame(residue_map_rows, columns=["id", "file"]).to_csv(args.residues_map_csv, index=False)
    pd.DataFrame(summary_rows, columns=SUMMARY_COLS).to_csv(args.summary_csv, index=False)
    pd.DataFrame(all_missing, columns=["id", "removed_by", "kind", "cause"]).to_csv(args.missing_csv, index=False)
    print(f"Pockets: {args.pockets_csv} ({len(pocket_rows)} rows)")
    print(f"Residues: {args.residues_csv} ({len(residue_rows)} rows)")
    print(f"Residue resi-csv files: {args.residues_map_csv} ({len(residue_map_rows)} rows)")
    print(f"Summary: {args.summary_csv} ({len(summary_rows)} rows)")
    print(f"Missing: {args.missing_csv} ({len(all_missing)} rows)")

    if missing_rows:
        print(f"Failed: {len(missing_rows)}/{len(missing_rows)+len(summary_rows)}", file=sys.stderr)
    if not summary_rows:
        sys.exit(1)


def _find_one(directory: str, suffix: str):
    """Return the single file in `directory` ending in `suffix`, or None."""
    if not os.path.isdir(directory):
        return None
    matches = [os.path.join(directory, f) for f in os.listdir(directory) if f.endswith(suffix)]
    return matches[0] if matches else None


if __name__ == "__main__":
    main()
