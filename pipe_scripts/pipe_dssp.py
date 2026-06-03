#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""DSSP runner. For each PDB, calls `mkdssp` (PDB-REDO) and parses the mmCIF
or legacy DSSP output into per-residue rows + a per-structure summary."""

import argparse
import os
import shutil
import subprocess
import sys

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from biopipelines.biopipelines_io import load_datastream, iterate_files, step_id_from_table_path  # noqa: E402
from biopipelines.pdb_parser import relative_accessibility  # noqa: E402


SEC_COLS = ["id", "chain", "resnum", "resi", "resname", "resn", "ss_code", "ss_simple"]
SUM_COLS = ["id", "n_residues", "helix_frac", "sheet_frac", "coil_frac"]
# Per-residue resi-csv stream consumable by the Selection tool. The numeric acc/rsa
# and 0/1 is_* columns let Selection's "column op value" filter threshold on
# accessibility or SS class (Selection casts the threshold to float).
SS_COLS = ["id", "chain", "resi", "resn", "ss_code", "ss_simple", "acc", "rsa", "is_helix", "is_sheet", "is_coil"]

# Simplified 3-state mapping from DSSP 8-state.
SS_MAP = {
    "H": "helix", "G": "helix", "I": "helix",
    "E": "sheet", "B": "sheet",
    "T": "coil", "S": "coil", " ": "coil", "-": "coil",
}


def find_dssp_binary() -> str:
    for name in ("mkdssp", "dssp"):
        if shutil.which(name):
            return name
    raise RuntimeError("Neither mkdssp nor dssp found on PATH")


def run_dssp(binary: str, pdb_path: str, out_path: str):
    """Try modern mkdssp invocation first; fall back to legacy positional form."""
    cmd1 = [binary, "--output-format", "dssp", pdb_path, out_path]
    res = subprocess.run(cmd1, capture_output=True, text=True)
    if res.returncode == 0 and os.path.exists(out_path) and os.path.getsize(out_path) > 0:
        return
    cmd2 = [binary, pdb_path, out_path]
    res = subprocess.run(cmd2, capture_output=True, text=True)
    if res.returncode == 0 and os.path.exists(out_path) and os.path.getsize(out_path) > 0:
        return
    raise RuntimeError(f"{binary} failed: {res.stderr.strip() or res.stdout.strip()}")


def parse_dssp_classic(out_path: str):
    """Parse the legacy fixed-width DSSP output. Returns list of dicts."""
    rows = []
    in_data = False
    with open(out_path) as f:
        for line in f:
            if line.startswith("  #  RESIDUE"):
                in_data = True
                continue
            if not in_data or len(line) < 17:
                continue
            try:
                resnum = int(line[5:10].strip())
            except ValueError:
                continue
            chain = line[11].strip()
            resname = line[13:14].strip() or "?"  # one-letter code
            # DSSP marks disulfide-bonded cysteines with a lowercase letter.
            if resname.islower():
                resname = "C"
            ss_code = line[16] if len(line) > 16 else " "
            try:
                acc = int(line[34:38].strip())
            except ValueError:
                acc = None
            rows.append({
                "chain": chain or "A",
                "resnum": resnum,
                "resname": resname,
                "ss_code": ss_code.strip() or "-",
                "acc": acc,
            })
    return rows


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--structures-json", required=True)
    p.add_argument("--dssp-dir", required=True)
    p.add_argument("--dssp-map-csv", required=True)
    p.add_argument("--ss-dir", required=True)
    p.add_argument("--ss-map-csv", required=True)
    p.add_argument("--secondary-csv", required=True)
    p.add_argument("--summary-csv", required=True)
    p.add_argument("--missing-csv", required=True)
    p.add_argument("--upstream-missing", default=None)
    args = p.parse_args()

    os.makedirs(args.dssp_dir, exist_ok=True)
    os.makedirs(args.ss_dir, exist_ok=True)
    binary = find_dssp_binary()
    ds = load_datastream(args.structures_json)

    sec_rows, sum_rows, dssp_rows, ss_rows, missing_rows = [], [], [], [], []
    step_id = step_id_from_table_path(args.missing_csv)
    for sid, pdb_path in iterate_files(ds):
        try:
            out_path = os.path.join(args.dssp_dir, f"{sid}.dssp")
            run_dssp(binary, pdb_path, out_path)
            residues = parse_dssp_classic(out_path)
            if not residues:
                raise RuntimeError("DSSP output produced no residue rows")

            counts = {"helix": 0, "sheet": 0, "coil": 0}
            per_id_ss = []
            for r in residues:
                simple = SS_MAP.get(r["ss_code"], "coil")
                counts[simple] += 1
                sec_rows.append({
                    "id": sid,
                    "chain": r["chain"],
                    "resnum": r["resnum"],
                    "resi": r["resnum"],
                    "resname": r["resname"],
                    "resn": r["resname"],
                    "ss_code": r["ss_code"],
                    "ss_simple": simple,
                })
                acc = r.get("acc")
                rsa = relative_accessibility(acc, r["resname"]) if acc is not None else None
                per_id_ss.append({
                    "id": sid,
                    "chain": r["chain"],
                    "resi": r["resnum"],
                    "resn": r["resname"],
                    "ss_code": r["ss_code"],
                    "ss_simple": simple,
                    "acc": "" if acc is None else acc,
                    "rsa": "" if rsa is None else rsa,
                    "is_helix": int(simple == "helix"),
                    "is_sheet": int(simple == "sheet"),
                    "is_coil": int(simple == "coil"),
                })

            # One resi-csv file per structure for the `ss` stream.
            ss_path = os.path.join(args.ss_dir, f"{sid}.csv")
            pd.DataFrame(per_id_ss, columns=SS_COLS).to_csv(ss_path, index=False)
            ss_rows.append({"id": sid, "file": ss_path})

            n = len(residues)
            sum_rows.append({
                "id": sid,
                "n_residues": n,
                "helix_frac": round(counts["helix"] / n, 3),
                "sheet_frac": round(counts["sheet"] / n, 3),
                "coil_frac": round(counts["coil"] / n, 3),
            })
            dssp_rows.append({"id": sid, "file": out_path})
            print(f"  {sid}: {n} residues  H={counts['helix']/n:.0%}  E={counts['sheet']/n:.0%}")
        except Exception as e:
            print(f"WARNING: {sid} DSSP failed: {e}", file=sys.stderr)
            missing_rows.append({"id": sid, "removed_by": step_id, "kind": "failure", "cause": str(e)[:200]})

    upstream_rows = []
    if args.upstream_missing and os.path.exists(args.upstream_missing):
        try:
            up_df = pd.read_csv(args.upstream_missing)
            if not up_df.empty:
                upstream_rows = up_df.to_dict("records")
        except Exception as e:
            print(f"Warning: could not read upstream missing.csv: {e}", file=sys.stderr)

    all_missing = upstream_rows + missing_rows

    for d in (args.secondary_csv, args.summary_csv, args.dssp_map_csv, args.ss_map_csv, args.missing_csv):
        os.makedirs(os.path.dirname(d), exist_ok=True)
    pd.DataFrame(sec_rows, columns=SEC_COLS).to_csv(args.secondary_csv, index=False)
    pd.DataFrame(sum_rows, columns=SUM_COLS).to_csv(args.summary_csv, index=False)
    pd.DataFrame(dssp_rows, columns=["id", "file"]).to_csv(args.dssp_map_csv, index=False)
    pd.DataFrame(ss_rows, columns=["id", "file"]).to_csv(args.ss_map_csv, index=False)
    pd.DataFrame(all_missing, columns=["id", "removed_by", "kind", "cause"]).to_csv(args.missing_csv, index=False)
    print(f"Secondary structure: {args.secondary_csv} ({len(sec_rows)} rows)")
    print(f"Summary: {args.summary_csv} ({len(sum_rows)} rows)")
    print(f"DSSP files: {args.dssp_map_csv} ({len(dssp_rows)} rows)")
    print(f"SS resi-csv files: {args.ss_map_csv} ({len(ss_rows)} rows)")
    print(f"Missing: {args.missing_csv} ({len(all_missing)} rows)")

    if missing_rows:
        print(f"Failed: {len(missing_rows)}/{len(missing_rows)+len(sum_rows)}", file=sys.stderr)
    if not sum_rows:
        sys.exit(1)


if __name__ == "__main__":
    main()
