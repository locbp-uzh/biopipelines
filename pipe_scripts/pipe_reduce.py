#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Reduce runner. For each input PDB, runs `reduce -BUILD` and writes the
protonated structure (stdout) to <output_dir>/<sid>.pdb plus a map CSV.

The reduce binary exits 1 when atoms are added (its normal success path);
exit 2 means a real error. The runner therefore treats both 0 and 1 as
success."""

import argparse
import os
import subprocess
import sys

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from biopipelines.biopipelines_io import load_datastream, iterate_files  # noqa: E402


def run_reduce(pdb_in: str, pdb_out: str) -> None:
    cmd = ["reduce", "-BUILD", "-NUClear", "-Quiet", pdb_in]
    with open(pdb_out, "w") as f:
        res = subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE, text=True)
    if res.returncode not in (0, 1) or not os.path.getsize(pdb_out):
        raise RuntimeError(f"reduce failed (exit {res.returncode}): {res.stderr.strip()}")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--structures-json", required=True)
    p.add_argument("--output-dir", required=True)
    p.add_argument("--map-csv", required=True)
    args = p.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    ds = load_datastream(args.structures_json)

    rows, failed = [], []
    for sid, pdb_path in iterate_files(ds):
        try:
            out_path = os.path.join(args.output_dir, f"{sid}.pdb")
            run_reduce(pdb_path, out_path)
            rows.append({"id": sid, "file": out_path})
            print(f"  {sid}: -> {out_path}")
        except Exception as e:
            print(f"WARNING: {sid} Reduce failed: {e}", file=sys.stderr)
            failed.append(sid)

    os.makedirs(os.path.dirname(args.map_csv), exist_ok=True)
    pd.DataFrame(rows, columns=["id", "file"]).to_csv(args.map_csv, index=False)
    print(f"Map: {args.map_csv} ({len(rows)} rows)")

    if failed:
        print(f"Failed: {len(failed)}/{len(failed)+len(rows)}: {failed}", file=sys.stderr)
    if not rows:
        sys.exit(1)


if __name__ == "__main__":
    main()
