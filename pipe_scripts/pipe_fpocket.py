#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""FPocket runner. For each PDB, runs `fpocket -f <pdb>`, parses the resulting
<basename>_out/<basename>_info.txt + pocket residue files, writes per-pocket
and summary CSVs."""

import argparse
import os
import shutil
import subprocess
import sys
from collections import defaultdict

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from biopipelines.biopipelines_io import load_datastream, iterate_files, step_id_from_table_path  # noqa: E402


P_COLS = ["id", "pocket_idx", "druggability", "volume",
          "n_alpha_spheres", "n_residues", "residues", "pocket_file"]
S_COLS = ["id", "n_pockets", "top_druggability", "top_volume", "pymol_script"]


def parse_info(info_path: str):
    """Parse fpocket *_info.txt; return list of dicts (one per pocket)."""
    pockets = []
    if not os.path.exists(info_path):
        return pockets
    current = None
    with open(info_path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith("Pocket"):
                if current is not None:
                    pockets.append(current)
                current = {"pocket_idx": int(line.split()[1].rstrip(":"))}
            elif current is not None and "\t" in line:
                key, val = line.split("\t", 1)
                key = key.strip().lower().replace(" ", "_").rstrip(":")
                try:
                    current[key] = float(val.strip())
                except ValueError:
                    current[key] = val.strip()
    if current is not None:
        pockets.append(current)
    return pockets


def parse_residues(pocket_pdb: str):
    """Read a pocket_atm.pdb file and return a set of (chain, resnum, resname)."""
    if not os.path.exists(pocket_pdb):
        return []
    res = set()
    with open(pocket_pdb) as f:
        for line in f:
            if line.startswith("ATOM"):
                resname = line[17:20].strip()
                chain = line[21].strip()
                resnum = line[22:26].strip()
                res.add((chain, int(resnum), resname))
    return sorted(res, key=lambda r: (r[0], r[1]))


def run_fpocket(pdb_path: str, work_dir: str, min_alpha_spheres: int,
                min_radius: float, max_radius: float, clustering_distance: float) -> str:
    """Run fpocket and return the path to the *_out directory."""
    os.makedirs(work_dir, exist_ok=True)
    work_pdb = os.path.join(work_dir, os.path.basename(pdb_path))
    shutil.copyfile(pdb_path, work_pdb)

    cmd = ["fpocket", "-f", work_pdb,
           "-i", str(min_alpha_spheres),
           "-m", str(min_radius),
           "-M", str(max_radius),
           "-D", str(clustering_distance)]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"fpocket failed: {result.stderr.strip() or result.stdout.strip()}")
    base = os.path.splitext(os.path.basename(work_pdb))[0]
    return os.path.join(work_dir, f"{base}_out")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--structures-json", required=True)
    p.add_argument("--output-structures-dir", required=True)
    p.add_argument("--structures-map-csv", required=True)
    p.add_argument("--scratch-dir", required=True)
    p.add_argument("--min-alpha-spheres", type=int, default=35)
    p.add_argument("--min-radius", type=float, default=3.0)
    p.add_argument("--max-radius", type=float, default=6.0)
    p.add_argument("--clustering-distance", type=float, default=2.4)
    p.add_argument("--pockets-csv", required=True)
    p.add_argument("--summary-csv", required=True)
    p.add_argument("--missing-csv", required=True)
    p.add_argument("--upstream-missing", default=None)
    args = p.parse_args()

    os.makedirs(args.output_structures_dir, exist_ok=True)
    os.makedirs(args.scratch_dir, exist_ok=True)
    ds = load_datastream(args.structures_json)
    pocket_rows, sum_rows, structures_rows, missing_rows = [], [], [], []
    step_id = step_id_from_table_path(args.missing_csv)

    for sid, pdb_path in iterate_files(ds):
        try:
            work_dir = os.path.join(args.scratch_dir, sid)
            out_dir = run_fpocket(pdb_path, work_dir, args.min_alpha_spheres,
                                  args.min_radius, args.max_radius, args.clustering_distance)
        except Exception as e:
            print(f"WARNING: {sid} FPocket failed: {e}", file=sys.stderr)
            missing_rows.append({"id": sid, "removed_by": step_id, "kind": "failure", "cause": str(e)[:200]})
            continue

        base = os.path.splitext(os.path.basename(pdb_path))[0]
        annotated_src = os.path.join(out_dir, f"{base}_out.pdb")
        annotated_dst = os.path.join(args.output_structures_dir, f"{sid}.pdb")
        if not os.path.exists(annotated_src):
            missing_rows.append({"id": sid, "removed_by": step_id, "kind": "failure", "cause": "fpocket produced no annotated PDB"})
            continue
        shutil.copyfile(annotated_src, annotated_dst)
        structures_rows.append({"id": sid, "file": annotated_dst})

        info = parse_info(os.path.join(out_dir, f"{base}_info.txt"))
        pml_src = os.path.join(out_dir, f"{base}.pml")
        pml_path = pml_src if os.path.exists(pml_src) else ""

        for pocket in info:
            idx = pocket["pocket_idx"]
            pocket_pdb = os.path.join(out_dir, "pockets", f"pocket{idx}_atm.pdb")
            residues = parse_residues(pocket_pdb)
            res_str = "+".join(f"{c}{n}" for c, n, _ in residues) if residues else ""
            pocket_rows.append({
                "id": sid,
                "pocket_idx": idx,
                "druggability": pocket.get("druggability_score", ""),
                "volume": pocket.get("real_volume_(approximation)", pocket.get("volume", "")),
                "n_alpha_spheres": pocket.get("number_of_alpha_spheres", ""),
                "n_residues": len(residues),
                "residues": res_str,
                "pocket_file": pocket_pdb if os.path.exists(pocket_pdb) else "",
            })

        top_drug = max((p.get("druggability_score", 0) for p in info), default=0)
        top_vol = max((p.get("real_volume_(approximation)", p.get("volume", 0)) for p in info), default=0)
        sum_rows.append({
            "id": sid,
            "n_pockets": len(info),
            "top_druggability": round(top_drug, 3) if top_drug else 0,
            "top_volume": round(top_vol, 1) if top_vol else 0,
            "pymol_script": pml_path,
        })
        print(f"  {sid}: {len(info)} pockets, top druggability {top_drug:.3f}")

    upstream_rows = []
    if args.upstream_missing and os.path.exists(args.upstream_missing):
        try:
            up_df = pd.read_csv(args.upstream_missing)
            if not up_df.empty:
                upstream_rows = up_df.to_dict("records")
        except Exception as e:
            print(f"Warning: could not read upstream missing.csv: {e}", file=sys.stderr)

    all_missing = upstream_rows + missing_rows

    for d in (args.pockets_csv, args.summary_csv, args.structures_map_csv, args.missing_csv):
        os.makedirs(os.path.dirname(d), exist_ok=True)
    pd.DataFrame(pocket_rows, columns=P_COLS).to_csv(args.pockets_csv, index=False)
    pd.DataFrame(sum_rows, columns=S_COLS).to_csv(args.summary_csv, index=False)
    pd.DataFrame(structures_rows, columns=["id", "file"]).to_csv(args.structures_map_csv, index=False)
    pd.DataFrame(all_missing, columns=["id", "removed_by", "kind", "cause"]).to_csv(args.missing_csv, index=False)
    print(f"Pockets: {args.pockets_csv} ({len(pocket_rows)} rows)")
    print(f"Summary: {args.summary_csv} ({len(sum_rows)} rows)")
    print(f"Structures map: {args.structures_map_csv} ({len(structures_rows)} rows)")
    print(f"Missing: {args.missing_csv} ({len(all_missing)} rows)")

    if missing_rows:
        print(f"Failed: {len(missing_rows)}/{len(missing_rows)+len(sum_rows)}", file=sys.stderr)
    if not sum_rows:
        sys.exit(1)


if __name__ == "__main__":
    main()
