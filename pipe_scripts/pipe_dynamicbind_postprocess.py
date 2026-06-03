#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
DynamicBind post-processing.

After the per-protein inference loop, DynamicBind writes:
    {results_root}/{prot_id}/index{N}_idx_{M}/
        affinity_prediction.csv
        rank1_ligand_lddt0.XX_affinity5.XX_relaxed.sdf
        rank2_ligand_lddt0.XX_affinity4.XX_relaxed.sdf
        ...
        rank1_receptor_relaxed.pdb
        ...

We flatten this to BioPipelines IDs `<prot>+<lig>_rank<N>`, copy the pose SDFs
into the structures stream folder, and build the affinity table by combining
the parsed lDDT/affinity from each SDF filename with the ligand mapping in
each protein's input CSV.
"""

import argparse
import csv
import os
import re
import shutil
import sys
from typing import Dict, List, Tuple

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from biopipelines.biopipelines_io import step_id_from_table_path  # noqa: E402


POSE_RE = re.compile(
    r"^rank(?P<rank>\d+)_ligand_lddt(?P<lddt>-?\d+(?:\.\d+)?)_affinity(?P<affinity>-?\d+(?:\.\d+)?)_relaxed\.sdf$"
)


def load_ligand_csv(path: str) -> List[str]:
    """Return ligand IDs in the same order as the per-protein CSV rows."""
    with open(path, newline="") as f:
        return [row["name"] for row in csv.DictReader(f)]


RECEPTOR_RE = re.compile(r"^rank(?P<rank>\d+)_receptor_relaxed\.pdb$")


def collect_protein_poses(prot_dir: str) -> Dict[int, List[Tuple[int, float, float, str, str]]]:
    """Return {ligand_index: [(rank, lddt, affinity, sdf_path, receptor_pdb), ...]}
    for a given protein output folder. ligand_index is the integer `M` from the
    `index{N}_idx_{M}/` folder name (upstream uses idx_M to identify which
    row of the ligand CSV the prediction came from). receptor_pdb is the
    matching `rank{N}_receptor_relaxed.pdb` (empty string if absent).
    """
    out: Dict[int, List[Tuple[int, float, float, str, str]]] = {}
    if not os.path.isdir(prot_dir):
        return out

    for entry in os.listdir(prot_dir):
        sub = os.path.join(prot_dir, entry)
        if not os.path.isdir(sub):
            continue
        m = re.match(r"^index(\d+)_idx_(\d+)$", entry)
        if not m:
            continue
        lig_idx = int(m.group(2))
        out.setdefault(lig_idx, [])

        # Index receptor PDBs by rank so each pose can carry its receptor.
        receptors: Dict[int, str] = {}
        for fname in os.listdir(sub):
            rm = RECEPTOR_RE.match(fname)
            if rm:
                receptors[int(rm.group("rank"))] = os.path.join(sub, fname)

        for fname in os.listdir(sub):
            pm = POSE_RE.match(fname)
            if not pm:
                continue
            rank = int(pm.group("rank"))
            out[lig_idx].append((
                rank,
                float(pm.group("lddt")),
                float(pm.group("affinity")),
                os.path.join(sub, fname),
                receptors.get(rank, ""),
            ))
        out[lig_idx].sort(key=lambda t: t[0])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results-root", required=True)
    ap.add_argument("--ligand-csvs-dir", required=True)
    ap.add_argument("--structures-folder", required=True)
    ap.add_argument("--structures-map", required=True)
    ap.add_argument("--affinity-csv", required=True)
    ap.add_argument("--missing-csv", required=True)
    ap.add_argument("--num-saved", type=int, required=True,
                    help="poses kept per pair (upstream --savings_per_complex)")
    ap.add_argument("--movie-jobs-csv", default="",
                    help="if set, write (id, receptor_pdb, ligand_sdf) rows for the movie renderer")
    args = ap.parse_args()

    os.makedirs(args.structures_folder, exist_ok=True)
    os.makedirs(os.path.dirname(args.structures_map), exist_ok=True)
    os.makedirs(os.path.dirname(args.affinity_csv), exist_ok=True)
    os.makedirs(os.path.dirname(args.missing_csv), exist_ok=True)

    map_rows: List[Dict[str, str]] = []
    aff_rows: List[Dict[str, str]] = []
    missing_rows: List[Dict[str, str]] = []
    movie_jobs: List[Dict[str, str]] = []
    step_id = step_id_from_table_path(args.missing_csv)

    # Discover proteins from the ligand_csvs directory (1 csv per protein).
    if not os.path.isdir(args.ligand_csvs_dir):
        print(f"ERROR: ligand csvs dir missing: {args.ligand_csvs_dir}", file=sys.stderr)
        sys.exit(1)

    for csv_name in sorted(os.listdir(args.ligand_csvs_dir)):
        if not csv_name.endswith(".csv"):
            continue
        prot_id = csv_name[:-4]
        ligand_ids = load_ligand_csv(os.path.join(args.ligand_csvs_dir, csv_name))

        prot_dir = os.path.join(args.results_root, prot_id)
        poses_by_idx = collect_protein_poses(prot_dir)

        for lig_idx, lig_id in enumerate(ligand_ids):
            pair_id = f"{prot_id}+{lig_id}"
            poses = poses_by_idx.get(lig_idx, [])
            if not poses:
                missing_rows.append({
                    "id": pair_id,
                    "removed_by": step_id,
                    "kind": "failure",
                    "cause": "no poses produced",
                })
                continue
            if len(poses) < args.num_saved:
                missing_rows.append({
                    "id": pair_id,
                    "removed_by": step_id,
                    "kind": "failure",
                    "cause": f"only {len(poses)}/{args.num_saved} ranks produced",
                })

            for rank, lddt, affinity, src, receptor_pdb in poses:
                out_id = f"{pair_id}_rank{rank}"
                dst = os.path.join(args.structures_folder, f"{out_id}.sdf")
                shutil.copyfile(src, dst)

                map_rows.append({
                    "id": out_id,
                    "file": dst,
                    "structures.id": prot_id,
                    "compounds.id": lig_id,
                })
                if args.movie_jobs_csv and receptor_pdb:
                    movie_jobs.append({
                        "id": out_id,
                        "receptor_pdb": receptor_pdb,
                        "ligand_sdf": dst,
                    })
                aff_rows.append({
                    "id": out_id,
                    "structures.id": prot_id,
                    "compounds.id": lig_id,
                    "rank": str(rank),
                    "lddt": f"{lddt:.4f}",
                    "affinity": f"{affinity:.4f}",
                })

    with open(args.structures_map, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["id", "file", "structures.id", "compounds.id"])
        w.writeheader()
        w.writerows(map_rows)

    with open(args.affinity_csv, "w", newline="") as f:
        w = csv.DictWriter(
            f, fieldnames=["id", "structures.id", "compounds.id", "rank", "lddt", "affinity"]
        )
        w.writeheader()
        w.writerows(aff_rows)

    with open(args.missing_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["id", "removed_by", "kind", "cause"])
        w.writeheader()
        w.writerows(missing_rows)

    if args.movie_jobs_csv:
        os.makedirs(os.path.dirname(args.movie_jobs_csv), exist_ok=True)
        with open(args.movie_jobs_csv, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["id", "receptor_pdb", "ligand_sdf"])
            w.writeheader()
            w.writerows(movie_jobs)
        print(f"Movie jobs: {args.movie_jobs_csv} ({len(movie_jobs)} pose(s))")

    print(
        f"DynamicBind post-process: {len(map_rows)} poses kept, "
        f"{len(missing_rows)} pair(s) failed"
    )

    if not map_rows:
        print("ERROR: DynamicBind produced no usable poses", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
