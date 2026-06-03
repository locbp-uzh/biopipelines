#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
DiffDock post-processing.

DiffDock at commit a6c5275 (the pinned version) writes results to --out_dir as
folders named `index{N}_<sanitised_protein_path>____<ligand_smiles>`, NOT as
`<complex_name>/` like the wrapper's CSV suggests. It also writes
`complex_names.npy` listing the CSV-row complex names in order; row i lines up
with the folder whose name starts `index{i}_`.

Inside each folder:
  - `rank1.sdf`                          — top-ranked pose (no confidence in name)
  - `rank{N}_confidence{score}.sdf`      — for N >= 2; score is float, signed

This script:
  1. Maps each `index{i}_*` folder back to the CSV's complex_name via
     complex_names.npy (or by sort order if numpy/the file is missing).
  2. Parses each `rank{N}[_confidence{X}].sdf` filename and copies it to
     <structures-folder>/<complex_name>_rank{N}.sdf so the names line up
     with the BioPipelines IDs the wrapper predicted.
  3. Writes the structures map_table CSV (id, file, structures.id,
     compounds.id) and the confidence table CSV (id, structures.id,
     compounds.id, rank, confidence).
  4. Writes missing.csv for any pair whose folder is missing or contains
     fewer ranks than requested.
"""

import argparse
import csv
import os
import re
import shutil
import sys
from typing import Dict, List, Optional, Tuple

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from biopipelines.biopipelines_io import step_id_from_table_path  # noqa: E402


RANK_RE = re.compile(r"^rank(?P<rank>\d+)(?:_confidence(?P<conf>-?\d+(?:\.\d+)?))?\.sdf$")
INDEX_DIR_RE = re.compile(r"^index(?P<idx>\d+)_")


def parse_input_csv(path: str) -> List[Tuple[str, str, str]]:
    """Return [(complex_name, protein_id, ligand_id), ...] from the DiffDock
    input CSV (complex_name is `<protein>+<ligand>` by construction).
    """
    out = []
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            complex_name = row["complex_name"]
            if "+" not in complex_name:
                raise ValueError(
                    f"complex_name {complex_name!r} does not contain '+'; the input CSV "
                    "must be the one written by pipe_diffdock_build_csv.py."
                )
            prot, lig = complex_name.split("+", 1)
            out.append((complex_name, prot, lig))
    return out


def collect_ranks(complex_dir: Optional[str]) -> List[Tuple[int, float, str]]:
    """Return sorted [(rank, confidence, abs_path), ...] inside `complex_dir`.

    DiffDock at a6c5275 writes rank 1 TWICE: once as `rank1.sdf` (no
    confidence in the name) and once as `rank1_confidence<X>.sdf`. Deduplicate
    by rank, keeping the entry that carries a confidence score so the table
    has the score for rank 1 too.
    """
    by_rank: Dict[int, Tuple[float, str]] = {}
    if not complex_dir or not os.path.isdir(complex_dir):
        return []
    for fname in os.listdir(complex_dir):
        m = RANK_RE.match(fname)
        if not m:
            continue
        rank = int(m.group("rank"))
        conf_str = m.group("conf")
        conf = float(conf_str) if conf_str is not None else float("nan")
        path = os.path.join(complex_dir, fname)
        # Prefer entries that carry a confidence value over the ones that don't.
        if rank not in by_rank or (conf == conf and by_rank[rank][0] != by_rank[rank][0]):
            by_rank[rank] = (conf, path)
    return [(rank, conf, path) for rank, (conf, path) in sorted(by_rank.items())]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw-out", required=True, help="DiffDock --out_dir (execution folder)")
    ap.add_argument("--input-csv", required=True, help="protein_ligand.csv from build step")
    ap.add_argument("--structures-folder", required=True, help="Destination for rank<N>.sdf files")
    ap.add_argument("--structures-map", required=True)
    ap.add_argument("--confidence-csv", required=True)
    ap.add_argument("--missing-csv", required=True)
    ap.add_argument("--samples-per-complex", type=int, required=True)
    args = ap.parse_args()

    pairs = parse_input_csv(args.input_csv)

    os.makedirs(args.structures_folder, exist_ok=True)
    os.makedirs(os.path.dirname(args.structures_map), exist_ok=True)
    os.makedirs(os.path.dirname(args.confidence_csv), exist_ok=True)
    os.makedirs(os.path.dirname(args.missing_csv), exist_ok=True)

    # Build index -> folder map from DiffDock's actual output layout.
    index_to_dir: Dict[int, str] = {}
    if os.path.isdir(args.raw_out):
        for entry in os.listdir(args.raw_out):
            m = INDEX_DIR_RE.match(entry)
            if not m:
                continue
            full = os.path.join(args.raw_out, entry)
            if os.path.isdir(full):
                index_to_dir[int(m.group("idx"))] = full

    map_rows: List[Dict[str, str]] = []
    conf_rows: List[Dict[str, str]] = []
    missing_rows: List[Dict[str, str]] = []
    step_id = step_id_from_table_path(args.missing_csv)

    for i, (complex_name, prot_id, lig_id) in enumerate(pairs):
        complex_dir: Optional[str] = index_to_dir.get(i)
        ranks = collect_ranks(complex_dir) if complex_dir else []
        if not ranks:
            missing_rows.append({
                "id": complex_name,
                "removed_by": step_id,
                "kind": "failure",
                "cause": "no rank*.sdf produced",
            })
            continue
        if len(ranks) < args.samples_per_complex:
            missing_rows.append({
                "id": complex_name,
                "removed_by": step_id,
                "kind": "failure",
                "cause": f"only {len(ranks)}/{args.samples_per_complex} ranks produced",
            })

        for rank, conf, src in ranks:
            out_id = f"{complex_name}_rank{rank}"
            dst = os.path.join(args.structures_folder, f"{out_id}.sdf")
            shutil.copyfile(src, dst)
            map_rows.append({
                "id": out_id,
                "file": dst,
                "structures.id": prot_id,
                "compounds.id": lig_id,
            })
            conf_rows.append({
                "id": out_id,
                "structures.id": prot_id,
                "compounds.id": lig_id,
                "rank": str(rank),
                "confidence": "" if conf != conf else f"{conf:.4f}",  # NaN-safe
            })

    # Write structures map_table (overwrites the config-time map with the
    # actual rows that DiffDock produced — failed pairs are dropped).
    with open(args.structures_map, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["id", "file", "structures.id", "compounds.id"])
        w.writeheader()
        w.writerows(map_rows)

    with open(args.confidence_csv, "w", newline="") as f:
        w = csv.DictWriter(
            f, fieldnames=["id", "structures.id", "compounds.id", "rank", "confidence"]
        )
        w.writeheader()
        w.writerows(conf_rows)

    with open(args.missing_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["id", "removed_by", "kind", "cause"])
        w.writeheader()
        w.writerows(missing_rows)

    print(
        f"DiffDock post-process: {len(map_rows)} poses written, "
        f"{len(missing_rows)} pair(s) flagged as missing/partial"
    )

    if not map_rows:
        print("ERROR: DiffDock produced no usable poses", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
