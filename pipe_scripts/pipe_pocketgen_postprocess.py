#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
PocketGen post-processing.

Upstream PocketGen (`generate_new.py` / `models.PD.Pocket_Design_new.generate`)
hardcodes 8 samples per pair (line ~227: `datalist = [data for _ in range(8)]`).
For each staged pair folder, the model writes per-sample files indexed 0..7:

    {staging}/{pair}/{i}.pdb            — designed pocket only
    {staging}/{pair}/{i}.sdf            — ligand pose (copy of input SDF)
    {staging}/{pair}/{i}_relaxed.pdb    — designed pocket relaxed (OpenMM)
    {staging}/{pair}/{i}_whole.pdb      — full protein with designed pocket grafted
    {staging}/{pair}/{i}_whole_relaxed.pdb — same, relaxed (CANONICAL output)
    {staging}/{pair}/{i}_orig.pdb       — original pocket (pre-design)
    {staging}/{pair}/{i}_orig_relaxed.pdb — original, relaxed

This script flattens those into BioPipelines IDs `<pair>_<i>` and copies
`{i}_whole_relaxed.pdb` into the structures stream folder. Other artifacts
stay in the staging tree under _configuration/.

Outputs the content-bearing sequences CSV (id | structures.id | compounds.id
| sample | sequence) parsed from the per-sample whole_relaxed PDB.
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
from biopipelines.pdb_parser import field_res_name, field_chain, field_res_seq  # noqa: E402


# Three-to-one letter code map (standard residues + a few common variants).
THREE_TO_ONE = {
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F",
    "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L",
    "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R",
    "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y",
    # Common protonation variants
    "HID": "H", "HIE": "H", "HIP": "H",
    "CYX": "C", "CYM": "C",
}

# Per-sample file template: an integer prefix (sample index) followed by
# the canonical "whole_relaxed.pdb" suffix.
SAMPLE_WHOLE_RELAXED_RE = re.compile(r"^(\d+)_whole_relaxed\.pdb$")


def parse_sequence_from_pdb(pdb_path: str) -> str:
    """Extract a single-letter sequence from a PDB. Walks ATOM records,
    keeps one residue per (chain, resnum, icode) triple in file order. Non-
    standard residues become 'X'.
    """
    seen = set()
    seq: List[str] = []
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            try:
                resname = field_res_name(line)
                chain = field_chain(line)
                resnum = field_res_seq(line)
                icode = line[26:27].strip()
            except IndexError:
                continue
            key = (chain, resnum, icode)
            if key in seen:
                continue
            seen.add(key)
            seq.append(THREE_TO_ONE.get(resname, "X"))
    return "".join(seq)


def parse_pair_id(pair_dir_name: str) -> Tuple[str, str]:
    if "+" not in pair_dir_name:
        raise ValueError(f"staged folder {pair_dir_name!r} is not a '<prot>+<lig>' pair")
    prot, lig = pair_dir_name.split("+", 1)
    return prot, lig


def collect_samples(pair_dir: str) -> List[Tuple[int, str]]:
    """Return sorted [(sample_index, whole_relaxed_pdb_path), ...] for a pair."""
    out = []
    if not os.path.isdir(pair_dir):
        return out
    for fname in os.listdir(pair_dir):
        m = SAMPLE_WHOLE_RELAXED_RE.match(fname)
        if not m:
            continue
        out.append((int(m.group(1)), os.path.join(pair_dir, fname)))
    out.sort(key=lambda t: t[0])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--staging-folder", required=True)
    ap.add_argument("--structures-folder", required=True)
    ap.add_argument("--sequences-folder", required=True)
    ap.add_argument("--structures-map", required=True)
    ap.add_argument("--sequences-csv", required=True)
    ap.add_argument("--missing-csv", required=True)
    args = ap.parse_args()

    os.makedirs(args.structures_folder, exist_ok=True)
    os.makedirs(args.sequences_folder, exist_ok=True)
    os.makedirs(os.path.dirname(args.structures_map), exist_ok=True)
    os.makedirs(os.path.dirname(args.sequences_csv), exist_ok=True)
    os.makedirs(os.path.dirname(args.missing_csv), exist_ok=True)

    map_rows: List[Dict[str, str]] = []
    seq_rows: List[Dict[str, str]] = []
    missing_rows: List[Dict[str, str]] = []
    step_id = step_id_from_table_path(args.missing_csv)

    if not os.path.isdir(args.staging_folder):
        print(f"ERROR: staging folder missing: {args.staging_folder}", file=sys.stderr)
        sys.exit(1)

    pair_names = sorted(
        d for d in os.listdir(args.staging_folder)
        if os.path.isdir(os.path.join(args.staging_folder, d))
    )

    for pair_id in pair_names:
        pair_dir = os.path.join(args.staging_folder, pair_id)
        samples = collect_samples(pair_dir)
        if not samples:
            missing_rows.append({
                "id": pair_id,
                "removed_by": step_id,
                "kind": "failure",
                "cause": "no <i>_whole_relaxed.pdb produced",
            })
            continue

        prot_id, lig_id = parse_pair_id(pair_id)
        for sample_idx, src_pdb in samples:
            out_id = f"{pair_id}_{sample_idx}"
            dst_pdb = os.path.join(args.structures_folder, f"{out_id}.pdb")
            shutil.copyfile(src_pdb, dst_pdb)

            seq = parse_sequence_from_pdb(dst_pdb)

            map_rows.append({
                "id": out_id,
                "file": dst_pdb,
                "structures.id": prot_id,
                "compounds.id": lig_id,
            })
            seq_rows.append({
                "id": out_id,
                "structures.id": prot_id,
                "compounds.id": lig_id,
                "sample": str(sample_idx),
                "sequence": seq,
            })

    with open(args.structures_map, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["id", "file", "structures.id", "compounds.id"])
        w.writeheader()
        w.writerows(map_rows)

    with open(args.sequences_csv, "w", newline="") as f:
        w = csv.DictWriter(
            f, fieldnames=["id", "structures.id", "compounds.id", "sample", "sequence"]
        )
        w.writeheader()
        w.writerows(seq_rows)

    with open(args.missing_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["id", "removed_by", "kind", "cause"])
        w.writeheader()
        w.writerows(missing_rows)

    print(
        f"PocketGen post-process: {len(map_rows)} sample(s) written from "
        f"{len(pair_names) - len(missing_rows)}/{len(pair_names)} pair(s); "
        f"{len(missing_rows)} pair(s) failed"
    )

    if not map_rows:
        print("ERROR: PocketGen produced no usable outputs", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
