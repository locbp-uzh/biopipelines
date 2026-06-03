#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
NeuralPLexer post-processing.

For each pair folder under --raw-out, NeuralPLexer (with --separate-pdb +
--rank-outputs-by-confidence) writes:
    {pair_id}/prot_0.pdb, prot_1.pdb, ...
    {pair_id}/lig_0.sdf, lig_1.sdf, ...
    {pair_id}/prot_all.pdb, lig_all.sdf
    {pair_id}/confidence.csv         (when --rank-outputs-by-confidence is set)

We merge prot_<k>.pdb + lig_<k>.sdf into a single PDB per rank, copy to
<structures_folder>/<pair_id>_rank<N>.pdb, and parse the confidence CSV
into our standard table.
"""

import argparse
import csv
import os
import re
import sys
from typing import Dict, List

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from biopipelines.biopipelines_io import step_id_from_table_path  # noqa: E402


def merge_protein_ligand_pdb(prot_pdb: str, lig_sdf: str, out_pdb: str) -> None:
    """Concatenate the protein ATOM records and the ligand atoms (converted
    from SDF to HETATM PDB records via RDKit) into a single PDB file.
    """
    from rdkit import Chem

    with open(out_pdb, "w") as out:
        with open(prot_pdb) as f:
            for line in f:
                if line.startswith(("ATOM", "TER")):
                    out.write(line)

        suppl = Chem.SDMolSupplier(lig_sdf, removeHs=False)
        mol = next(iter(suppl), None)
        if mol is not None:
            # RDKit's MolToPDBBlock yields HETATM lines for non-polymer mols.
            block = Chem.MolToPDBBlock(mol)
            for line in block.splitlines():
                if line.startswith(("HETATM", "CONECT")):
                    out.write(line + "\n")
        out.write("END\n")


def parse_confidence_csv(path: str) -> Dict[int, float]:
    """Read NeuralPLexer's per-rank confidence CSV if present. Returns
    {rank_index_0based: confidence}. Tolerant of the exact column layout —
    we look for an `index`/`rank` column plus a `confidence`/`score` column.
    """
    if not os.path.isfile(path):
        return {}
    out: Dict[int, float] = {}
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        rank_col = next(
            (c for c in (reader.fieldnames or []) if c.lower() in ("index", "rank", "idx")),
            None,
        )
        conf_col = next(
            (c for c in (reader.fieldnames or []) if c.lower() in ("confidence", "score", "conf")),
            None,
        )
        if rank_col is None or conf_col is None:
            return {}
        for row in reader:
            try:
                out[int(row[rank_col])] = float(row[conf_col])
            except (ValueError, KeyError):
                continue
    return out


def parse_pair_id(name: str):
    if "+" not in name:
        raise ValueError(f"folder {name!r} is not a '<prot>+<lig>' pair")
    return name.split("+", 1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw-out", required=True)
    ap.add_argument("--structures-folder", required=True)
    ap.add_argument("--structures-map", required=True)
    ap.add_argument("--confidence-csv", required=True)
    ap.add_argument("--missing-csv", required=True)
    ap.add_argument("--n-samples", type=int, required=True)
    args = ap.parse_args()

    os.makedirs(args.structures_folder, exist_ok=True)
    os.makedirs(os.path.dirname(args.structures_map), exist_ok=True)
    os.makedirs(os.path.dirname(args.confidence_csv), exist_ok=True)
    os.makedirs(os.path.dirname(args.missing_csv), exist_ok=True)

    map_rows: List[Dict[str, str]] = []
    conf_rows: List[Dict[str, str]] = []
    missing_rows: List[Dict[str, str]] = []
    step_id = step_id_from_table_path(args.missing_csv)

    if not os.path.isdir(args.raw_out):
        print(f"ERROR: raw output dir missing: {args.raw_out}", file=sys.stderr)
        sys.exit(1)

    pair_names = sorted(
        d for d in os.listdir(args.raw_out)
        if os.path.isdir(os.path.join(args.raw_out, d))
    )

    for pair_id in pair_names:
        pair_dir = os.path.join(args.raw_out, pair_id)
        prot_id, lig_id = parse_pair_id(pair_id)

        # Discover per-rank files. With --rank-outputs-by-confidence (which the
        # wrapper passes) upstream writes prot_rank<N>_plddt<score>.pdb +
        # lig_rank<N>_plddt<score>.sdf, with rank/confidence in the filename.
        # Fall back to the un-ranked prot_<k>.pdb / lig_<k>.sdf naming otherwise.
        ranked = []
        for f in os.listdir(pair_dir):
            m = re.match(r"^prot_rank(\d+)_plddt([0-9.]+)\.pdb$", f)
            if m:
                rank = int(m.group(1))
                score = float(m.group(2))
                lig_sdf = os.path.join(pair_dir, f"lig_rank{rank}_plddt{m.group(2)}.sdf")
                ranked.append((rank, score, os.path.join(pair_dir, f), lig_sdf))
        if ranked:
            # (rank, confidence, prot_pdb, lig_sdf), already 1-indexed by upstream.
            prot_entries = sorted(ranked)
        else:
            # Legacy un-ranked naming: prot_<k>.pdb / lig_<k>.sdf (k 0-based),
            # with confidence from the separate confidence.csv if present.
            confidence_lookup = parse_confidence_csv(os.path.join(pair_dir, "confidence.csv"))
            prot_entries = sorted(
                (int(m.group(1)) + 1, confidence_lookup.get(int(m.group(1))),
                 os.path.join(pair_dir, f), os.path.join(pair_dir, f"lig_{m.group(1)}.sdf"))
                for f in os.listdir(pair_dir)
                for m in [re.match(r"^prot_(\d+)\.pdb$", f)] if m
            )

        if not prot_entries:
            missing_rows.append({
                "id": pair_id,
                "removed_by": step_id,
                "kind": "failure",
                "cause": "no prot_*.pdb produced",
            })
            continue
        if len(prot_entries) < args.n_samples:
            missing_rows.append({
                "id": pair_id,
                "removed_by": step_id,
                "kind": "failure",
                "cause": f"only {len(prot_entries)}/{args.n_samples} samples produced",
            })

        for rank, conf_val, prot_pdb, lig_sdf in prot_entries:
            out_id = f"{pair_id}_rank{rank}"
            dst = os.path.join(args.structures_folder, f"{out_id}.pdb")
            try:
                merge_protein_ligand_pdb(prot_pdb, lig_sdf, dst)
            except Exception as exc:
                # Don't silently emit a protein-only PDB labelled as a
                # successful complex — flag the rank as missing instead.
                cause = f"merge failed: {type(exc).__name__}: {exc}".replace(",", ";").replace("\n", " ")
                print(f"  [fail] {out_id}: {cause}", file=sys.stderr)
                missing_rows.append({
                    "id": out_id,
                    "removed_by": step_id,
                    "kind": "failure",
                    "cause": cause,
                })
                continue

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
                "confidence": "" if conf_val is None else f"{conf_val:.4f}",
            })

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
        f"NeuralPLexer post-process: {len(map_rows)} complexes kept, "
        f"{len(missing_rows)} pair(s) failed/partial"
    )

    if not map_rows:
        print("ERROR: NeuralPLexer produced no usable outputs", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
