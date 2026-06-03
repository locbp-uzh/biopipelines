#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Distance-based residue selection.

Computes per-residue distance from each protein residue to a reference (any
selection expression accepted by ``pdb_parser.resolve_selection``), then
partitions residues into within/beyond by distance cutoff, top-K cap, or both.

Outputs:
  - selections CSV  (id, pdb, within, beyond, distance_cutoff, top_k, mode, reference)
  - per-id resi-csv files (id, chain, resi, distance) under <distances-dir>/<id>.csv
  - distances map CSV (id, file)
"""

import argparse
import math
import os
import re
import sys
from typing import Dict, List, Optional, Tuple

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from biopipelines.pdb_parser import (
    Atom,
    STANDARD_RESIDUES,
    parse_pdb_file,
    resolve_selection,
    select_atoms_by_ligand,
    select_atoms_by_residue_number,
)
from biopipelines.sele_utils import sele_to_list, chain_aware_sele
from biopipelines.biopipelines_io import (
    load_datastream,
    iterate_files,
    load_table,
    lookup_table_value,
)


BACKBONE_ATOMS = {"N", "CA", "C", "O"}


def is_placeholder_atom(atom: Atom) -> bool:
    return atom.x == 0.0 and atom.y == 0.0 and atom.z == 0.0


def filter_by_atom_class(atoms: List[Atom], atom_class: str) -> List[Atom]:
    """Apply the protein/reference atom-class filter."""
    if atom_class == "all":
        return atoms
    if atom_class == "CA":
        return [a for a in atoms if a.atom_name.strip() == "CA"]
    if atom_class == "backbone":
        return [a for a in atoms if a.atom_name.strip() in BACKBONE_ATOMS]
    if atom_class == "sidechain":
        return [a for a in atoms if a.atom_name.strip() not in BACKBONE_ATOMS]
    raise ValueError(f"Unknown atom class: {atom_class!r}")


def get_protein_residues(atoms: List[Atom]) -> Dict[Tuple[str, int], List[Atom]]:
    by_residue: Dict[Tuple[str, int], List[Atom]] = {}
    for atom in atoms:
        if atom.res_name not in STANDARD_RESIDUES:
            continue
        key = (atom.chain, atom.res_num)
        by_residue.setdefault(key, []).append(atom)
    if not by_residue:
        raise ValueError("No protein residues found")
    return by_residue


def parse_restrict_spec(restrict_spec: str, structure_id: str, available_chains: set) -> Tuple[Optional[set], Optional[set], Optional[set]]:
    """Resolve restrict-to spec into (chained_residues, chainless_residues, chains).

    Returns three optional sets — chained ``{(chain, resnum)}``, chainless
    ``{resnum}``, and chain-only ``{chain}`` — describing the restriction.
    A return of ``(None, None, None)`` means no restriction.
    """
    if not restrict_spec:
        return None, None, None

    if restrict_spec.startswith("TABLE_REFERENCE:"):
        try:
            table, column = load_table(restrict_spec)
        except FileNotFoundError as e:
            print(f"Warning: Restriction table not found: {e}", file=sys.stderr)
            return None, None, None
        try:
            value = lookup_table_value(table, structure_id, column)
        except KeyError as e:
            print(f"ERROR: No restriction-table entry for ID '{structure_id}': {e}", file=sys.stderr)
            return set(), set(), set()  # empty restriction = nothing passes
        residues = sele_to_list(value)
        chained = {(c, r) for c, r in residues if c}
        chainless = {r for c, r in residues if not c}
        return chained or None, chainless or None, None

    # Bare chain spec: "B", "chain B"
    m = re.match(r"^\s*(?:chain\s+)?([A-Za-z])\s*$", restrict_spec)
    if m and m.group(1).upper() in {c.upper() for c in available_chains if c}:
        return None, None, {m.group(1)}

    # Residue selection
    residues = sele_to_list(restrict_spec)
    chained = {(c, r) for c, r in residues if c}
    chainless = {r for c, r in residues if not c}
    if not chained and not chainless:
        return set(), set(), set()
    return chained or None, chainless or None, None


def passes_restriction(chain: str, resnum: int,
                       chained: Optional[set], chainless: Optional[set], chains: Optional[set]) -> bool:
    if chained is None and chainless is None and chains is None:
        return True
    if chains is not None and chain in chains:
        return True
    if chained is not None and (chain, resnum) in chained:
        return True
    if chainless is not None and resnum in chainless:
        return True
    return False


def centroid(atoms: List[Atom]) -> Tuple[float, float, float]:
    n = len(atoms)
    if n == 0:
        raise ValueError("Cannot take centroid of empty atom list")
    sx = sum(a.x for a in atoms) / n
    sy = sum(a.y for a in atoms) / n
    sz = sum(a.z for a in atoms) / n
    return sx, sy, sz


def min_distance_atoms(group: List[Atom], ref_atoms: List[Atom]) -> float:
    best = math.inf
    for a in group:
        if is_placeholder_atom(a):
            continue
        for b in ref_atoms:
            if is_placeholder_atom(b):
                continue
            dx = a.x - b.x
            dy = a.y - b.y
            dz = a.z - b.z
            d = math.sqrt(dx * dx + dy * dy + dz * dz)
            if d < best:
                best = d
    return best


def min_distance_to_centroid(group: List[Atom], ref_centroid: Tuple[float, float, float]) -> float:
    best = math.inf
    cx, cy, cz = ref_centroid
    for a in group:
        if is_placeholder_atom(a):
            continue
        dx = a.x - cx
        dy = a.y - cy
        dz = a.z - cz
        d = math.sqrt(dx * dx + dy * dy + dz * dz)
        if d < best:
            best = d
    return best


def resolve_reference_atoms(reference: str, atoms: List[Atom]) -> List[Atom]:
    """Resolve a reference string to a list of atoms.

    Layered resolution:
      1. ``sele_to_list`` first — handles bare residue numbers, chain-prefixed
         ranges (``A300-310``), and multi-residue selections (``A370+B372``).
         Hits the protein-residue path with full chain awareness.
      2. ``resolve_selection`` for everything else — atom-level
         (``LIG.O5``, ``87.CA``, ``A141.CB``), sequence context, keywords.
      3. ``select_atoms_by_ligand`` as a bare-name fallback for ligand codes
         (``GIV``, ``STI``, ``LIG``) that ``resolve_selection``'s atom-name
         fallback wouldn't recognise.
    """
    # Only try sele_to_list when the string looks like a residue selection
    # (digits, optional chain letters, '+' / '-' separators). Skip when it
    # contains '.' (atom-level), ' in ' (sequence context), or any other
    # token sele_to_list would warn about.
    looks_like_residue = bool(re.fullmatch(r"[A-Za-z]?-?\d+(?:[-+][A-Za-z]?-?\d+)*", reference.strip()))
    parsed = sele_to_list(reference) if looks_like_residue else []
    if parsed:
        # Group by chain so a multi-chain selection like A300-310+B400 is
        # honoured. Empty chain ('') is the chain-agnostic case.
        out: List[Atom] = []
        seen_atom_ids = set()
        by_chain: Dict[str, List[int]] = {}
        for chain, resnum in parsed:
            by_chain.setdefault(chain, []).append(resnum)
        for chain, nums in by_chain.items():
            picked = select_atoms_by_residue_number(atoms, nums, chain=chain or None)
            for a in picked:
                key = id(a)
                if key not in seen_atom_ids:
                    seen_atom_ids.add(key)
                    out.append(a)
        if out:
            return out

    selected = resolve_selection(reference, atoms)
    if selected:
        return selected

    return select_atoms_by_ligand(atoms, reference)


def reference_is_protein(ref_atoms: List[Atom]) -> bool:
    """True if any reference atom belongs to a standard protein residue."""
    return any(a.res_name in STANDARD_RESIDUES for a in ref_atoms)


def reference_residue_keys(reference: str, ref_atoms: List[Atom]) -> set:
    """Set of (chain, resnum) for any protein residues the reference covers.

    Used to exclude the reference itself from `within` when
    `include_reference=False`. Ligand/HETATM references contribute no protein
    residues and naturally yield an empty set.
    """
    keys = set()
    for a in ref_atoms:
        if a.res_name in STANDARD_RESIDUES:
            keys.add((a.chain, a.res_num))
    return keys


def analyze_structure(structure_id: str,
                      pdb_file: str,
                      reference: str,
                      distance_cutoff: Optional[float],
                      top_k: Optional[int],
                      mode: str,
                      atom_class: str,
                      restrict_spec: str,
                      include_reference: bool,
                      distances_dir: str) -> Tuple[Dict[str, str], str]:
    """Analyze one structure. Returns (selections_row, distances_file_path)."""
    atoms = parse_pdb_file(pdb_file)
    if not atoms:
        raise ValueError(f"No atoms in {pdb_file}")

    ref_atoms_raw = resolve_reference_atoms(reference, atoms)
    if not ref_atoms_raw:
        raise ValueError(f"Reference '{reference}' selected no atoms in {pdb_file}")

    # atoms= is a protein-side concept (backbone / CA / sidechain only make
    # sense for standard residues). Apply it to the reference only when the
    # reference itself covers protein residues; for ligand / HETATM
    # references the full atom set is kept.
    if reference_is_protein(ref_atoms_raw):
        ref_atoms = filter_by_atom_class(ref_atoms_raw, atom_class)
        if not ref_atoms:
            raise ValueError(
                f"Reference '{reference}' has no atoms left after applying atoms='{atom_class}'"
            )
    else:
        ref_atoms = ref_atoms_raw

    ref_keys = reference_residue_keys(reference, ref_atoms_raw)

    available_chains = {a.chain for a in atoms if a.res_name in STANDARD_RESIDUES}
    r_chained, r_chainless, r_chains = parse_restrict_spec(
        restrict_spec, structure_id, available_chains
    )

    protein_residues = get_protein_residues(atoms)

    ref_cent = centroid(ref_atoms) if mode == "centroid" else None

    per_residue_rows = []  # for the resi-csv stream
    candidates: List[Tuple[float, str, int]] = []  # (distance, chain, resnum) passing restriction

    for (chain, resnum), residue_atoms in protein_residues.items():
        if not passes_restriction(chain, resnum, r_chained, r_chainless, r_chains):
            continue

        group = filter_by_atom_class(residue_atoms, atom_class)
        if not group:
            continue

        if mode == "centroid":
            d = min_distance_to_centroid(group, ref_cent)
        else:
            d = min_distance_atoms(group, ref_atoms)

        if not math.isfinite(d):
            continue

        per_residue_rows.append({
            "id": structure_id,
            "chain": chain if chain else "A",
            "resi": resnum,
            "distance": round(d, 4),
        })
        candidates.append((d, chain, resnum))

    # Apply cutoff and/or top-K
    if distance_cutoff is not None:
        passing = [c for c in candidates if c[0] <= distance_cutoff]
    else:
        passing = list(candidates)

    passing.sort(key=lambda x: x[0])
    if top_k is not None:
        passing = passing[:top_k]

    passing_set = {(c, r) for _, c, r in passing}

    within_residues = []
    beyond_residues = []
    for d, chain, resnum in candidates:
        chain_id = chain if chain else "A"
        if (chain, resnum) in passing_set:
            if (not include_reference) and (chain, resnum) in ref_keys:
                beyond_residues.append((chain_id, resnum))
            else:
                within_residues.append((chain_id, resnum))
        else:
            beyond_residues.append((chain_id, resnum))

    within_sel = chain_aware_sele(within_residues)
    beyond_sel = chain_aware_sele(beyond_residues)

    # Write per-id resi-csv
    os.makedirs(distances_dir, exist_ok=True)
    distances_file = os.path.join(distances_dir, f"{structure_id}.csv")
    pd.DataFrame(per_residue_rows, columns=["id", "chain", "resi", "distance"]).to_csv(
        distances_file, index=False
    )

    row = {
        "id": structure_id,
        "pdb": pdb_file,
        "within": within_sel,
        "beyond": beyond_sel,
        "distance_cutoff": "" if distance_cutoff is None else distance_cutoff,
        "top_k": "" if top_k is None else top_k,
        "mode": mode,
        "reference": reference,
    }
    return row, distances_file


def main():
    p = argparse.ArgumentParser(description="DistanceSelector pipe script")
    p.add_argument("--structures-json", required=True)
    p.add_argument("--reference", required=True)
    p.add_argument("--distance", default="")
    p.add_argument("--top-k", default="")
    p.add_argument("--mode", default="min", choices=["min", "centroid"])
    p.add_argument("--atoms", default="all", choices=["all", "backbone", "CA", "sidechain"])
    p.add_argument("--restrict-to", default="")
    p.add_argument("--include-reference", default="true")
    p.add_argument("--selections-csv", required=True)
    p.add_argument("--distances-dir", required=True)
    p.add_argument("--distances-map-csv", required=True)
    args = p.parse_args()

    distance_cutoff = float(args.distance) if args.distance not in ("", "None") else None
    top_k = int(args.top_k) if args.top_k not in ("", "None") else None
    include_reference = args.include_reference.lower() in ("true", "1", "yes")

    if distance_cutoff is None and top_k is None:
        print("ERROR: at least one of --distance or --top-k must be set", file=sys.stderr)
        sys.exit(1)

    ds = load_datastream(args.structures_json)

    print(f"Analyzing {len(ds.ids_expanded)} structures")
    print(f"Reference: {args.reference}")
    print(f"Distance cutoff: {distance_cutoff}")
    print(f"Top-K: {top_k}")
    print(f"Mode: {args.mode}, Atoms: {args.atoms}")
    if args.restrict_to:
        print(f"Restriction: {args.restrict_to}")

    sel_rows = []
    map_rows = []
    failed = []

    for structure_id, pdb_file in iterate_files(ds):
        if not os.path.exists(pdb_file):
            print(f"WARNING: PDB file not found: {pdb_file}", file=sys.stderr)
            failed.append(structure_id)
            continue
        try:
            print(f"\nAnalyzing: {structure_id} ({os.path.basename(pdb_file)})")
            row, dist_file = analyze_structure(
                structure_id=structure_id,
                pdb_file=pdb_file,
                reference=args.reference,
                distance_cutoff=distance_cutoff,
                top_k=top_k,
                mode=args.mode,
                atom_class=args.atoms,
                restrict_spec=args.restrict_to,
                include_reference=include_reference,
                distances_dir=args.distances_dir,
            )
            sel_rows.append(row)
            map_rows.append({"id": structure_id, "file": dist_file})
            within_n = len(row["within"].split("+")) if row["within"] else 0
            beyond_n = len(row["beyond"].split("+")) if row["beyond"] else 0
            print(f"  within: {within_n} segments, beyond: {beyond_n} segments")
        except Exception as e:
            print(f"ERROR analyzing {structure_id}: {e}", file=sys.stderr)
            failed.append(structure_id)
            continue

    os.makedirs(os.path.dirname(args.selections_csv), exist_ok=True)
    os.makedirs(os.path.dirname(args.distances_map_csv), exist_ok=True)

    pd.DataFrame(
        sel_rows,
        columns=["id", "pdb", "within", "beyond", "distance_cutoff", "top_k", "mode", "reference"],
    ).to_csv(args.selections_csv, index=False)
    pd.DataFrame(map_rows, columns=["id", "file"]).to_csv(args.distances_map_csv, index=False)

    print(f"\nSelections: {args.selections_csv} ({len(sel_rows)} rows)")
    print(f"Distances resi-csv map: {args.distances_map_csv} ({len(map_rows)} rows)")
    if failed:
        print(f"Failed {len(failed)}/{len(failed)+len(sel_rows)}: {failed}", file=sys.stderr)
    if not sel_rows:
        sys.exit(1)


if __name__ == "__main__":
    main()
