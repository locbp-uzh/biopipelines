#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Distance-based residue selection referenced to a SUBSET of atoms within a
named ligand residue.

Use case: a chimeric ligand (dye + glutathione conjugate) sits in a single
ligand residue; you want protein residues within X Å of only the glutathione
atoms, not the dye atoms.

Usage:
    python pipe_ligand_atom_selector.py <structures_json> <ligand> <atoms> <distance> <restrict_spec> <output_csv> <include_reference>

Arguments:
    structures_json: JSON file containing DataStream with ids and files
    ligand: Ligand residue name (e.g. "LIG")
    atoms: ``+``-joined list of atom names within the ligand residue
           (e.g. "C61+C62+S57+O49")
    distance: Distance cutoff in Angstroms
    restrict_spec: Restriction specification (table reference, direct selection, or "")
    output_csv: Path to output CSV file
    include_reference: 'true' or 'false' (no effect; ligand is not a protein residue)

Output:
    CSV file with columns: id, pdb, within, beyond, distance_cutoff, reference_ligand
"""

import sys
import os
import math
import pandas as pd
from typing import List, Dict, Tuple

# Reuse helpers from pipe_distance_selector and biopipelines
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from biopipelines.pdb_parser import parse_pdb_file, select_atoms_by_ligand, Atom, STANDARD_RESIDUES
from biopipelines.sele_utils import sele_to_list, list_to_sele
from biopipelines.biopipelines_io import load_datastream, iterate_files, load_table, lookup_table_value


def get_ligand_atom_subset(atoms: List[Atom], ligand_name: str, atom_names: List[str]) -> List[Atom]:
    """
    Select atoms in the named ligand residue whose atom name is in atom_names.

    Args:
        atoms: All atoms parsed from PDB
        ligand_name: Ligand residue name (e.g. "LIG")
        atom_names: List of atom names to keep (e.g. ["C61", "S57", "O49"])

    Returns:
        List of Atom objects matching ligand + atom-name filter
    """
    ligand_atoms = select_atoms_by_ligand(atoms, ligand_name)
    if not ligand_atoms:
        residue_names = set(atom.res_name for atom in atoms)
        raise ValueError(
            f"Could not find ligand '{ligand_name}' in structure. "
            f"Available residue names: {residue_names}"
        )

    wanted = set(atom_names)
    subset = [a for a in ligand_atoms if a.atom_name in wanted]

    if not subset:
        available = sorted(set(a.atom_name for a in ligand_atoms))
        raise ValueError(
            f"None of the requested atoms {sorted(wanted)} found in ligand "
            f"'{ligand_name}'. Available atoms in this ligand: {available}"
        )

    missing = wanted - set(a.atom_name for a in subset)
    if missing:
        # Warn but don't fail — caller may have given a superset (e.g., GSH atom
        # list that includes atoms only present in some standardization variants)
        print(f"  WARNING: {len(missing)} requested atom(s) not found in {ligand_name}: {sorted(missing)}")

    return subset


def get_protein_residues(atoms: List[Atom]) -> Dict[Tuple[str, int], List[Atom]]:
    protein_residues: Dict[Tuple[str, int], List[Atom]] = {}
    for atom in atoms:
        if atom.res_name in STANDARD_RESIDUES:
            res_key = (atom.chain, atom.res_num)
            protein_residues.setdefault(res_key, []).append(atom)
    if not protein_residues:
        raise ValueError("No protein residues found")
    return protein_residues


def is_placeholder_atom(atom: Atom) -> bool:
    return atom.x == 0.0 and atom.y == 0.0 and atom.z == 0.0


def calculate_distance(atom1: Atom, atom2: Atom) -> float:
    dx = atom1.x - atom2.x
    dy = atom1.y - atom2.y
    dz = atom1.z - atom2.z
    return math.sqrt(dx * dx + dy * dy + dz * dz)


def resolve_restriction_spec(restrict_spec: str, structure_id: str) -> List[Tuple[str, int]]:
    if not restrict_spec or restrict_spec == "":
        return []
    if restrict_spec.startswith("TABLE_REFERENCE:"):
        try:
            table, column_name = load_table(restrict_spec)
        except FileNotFoundError as e:
            print(f"Warning: Restriction table not found: {e}")
            return []
        try:
            selection_str = lookup_table_value(table, structure_id, column_name)
            residues = sele_to_list(selection_str)
            print(f"Restricting to {len(residues)} residues: {list_to_sele(residues)}")
            return residues
        except KeyError as e:
            print(f"ERROR: No restriction table entry found for ID '{structure_id}'. {e}")
            return []
    else:
        return sele_to_list(restrict_spec)


def calculate_residue_distances(
    atoms: List[Atom],
    reference_atoms: List[Atom],
    distance_cutoff: float,
    restrict_to_residues: List[Tuple[str, int]] = None,
) -> Tuple[List[str], List[str]]:
    protein_residues = get_protein_residues(atoms)

    if restrict_to_residues:
        restrict_chained = set((c, r) for c, r in restrict_to_residues if c)
        restrict_chainless = set(r for c, r in restrict_to_residues if not c)
    else:
        restrict_chained = set()
        restrict_chainless = set()

    within_residues: List[str] = []
    beyond_residues: List[str] = []

    for res_key, residue_atoms in protein_residues.items():
        chain, res_num = res_key

        if restrict_to_residues:
            if (chain, res_num) not in restrict_chained and res_num not in restrict_chainless:
                continue

        min_distance = float("inf")
        for res_atom in residue_atoms:
            if is_placeholder_atom(res_atom):
                continue
            for ref_atom in reference_atoms:
                if is_placeholder_atom(ref_atom):
                    continue
                d = calculate_distance(res_atom, ref_atom)
                if d < min_distance:
                    min_distance = d

        chain_id = chain if chain else "A"
        residue_id = f"{chain_id}{res_num}"

        if min_distance <= distance_cutoff:
            within_residues.append(residue_id)
        else:
            beyond_residues.append(residue_id)

    return within_residues, beyond_residues


def format_pymol_selection(residue_list: List[str]) -> str:
    if not residue_list:
        return ""

    chain_residues: Dict[str, List[int]] = {}
    for res_id in residue_list:
        chain = res_id[0]
        res_num = int(res_id[1:])
        chain_residues.setdefault(chain, []).append(res_num)

    range_parts: List[str] = []
    for chain, res_nums in chain_residues.items():
        res_nums.sort()
        ranges: List[str] = []
        start = res_nums[0]
        end = res_nums[0]

        for i in range(1, len(res_nums)):
            if res_nums[i] == end + 1:
                end = res_nums[i]
            else:
                if start == end:
                    ranges.append(f"{chain}{start}")
                else:
                    ranges.append(f"{chain}{start}-{end}")
                start = end = res_nums[i]

        if start == end:
            ranges.append(f"{chain}{start}")
        else:
            ranges.append(f"{chain}{start}-{end}")

        range_parts.extend(ranges)

    return "+".join(range_parts)


def analyze_structure(
    structure_id: str,
    pdb_file: str,
    ligand_name: str,
    atom_names: List[str],
    distance_cutoff: float,
    restrict_spec: str = "",
) -> Dict[str, any]:
    try:
        atoms = parse_pdb_file(pdb_file)
    except Exception as e:
        raise ValueError(f"Could not load PDB file {pdb_file}: {e}")
    if not atoms:
        raise ValueError(f"No atoms found in PDB file: {pdb_file}")

    reference_atoms = get_ligand_atom_subset(atoms, ligand_name, atom_names)
    print(f"  Reference: {ligand_name} atoms {sorted(a.atom_name for a in reference_atoms)} ({len(reference_atoms)} atoms)")

    restrict_to_residues = resolve_restriction_spec(restrict_spec, structure_id)
    if restrict_to_residues:
        print(f"  Restricting to {len(restrict_to_residues)} residues: {list_to_sele(restrict_to_residues)}")

    within_residues, beyond_residues = calculate_residue_distances(
        atoms, reference_atoms, distance_cutoff, restrict_to_residues
    )

    within_selection = format_pymol_selection(within_residues)
    beyond_selection = format_pymol_selection(beyond_residues)

    return {
        "id": structure_id,
        "pdb": pdb_file,
        "within": within_selection,
        "beyond": beyond_selection,
        "distance_cutoff": distance_cutoff,
        "reference_ligand": f"{ligand_name}:{'+'.join(atom_names)}",
    }


def main():
    if len(sys.argv) != 8:
        print("Usage: python pipe_ligand_atom_selector.py "
              "<structures_json> <ligand> <atoms> <distance> "
              "<restrict_spec> <output_csv> <include_reference>")
        sys.exit(1)

    structures_json = sys.argv[1]
    ligand_name = sys.argv[2]
    atom_names_str = sys.argv[3]
    distance_cutoff = float(sys.argv[4])
    restrict_spec = sys.argv[5]
    output_csv = sys.argv[6]
    # include_reference is accepted for symmetry with DistanceSelector but is
    # meaningless here (the ligand isn't a protein residue).
    _ = sys.argv[7]

    atom_names = [a for a in atom_names_str.split("+") if a]
    if not atom_names:
        raise ValueError(f"No atom names parsed from '{atom_names_str}'")

    structures_ds = load_datastream(structures_json)
    print(f"Analyzing {len(structures_ds.ids_expanded)} structures with distance cutoff {distance_cutoff}Å")
    print(f"Ligand: {ligand_name}")
    print(f"Reference atoms ({len(atom_names)}): {'+'.join(atom_names)}")
    if restrict_spec:
        print(f"Restriction: {restrict_spec}")

    results: List[Dict[str, any]] = []
    for structure_id, pdb_file in iterate_files(structures_ds):
        if not os.path.exists(pdb_file):
            print(f"Warning: PDB file not found: {pdb_file}")
            continue

        try:
            print(f"\nAnalyzing: {structure_id} ({os.path.basename(pdb_file)})")
            result = analyze_structure(
                structure_id, pdb_file, ligand_name, atom_names,
                distance_cutoff, restrict_spec
            )
            results.append(result)
            within_n = len(result["within"].split("+")) if result["within"] else 0
            beyond_n = len(result["beyond"].split("+")) if result["beyond"] else 0
            print(f"  Within {distance_cutoff}Å: {within_n} segments")
            print(f"  Beyond {distance_cutoff}Å: {beyond_n} segments")
        except Exception as e:
            print(f"Error analyzing {structure_id}: {e}")
            continue

    if not results:
        raise ValueError("No structures could be analyzed successfully")

    df = pd.DataFrame(results)
    os.makedirs(os.path.dirname(output_csv), exist_ok=True)
    df.to_csv(output_csv, index=False)

    print(f"\nDistance analysis completed!")
    print(f"Results saved to: {output_csv}")
    print(f"Analyzed {len(results)} structures successfully")


if __name__ == "__main__":
    main()
