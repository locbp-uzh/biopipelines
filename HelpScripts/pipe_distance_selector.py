#!/usr/bin/env python3
"""
Distance-based residue selection for protein structures.

This script analyzes protein structures to identify residues within or beyond
a specified distance from a reference (ligand, atoms, or residues), generating
PyMOL-formatted selections for use in downstream tools like LigandMPNN.

Usage:
    python pipe_distance_selector.py <structures_list_file> <reference_spec> <distance> <restrict_spec> <output_csv> <id_map>

Arguments:
    structures_list_file: File containing list of PDB file paths (one per line)
    reference_spec: Reference specification in format "type:selection"
                   - "ligand:LIG" - Use ligand with name LIG
                   - "atoms:resname ATP" - Use PyMOL atom selection
                   - "residues:resi 100-105" - Use PyMOL residue selection
    distance: Distance cutoff in Angstroms
    restrict_spec: Restriction specification (table reference, direct selection, or "" for no restriction)
    output_csv: Path to output CSV file with selections
    id_map: JSON string with ID mapping pattern (e.g., '{"*": "*_<N>"}') for matching structure IDs to table IDs

Output:
    CSV file with columns: id, pdb, within, beyond, distance_cutoff, reference_ligand
"""

import sys
import os
import pandas as pd
import math
from pathlib import Path

# Import our native PDB parser
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from pdb_parser import parse_pdb_file, select_atoms_by_ligand, Atom
from typing import List, Dict, Tuple

# Import unified ID mapping utilities
from id_map_utils import map_table_ids_to_ids


def parse_reference_spec(reference_spec):
    """
    Parse reference specification into type and selection.

    Args:
        reference_spec: String in format "type:selection"

    Returns:
        Tuple of (reference_type, selection)
    """
    if ":" not in reference_spec:
        raise ValueError(f"Invalid reference specification: {reference_spec}. Expected format 'type:selection'")

    reference_type, selection = reference_spec.split(":", 1)

    if reference_type not in ["ligand", "atoms", "residues"]:
        raise ValueError(f"Invalid reference type: {reference_type}. Must be 'ligand', 'atoms', or 'residues'")

    return reference_type, selection


def get_ligand_atoms(atoms: List[Atom], ligand_name: str) -> List[Atom]:
    """
    Get atoms for a ligand by name using native PDB parser.

    Args:
        atoms: List of all atoms from PDB file
        ligand_name: Name of the ligand (e.g., "LIG", "ATP")

    Returns:
        List of Atom objects for the ligand
    """
    ligand_atoms = select_atoms_by_ligand(atoms, ligand_name)

    if not ligand_atoms:
        # Get all unique residue names for error message
        residue_names = set(atom.res_name for atom in atoms)
        raise ValueError(f"Could not find ligand '{ligand_name}' in structure. Available residue names: {residue_names}")

    print(f"Found ligand '{ligand_name}': {len(ligand_atoms)} atoms")
    return ligand_atoms


def get_residue_atoms(atoms: List[Atom], residue_selection: str) -> List[Atom]:
    """
    Get atoms for specific residues using PyMOL-style selection.

    Args:
        atoms: List of all atoms from PDB file
        residue_selection: PyMOL-style residue selection (e.g., "87-100", "10+15+20-25")

    Returns:
        List of Atom objects for the selected residues
    """
    # Parse selection to get list of residue numbers
    selected_residues = set(sele_to_list(residue_selection))

    if not selected_residues:
        raise ValueError(f"Could not parse residue selection: {residue_selection}")

    # Standard amino acid names
    standard_residues = {
        'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
        'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'
    }

    # Select atoms from specified residues
    residue_atoms = []
    for atom in atoms:
        if atom.res_name in standard_residues and atom.res_num in selected_residues:
            residue_atoms.append(atom)

    if not residue_atoms:
        # Get all available residue numbers for error message
        available_residues = sorted(set(atom.res_num for atom in atoms if atom.res_name in standard_residues))
        raise ValueError(f"Could not find residues matching selection '{residue_selection}'. Available protein residues: {format_ligandmpnn_selection_from_list(available_residues)}")

    print(f"Found residues '{residue_selection}': {len(residue_atoms)} atoms")
    return residue_atoms


def get_protein_residues(atoms: List[Atom]) -> Dict[int, List[Atom]]:
    """
    Get protein residues organized by residue number.

    Args:
        atoms: List of all atoms from PDB file

    Returns:
        Dictionary mapping residue number to list of atoms in that residue
    """
    # Standard amino acid names
    standard_residues = {
        'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
        'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'
    }

    protein_residues = {}
    for atom in atoms:
        if atom.res_name in standard_residues:
            res_key = (atom.chain, atom.res_num)  # Use chain+resnum as key to handle multiple chains
            if res_key not in protein_residues:
                protein_residues[res_key] = []
            protein_residues[res_key].append(atom)

    if not protein_residues:
        raise ValueError("No protein residues found")

    return protein_residues


def calculate_distance(atom1: Atom, atom2: Atom) -> float:
    """
    Calculate Euclidean distance between two atoms.
    """
    dx = atom1.x - atom2.x
    dy = atom1.y - atom2.y
    dz = atom1.z - atom2.z
    return math.sqrt(dx*dx + dy*dy + dz*dz)


def sele_to_list(sele_str: str) -> List[int]:
    """
    Convert selection string to list of residue numbers.

    Args:
        sele_str: PyMOL-style selection (e.g., "10-20+30-40")

    Returns:
        List of residue numbers
    """
    if not sele_str or sele_str == "-" or pd.isna(sele_str):
        return []

    residues = []
    parts = str(sele_str).split('+')

    for part in parts:
        part = part.strip()
        if '-' in part and not part.startswith('-'):
            # Range like "10-20"
            start, end = map(int, part.split('-', 1))
            residues.extend(range(start, end + 1))
        elif part:
            # Single residue
            residues.append(int(part))

    return sorted(set(residues))  # Remove duplicates and sort


def format_ligandmpnn_selection_from_list(res_nums: List[int]) -> str:
    """
    Format residue number list as LigandMPNN selection string.

    Args:
        res_nums: List of residue numbers

    Returns:
        Selection string with ranges (e.g., "10-11+15")
    """
    if not res_nums:
        return ""

    # Create range strings
    res_nums = sorted(set(res_nums))
    ranges = []
    start = res_nums[0]
    end = res_nums[0]

    for i in range(1, len(res_nums)):
        if res_nums[i] == end + 1:
            end = res_nums[i]
        else:
            # Add completed range
            if start == end:
                ranges.append(f"{start}")
            else:
                ranges.append(f"{start}-{end}")
            start = end = res_nums[i]

    # Add final range
    if start == end:
        ranges.append(f"{start}")
    else:
        ranges.append(f"{start}-{end}")

    return "+".join(ranges)


# Note: map_structure_id_to_table_id is now imported from id_map_utils


def resolve_restriction_spec(restrict_spec: str, pdb_file: str, id_map: Dict[str, str] = None) -> List[int]:
    """
    Resolve restriction specification to list of residue numbers.

    Args:
        restrict_spec: Either:
                      - "" (empty): No restriction
                      - "DATASHEET_REFERENCE:path:column": Table reference
                      - "10-20+30-40": Direct PyMOL selection
        pdb_file: PDB file being analyzed
        id_map: ID mapping dictionary for matching structure IDs to table IDs

    Returns:
        List of residue numbers to restrict to (empty list = no restriction)
    """
    if not restrict_spec or restrict_spec == "":
        return []  # No restriction

    # Handle table reference
    if restrict_spec.startswith("DATASHEET_REFERENCE:"):
        _, table_path, column_name = restrict_spec.split(":", 2)

        if not os.path.exists(table_path):
            print(f"Warning: Restriction table not found: {table_path}")
            return []

        df = pd.read_csv(table_path)
        pdb_name = os.path.basename(pdb_file)
        pdb_base = os.path.splitext(pdb_name)[0]

        # Track attempted IDs for error reporting
        attempted_ids = [pdb_name]

        # Find matching row - try multiple lookup strategies in order:
        # 1. Try pdb filename match
        matching_rows = df[df['pdb'] == pdb_name]

        # 2. Apply id_map to get candidate IDs and try each in priority order
        if matching_rows.empty and id_map:
            candidate_ids = map_table_ids_to_ids(pdb_base, id_map)
            attempted_ids.extend(candidate_ids)
            for candidate_id in candidate_ids:
                matching_rows = df[df['id'] == candidate_id]
                if not matching_rows.empty:
                    if candidate_id != pdb_base:
                        print(f"Mapped structure ID '{pdb_base}' -> table ID '{candidate_id}'")
                    break

        # 3. If no id_map or no match yet, try original structure ID
        if matching_rows.empty:
            matching_rows = df[df['id'] == pdb_base]
            if pdb_base not in attempted_ids:
                attempted_ids.append(pdb_base)

        if not matching_rows.empty:
            row = matching_rows.iloc[0]
            selection_str = row.get(column_name, '')
            residues = sele_to_list(selection_str)
            print(f"Restricting to {len(residues)} residues: {format_ligandmpnn_selection_from_list(residues)}")
            return residues
        else:
            print(f"ERROR: No restriction table entry found. Tried: {', '.join(attempted_ids)}")
            print(f"ERROR: This means the restriction cannot be applied!")
            return []

    # Handle direct selection
    else:
        return sele_to_list(restrict_spec)

def calculate_residue_distances(atoms: List[Atom], reference_atoms: List[Atom], distance_cutoff: float, restrict_to_residues: List[int] = None, exclude_reference_residues: List[int] = None) -> Tuple[List[str], List[str]]:
    """
    Calculate distances from protein residues to reference atoms using native PDB parser.

    Args:
        atoms: List of all atoms from PDB file
        reference_atoms: List of reference atoms
        distance_cutoff: Distance cutoff in Angstroms
        restrict_to_residues: Optional list of residue numbers to restrict search to
        exclude_reference_residues: Optional list of reference residue numbers to exclude from "within"

    Returns:
        Tuple of (within_residues, beyond_residues) as lists of residue identifiers
    """
    protein_residues = get_protein_residues(atoms)

    within_residues = []
    beyond_residues = []

    for res_key, residue_atoms in protein_residues.items():
        chain, res_num = res_key

        # If restriction is specified, skip residues not in restriction
        if restrict_to_residues is not None and len(restrict_to_residues) > 0:
            if res_num not in restrict_to_residues:
                continue  # Skip this residue - not in restriction set

        # Calculate minimum distance from any atom in this residue to any reference atom
        min_distance = float('inf')

        for res_atom in residue_atoms:
            for ref_atom in reference_atoms:
                dist = calculate_distance(res_atom, ref_atom)
                min_distance = min(min_distance, dist)

        # Format residue identifier (chain + residue number)
        chain_id = chain if chain else "A"
        residue_id = f"{chain_id}{res_num}"

        if min_distance <= distance_cutoff:
            # Exclude reference residues if requested
            if exclude_reference_residues and res_num in exclude_reference_residues:
                beyond_residues.append(residue_id)
            else:
                within_residues.append(residue_id)
        else:
            beyond_residues.append(residue_id)

    return within_residues, beyond_residues


def format_ligandmpnn_selection(residue_list):
    """
    Format residue list as LigandMPNN selection string (numbers only).

    Args:
        residue_list: List of residue identifiers (e.g., ["A10", "A11", "B15"])

    Returns:
        LigandMPNN selection string with ranges (e.g., "10-11+15")
    """
    if not residue_list:
        return ""

    # Extract just the numbers, assuming all residues are from the same chain
    res_nums = []
    for res_id in residue_list:
        res_num = int(res_id[1:])  # Remove chain prefix
        res_nums.append(res_num)

    # Create range strings
    res_nums.sort()
    ranges = []
    start = res_nums[0]
    end = res_nums[0]

    for i in range(1, len(res_nums)):
        if res_nums[i] == end + 1:
            end = res_nums[i]
        else:
            # Add completed range
            if start == end:
                ranges.append(f"{start}")
            else:
                ranges.append(f"{start}-{end}")
            start = end = res_nums[i]

    # Add final range
    if start == end:
        ranges.append(f"{start}")
    else:
        ranges.append(f"{start}-{end}")

    return "+".join(ranges)


def format_pymol_selection(residue_list):
    """
    Format residue list as PyMOL range selection string.

    Args:
        residue_list: List of residue identifiers (e.g., ["A10", "A11", "B15"])

    Returns:
        PyMOL selection string with ranges (e.g., "A10-11+B15")
    """
    if not residue_list:
        return ""

    # Group by chain and create ranges
    chain_residues = {}
    for res_id in residue_list:
        chain = res_id[0]  # First character is chain
        res_num = int(res_id[1:])  # Rest is residue number
        if chain not in chain_residues:
            chain_residues[chain] = []
        chain_residues[chain].append(res_num)

    # Create range strings for each chain
    range_parts = []
    for chain, res_nums in chain_residues.items():
        res_nums.sort()
        ranges = []
        start = res_nums[0]
        end = res_nums[0]

        for i in range(1, len(res_nums)):
            if res_nums[i] == end + 1:
                end = res_nums[i]
            else:
                # Add completed range
                if start == end:
                    ranges.append(f"{chain}{start}")
                else:
                    ranges.append(f"{chain}{start}-{end}")
                start = end = res_nums[i]

        # Add final range
        if start == end:
            ranges.append(f"{chain}{start}")
        else:
            ranges.append(f"{chain}{start}-{end}")

        range_parts.extend(ranges)

    return "+".join(range_parts)


def analyze_structure_distance(pdb_file: str, reference_spec: str, distance_cutoff: float, restrict_spec: str = "", id_map: Dict[str, str] = None, include_reference: bool = True) -> Dict[str, any]:
    """
    Analyze a single structure for distance-based residue selection using native PDB parser.

    Args:
        pdb_file: Path to PDB file
        reference_spec: Reference specification string
        distance_cutoff: Distance cutoff in Angstroms
        restrict_spec: Restriction specification (table reference or direct selection)
        id_map: ID mapping dictionary for matching structure IDs to table IDs
        include_reference: Whether to include reference residues in "within" selection

    Returns:
        Dictionary with analysis results
    """
    # Load structure using native PDB parser
    try:
        atoms = parse_pdb_file(pdb_file)
    except Exception as e:
        raise ValueError(f"Could not load PDB file {pdb_file}: {e}")

    if not atoms:
        raise ValueError(f"No atoms found in PDB file: {pdb_file}")

    # Parse reference specification
    reference_type, selection = parse_reference_spec(reference_spec)

    # Get reference atoms based on type
    if reference_type == "ligand":
        reference_atoms = get_ligand_atoms(atoms, selection)
        reference_description = selection
        exclude_reference_residues = None  # Ligands are not protein residues
    elif reference_type == "residues":
        reference_atoms = get_residue_atoms(atoms, selection)
        reference_description = selection
        # If include_reference is False, exclude the reference residues from "within"
        if not include_reference:
            exclude_reference_residues = sele_to_list(selection)
        else:
            exclude_reference_residues = None
    else:
        # For atoms, we'll need to implement native selection parsing
        raise ValueError(f"Reference type '{reference_type}' not yet supported with native PDB parser. Supported types: 'ligand', 'residues'")

    print(f"Using reference: {reference_description} ({len(reference_atoms)} atoms)")

    # Resolve restriction with id_map
    restrict_to_residues = resolve_restriction_spec(restrict_spec, pdb_file, id_map)
    if restrict_to_residues:
        print(f"Restricting to {len(restrict_to_residues)} residues: {format_ligandmpnn_selection_from_list(restrict_to_residues)}")

    # Calculate distances with restriction
    within_residues, beyond_residues = calculate_residue_distances(
        atoms, reference_atoms, distance_cutoff, restrict_to_residues, exclude_reference_residues
    )

    # Format as universal selections (numbers only) - works for all tools
    within_selection = format_ligandmpnn_selection(within_residues)
    beyond_selection = format_ligandmpnn_selection(beyond_residues)

    # Extract structure ID from filename
    structure_id = os.path.splitext(os.path.basename(pdb_file))[0]

    return {
        "id": structure_id,
        "pdb": pdb_file,
        "within": within_selection,
        "beyond": beyond_selection,
        "distance_cutoff": distance_cutoff,
        "reference_ligand": reference_description
    }


def main():
    if len(sys.argv) != 8:
        print("Usage: python pipe_distance_selector.py <structures_list_file> <reference_spec> <distance> <restrict_spec> <output_csv> <id_map> <include_reference>")
        print("")
        print("Arguments:")
        print("  structures_list_file: File containing list of PDB file paths (one per line)")
        print("  reference_spec: Reference specification (e.g., 'ligand:LIG', 'residues:87-100')")
        print("  distance: Distance cutoff in Angstroms")
        print("  restrict_spec: Restriction specification (table reference, direct selection, or empty string)")
        print("  output_csv: Path to output CSV file")
        print("  id_map: JSON string with ID mapping pattern (e.g., '{\"*\": \"*_<N>\"}')")
        print("  include_reference: 'true' or 'false' - whether to include reference residues in 'within'")
        sys.exit(1)

    structures_list_file = sys.argv[1]
    reference_spec = sys.argv[2]
    distance_cutoff = float(sys.argv[3])
    restrict_spec = sys.argv[4]
    output_csv = sys.argv[5]
    id_map_str = sys.argv[6]
    include_reference_str = sys.argv[7]

    # Read structure files from list file (one path per line)
    with open(structures_list_file, 'r') as f:
        structure_files = [line.strip() for line in f if line.strip()]
    if not structure_files:
        raise ValueError(f"No structure files found in: {structures_list_file}")

    # Parse id_map JSON string
    import json
    try:
        id_map = json.loads(id_map_str)
    except json.JSONDecodeError as e:
        raise ValueError(f"Invalid id_map JSON: {e}")

    # Parse include_reference boolean
    include_reference = include_reference_str.lower() in ['true', '1', 'yes']

    print(f"Analyzing {len(structure_files)} structures with distance cutoff {distance_cutoff}Å")
    print(f"Reference: {reference_spec}")
    if restrict_spec:
        print(f"Restriction: {restrict_spec}")
    print(f"ID mapping: {id_map}")
    print(f"Include reference: {include_reference}")

    # Analyze each structure
    results = []
    for pdb_file in structure_files:
        if not os.path.exists(pdb_file):
            print(f"Warning: PDB file not found: {pdb_file}")
            continue

        try:
            print(f"\nAnalyzing: {os.path.basename(pdb_file)}")
            result = analyze_structure_distance(pdb_file, reference_spec, distance_cutoff, restrict_spec, id_map, include_reference)
            results.append(result)

            print(f"  Within {distance_cutoff}Å: {len(result['within'].split('+')) if result['within'] else 0} segments")
            print(f"  Beyond {distance_cutoff}Å: {len(result['beyond'].split('+')) if result['beyond'] else 0} segments")

        except Exception as e:
            print(f"Error analyzing {pdb_file}: {e}")
            continue

    if not results:
        raise ValueError("No structures could be analyzed successfully")

    # Create output DataFrame
    df = pd.DataFrame(results)

    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_csv), exist_ok=True)

    # Save to CSV
    df.to_csv(output_csv, index=False)

    print(f"\nDistance analysis completed!")
    print(f"Results saved to: {output_csv}")
    print(f"Analyzed {len(results)} structures successfully")


if __name__ == "__main__":
    main()