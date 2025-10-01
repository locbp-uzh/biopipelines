#!/usr/bin/env python3
"""
Distance-based residue selection for protein structures.

This script analyzes protein structures to identify residues within or beyond
a specified distance from a reference (ligand, atoms, or residues), generating
PyMOL-formatted selections for use in downstream tools like LigandMPNN.

Usage:
    python pipe_distance_selector.py <structure_files> <reference_spec> <distance> <output_csv>

Arguments:
    structure_files: Comma-separated list of PDB file paths
    reference_spec: Reference specification in format "type:selection"
                   - "ligand:LIG" - Use ligand with name LIG
                   - "atoms:resname ATP" - Use PyMOL atom selection
                   - "residues:resi 100-105" - Use PyMOL residue selection
    distance: Distance cutoff in Angstroms
    output_csv: Path to output CSV file with selections

Output:
    CSV file with columns: id, within, beyond, distance_cutoff, reference_ligand
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

def calculate_residue_distances(atoms: List[Atom], reference_atoms: List[Atom], distance_cutoff: float) -> Tuple[List[str], List[str]]:
    """
    Calculate distances from protein residues to reference atoms using native PDB parser.

    Args:
        atoms: List of all atoms from PDB file
        reference_atoms: List of reference atoms
        distance_cutoff: Distance cutoff in Angstroms

    Returns:
        Tuple of (within_residues, beyond_residues) as lists of residue identifiers
    """
    protein_residues = get_protein_residues(atoms)

    within_residues = []
    beyond_residues = []

    for res_key, residue_atoms in protein_residues.items():
        chain, res_num = res_key

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
            within_residues.append(residue_id)
        else:
            beyond_residues.append(residue_id)

    return within_residues, beyond_residues


def format_pymol_selection(residue_list):
    """
    Format residue list as PyMOL selection string.

    Args:
        residue_list: List of residue identifiers (e.g., ["A10", "A11", "B15"])

    Returns:
        PyMOL selection string (e.g., "A10 A11 B15")
    """
    if not residue_list:
        return ""

    return " ".join(residue_list)


def analyze_structure_distance(pdb_file: str, reference_spec: str, distance_cutoff: float) -> Dict[str, any]:
    """
    Analyze a single structure for distance-based residue selection using native PDB parser.

    Args:
        pdb_file: Path to PDB file
        reference_spec: Reference specification string
        distance_cutoff: Distance cutoff in Angstroms

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
        reference_description = f"ligand {selection}"
    else:
        # For atoms/residues, we'll need to implement native selection parsing
        # For now, only ligand selection is supported
        raise ValueError(f"Reference type '{reference_type}' not yet supported with native PDB parser. Only 'ligand' is supported.")

    print(f"Using reference: {reference_description} ({len(reference_atoms)} atoms)")

    # Calculate distances
    within_residues, beyond_residues = calculate_residue_distances(
        atoms, reference_atoms, distance_cutoff
    )

    # Format as PyMOL selections
    within_selection = format_pymol_selection(within_residues)
    beyond_selection = format_pymol_selection(beyond_residues)

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
    if len(sys.argv) != 5:
        print("Usage: python pipe_distance_selector.py <structure_files> <reference_spec> <distance> <output_csv>")
        print("")
        print("Arguments:")
        print("  structure_files: Comma-separated list of PDB file paths")
        print("  reference_spec: Reference specification (e.g., 'ligand:LIG', 'atoms:resname ATP')")
        print("  distance: Distance cutoff in Angstroms")
        print("  output_csv: Path to output CSV file")
        sys.exit(1)

    structure_files_str = sys.argv[1]
    reference_spec = sys.argv[2]
    distance_cutoff = float(sys.argv[3])
    output_csv = sys.argv[4]

    # Parse comma-separated list of structure files
    structure_files = [f.strip() for f in structure_files_str.split(",") if f.strip()]
    if not structure_files:
        raise ValueError(f"No structure files provided: {structure_files_str}")

    print(f"Analyzing {len(structure_files)} structures with distance cutoff {distance_cutoff}Å")
    print(f"Reference: {reference_spec}")

    # Analyze each structure
    results = []
    for pdb_file in structure_files:
        if not os.path.exists(pdb_file):
            print(f"Warning: PDB file not found: {pdb_file}")
            continue

        try:
            print(f"\nAnalyzing: {os.path.basename(pdb_file)}")
            result = analyze_structure_distance(pdb_file, reference_spec, distance_cutoff)
            results.append(result)

            print(f"  Within {distance_cutoff}Å: {len(result['within'].split()) if result['within'] else 0} residues")
            print(f"  Beyond {distance_cutoff}Å: {len(result['beyond'].split()) if result['beyond'] else 0} residues")

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