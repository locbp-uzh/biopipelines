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
import numpy as np
from pathlib import Path

try:
    import MDAnalysis as mda
    from MDAnalysis.analysis import distances
    HAS_MDANALYSIS = True
except ImportError:
    HAS_MDANALYSIS = False


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


def get_ligand_selection(universe, ligand_name):
    """
    Get MDAnalysis selection for a ligand by name.

    Args:
        universe: MDAnalysis Universe object
        ligand_name: Name of the ligand (e.g., "LIG", "ATP")

    Returns:
        AtomGroup for the ligand
    """
    # Try different common ligand selection patterns
    selection_patterns = [
        f"resname {ligand_name}",
        f"segid {ligand_name}",
        f"resname {ligand_name.upper()}",
        f"resname {ligand_name.lower()}"
    ]

    for pattern in selection_patterns:
        try:
            ligand_atoms = universe.select_atoms(pattern)
            if len(ligand_atoms) > 0:
                print(f"Found ligand using selection: {pattern} ({len(ligand_atoms)} atoms)")
                return ligand_atoms
        except Exception as e:
            continue

    raise ValueError(f"Could not find ligand '{ligand_name}' in structure. Available residue names: {set(universe.residues.resnames)}")


def get_protein_residues(universe):
    """
    Get protein residues from the universe.

    Args:
        universe: MDAnalysis Universe object

    Returns:
        ResidueGroup for protein residues
    """
    # Standard protein residue selection
    protein_selection = "protein"

    try:
        protein_residues = universe.select_atoms(protein_selection).residues
        if len(protein_residues) == 0:
            raise ValueError("No protein residues found")
        return protein_residues
    except Exception as e:
        raise ValueError(f"Could not select protein residues: {e}")


def calculate_residue_distances(universe, reference_atoms, distance_cutoff):
    """
    Calculate distances from protein residues to reference atoms.

    Args:
        universe: MDAnalysis Universe object
        reference_atoms: AtomGroup for reference atoms
        distance_cutoff: Distance cutoff in Angstroms

    Returns:
        Tuple of (within_residues, beyond_residues) as lists of residue numbers
    """
    protein_residues = get_protein_residues(universe)

    within_residues = []
    beyond_residues = []

    for residue in protein_residues:
        # Get all atoms in this residue
        residue_atoms = residue.atoms

        # Calculate minimum distance from any atom in this residue to any reference atom
        distances_matrix = distances.distance_array(
            residue_atoms.positions,
            reference_atoms.positions,
            box=universe.dimensions
        )

        min_distance = np.min(distances_matrix)

        # Format residue identifier (chain + residue number)
        chain_id = residue.segment.segid if residue.segment.segid else "A"
        residue_id = f"{chain_id}{residue.resid}"

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


def analyze_structure_distance(pdb_file, reference_spec, distance_cutoff):
    """
    Analyze a single structure for distance-based residue selection.

    Args:
        pdb_file: Path to PDB file
        reference_spec: Reference specification string
        distance_cutoff: Distance cutoff in Angstroms

    Returns:
        Dictionary with analysis results
    """
    if not HAS_MDANALYSIS:
        raise ImportError("MDAnalysis is required for distance calculations. Install with: pip install MDAnalysis")

    # Load structure
    try:
        universe = mda.Universe(pdb_file)
    except Exception as e:
        raise ValueError(f"Could not load PDB file {pdb_file}: {e}")

    # Parse reference specification
    reference_type, selection = parse_reference_spec(reference_spec)

    # Get reference atoms based on type
    if reference_type == "ligand":
        reference_atoms = get_ligand_selection(universe, selection)
        reference_description = f"ligand {selection}"
    else:
        # For "atoms" or "residues", use the selection directly
        try:
            reference_atoms = universe.select_atoms(selection)
            if len(reference_atoms) == 0:
                raise ValueError(f"No atoms found for selection: {selection}")
            reference_description = f"{reference_type} {selection}"
        except Exception as e:
            raise ValueError(f"Invalid {reference_type} selection '{selection}': {e}")

    print(f"Using reference: {reference_description} ({len(reference_atoms)} atoms)")

    # Calculate distances
    within_residues, beyond_residues = calculate_residue_distances(
        universe, reference_atoms, distance_cutoff
    )

    # Format as PyMOL selections
    within_selection = format_pymol_selection(within_residues)
    beyond_selection = format_pymol_selection(beyond_residues)

    # Extract structure ID from filename
    structure_id = os.path.splitext(os.path.basename(pdb_file))[0]

    return {
        "id": structure_id,
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