#!/usr/bin/env python3
"""
Runtime helper script for ProteinLigandContacts analysis.

This script analyzes protein structures to calculate minimum distances between
selected protein residues and ligands, returning contact count and normalized distance sum.
Uses native PDB file parsing instead of external libraries.
"""

import os
import sys
import argparse
import json
import pandas as pd
import numpy as np
import math
from typing import Dict, List, Any, Optional, Tuple, NamedTuple

# Import unified I/O utilities
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from pipe_biopipelines_io import load_datastream, iterate_files

# Import unified ID mapping utilities
from id_map_utils import map_table_ids_to_ids


class Atom(NamedTuple):
    """Simple atom representation."""
    atom_name: str
    residue_name: str
    residue_number: int
    chain_id: str
    x: float
    y: float
    z: float


def parse_pdb_file(pdb_path: str) -> List[Atom]:
    """
    Parse PDB file and extract atom information.

    Args:
        pdb_path: Path to PDB file

    Returns:
        List of Atom objects (excludes atoms at origin which are likely missing/placeholder atoms)
    """
    atoms = []
    skipped_origin_atoms = 0

    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                try:
                    atom_name = line[12:16].strip()
                    residue_name = line[17:20].strip()
                    chain_id = line[21:22].strip()
                    residue_number = int(line[22:26].strip())
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())

                    # Skip atoms at or very near origin (0,0,0) - these are likely missing/placeholder atoms
                    if abs(x) < 0.001 and abs(y) < 0.001 and abs(z) < 0.001:
                        skipped_origin_atoms += 1
                        continue

                    atom = Atom(
                        atom_name=atom_name,
                        residue_name=residue_name,
                        residue_number=residue_number,
                        chain_id=chain_id,
                        x=x,
                        y=y,
                        z=z
                    )
                    atoms.append(atom)

                except (ValueError, IndexError) as e:
                    print(f"Warning: Could not parse line: {line.strip()}")
                    continue

    if skipped_origin_atoms > 0:
        print(f"  - Skipped {skipped_origin_atoms} atoms at origin (0,0,0) - likely missing/placeholder atoms")

    return atoms


def parse_residue_selection(selection_str: str) -> List[int]:
    """
    Parse residue selection string into list of residue numbers.

    Args:
        selection_str: Selection string like '10-20+30-40+145'

    Returns:
        List of residue numbers
    """
    residue_numbers = []

    if not selection_str or selection_str.lower() in ['all', 'none', 'null']:
        return residue_numbers

    # Split by '+' to get individual parts
    parts = selection_str.split('+')

    for part in parts:
        part = part.strip()
        if '-' in part:
            # Range specification like '10-20'
            try:
                start, end = part.split('-')
                start_num = int(start.strip())
                end_num = int(end.strip())
                residue_numbers.extend(range(start_num, end_num + 1))
            except ValueError:
                print(f"Warning: Could not parse range: {part}")
        else:
            # Single residue number
            try:
                residue_numbers.append(int(part))
            except ValueError:
                print(f"Warning: Could not parse residue number: {part}")

    return sorted(list(set(residue_numbers)))  # Remove duplicates and sort


def get_protein_atoms(atoms: List[Atom], ligand_name: str) -> List[Atom]:
    """
    Filter atoms to get only protein atoms (excluding specified ligands).

    Args:
        atoms: List of all atoms
        ligand_name: Ligand residue name to exclude

    Returns:
        List of protein atoms
    """
    protein_atoms = []

    # Standard amino acid names
    standard_aa = {
        'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE',
        'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'
    }

    for atom in atoms:
        # Include standard amino acids, exclude specified ligand
        if atom.residue_name in standard_aa and atom.residue_name != ligand_name:
            protein_atoms.append(atom)

    return protein_atoms


def get_ligand_atoms(atoms: List[Atom], ligand_name: str) -> List[Atom]:
    """
    Filter atoms to get only ligand atoms.

    Args:
        atoms: List of all atoms
        ligand_name: Ligand residue name to include

    Returns:
        List of ligand atoms
    """
    return [atom for atom in atoms if atom.residue_name == ligand_name]


def filter_atoms_by_selection(atoms: List[Atom], residue_numbers: List[int]) -> List[Atom]:
    """
    Filter atoms to include only specified residue numbers.

    Args:
        atoms: List of atoms
        residue_numbers: List of residue numbers to include

    Returns:
        Filtered list of atoms
    """
    if not residue_numbers:  # Empty list means no filtering
        return atoms

    return [atom for atom in atoms if atom.residue_number in residue_numbers]


def calculate_distance(atom1: Atom, atom2: Atom) -> float:
    """
    Calculate Euclidean distance between two atoms.

    Args:
        atom1: First atom
        atom2: Second atom

    Returns:
        Distance in Angstroms
    """
    dx = atom1.x - atom2.x
    dy = atom1.y - atom2.y
    dz = atom1.z - atom2.z
    return math.sqrt(dx*dx + dy*dy + dz*dz)


def calculate_protein_ligand_contacts(atoms: List[Atom], protein_selections: str, ligand_name: str, contact_threshold: float) -> Tuple[Optional[int], Optional[float], Optional[float], Optional[float], Optional[float]]:
    """
    Calculate protein-ligand contact metrics for a structure.

    Args:
        atoms: List of atoms from PDB file
        protein_selections: Protein region selections (e.g., '10-20+30-40' or None for all)
        ligand_name: Ligand residue name (3-letter code)
        contact_threshold: Distance threshold for counting contacts

    Returns:
        Tuple of (contact_count, min_distance, max_distance, mean_distance, sum_distances_sqrt_normalized)
        or (None, None, None, None, None) if calculation failed
    """
    try:
        # Get protein and ligand atoms
        protein_atoms = get_protein_atoms(atoms, ligand_name)
        ligand_atoms = get_ligand_atoms(atoms, ligand_name)

        print(f"[DEBUG] - Found {len(protein_atoms)} protein atoms")
        print(f"[DEBUG] - Found {len(ligand_atoms)} ligand atoms")

        if not ligand_atoms:
            print(f"[WARNING] - No ligand atoms found for '{ligand_name}'")
            return None, None, None, None, None

        # Filter protein atoms by selection if specified
        if protein_selections and protein_selections.lower() not in ['all', 'none', 'null']:
            residue_numbers = parse_residue_selection(protein_selections)
            if residue_numbers:
                protein_atoms = filter_atoms_by_selection(protein_atoms, residue_numbers)
                print(f"[DEBUG] - After selection filtering: {len(protein_atoms)} protein atoms")

        if not protein_atoms:
            print(f"[WARNING] - No protein atoms found after selection filtering")
            return None, None, None, None, None

        # Calculate minimum distances per residue
        distances = {}
        closest_ligand_atoms = {}  # Track which ligand atom is closest to each residue
        closest_protein_atoms = {}  # Track which protein atom is closest for each residue
        closest_coords = {}  # Track coordinates for debugging

        for protein_atom in protein_atoms:
            residue_key = (protein_atom.chain_id, protein_atom.residue_number)

            if residue_key not in distances:
                distances[residue_key] = float('inf')
                closest_ligand_atoms[residue_key] = None
                closest_protein_atoms[residue_key] = None

            # Find minimum distance to any ligand atom
            for ligand_atom in ligand_atoms:
                dist = calculate_distance(protein_atom, ligand_atom)
                if dist < distances[residue_key]:
                    distances[residue_key] = dist
                    closest_ligand_atoms[residue_key] = ligand_atom.atom_name
                    closest_protein_atoms[residue_key] = protein_atom.atom_name
                    closest_coords[residue_key] = {
                        'protein': (protein_atom.x, protein_atom.y, protein_atom.z),
                        'ligand': (ligand_atom.x, ligand_atom.y, ligand_atom.z),
                        'distance': dist
                    }

        # Debug: Print first few residue details and identify which ligand atom is causing issues
        if len(distances) > 0:
            first_residues = list(distances.items())[:3]
            print(f"[DEBUG] - First 3 residues: {first_residues}")
            # Also print some sample protein and ligand atom coords
            if len(protein_atoms) > 0:
                print(f"[DEBUG] - Sample protein atom: {protein_atoms[0]}")
            if len(ligand_atoms) > 0:
                print(f"[DEBUG] - Sample ligand atom: {ligand_atoms[0]}")

            # Check which ligand atoms are closest to residues with distance 1.419
            suspicious_distance = 1.4191951944676249
            suspicious_residues = {k: v for k, v in distances.items() if abs(v - suspicious_distance) < 0.001}
            if len(suspicious_residues) > 5:
                # Get the ligand atoms responsible
                ligand_atom_counts = {}
                for res_key in list(suspicious_residues.keys())[:10]:  # Sample first 10
                    lig_atom = closest_ligand_atoms.get(res_key)
                    if lig_atom:
                        ligand_atom_counts[lig_atom] = ligand_atom_counts.get(lig_atom, 0) + 1
                print(f"[DEBUG] - {len(suspicious_residues)} residues have distance ~1.419 Å")
                print(f"[DEBUG] - Closest ligand atoms causing this: {ligand_atom_counts}")

                # Print detailed coordinates for first suspicious case
                first_suspicious = list(suspicious_residues.keys())[0]
                if first_suspicious in closest_coords:
                    coords = closest_coords[first_suspicious]
                    print(f"[DEBUG] - Detailed coords for {first_suspicious}:")
                    print(f"[DEBUG] -   Protein atom {closest_protein_atoms[first_suspicious]}: {coords['protein']}")
                    print(f"[DEBUG] -   Ligand atom {closest_ligand_atoms[first_suspicious]}: {coords['ligand']}")
                    print(f"[DEBUG] -   Calculated distance: {coords['distance']:.4f} Å")
                    # Verify the distance calculation
                    p = coords['protein']
                    l = coords['ligand']
                    manual_dist = math.sqrt((p[0]-l[0])**2 + (p[1]-l[1])**2 + (p[2]-l[2])**2)
                    print(f"[DEBUG] -   Manual verification: {manual_dist:.4f} Å")

        if not distances:
            print(f"[WARNING] - No distances calculated")
            return None, None, None, None, None

        # Count contacts below threshold
        contact_count = sum(1 for dist in distances.values() if dist < contact_threshold)

        # Calculate distance statistics
        distance_values = list(distances.values())
        min_distance = min(distance_values) if distance_values else None
        max_distance = max(distance_values) if distance_values else None
        mean_distance = sum(distance_values) / len(distance_values) if distance_values else None

        # Calculate sum of distances normalized by sqrt(num_residues)
        # This metric grows with sum of distances but is dampened by sqrt of number of residues
        num_residues = len(distances)
        sum_distances_sqrt_normalized = sum(distance_values) / math.sqrt(num_residues) if num_residues > 0 else None

        print(f"[DEBUG] - Residue distances: {dict(distances)}")
        print(f"[DEBUG] - Contact count: {contact_count}")
        print(f"[DEBUG] - Min distance: {min_distance:.3f} Å")
        print(f"[DEBUG] - Max distance: {max_distance:.3f} Å")
        print(f"[DEBUG] - Mean distance: {mean_distance:.3f} Å")
        print(f"[DEBUG] - Sum/sqrt(N): {sum_distances_sqrt_normalized:.3f}")

        return contact_count, min_distance, max_distance, mean_distance, sum_distances_sqrt_normalized

    except Exception as e:
        print(f"[ERROR] - Error calculating contact metrics: {e}")
        return None, None, None, None, None


# Note: map_structure_id_to_table_id is now imported from id_map_utils


def load_selections_from_table(table_path: str, column_name: str, id_map: Dict[str, str] = None) -> Dict[str, str]:
    """
    Load protein selections from table CSV file.

    Args:
        table_path: Path to CSV file
        column_name: Column containing selection specifications
        id_map: ID mapping pattern for matching structure IDs to table IDs

    Returns:
        Dictionary mapping structure IDs to selection strings
    """
    if not os.path.exists(table_path):
        raise FileNotFoundError(f"Table file not found: {table_path}")

    df = pd.read_csv(table_path)
    if column_name not in df.columns:
        raise ValueError(f"Column '{column_name}' not found in table. Available columns: {list(df.columns)}")

    # Assuming the first column contains IDs
    id_column = df.columns[0]
    table_selections = {}

    # Load table selections with table IDs
    for _, row in df.iterrows():
        table_id = row[id_column]
        selection_value = row[column_name]
        table_selections[str(table_id)] = str(selection_value)

    print(f"Loaded protein selections for {len(table_selections)} table IDs from {table_path}")

    # Return table selections - they will be mapped to structure IDs in the caller
    return table_selections


def calculate_contact_metrics(structure_path: str, selections: str, ligand: str, threshold: float) -> Tuple[Optional[int], Optional[float], Optional[float], Optional[float], Optional[float]]:
    """
    Calculate protein-ligand contact metrics for a structure.

    Args:
        structure_path: Path to structure file
        selections: Protein region selections (e.g., '10-20+30-40' or None for all protein)
        ligand: Ligand residue name (3-letter code)
        threshold: Distance threshold for counting contacts

    Returns:
        Tuple of (contact_count, min_distance, max_distance, mean_distance, sum_distances_sqrt_normalized)
        or (None, None, None, None, None) if calculation failed
    """
    try:
        print(f"  - Loading structure: {structure_path}")
        print(f"  - Protein selections: {selections if selections else 'All protein residues'}")
        print(f"  - Ligand: {ligand}")
        print(f"  - Threshold: {threshold} Å")

        # Parse PDB file
        atoms = parse_pdb_file(structure_path)
        if not atoms:
            print(f"  - Error: No atoms found in structure")
            return None, None, None, None, None

        print(f"  - Parsed {len(atoms)} atoms from structure")

        # Calculate contact metrics
        contact_count, min_dist, max_dist, mean_dist, sum_sqrt_norm = calculate_protein_ligand_contacts(
            atoms, selections, ligand, threshold
        )

        print(f"  - Contact count: {contact_count}")
        print(f"  - Min distance: {min_dist:.3f} Å" if min_dist is not None else "  - Min distance: None")
        print(f"  - Max distance: {max_dist:.3f} Å" if max_dist is not None else "  - Max distance: None")
        print(f"  - Mean distance: {mean_dist:.3f} Å" if mean_dist is not None else "  - Mean distance: None")
        print(f"  - Sum/sqrt(N): {sum_sqrt_norm:.3f}" if sum_sqrt_norm is not None else "  - Sum/sqrt(N): None")

        return contact_count, min_dist, max_dist, mean_dist, sum_sqrt_norm

    except Exception as e:
        print(f"  - Error calculating contact metrics: {e}")
        return None, None, None, None, None


def analyze_protein_ligand_contacts(config_data: Dict[str, Any]) -> None:
    """
    Analyze protein-ligand contacts in structures.

    Args:
        config_data: Configuration dictionary with analysis parameters
    """
    # Load structures DataStream using pipe_biopipelines_io
    structures_ds = load_datastream(config_data['structures_json'])

    selections_config = config_data['protein_selections']
    ligand_name = config_data['ligand_name']
    contact_threshold = config_data['contact_threshold']
    contact_metric_name = config_data['contact_metric_name']
    id_map = config_data.get('id_map', {"*": "*_<N>"})  # Default to standard pattern
    output_csv = config_data['output_csv']

    print(f"Analyzing protein-ligand contacts")
    print(f"Structures: {len(structures_ds.ids)}")
    print(f"Protein selections: {selections_config}")
    print(f"Ligand: {ligand_name}")
    print(f"Contact threshold: {contact_threshold} Å")
    print(f"Contact metric: {contact_metric_name}")
    print(f"ID mapping pattern: {id_map}")

    # Handle protein selections - build map using IDs from DataStream
    selections_map = {}
    if selections_config['type'] == 'all_protein':
        # All protein residues for all structures
        print(f"Using all protein residues for all structures")
        for structure_id in structures_ds.ids:
            selections_map[structure_id] = None  # None means all protein
    elif selections_config['type'] == 'fixed':
        # Fixed selection for all structures
        fixed_selection = selections_config['value']
        print(f"Using fixed protein selection: {fixed_selection}")
        for structure_id in structures_ds.ids:
            selections_map[structure_id] = fixed_selection
    else:
        # Load from table with ID mapping
        table_path = selections_config['table_path']
        column_name = selections_config['column_name']
        table_selections = load_selections_from_table(table_path, column_name, id_map)

        # Map structure IDs to table IDs and populate selections_map
        for structure_id in structures_ds.ids:
            candidate_ids = map_table_ids_to_ids(structure_id, id_map)

            # Try all candidate IDs in priority order (most specific to least specific)
            found = False
            for candidate_id in candidate_ids:
                if candidate_id in table_selections:
                    selections_map[structure_id] = table_selections[candidate_id]
                    if candidate_id != structure_id:
                        print(f"Mapped structure ID '{structure_id}' -> table ID '{candidate_id}'")
                    found = True
                    break

            if not found:
                print(f"Warning: No table entry for structure ID '{structure_id}'. Tried: {', '.join(candidate_ids)}")

    # Process structures using iterate_files for proper ID-file matching
    results = []
    structure_items = list(iterate_files(structures_ds))
    total = len(structure_items)

    for i, (structure_id, structure_path) in enumerate(structure_items):
        if not os.path.exists(structure_path):
            print(f"Warning: Structure file not found: {structure_path}")
            continue

        print(f"\nProcessing structure {i+1}/{total}: {structure_path}")
        print(f"  - ID: {structure_id}")

        # Get selection for this structure
        selections = selections_map.get(structure_id)
        if selections_config['type'] == 'table' and selections is None:
            print(f"  - Warning: No selection found for structure ID: {structure_id}")
            continue

        # Calculate contact metrics
        contact_count, min_dist, max_dist, mean_dist, sum_sqrt_norm = calculate_contact_metrics(
            structure_path, selections, ligand_name, contact_threshold
        )

        # Store result using proper ID from DataStream
        result = {
            'id': structure_id,
            'source_structure': structure_path,
            'selections': selections if selections else 'all_protein',
            'ligand': ligand_name,
            contact_metric_name: contact_count,
            'min_distance': min_dist,
            'max_distance': max_dist,
            'mean_distance': mean_dist,
            'sum_distances_sqrt_normalized': sum_sqrt_norm
        }
        results.append(result)

        print(f"  - Result: {contact_metric_name} = {contact_count}, min_distance = {min_dist}, max_distance = {max_dist}, mean_distance = {mean_dist}")

    # Create DataFrame and save
    if results:
        df = pd.DataFrame(results)

        # Create output directory
        output_dir = os.path.dirname(output_csv)
        print(f"Creating output directory: {output_dir}")
        os.makedirs(output_dir, exist_ok=True)

        # Save results
        print(f"Writing results to: {output_csv}")
        df.to_csv(output_csv, index=False)

        print(f"\nProtein-ligand contact analysis completed successfully!")
        print(f"Analyzed {len(results)} structures")
        print(f"Results saved to: {output_csv}")
        print(f"\nResults summary:")
        print(df)

        # Statistics
        contacts = [r[contact_metric_name] for r in results if r[contact_metric_name] is not None]
        min_dists = [r['min_distance'] for r in results if r['min_distance'] is not None]
        max_dists = [r['max_distance'] for r in results if r['max_distance'] is not None]
        mean_dists = [r['mean_distance'] for r in results if r['mean_distance'] is not None]

        if contacts:
            print(f"\nContact count statistics:")
            print(f"  Min: {min(contacts)}")
            print(f"  Max: {max(contacts)}")
            print(f"  Mean: {np.mean(contacts):.1f}")
            print(f"  Std: {np.std(contacts):.1f}")

        if min_dists:
            print(f"\nMinimum distance statistics:")
            print(f"  Min: {min(min_dists):.3f} Å")
            print(f"  Max: {max(min_dists):.3f} Å")
            print(f"  Mean: {np.mean(min_dists):.3f} Å")
            print(f"  Std: {np.std(min_dists):.3f} Å")

        if mean_dists:
            print(f"\nMean distance statistics:")
            print(f"  Min: {min(mean_dists):.3f} Å")
            print(f"  Max: {max(mean_dists):.3f} Å")
            print(f"  Mean: {np.mean(mean_dists):.3f} Å")
            print(f"  Std: {np.std(mean_dists):.3f} Å")

        if not contacts:
            print(f"\nError: No valid contact measurements found!")
            raise ValueError("No valid contact measurements - all results were None")
    else:
        raise ValueError("No valid results generated - check structure files and protein selections")


def main():
    parser = argparse.ArgumentParser(description='Analyze protein-ligand contacts in structures')
    parser.add_argument('--config', required=True, help='JSON config file with analysis parameters')

    args = parser.parse_args()

    # Load configuration
    if not os.path.exists(args.config):
        print(f"Error: Config file not found: {args.config}")
        sys.exit(1)

    try:
        with open(args.config, 'r') as f:
            config_data = json.load(f)
    except Exception as e:
        print(f"Error loading config: {e}")
        sys.exit(1)

    # Validate required parameters
    required_params = ['structures_json', 'protein_selections', 'ligand_name',
                       'contact_threshold', 'contact_metric_name', 'output_csv']
    for param in required_params:
        if param not in config_data:
            print(f"Error: Missing required parameter: {param}")
            sys.exit(1)

    try:
        analyze_protein_ligand_contacts(config_data)

    except Exception as e:
        print(f"Error analyzing protein-ligand contacts: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()