#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Runtime helper script for Contacts analysis.

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
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from biopipelines.biopipelines_io import load_datastream, iterate_files, lookup_table_value
from biopipelines.id_map_utils import get_mapped_ids

# Import PDB parser and selection utilities
from biopipelines.pdb_parser import Atom as _PdbAtom, parse_pdb_file as _pdb_parse_pdb_file, STANDARD_RESIDUES
from biopipelines.sele_utils import sele_to_list as _sele_to_list


def parse_pdb_file(pdb_path: str) -> List[_PdbAtom]:
    """
    Parse PDB file, excluding atoms at origin (placeholder atoms).

    Wraps :func:`pdb_parser.parse_pdb_file` and filters out (0,0,0) atoms.
    """
    all_atoms = _pdb_parse_pdb_file(pdb_path)
    skipped = 0
    atoms = []
    for atom in all_atoms:
        if abs(atom.x) < 0.001 and abs(atom.y) < 0.001 and abs(atom.z) < 0.001:
            skipped += 1
            continue
        atoms.append(atom)
    if skipped > 0:
        print(f"  - Skipped {skipped} atoms at origin (0,0,0) - likely missing/placeholder atoms")
    return atoms


# Use pdb_parser.Atom everywhere — alias for convenience in this module
Atom = _PdbAtom


def parse_residue_selection(selection_str: str) -> List[Tuple[str, int]]:
    """
    Parse residue selection string into list of (chain, resnum) tuples.

    Delegates to :func:`sele_utils.sele_to_list` which handles both
    chain-aware (``A1-117+B10``) and chainless (``1-50+60``) formats.
    """
    if not selection_str or selection_str.lower() in ['all', 'none', 'null']:
        return []

    return _sele_to_list(str(selection_str))


def load_selection_from_table(table_path: str, column_name: str) -> Dict[str, str]:
    """
    Load per-structure selection specifications from a table CSV file.

    The first column is taken as the ID column. Returns a dict mapping
    structure ID -> selection string.
    """
    if not os.path.exists(table_path):
        raise FileNotFoundError(f"Table file not found: {table_path}")

    df = pd.read_csv(table_path)
    if column_name not in df.columns:
        raise ValueError(f"Column '{column_name}' not found in table. Available columns: {list(df.columns)}")

    id_column = df.columns[0]
    selection_map = {}
    for _, row in df.iterrows():
        selection_map[str(row[id_column])] = str(row[column_name])

    print(f"Loaded selections for {len(selection_map)} structures from {table_path}")
    return selection_map


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

    for atom in atoms:
        # Include standard amino acids, exclude specified ligand
        if atom.res_name in STANDARD_RESIDUES and atom.res_name != ligand_name:
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
    return [atom for atom in atoms if atom.res_name == ligand_name]


def filter_atoms_by_selection(atoms: List[Atom], residue_tuples: List[Tuple[str, int]]) -> List[Atom]:
    """
    Filter atoms to include only specified (chain, resnum) pairs.

    Args:
        atoms: List of atoms
        residue_tuples: List of (chain, resnum) tuples to include.
            If chain is '' (chainless selection), matches any chain.

    Returns:
        Filtered list of atoms
    """
    if not residue_tuples:  # Empty list means no filtering
        return atoms

    # Build a lookup set; separate chainless entries (match any chain)
    chained = set()
    chainless_nums = set()
    for chain, resnum in residue_tuples:
        if chain:
            chained.add((chain, resnum))
        else:
            chainless_nums.add(resnum)

    return [atom for atom in atoms
            if (atom.chain, atom.res_num) in chained or atom.res_num in chainless_nums]


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


def calculate_contacts(atoms: List[Atom], protein_selections: str, ligand_name: str, contact_threshold: float, reference: str = "ligand") -> Tuple[Optional[int], Optional[float], Optional[float], Optional[float], Optional[float]]:
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
        # Protein atoms always exclude HETATM; the reference atoms are either the
        # named ligand or a residue selection (residues-vs-residues contacts).
        protein_atoms = get_protein_atoms(atoms, ligand_name)
        if reference == "ligand":
            ligand_atoms = get_ligand_atoms(atoms, ligand_name)
        else:
            ref_residues = parse_residue_selection(reference)
            ligand_atoms = filter_atoms_by_selection(protein_atoms, ref_residues) if ref_residues else []

        print(f"[DEBUG] - Found {len(protein_atoms)} protein atoms")
        print(f"[DEBUG] - Found {len(ligand_atoms)} reference atoms ({reference})")

        if not ligand_atoms:
            print(f"[WARNING] - No reference atoms found for '{reference}'")
            return None, None, None, None, None

        # Filter protein atoms by selection if specified
        if protein_selections and protein_selections.lower() not in ['all', 'none', 'null']:
            residue_numbers = parse_residue_selection(protein_selections)
            if residue_numbers:
                protein_atoms = filter_atoms_by_selection(protein_atoms, residue_numbers)
                print(f"[DEBUG] - After selection filtering: {len(protein_atoms)} protein atoms")

        # In residue-reference mode, drop reference residues from the protein side
        # so a residue is not counted as contacting itself. Chainless ref entries
        # ('', resnum) match any chain, mirroring filter_atoms_by_selection.
        if reference != "ligand" and ref_residues:
            ref_with_chain = {(c, r) for c, r in ref_residues if c}
            ref_any_chain = {r for c, r in ref_residues if not c}
            protein_atoms = [a for a in protein_atoms
                             if a.res_num not in ref_any_chain
                             and (a.chain, a.res_num) not in ref_with_chain]
            print(f"[DEBUG] - After reference-exclusion: {len(protein_atoms)} protein atoms")

        if not protein_atoms:
            print(f"[WARNING] - No protein atoms found after selection filtering")
            return None, None, None, None, None

        # Calculate minimum distances per residue
        distances = {}
        closest_ligand_atoms = {}  # Track which ligand atom is closest to each residue
        closest_protein_atoms = {}  # Track which protein atom is closest for each residue
        closest_coords = {}  # Track coordinates for debugging

        for protein_atom in protein_atoms:
            residue_key = (protein_atom.chain, protein_atom.res_num)

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


def calculate_contact_metrics(structure_path: str, selections: str, ligand: str, threshold: float, reference: str = "ligand") -> Tuple[Optional[int], Optional[float], Optional[float], Optional[float], Optional[float]]:
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
        contact_count, min_dist, max_dist, mean_dist, sum_sqrt_norm = calculate_contacts(
            atoms, selections, ligand, threshold, reference
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


def analyze_contacts(config_data: Dict[str, Any]) -> None:
    """
    Analyze protein-ligand contacts in structures.

    Args:
        config_data: Configuration dictionary with analysis parameters
    """
    # Load structures DataStream using pipe_biopipelines_io
    structures_ds = load_datastream(config_data['structures_json'])

    selections_config = config_data['protein_selections']
    ligand_name = config_data['ligand_name']
    reference = config_data.get('reference', 'ligand')
    contact_threshold = config_data['contact_threshold']
    contact_metric_name = config_data['contact_metric_name']
    output_csv = config_data['output_csv']

    print(f"Analyzing contacts (reference: {reference})")
    print(f"Structures: {len(structures_ds.ids_expanded)}")
    print(f"Protein selections: {selections_config}")
    if reference == "ligand":
        print(f"Ligand: {ligand_name}")
    print(f"Contact threshold: {contact_threshold} Å")
    print(f"Contact metric: {contact_metric_name}")

    # Handle protein selections - build map using IDs from DataStream.
    # For all_protein / fixed the same value applies to every structure;
    # for table_column each structure's selection is resolved by ID match
    # (single-row table broadcasts to all; otherwise unmatched IDs are skipped).
    structure_ids = list(structures_ds.ids_expanded)
    selections_map = {}
    structure_to_sele_id = None
    if selections_config['type'] == 'all_protein':
        print(f"Using all protein residues for all structures")
        for structure_id in structure_ids:
            selections_map[structure_id] = None
    elif selections_config['type'] == 'fixed':
        fixed_selection = selections_config['value']
        print(f"Using fixed protein selection: {fixed_selection}")
        for structure_id in structure_ids:
            selections_map[structure_id] = fixed_selection
    elif selections_config['type'] == 'table_column':
        table_path = selections_config['table_path']
        column_name = selections_config['column_name']
        selections_map = load_selection_from_table(table_path, column_name)
        if len(selections_map) == 1:
            single_value = next(iter(selections_map.values()))
            print(f"Using single table selection for all structures: {single_value}")
            for structure_id in structure_ids:
                selections_map[structure_id] = single_value
        else:
            structure_to_sele_id = get_mapped_ids(
                source_ids=structure_ids,
                target_ids=list(selections_map.keys()),
                unique=True,
            )
    else:
        raise ValueError(f"Unsupported protein_selections type: {selections_config['type']}")

    # Handle a per-structure reference table the same way as selections:
    # `reference == "table_column"` means the contact-target residues come from
    # a table column resolved by ID match (single-row broadcasts; unmatched
    # IDs are skipped). For "ligand" / static-string references this is a no-op.
    reference_map = {}
    structure_to_ref_id = None
    if reference == "table_column":
        ref_table = config_data['reference_table']
        reference_map = load_selection_from_table(ref_table['table_path'], ref_table['column_name'])
        if len(reference_map) == 1:
            single_ref = next(iter(reference_map.values()))
            print(f"Using single table reference for all structures: {single_ref}")
            for structure_id in structure_ids:
                reference_map[structure_id] = single_ref
        else:
            structure_to_ref_id = get_mapped_ids(
                source_ids=structure_ids,
                target_ids=list(reference_map.keys()),
                unique=True,
            )

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

        if structure_to_sele_id is not None:
            matched_sele_id = structure_to_sele_id.get(structure_id)
            if matched_sele_id is None:
                print(f"Warning: No selection-table entry for structure ID '{structure_id}', skipping")
                continue
            if matched_sele_id != structure_id:
                print(f"  - Matched selection ID: {structure_id} -> {matched_sele_id}")
            selections = selections_map[matched_sele_id]
        else:
            selections = selections_map.get(structure_id)

        # Resolve the per-structure reference (table_column mode); otherwise the
        # config reference string applies unchanged.
        if reference == "table_column":
            if structure_to_ref_id is not None:
                matched_ref_id = structure_to_ref_id.get(structure_id)
                if matched_ref_id is None:
                    print(f"Warning: No reference-table entry for structure ID '{structure_id}', skipping")
                    continue
                if matched_ref_id != structure_id:
                    print(f"  - Matched reference ID: {structure_id} -> {matched_ref_id}")
                structure_reference = reference_map[matched_ref_id]
            else:
                structure_reference = reference_map[structure_id]
        else:
            structure_reference = reference

        # Calculate contact metrics
        contact_count, min_dist, max_dist, mean_dist, sum_sqrt_norm = calculate_contact_metrics(
            structure_path, selections, ligand_name, contact_threshold, structure_reference
        )

        # Store result using proper ID from DataStream
        result = {
            'id': structure_id,
            'source_structure': structure_path,
            'selections': selections if selections else 'all_protein',
            'ligand': ligand_name if structure_reference == "ligand" else structure_reference,
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
    parser.add_argument('--ligand', default=None,
                        help='Ligand residue code, overrides ligand_name in the config '
                             '(used when the code is resolved from a compounds stream at runtime)')

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

    if args.ligand:
        config_data['ligand_name'] = args.ligand

    # ligand_name is only required when reference="ligand"; residue-reference
    # mode (static string or table_column) resolves contact targets from the
    # residue selection instead.
    reference = config_data.get('reference', 'ligand')
    required_params = ['structures_json', 'protein_selections',
                       'contact_threshold', 'contact_metric_name', 'output_csv']
    if reference == 'ligand':
        required_params.append('ligand_name')
    elif reference == 'table_column':
        required_params.append('reference_table')
    for param in required_params:
        if param not in config_data:
            print(f"Error: Missing required parameter: {param}")
            sys.exit(1)

    if reference == 'ligand' and not config_data['ligand_name']:
        print("Error: ligand residue code is empty (no --ligand override and no ligand_name in config)")
        sys.exit(1)

    try:
        analyze_contacts(config_data)

    except Exception as e:
        print(f"Error analyzing protein-ligand contacts: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()