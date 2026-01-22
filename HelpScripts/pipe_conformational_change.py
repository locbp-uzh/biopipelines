#!/usr/bin/env python3
"""
Runtime helper script for ConformationalChange analysis.

This script analyzes protein structures to quantify conformational changes using multiple metrics:
- RMSD: Root Mean Square Deviation of Cα atoms
- max_distance: Maximum distance between any Cα pair
- mean_distance: Mean distance between Cα atoms
- sum_over_square_root: Sum of distances normalized by √(number of residues)

Supports structural alignment using PyMOL's align, super, or cealign methods.
"""

import os
import sys
import argparse
import json
import pandas as pd
import numpy as np
import math
from typing import Dict, List, Any, Optional, Tuple

# Import PyMOL for structure analysis
import pymol
from pymol import cmd


def align_structures(ref_obj: str, target_obj: str, selection: str, alignment_method: str) -> None:
    """
    Align target structure to reference using specified method.

    Args:
        ref_obj: Reference PyMOL object name
        target_obj: Target PyMOL object name (will be aligned to reference)
        selection: Selection specification for alignment, or None/"all" for whole structure
        alignment_method: "align", "super", or "cealign"
    """
    # Create selection strings
    if selection is None or selection == "all":
        # Use all CA atoms
        ref_sel = f"{ref_obj} and name CA"
        target_sel = f"{target_obj} and name CA"
    else:
        ref_sel = f"{ref_obj} and resi {selection} and name CA"
        target_sel = f"{target_obj} and resi {selection} and name CA"

    if alignment_method == "align":
        # Sequence-dependent alignment
        cmd.align(target_sel, ref_sel)
    elif alignment_method == "super":
        # Structure-based superposition
        cmd.super(target_sel, ref_sel)
    elif alignment_method == "cealign":
        # Combinatorial extension alignment
        cmd.cealign(ref_sel, target_sel)
    else:
        raise ValueError(f"Unknown alignment method: {alignment_method}")

    print(f"  - Aligned using: {alignment_method}")


def calculate_metrics(ref_obj: str, target_obj: str, selection: str) -> Dict[str, float]:
    """
    Calculate multiple conformational change metrics between reference and target structures.

    Args:
        ref_obj: Reference PyMOL object name
        target_obj: Target PyMOL object name
        selection: Selection specification for analysis, or None/"all" for whole structure

    Returns:
        Dictionary with metrics: num_residues, RMSD, max_distance, mean_distance, sum_over_square_root
    """
    # Fetch Cα atoms for selection in each object
    if selection is None or selection == "all":
        ref_model = cmd.get_model(f"{ref_obj} and name CA")
        target_model = cmd.get_model(f"{target_obj} and name CA")
    else:
        ref_model = cmd.get_model(f"{ref_obj} and resi {selection} and name CA")
        target_model = cmd.get_model(f"{target_obj} and resi {selection} and name CA")

    # Index target atoms by (chain, resi) for quick lookup
    target_coords = {(a.chain, a.resi): np.array([a.coord[0], a.coord[1], a.coord[2]])
                     for a in target_model.atom}

    # Collect matched coordinates
    ref_coords_list = []
    target_coords_list = []
    distances = []

    for ref_atom in ref_model.atom:
        key = (ref_atom.chain, ref_atom.resi)
        if key in target_coords:
            ref_coord = np.array([ref_atom.coord[0], ref_atom.coord[1], ref_atom.coord[2]])
            target_coord = target_coords[key]

            ref_coords_list.append(ref_coord)
            target_coords_list.append(target_coord)

            # Calculate distance
            dist = np.linalg.norm(ref_coord - target_coord)
            distances.append(dist)

    num_residues = len(distances)

    if num_residues == 0:
        return {
            'num_residues': 0,
            'RMSD': 0.0,
            'max_distance': 0.0,
            'mean_distance': 0.0,
            'sum_over_square_root': 0.0
        }

    # Convert to numpy arrays
    ref_coords_array = np.array(ref_coords_list)
    target_coords_array = np.array(target_coords_list)
    distances_array = np.array(distances)

    # Calculate RMSD: sqrt(mean(distances²))
    rmsd = np.sqrt(np.mean(distances_array ** 2))

    # Calculate max distance
    max_distance = np.max(distances_array)

    # Calculate mean distance
    mean_distance = np.mean(distances_array)

    # Calculate sum over square root: sum(distances) / sqrt(n)
    sum_over_square_root = np.sum(distances_array) / np.sqrt(num_residues)

    return {
        'num_residues': num_residues,
        'RMSD': rmsd,
        'max_distance': max_distance,
        'mean_distance': mean_distance,
        'sum_over_square_root': sum_over_square_root
    }


def load_selection_from_table(table_path: str, column_name: str) -> Dict[str, str]:
    """
    Load selection specifications from table CSV file.

    Args:
        table_path: Path to CSV file
        column_name: Column containing selection specifications

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
    selection_map = {}

    for _, row in df.iterrows():
        structure_id = row[id_column]
        selection_value = row[column_name]
        selection_map[str(structure_id)] = str(selection_value)

    print(f"Loaded selections for {len(selection_map)} structures from {table_path}")

    return selection_map


def analyze_conformational_change(ref_path: str, target_path: str, selection: str,
                                  alignment_method: str) -> Optional[Dict[str, Any]]:
    """
    Analyze conformational change between reference and target structures.

    Args:
        ref_path: Path to reference structure file
        target_path: Path to target structure file
        selection: Selection specification (e.g., '10-20+30-40')
        alignment_method: Alignment method ("align", "super", or "cealign")

    Returns:
        Dictionary with metrics or None if calculation failed
    """
    try:
        # Extract structure IDs from filenames for PyMOL object names
        ref_id = os.path.splitext(os.path.basename(ref_path))[0]
        target_id = os.path.splitext(os.path.basename(target_path))[0]

        ref_obj = f"ref_{ref_id}"
        target_obj = f"target_{target_id}"

        # Load structures into PyMOL
        cmd.load(ref_path, ref_obj)
        cmd.load(target_path, target_obj)

        print(f"  - Loaded reference: {ref_obj}")
        print(f"  - Loaded target: {target_obj}")
        print(f"  - Selection: {selection}")

        # Align structures
        align_structures(ref_obj, target_obj, selection, alignment_method)

        # Calculate metrics
        metrics = calculate_metrics(ref_obj, target_obj, selection)

        print(f"  - Num residues: {metrics['num_residues']}")
        print(f"  - RMSD: {metrics['RMSD']:.3f}")
        print(f"  - Max distance: {metrics['max_distance']:.3f}")
        print(f"  - Mean distance: {metrics['mean_distance']:.3f}")
        print(f"  - Sum over sqrt: {metrics['sum_over_square_root']:.3f}")

        # Clean up PyMOL objects
        cmd.delete(ref_obj)
        cmd.delete(target_obj)

        return metrics

    except Exception as e:
        print(f"  - Error analyzing conformational change: {e}")
        import traceback
        traceback.print_exc()
        return None


def analyze_all_conformational_changes(config_data: Dict[str, Any]) -> None:
    """
    Analyze conformational changes for all structure pairs.

    Args:
        config_data: Configuration dictionary with analysis parameters
    """
    reference_structures = config_data['reference_structures']
    target_structures = config_data['target_structures']
    selection_config = config_data['selection']
    alignment_method = config_data['alignment_method']
    output_csv = config_data['output_csv']

    print(f"Analyzing conformational changes")
    print(f"Reference structures: {len(reference_structures)}")
    print(f"Target structures: {len(target_structures)}")
    print(f"Selection: {selection_config}")
    print(f"Alignment method: {alignment_method}")

    # Initialize PyMOL in headless mode
    pymol.pymol_argv = ['pymol', '-c']
    pymol.finish_launching()
    cmd.set("cartoon_gap_cutoff", 0)

    # Handle selection
    selection_map = {}
    use_all_selection = False
    if selection_config['type'] == 'all':
        # Use all CA atoms (whole structure RMSD)
        print("Using all CA atoms (whole structure RMSD)")
        use_all_selection = True
        for target_path in target_structures:
            structure_id = os.path.splitext(os.path.basename(target_path))[0]
            selection_map[structure_id] = "all"
    elif selection_config['type'] == 'fixed':
        # Fixed selection for all structures
        fixed_selection = selection_config['value']
        print(f"Using fixed selection: {fixed_selection}")
        # Create mapping for all structures (we'll use structure IDs as keys)
        for target_path in target_structures:
            structure_id = os.path.splitext(os.path.basename(target_path))[0]
            selection_map[structure_id] = fixed_selection
    else:
        # Load from table
        table_path = selection_config['table_path']
        column_name = selection_config['column_name']
        selection_map = load_selection_from_table(table_path, column_name)

    # Determine if reference is single or multiple
    use_single_reference = len(reference_structures) == 1
    if use_single_reference:
        print(f"Using single reference structure: {reference_structures[0]}")
    else:
        print(f"Using paired reference structures")

    # Ensure we have compatible number of structures
    if not use_single_reference and len(reference_structures) != len(target_structures):
        print(f"Warning: Reference structures ({len(reference_structures)}) and target structures ({len(target_structures)}) count mismatch")

    # Process structure pairs
    results = []

    for i, target_path in enumerate(target_structures):
        if not os.path.exists(target_path):
            print(f"Warning: Target structure file not found: {target_path}")
            continue

        # Get reference structure (single or paired)
        if use_single_reference:
            ref_path = reference_structures[0]
        else:
            if i >= len(reference_structures):
                print(f"Warning: No matching reference for target {i+1}")
                continue
            ref_path = reference_structures[i]

        if not os.path.exists(ref_path):
            print(f"Warning: Reference structure file not found: {ref_path}")
            continue

        print(f"\nProcessing structure pair {i+1}/{len(target_structures)}")
        print(f"Reference: {ref_path}")
        print(f"Target: {target_path}")

        # Extract structure ID from target filename
        structure_id = os.path.splitext(os.path.basename(target_path))[0]

        # Get selection for this structure
        selection = selection_map.get(structure_id)
        if not selection:
            print(f"  - Warning: No selection found for structure ID: {structure_id}")
            continue

        # Analyze conformational change
        metrics = analyze_conformational_change(ref_path, target_path, selection, alignment_method)

        if metrics is None:
            continue

        # Store result
        result = {
            'id': structure_id,
            'reference_structure': ref_path,
            'target_structure': target_path,
            'selection': selection,
            'num_residues': metrics['num_residues'],
            'RMSD': metrics['RMSD'],
            'max_distance': metrics['max_distance'],
            'mean_distance': metrics['mean_distance'],
            'sum_over_square_root': metrics['sum_over_square_root']
        }
        results.append(result)

        print(f"  - Result stored")

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

        print(f"\nConformational change analysis completed successfully!")
        print(f"Analyzed {len(results)} structure pairs")
        print(f"Results saved to: {output_csv}")
        print(f"\nResults summary:")
        print(df)

        # Statistics for each metric
        print(f"\nMetrics statistics:")
        for metric in ['RMSD', 'max_distance', 'mean_distance', 'sum_over_square_root']:
            values = df[metric].dropna()
            if len(values) > 0:
                print(f"\n{metric}:")
                print(f"  Min: {values.min():.3f}")
                print(f"  Max: {values.max():.3f}")
                print(f"  Mean: {values.mean():.3f}")
                print(f"  Std: {values.std():.3f}")
    else:
        raise ValueError("No valid results generated - check structure files and selections")

    # Clean up PyMOL
    cmd.quit()


def main():
    parser = argparse.ArgumentParser(description='Analyze conformational changes between reference and target structures')
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
    required_params = ['reference_structures', 'target_structures', 'selection', 'alignment_method', 'output_csv']
    for param in required_params:
        if param not in config_data:
            print(f"Error: Missing required parameter: {param}")
            sys.exit(1)

    try:
        analyze_all_conformational_changes(config_data)

    except Exception as e:
        print(f"Error analyzing conformational changes: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
