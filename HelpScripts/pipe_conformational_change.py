#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Runtime helper script for ConformationalChange analysis.

This script analyzes protein structures to quantify conformational changes by computing
RMSD using PyMOL's align, super, or cealign methods.
"""

import os
import sys
import argparse
import json
import pandas as pd
from typing import Dict, Any, Optional

# Import PyMOL for structure analysis
import pymol
from pymol import cmd

# Import unified I/O utilities
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from biopipelines_io import load_datastream, iterate_files


def align_and_compute_rmsd(ref_obj: str, target_obj: str, selection: str,
                           alignment_method: str) -> Dict[str, Any]:
    """
    Align target structure to reference and return RMSD from PyMOL.

    Args:
        ref_obj: Reference PyMOL object name
        target_obj: Target PyMOL object name (will be aligned to reference)
        selection: Selection specification, or None/"all" for whole structure
        alignment_method: "align", "super", or "cealign"

    Returns:
        Dictionary with RMSD and num_aligned_atoms from PyMOL's alignment
    """
    # Create selection strings
    if selection is None or selection == "all":
        ref_sel = f"{ref_obj}"
        target_sel = f"{target_obj}"
    else:
        ref_sel = f"{ref_obj} and resi {selection}"
        target_sel = f"{target_obj} and resi {selection}"

    if alignment_method == "align":
        # Returns: [RMSD_after, atoms_after, cycles, RMSD_before, atoms_before, score, residues_aligned]
        result = cmd.align(target_sel, ref_sel)
        rmsd = result[0]
        num_atoms = result[1]
    elif alignment_method == "super":
        # Same return format as align
        result = cmd.super(target_sel, ref_sel)
        rmsd = result[0]
        num_atoms = result[1]
    elif alignment_method == "cealign":
        # Returns: {'RMSD': float, 'alignment_length': int, ...}
        result = cmd.cealign(ref_sel, target_sel)
        rmsd = result['RMSD']
        num_atoms = result['alignment_length']
    else:
        raise ValueError(f"Unknown alignment method: {alignment_method}")

    print(f"  - Aligned using: {alignment_method}, RMSD: {rmsd:.3f}, atoms: {num_atoms}")

    return {
        'RMSD': rmsd,
        'num_aligned_atoms': num_atoms
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


def lookup_selection_with_base_id(structure_id: str, selection_map: Dict[str, str]) -> Optional[str]:
    """
    Look up selection for a structure ID, recursively stripping _N suffixes if not found.

    This handles cases where target structures have IDs with additional suffixes
    (e.g., RFD3_5_Test_007_1_1_1) but the selection table has base IDs
    (e.g., RFD3_5_Test_007_1).

    Args:
        structure_id: Structure ID to look up (e.g., "RFD3_5_Test_007_1_1_1")
        selection_map: Dictionary mapping IDs to selection strings

    Returns:
        Selection string if found, None otherwise
    """
    import re

    # Try exact match first
    if structure_id in selection_map:
        return selection_map[structure_id]

    # Recursively strip trailing _N suffixes and try again
    current_id = structure_id
    while True:
        # Try to strip trailing _N suffix (where N is one or more digits)
        match = re.match(r'^(.+)_\d+$', current_id)
        if not match:
            # No more suffixes to strip
            break

        current_id = match.group(1)
        if current_id in selection_map:
            print(f"  - Matched {structure_id} to base ID {current_id}")
            return selection_map[current_id]

    return None


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
        Dictionary with RMSD and num_aligned_atoms, or None if failed
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

        # Align and get RMSD from PyMOL
        metrics = align_and_compute_rmsd(ref_obj, target_obj, selection, alignment_method)

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
    # Load DataStreams using pipe_biopipelines_io
    reference_ds = load_datastream(config_data['reference_structures_json'])
    target_ds = load_datastream(config_data['target_structures_json'])

    selection_config = config_data['selection']
    alignment_method = config_data['alignment_method']
    output_csv = config_data['output_csv']

    print(f"Analyzing conformational changes")
    print(f"Reference structures: {len(reference_ds.ids)}")
    print(f"Target structures: {len(target_ds.ids)}")
    print(f"Selection: {selection_config}")
    print(f"Alignment method: {alignment_method}")

    # Initialize PyMOL in headless mode
    pymol.pymol_argv = ['pymol', '-c']
    pymol.finish_launching()
    cmd.set("cartoon_gap_cutoff", 0)

    # Build reference lookup by ID for efficient matching
    reference_files_by_id = {}
    for ref_id, ref_file in iterate_files(reference_ds):
        reference_files_by_id[ref_id] = ref_file

    # Handle selection
    selection_map = {}
    if selection_config['type'] == 'all':
        # Use all CA atoms (whole structure RMSD)
        print("Using all atoms (whole structure RMSD)")
        for target_id in target_ds.ids:
            selection_map[target_id] = "all"
    elif selection_config['type'] == 'fixed':
        # Fixed selection for all structures
        fixed_selection = selection_config['value']
        print(f"Using fixed selection: {fixed_selection}")
        for target_id in target_ds.ids:
            selection_map[target_id] = fixed_selection
    else:
        # Load from table
        table_path = selection_config['table_path']
        column_name = selection_config['column_name']
        selection_map = load_selection_from_table(table_path, column_name)

    # Determine if reference is single or multiple
    use_single_reference = len(reference_ds.ids) == 1
    if use_single_reference:
        single_ref_id, single_ref_path = next(iterate_files(reference_ds))
        print(f"Using single reference structure: {single_ref_path}")
    else:
        print(f"Using paired reference structures")

    # Ensure we have compatible number of structures
    if not use_single_reference and len(reference_ds.ids) != len(target_ds.ids):
        print(f"Warning: Reference structures ({len(reference_ds.ids)}) and target structures ({len(target_ds.ids)}) count mismatch")

    # Process structure pairs using iterate_files for proper ID-file matching
    results = []
    target_items = list(iterate_files(target_ds))

    for i, (target_id, target_path) in enumerate(target_items):
        if not os.path.exists(target_path):
            print(f"Warning: Target structure file not found: {target_path}")
            continue

        # Get reference structure (single or paired by ID)
        if use_single_reference:
            ref_path = single_ref_path
        else:
            # Match by ID - reference and target should have same IDs
            if target_id in reference_files_by_id:
                ref_path = reference_files_by_id[target_id]
            else:
                print(f"Warning: No matching reference for target ID: {target_id}")
                continue

        if not os.path.exists(ref_path):
            print(f"Warning: Reference structure file not found: {ref_path}")
            continue

        print(f"\nProcessing structure pair {i+1}/{len(target_items)}")
        print(f"Reference: {ref_path}")
        print(f"Target: {target_path}")
        print(f"ID: {target_id}")

        # Get selection for this structure (with recursive base ID lookup)
        selection = lookup_selection_with_base_id(target_id, selection_map)
        if not selection:
            print(f"  - Warning: No selection found for structure ID: {target_id}")
            print(f"    Available IDs in selection map: {list(selection_map.keys())[:5]}...")
            continue

        # Analyze conformational change
        metrics = analyze_conformational_change(ref_path, target_path, selection, alignment_method)

        if metrics is None:
            continue

        # Store result using the proper ID from DataStream
        result = {
            'id': target_id,
            'reference_structure': ref_path,
            'target_structure': target_path,
            'selection': selection,
            'num_aligned_atoms': metrics['num_aligned_atoms'],
            'RMSD': metrics['RMSD']
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

        # Statistics
        values = df['RMSD'].dropna()
        if len(values) > 0:
            print(f"\nRMSD statistics:")
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
    required_params = ['reference_structures_json', 'target_structures_json', 'selection', 'alignment_method', 'output_csv']
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
