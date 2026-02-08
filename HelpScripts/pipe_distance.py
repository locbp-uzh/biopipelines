#!/usr/bin/env python3
# Copyright (c) 2026 Gianluca Quargnali @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# This software is freely available for use, modification, and redistribution.
# If you use this software or any derivative of it in published work,
# you must cite the original author and this repository.
"""
Runtime helper script for Distance analysis.

This script analyzes protein structures to calculate distances between specific atoms and residues,
outputting CSV with distance metrics for all structures. Parses PDB files directly as text.
"""

import os
import sys
import argparse
import json
import pandas as pd
import numpy as np
from typing import Dict, List, Any, Optional, Tuple

# Import unified I/O utilities
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from biopipelines_io import load_datastream, iterate_files
from pdb_parser import parse_pdb_file, parse_selection, calculate_distances, debug_ligand_atoms



def calculate_distance(structure_path: str, atom_selection, residue_selection,
                      distance_metric: str) -> Optional[float]:
    """
    Calculate distance between selections by parsing PDB file.
    Supports three modes:
    - Atom-Residue: atom=str, residue=str
    - Atom-Atom: atom=[str, str], residue=None
    - Residue-Residue: atom=None, residue=[str, str]

    Args:
        structure_path: Path to structure file
        atom_selection: Atom selection string, list of two strings, or None
        residue_selection: Residue selection string, list of two strings, or None
        distance_metric: Type of distance calculation (min, max, mean)

    Returns:
        Calculated distance or None if failed
    """
    try:
        # Parse PDB file
        atoms = parse_pdb_file(structure_path)

        if not atoms:
            print(f"  - Warning: No atoms found in {structure_path}")
            return None

        print(f"  - Total atoms in structure: {len(atoms)}")

        # Determine mode and parse selections
        if isinstance(atom_selection, list) and len(atom_selection) == 2:
            # Atom-Atom mode
            selection1 = parse_selection(atom_selection[0], atoms)
            selection2 = parse_selection(atom_selection[1], atoms)
            print(f"  - Mode: Atom-Atom")
            print(f"  - Atom 1 selection '{atom_selection[0]}': {len(selection1)} atoms")
            print(f"  - Atom 2 selection '{atom_selection[1]}': {len(selection2)} atoms")

        elif isinstance(residue_selection, list) and len(residue_selection) == 2:
            # Residue-Residue mode
            selection1 = parse_selection(residue_selection[0], atoms)
            selection2 = parse_selection(residue_selection[1], atoms)
            print(f"  - Mode: Residue-Residue")
            print(f"  - Residue 1 selection '{residue_selection[0]}': {len(selection1)} atoms")
            print(f"  - Residue 2 selection '{residue_selection[1]}': {len(selection2)} atoms")

        else:
            # Atom-Residue mode (original behavior)
            selection1 = parse_selection(atom_selection, atoms)
            selection2 = parse_selection(residue_selection, atoms)
            print(f"  - Mode: Atom-Residue")
            print(f"  - Atom selection '{atom_selection}': {len(selection1)} atoms")
            print(f"  - Residue selection '{residue_selection}': {len(selection2)} atoms")

        if not selection1:
            print(f"  - Warning: No atoms found for first selection")
            return None

        if not selection2:
            print(f"  - Warning: No atoms found for second selection")
            return None

        # Calculate distances using pdb_parser
        distance = calculate_distances(selection1, selection2, distance_metric)

        if distance is not None:
            print(f"  - Distance ({distance_metric}): {distance:.3f} Å")

        return distance

    except Exception as e:
        print(f"  - Error calculating distance: {e}")
        import traceback
        traceback.print_exc()
        return None


def analyze_distances(config_data: Dict[str, Any]) -> None:
    """
    Analyze distances between selections in structures by parsing PDB files.
    Supports atom-residue, atom-atom, and residue-residue modes.

    Args:
        config_data: Configuration dictionary with analysis parameters
    """
    # Load structures DataStream using pipe_biopipelines_io
    structures_ds = load_datastream(config_data['structures_json'])

    atom_selection = config_data['atom_selection']
    residue_selection = config_data['residue_selection']
    distance_metric = config_data['distance_metric']
    metric_name = config_data['metric_name']
    output_csv = config_data['output_csv']

    print(f"Analyzing distances in {len(structures_ds.ids)} structures")

    # Determine and display mode
    if isinstance(atom_selection, list):
        print(f"Mode: Atom-Atom")
        print(f"Atom 1: {atom_selection[0]}")
        print(f"Atom 2: {atom_selection[1]}")
    elif isinstance(residue_selection, list):
        print(f"Mode: Residue-Residue")
        print(f"Residue 1: {residue_selection[0]}")
        print(f"Residue 2: {residue_selection[1]}")
    else:
        print(f"Mode: Atom-Residue")
        print(f"Atom selection: {atom_selection}")
        print(f"Residue selection: {residue_selection}")

    print(f"Distance metric: {distance_metric}")
    print(f"Output column: {metric_name}")
    
    # Process each structure using iterate_files for proper ID-file matching
    results = []
    structure_items = list(iterate_files(structures_ds))
    total = len(structure_items)

    for i, (structure_id, structure_path) in enumerate(structure_items):
        if not os.path.exists(structure_path):
            print(f"Warning: Structure file not found: {structure_path}")
            continue

        print(f"\nProcessing structure {i+1}/{total}: {structure_path}")
        print(f"  - ID: {structure_id}")

        # Calculate distance
        distance = calculate_distance(structure_path, atom_selection, residue_selection, distance_metric)

        # Store result using proper ID from DataStream
        result = {
            'id': structure_id,
            'source_structure': structure_path,
            metric_name: distance
        }
        results.append(result)

        print(f"  - Result: {metric_name} = {distance}")
    
    # Create DataFrame and save
    if results:
        df = pd.DataFrame(results)
        
        # Create output directory with debug logging
        output_dir = os.path.dirname(output_csv)
        print(f"Debug: Creating output directory: {output_dir}")
        try:
            os.makedirs(output_dir, exist_ok=True)
            print(f"Debug: Output directory created/verified: {output_dir}")
            print(f"Debug: Directory exists: {os.path.exists(output_dir)}")
            print(f"Debug: Directory is writable: {os.access(output_dir, os.W_OK)}")
        except Exception as e:
            print(f"Debug: Error creating output directory: {e}")
            raise
        
        # Save results with debug logging
        print(f"Debug: About to write CSV to: {output_csv}")
        try:
            df.to_csv(output_csv, index=False)
            print(f"Debug: CSV write call completed")
            
            # Verify file exists immediately after write
            if os.path.exists(output_csv):
                file_size = os.path.getsize(output_csv)
                print(f"Debug: File verified exists, size: {file_size} bytes")
                print(f"Debug: File is readable: {os.access(output_csv, os.R_OK)}")
            else:
                print(f"Debug: ERROR - File does not exist immediately after write!")
                # List directory contents to see what's there
                try:
                    dir_contents = os.listdir(output_dir)
                    print(f"Debug: Directory contents: {dir_contents}")
                except Exception as e:
                    print(f"Debug: Error listing directory: {e}")
        except Exception as e:
            print(f"Debug: Error writing CSV file: {e}")
            import traceback
            traceback.print_exc()
            raise
        
        print(f"\nDistance analysis completed successfully!")
        print(f"Analyzed {len(results)} structures")
        print(f"Results saved to: {output_csv}")
        print(f"\nResults summary:")
        print(df)
        
        # Statistics
        distances = [r[metric_name] for r in results if r[metric_name] is not None]
        if distances:
            print(f"\nDistance statistics:")
            print(f"  Min: {min(distances):.3f} Å")
            print(f"  Max: {max(distances):.3f} Å")
            print(f"  Mean: {np.mean(distances):.3f} Å")
            print(f"  Std: {np.std(distances):.3f} Å")
        else:
            print(f"\nError: No valid distance measurements found!")
            print(f"All {len(results)} structures returned None for distance calculations.")
            print(f"Check that:")
            print(f"  - Atom selection '{atom_selection}' matches atoms in structures")
            print(f"  - Residue selection '{residue_selection}' matches residues in structures")
            raise ValueError("No valid distance measurements - all results were None")
    else:
        raise ValueError("No valid results generated - check structure files and selections")


def main():
    parser = argparse.ArgumentParser(description='Analyze residue-atom distances in structures')
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
    required_params = ['structures_json', 'atom_selection', 'residue_selection',
                      'distance_metric', 'metric_name', 'output_csv']
    for param in required_params:
        if param not in config_data:
            print(f"Error: Missing required parameter: {param}")
            sys.exit(1)
    
    try:
        analyze_distances(config_data)
        
    except Exception as e:
        print(f"Error analyzing distances: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()