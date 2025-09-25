#!/usr/bin/env python3
"""
Runtime helper script for MobileDistance analysis.

This script analyzes protein structures to calculate sum of Cα distances between
APO and HOLO forms over mobile residues, normalized by √(number of residues).
Based on calculate_dist_mobiles function from rfdaa_results.py.
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


def calculate_dist_mobiles(p1, p2, mobile):
    """
    Sum Cα distances between p1 and p2 over residues in mobile,
    then normalize by √(number of residues).
    Based on the function from rfdaa_results.py
    """
    # Fetch Cα atoms for mobile region in each object
    m1 = cmd.get_model(f"{p1} and resi {mobile} and name CA")
    m2 = cmd.get_model(f"{p2} and resi {mobile} and name CA")

    # Index p2 atoms by (chain, resi) for quick lookup
    coords2 = {(a.chain, a.resi): (a.coord[0], a.coord[1], a.coord[2]) for a in m2.atom}

    total_dist = 0.0
    count = 0

    for a1 in m1.atom:
        key = (a1.chain, a1.resi)
        if key in coords2:
            x1, y1, z1 = a1.coord
            x2, y2, z2 = coords2[key]
            dx, dy, dz = x1 - x2, y1 - y2, z1 - z2
            total_dist += math.sqrt(dx*dx + dy*dy + dz*dz)
            count += 1

    if count == 0:
        return 0.0

    return total_dist / math.sqrt(count)


def load_mobile_from_datasheet(datasheet_path: str, column_name: str) -> Dict[str, str]:
    """
    Load mobile regions from datasheet CSV file.

    Args:
        datasheet_path: Path to CSV file
        column_name: Column containing mobile region specifications

    Returns:
        Dictionary mapping structure IDs to mobile region strings
    """
    if not os.path.exists(datasheet_path):
        raise FileNotFoundError(f"Datasheet file not found: {datasheet_path}")

    df = pd.read_csv(datasheet_path)
    if column_name not in df.columns:
        raise ValueError(f"Column '{column_name}' not found in datasheet. Available columns: {list(df.columns)}")

    # Assuming the first column contains IDs
    id_column = df.columns[0]
    mobile_map = {}

    for _, row in df.iterrows():
        structure_id = row[id_column]
        mobile_value = row[column_name]
        mobile_map[str(structure_id)] = str(mobile_value)

    print(f"Loaded mobile regions for {len(mobile_map)} structures from {datasheet_path}")

    return mobile_map


def calculate_mobile_distance(apo_path: str, holo_path: str, mobile: str) -> Optional[float]:
    """
    Calculate mobile distance between APO and HOLO structures.

    Args:
        apo_path: Path to APO structure file
        holo_path: Path to HOLO structure file
        mobile: Mobile region specification (e.g., '10-20+30-40')

    Returns:
        Mobile distance or None if calculation failed
    """
    try:
        # Extract structure IDs from filenames for PyMOL object names
        apo_id = os.path.splitext(os.path.basename(apo_path))[0]
        holo_id = os.path.splitext(os.path.basename(holo_path))[0]

        apo_obj = f"apo_{apo_id}"
        holo_obj = f"holo_{holo_id}"

        # Load structures into PyMOL
        cmd.load(apo_path, apo_obj)
        cmd.load(holo_path, holo_obj)

        print(f"  - Loaded APO: {apo_obj}")
        print(f"  - Loaded HOLO: {holo_obj}")
        print(f"  - Mobile region: {mobile}")

        # Calculate mobile distance
        distance = calculate_dist_mobiles(apo_obj, holo_obj, mobile)

        print(f"  - Mobile distance: {distance:.3f}")

        # Clean up PyMOL objects
        cmd.delete(apo_obj)
        cmd.delete(holo_obj)

        return distance

    except Exception as e:
        print(f"  - Error calculating mobile distance: {e}")
        return None


def analyze_mobile_distances(config_data: Dict[str, Any]) -> None:
    """
    Analyze mobile distances between APO and HOLO structures.

    Args:
        config_data: Configuration dictionary with analysis parameters
    """
    apo_structures = config_data['apo_structures']
    holo_structures = config_data['holo_structures']
    mobile_config = config_data['mobile_regions']
    metric_name = config_data['metric_name']
    output_csv = config_data['output_csv']

    print(f"Analyzing mobile distances")
    print(f"APO structures: {len(apo_structures)}")
    print(f"HOLO structures: {len(holo_structures)}")
    print(f"Mobile regions: {mobile_config}")
    print(f"Output column: {metric_name}")

    # Initialize PyMOL in headless mode
    pymol.pymol_argv = ['pymol', '-c']
    pymol.finish_launching()
    cmd.set("cartoon_gap_cutoff", 0)

    # Handle mobile regions
    mobile_map = {}
    if mobile_config['type'] == 'fixed':
        # Fixed mobile region for all structures
        fixed_mobile = mobile_config['value']
        print(f"Using fixed mobile region: {fixed_mobile}")
        # Create mapping for all structures (we'll use structure IDs as keys)
        for apo_path in apo_structures:
            structure_id = os.path.splitext(os.path.basename(apo_path))[0]
            mobile_map[structure_id] = fixed_mobile
    else:
        # Load from datasheet
        datasheet_path = mobile_config['datasheet_path']
        column_name = mobile_config['column_name']
        mobile_map = load_mobile_from_datasheet(datasheet_path, column_name)

    # Ensure we have the same number of APO and HOLO structures
    if len(apo_structures) != len(holo_structures):
        print(f"Warning: APO structures ({len(apo_structures)}) and HOLO structures ({len(holo_structures)}) count mismatch")

    # Process structure pairs
    results = []

    for i, (apo_path, holo_path) in enumerate(zip(apo_structures, holo_structures)):
        if not os.path.exists(apo_path):
            print(f"Warning: APO structure file not found: {apo_path}")
            continue

        if not os.path.exists(holo_path):
            print(f"Warning: HOLO structure file not found: {holo_path}")
            continue

        print(f"\nProcessing structure pair {i+1}/{len(apo_structures)}")
        print(f"APO: {apo_path}")
        print(f"HOLO: {holo_path}")

        # Extract structure ID from APO filename
        structure_id = os.path.splitext(os.path.basename(apo_path))[0]

        # Get mobile region for this structure
        mobile = mobile_map.get(structure_id)
        if not mobile:
            print(f"  - Warning: No mobile region found for structure ID: {structure_id}")
            continue

        # Calculate mobile distance
        distance = calculate_mobile_distance(apo_path, holo_path, mobile)

        # Store result
        result = {
            'id': structure_id,
            'apo_structure': apo_path,
            'holo_structure': holo_path,
            'mobile': mobile,
            metric_name: distance
        }
        results.append(result)

        print(f"  - Result: {metric_name} = {distance}")

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

        print(f"\nMobile distance analysis completed successfully!")
        print(f"Analyzed {len(results)} structure pairs")
        print(f"Results saved to: {output_csv}")
        print(f"\nResults summary:")
        print(df)

        # Statistics
        distances = [r[metric_name] for r in results if r[metric_name] is not None]
        if distances:
            print(f"\nDistance statistics:")
            print(f"  Min: {min(distances):.3f}")
            print(f"  Max: {max(distances):.3f}")
            print(f"  Mean: {np.mean(distances):.3f}")
            print(f"  Std: {np.std(distances):.3f}")
        else:
            print(f"\nError: No valid distance measurements found!")
            raise ValueError("No valid distance measurements - all results were None")
    else:
        raise ValueError("No valid results generated - check structure files and mobile regions")

    # Clean up PyMOL
    cmd.quit()


def main():
    parser = argparse.ArgumentParser(description='Analyze mobile distances between APO and HOLO structures')
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
    required_params = ['apo_structures', 'holo_structures', 'mobile_regions', 'metric_name', 'output_csv']
    for param in required_params:
        if param not in config_data:
            print(f"Error: Missing required parameter: {param}")
            sys.exit(1)

    try:
        analyze_mobile_distances(config_data)

    except Exception as e:
        print(f"Error analyzing mobile distances: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()