#!/usr/bin/env python3
"""
Helper script for BoltzGen postprocessing.

This script handles any additional processing of BoltzGen outputs
that may be needed for integration with BioPipelines.
"""

import os
import sys
import argparse
import pandas as pd
from pathlib import Path


def collect_final_structures(final_designs_folder, output_list_file):
    """
    Collect all structure files from final designs folder.

    Args:
        final_designs_folder: Path to final_ranked_designs/final_<budget>_designs/
        output_list_file: Path to output file listing all structures
    """
    structures = []

    if not os.path.exists(final_designs_folder):
        print(f"Warning: Final designs folder not found: {final_designs_folder}")
        return structures

    # Collect all .cif and .pdb files
    for ext in ['*.cif', '*.pdb']:
        structures.extend(Path(final_designs_folder).glob(ext))

    # Sort by name for consistent ordering
    structures.sort()

    # Write structure list to file
    with open(output_list_file, 'w') as f:
        for struct in structures:
            f.write(f"{struct.absolute()}\n")

    print(f"Found {len(structures)} final structures")
    return [str(s) for s in structures]


def validate_metrics_csv(metrics_csv):
    """
    Validate that metrics CSV was generated correctly.

    Args:
        metrics_csv: Path to metrics CSV file

    Returns:
        True if valid, False otherwise
    """
    if not os.path.exists(metrics_csv):
        print(f"Error: Metrics CSV not found: {metrics_csv}")
        return False

    try:
        df = pd.read_csv(metrics_csv)
        if len(df) == 0:
            print(f"Warning: Metrics CSV is empty: {metrics_csv}")
            return False

        print(f"Metrics CSV validated: {len(df)} entries")
        return True
    except Exception as e:
        print(f"Error reading metrics CSV: {e}")
        return False


def main():
    """Main entry point for BoltzGen helper script."""
    parser = argparse.ArgumentParser(
        description='Helper script for BoltzGen postprocessing'
    )
    parser.add_argument(
        '--output_folder',
        type=str,
        required=True,
        help='BoltzGen output folder'
    )
    parser.add_argument(
        '--budget',
        type=int,
        required=True,
        help='Budget parameter used in BoltzGen run'
    )
    parser.add_argument(
        '--collect_structures',
        action='store_true',
        help='Collect final structure files into a list'
    )
    parser.add_argument(
        '--validate',
        action='store_true',
        help='Validate output files were generated correctly'
    )

    args = parser.parse_args()

    # Define key paths
    final_designs_folder = os.path.join(
        args.output_folder,
        'final_ranked_designs',
        f'final_{args.budget}_designs'
    )
    metrics_csv = os.path.join(
        args.output_folder,
        'final_ranked_designs',
        f'final_designs_metrics_{args.budget}.csv'
    )
    structures_list = os.path.join(
        args.output_folder,
        'final_structures.txt'
    )

    success = True

    # Collect structures if requested
    if args.collect_structures:
        print("Collecting final structures...")
        structures = collect_final_structures(final_designs_folder, structures_list)
        if not structures:
            print("Warning: No structures found")
            success = False

    # Validate outputs if requested
    if args.validate:
        print("Validating outputs...")
        if not validate_metrics_csv(metrics_csv):
            success = False

        if not os.path.exists(final_designs_folder):
            print(f"Error: Final designs folder not found: {final_designs_folder}")
            success = False

    if success:
        print("BoltzGen postprocessing completed successfully")
        return 0
    else:
        print("BoltzGen postprocessing completed with warnings/errors")
        return 1


if __name__ == "__main__":
    sys.exit(main())
