#!/usr/bin/env python3
"""
Boltz2 completion check with missing sequence filtering.

Checks that all expected output structures exist, excluding sequences
that were marked as missing by upstream tools.
"""

import argparse
import os
import sys

import pandas as pd


def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Check Boltz2 completion with missing sequence filtering'
    )
    parser.add_argument(
        '--output-folder', required=True,
        help='Boltz2 output folder'
    )
    parser.add_argument(
        '--sequence-ids-file', required=True,
        help='CSV file with expected sequence IDs'
    )
    parser.add_argument(
        '--missing-file',
        help='Path to upstream missing.csv (optional)'
    )
    parser.add_argument(
        '--output-format', default='pdb',
        choices=['pdb', 'mmcif'],
        help='Output format (pdb or mmcif)'
    )
    return parser.parse_args()


def load_sequence_ids(sequence_ids_file: str) -> set:
    """Load expected sequence IDs from CSV file."""
    if not os.path.exists(sequence_ids_file):
        print(f"Error: Sequence IDs file not found: {sequence_ids_file}")
        sys.exit(1)

    try:
        df = pd.read_csv(sequence_ids_file)
        if 'id' not in df.columns:
            print(f"Error: No 'id' column in {sequence_ids_file}")
            sys.exit(1)
        return set(df['id'].astype(str).tolist())
    except Exception as e:
        print(f"Error reading sequence IDs file: {e}")
        sys.exit(1)


def load_missing_ids(missing_file: str) -> set:
    """Load missing sequence IDs from CSV file."""
    if not missing_file or not os.path.exists(missing_file):
        return set()

    try:
        df = pd.read_csv(missing_file)
        if 'id' not in df.columns:
            print(f"Warning: No 'id' column in {missing_file}")
            return set()
        return set(df['id'].astype(str).tolist())
    except Exception as e:
        print(f"Warning: Could not read missing file: {e}")
        return set()


def check_completion(
    output_folder: str,
    all_ids: set,
    missing_ids: set,
    output_format: str
) -> tuple:
    """
    Check if all expected output files exist.

    Args:
        output_folder: Boltz2 output folder
        all_ids: All expected sequence IDs
        missing_ids: IDs to exclude (from upstream missing table)
        output_format: Output format (pdb or mmcif)

    Returns:
        Tuple of (success: bool, missing_files: list)
    """
    expected_ids = all_ids - missing_ids

    print(f"Total sequence IDs: {len(all_ids)}")
    print(f"Missing from upstream: {len(missing_ids)}")
    print(f"Expected outputs: {len(expected_ids)}")

    if missing_ids:
        print(f"Excluded IDs: {sorted(missing_ids)}")

    ext = ".pdb" if output_format == "pdb" else ".cif"
    missing_files = []

    for seq_id in expected_ids:
        expected_file = os.path.join(output_folder, f"{seq_id}{ext}")
        if not os.path.exists(expected_file):
            missing_files.append(expected_file)

    if missing_files:
        print(f"\nMissing {len(missing_files)} output files:")
        for f in missing_files[:10]:
            print(f"  - {f}")
        if len(missing_files) > 10:
            print(f"  ... and {len(missing_files) - 10} more")
        return False, missing_files

    print(f"\nAll {len(expected_ids)} expected outputs found")
    return True, []


def main():
    args = parse_arguments()

    print("Boltz2 Completion Check")
    print("=" * 50)

    # Load sequence IDs
    all_ids = load_sequence_ids(args.sequence_ids_file)

    # Load missing IDs (if provided)
    missing_ids = load_missing_ids(args.missing_file)

    # Check completion
    success, missing_files = check_completion(
        args.output_folder,
        all_ids,
        missing_ids,
        args.output_format
    )

    if success:
        print("\nBoltz2 completed successfully")
        sys.exit(0)
    else:
        print("\nBoltz2 failed - some outputs missing")
        sys.exit(1)


if __name__ == "__main__":
    main()
