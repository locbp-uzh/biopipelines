#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Copy MSA files from a previous Boltz2 run for MSA recycling.

Reads an MSA table CSV and copies the referenced MSA files to the output folder.
"""

import argparse
import os
import shutil
import sys

import pandas as pd


def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Copy MSA files from a table to output folder'
    )
    parser.add_argument(
        '--msa-table', required=True,
        help='Path to MSA table CSV file (columns: id, sequences.id, msa_file)'
    )
    parser.add_argument(
        '--output-folder', required=True,
        help='Output folder to copy MSA files to'
    )
    return parser.parse_args()


def copy_msa_files(msa_table_path: str, output_folder: str) -> int:
    """
    Copy MSA files from table to output folder.

    Args:
        msa_table_path: Path to MSA table CSV
        output_folder: Destination folder for MSA files

    Returns:
        Number of MSA files successfully copied
    """
    if not os.path.exists(msa_table_path):
        print(f"Error: MSA table not found: {msa_table_path}")
        sys.exit(1)

    os.makedirs(output_folder, exist_ok=True)

    try:
        df = pd.read_csv(msa_table_path)
    except Exception as e:
        print(f"Error reading MSA table: {e}")
        sys.exit(1)

    if 'msa_file' not in df.columns:
        print(f"Error: MSA table missing 'msa_file' column")
        sys.exit(1)

    copied_count = 0
    for _, row in df.iterrows():
        msa_file = row.get('msa_file', '')
        seq_id = row.get('sequences.id', row.get('id', ''))

        if not msa_file or not os.path.exists(msa_file):
            print(f"Warning: MSA file not found: {msa_file}")
            continue

        # Determine destination filename
        ext = os.path.splitext(msa_file)[1]
        dest_file = os.path.join(output_folder, f"{seq_id}{ext}")

        try:
            shutil.copy2(msa_file, dest_file)
            print(f"Copied MSA: {os.path.basename(msa_file)} -> {os.path.basename(dest_file)}")
            copied_count += 1
        except Exception as e:
            print(f"Error copying {msa_file}: {e}")

    return copied_count


def main():
    args = parse_arguments()

    print(f"Copying MSA files from: {args.msa_table}")
    print(f"To folder: {args.output_folder}")

    copied = copy_msa_files(args.msa_table, args.output_folder)

    print(f"\nSuccessfully copied {copied} MSA files")

    if copied == 0:
        print("Warning: No MSA files were copied")


if __name__ == "__main__":
    main()
