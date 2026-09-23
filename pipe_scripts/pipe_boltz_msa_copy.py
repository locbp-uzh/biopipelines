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
        help='Path to MSA table CSV file (columns: id, sequences.id, file)'
    )
    parser.add_argument(
        '--output-folder', required=True,
        help='Output folder to copy MSA files to'
    )
    parser.add_argument(
        '--queries-csv', default=None,
        help='Optional queries CSV. When given, only MSAs whose destination name '
             'appears in its id column are copied.'
    )
    return parser.parse_args()


def read_wanted_ids(queries_csv: str) -> set:
    """Ids being predicted, from the queries CSV's id column."""
    if not os.path.exists(queries_csv):
        print(f"Error: queries CSV not found: {queries_csv}")
        sys.exit(1)

    try:
        df = pd.read_csv(queries_csv)
    except Exception as e:
        print(f"Error reading queries CSV: {e}")
        sys.exit(1)

    if 'id' not in df.columns:
        print(f"Error: queries CSV missing 'id' column: {queries_csv}")
        sys.exit(1)

    return {str(v) for v in df['id'].tolist()}


def copy_msa_files(msa_table_path: str, output_folder: str,
                   queries_csv: str = None) -> int:
    """
    Copy MSA files from table to output folder.

    Args:
        msa_table_path: Path to MSA table CSV
        output_folder: Destination folder for MSA files
        queries_csv: Optional queries CSV; restricts the copy to its ids

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

    if 'file' not in df.columns:
        print(f"Error: MSA table missing 'file' column")
        sys.exit(1)

    # An msas stream is often wider than the query set — a filtered stream, or a
    # whole previous run recycled for one sequence. Copying all of it is wasted
    # IO, and the extra a3m files sit in the same folder the predictor scans.
    wanted = read_wanted_ids(queries_csv) if queries_csv else None

    copied_count = 0
    skipped_count = 0
    for _, row in df.iterrows():
        msa_file = row.get('file', '')
        seq_id = row.get('sequences.id', row.get('id', ''))

        # Match on the destination name, which is what the predictor looks up.
        if wanted is not None and str(seq_id) not in wanted:
            skipped_count += 1
            continue

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

    if skipped_count:
        print(f"Skipped {skipped_count} MSA(s) not among the queried ids")

    return copied_count


def main():
    args = parse_arguments()

    print(f"Copying MSA files from: {args.msa_table}")
    print(f"To folder: {args.output_folder}")
    if args.queries_csv:
        print(f"Restricted to the ids in: {args.queries_csv}")

    copied = copy_msa_files(args.msa_table, args.output_folder, args.queries_csv)

    print(f"\nSuccessfully copied {copied} MSA files")

    if copied == 0:
        print("Warning: No MSA files were copied")


if __name__ == "__main__":
    main()
