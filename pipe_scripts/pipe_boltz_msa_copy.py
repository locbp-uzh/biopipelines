#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Copy MSA files from a previous run for MSA recycling.

Reads an msas map_table and copies the referenced MSA files to the output folder, each named after the id of the sequence it belongs to.
With --queries-csv only the queried ids are copied, each found through the shared MsaLookup, and an id with no alignment is recorded as a failure row of --missing-csv.
"""

import argparse
import os
import shutil
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from biopipelines.biopipelines_io import MISSING_COLUMNS, step_id_from_table_path  # noqa: E402
from biopipelines.id_map_utils import MsaLookup  # noqa: E402


def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Copy MSA files from a table to output folder'
    )
    parser.add_argument(
        '--msa-table', required=True,
        help='Path to the msas map_table CSV (columns: id, file)'
    )
    parser.add_argument(
        '--output-folder', required=True,
        help='Output folder to copy MSA files to'
    )
    parser.add_argument(
        '--queries-csv', default=None,
        help='Optional queries CSV. When given, only the alignments of the ids in its id column '
             'are copied, matched to them by the shared MsaLookup.'
    )
    parser.add_argument(
        '--missing-csv', default=None,
        help='With --queries-csv: the missing table that gets a failure row for each queried id without an alignment.'
    )
    return parser.parse_args()


def read_queries(queries_csv: str) -> list:
    """(id, sequence) of every query, from the queries CSV's id and sequence columns."""
    if not os.path.exists(queries_csv):
        print(f"Error: queries CSV not found: {queries_csv}")
        sys.exit(1)

    try:
        df = pd.read_csv(queries_csv, dtype=str, keep_default_na=False)
    except Exception as e:
        print(f"Error reading queries CSV: {e}")
        sys.exit(1)

    if 'id' not in df.columns:
        print(f"Error: queries CSV missing 'id' column: {queries_csv}")
        sys.exit(1)

    sequences = df['sequence'] if 'sequence' in df.columns else [''] * len(df)
    return list(zip(df['id'].tolist(), sequences))


def copy_one(msa_file: str, seq_id: str, output_folder: str) -> bool:
    """Copy one alignment to <output_folder>/<seq_id><ext>; False when it cannot be copied."""
    if not msa_file or not os.path.exists(msa_file):
        print(f"Warning: MSA file not found: {msa_file}")
        return False
    dest_file = os.path.join(output_folder, f"{seq_id}{os.path.splitext(msa_file)[1]}")
    try:
        shutil.copy2(msa_file, dest_file)
    except Exception as e:
        print(f"Error copying {msa_file}: {e}")
        return False
    print(f"Copied MSA: {os.path.basename(msa_file)} -> {os.path.basename(dest_file)}")
    return True


def copy_msa_files(msa_table_path: str, output_folder: str,
                   queries_csv: str = None, missing_csv: str = None) -> int:
    """
    Copy MSA files from table to output folder.

    Args:
        msa_table_path: Path to the msas map_table
        output_folder: Destination folder for MSA files
        queries_csv: Optional queries CSV; restricts the copy to its ids
        missing_csv: With queries_csv, where the ids without an alignment are written

    Returns:
        Number of MSA files successfully copied
    """
    if not os.path.exists(msa_table_path):
        print(f"Error: MSA table not found: {msa_table_path}")
        sys.exit(1)

    os.makedirs(output_folder, exist_ok=True)

    try:
        lookup = MsaLookup.load(msa_table_path)
    except ValueError as e:
        print(f"Error: {e}")
        sys.exit(1)

    copied_count = 0
    if queries_csv:
        # An msas stream is often wider than the query set — a filtered stream, or a
        # whole previous run recycled for one sequence. Copying all of it is wasted
        # IO, and the extra a3m files sit in the same folder the predictor scans.
        missing = {}
        queries = read_queries(queries_csv)
        for seq_id, sequence in queries:
            msa_file = lookup.path(seq_id, sequence)
            if msa_file is None:
                missing[seq_id] = lookup.cause(seq_id)
            elif copy_one(msa_file, seq_id, output_folder):
                copied_count += 1
            else:
                missing[seq_id] = f"the alignment {msa_file} of {seq_id} does not exist"
        for seq_id, cause in missing.items():
            print(f"Warning: {cause}")
        if missing_csv:
            os.makedirs(os.path.dirname(os.path.abspath(missing_csv)), exist_ok=True)
            step_id = step_id_from_table_path(missing_csv)
            pd.DataFrame([{"id": seq_id, "removed_by": step_id, "kind": "failure", "cause": cause}
                          for seq_id, cause in missing.items()],
                         columns=MISSING_COLUMNS).to_csv(missing_csv, index=False)
        print(f"{len(queries) - len(missing)} of {len(queries)} queried id(s) have an alignment")
    else:
        for seq_id, msa_file in lookup.entries:
            if copy_one(msa_file, seq_id, output_folder):
                copied_count += 1

    return copied_count


def main():
    args = parse_arguments()

    print(f"Copying MSA files from: {args.msa_table}")
    print(f"To folder: {args.output_folder}")
    if args.queries_csv:
        print(f"Restricted to the ids in: {args.queries_csv}")

    copied = copy_msa_files(args.msa_table, args.output_folder, args.queries_csv, args.missing_csv)

    print(f"\nSuccessfully copied {copied} MSA files")

    if copied == 0:
        print("Warning: No MSA files were copied")


if __name__ == "__main__":
    main()
