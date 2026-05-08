#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Restore original IDs in output files and CSVs after external tools sanitize them.

Some external tools (e.g., ColabFold) replace characters like '+' with '_' in
output filenames. This script builds a mapping from sanitized IDs back to the
original predicted IDs, then renames files and fixes CSV id columns accordingly.

Usage:
    python pipe_unsanitize_ids.py \\
        --predicted-ids /path/to/structures_map.csv \\
        --rename-files /path/to/output pdb \\
        --rename-files /path/to/msas a3m \\
        --fix-csv /path/to/confidence.csv \\
        --fix-csv /path/to/msas.csv
"""

import argparse
import os
import glob
import pandas as pd
from typing import Dict, List, Tuple


# Characters that external tools may replace
SANITIZE_MAP = {'+': '_'}


def sanitize_id(original_id: str) -> str:
    """Apply the same sanitization that external tools do."""
    result = original_id
    for char, replacement in SANITIZE_MAP.items():
        result = result.replace(char, replacement)
    return result


def build_unsanitize_map(predicted_ids: List[str]) -> Dict[str, str]:
    """
    Build a mapping from sanitized IDs back to original IDs.

    Only includes entries where sanitization actually changes the ID.

    Args:
        predicted_ids: List of original predicted IDs

    Returns:
        Dict mapping sanitized_id -> original_id
    """
    mapping = {}
    for original_id in predicted_ids:
        sanitized = sanitize_id(original_id)
        if sanitized != original_id:
            mapping[sanitized] = original_id
    return mapping


def load_predicted_ids(csv_path: str) -> List[str]:
    """Load predicted IDs from a CSV file with an 'id' column."""
    if not os.path.exists(csv_path):
        print(f"Warning: Predicted IDs file not found: {csv_path}")
        return []
    df = pd.read_csv(csv_path)
    if 'id' not in df.columns:
        print(f"Warning: No 'id' column in {csv_path}")
        return []
    return df['id'].astype(str).tolist()


def rename_files(folder: str, extension: str, unsanitize_map: Dict[str, str]) -> int:
    """
    Rename files in a folder from sanitized IDs back to original IDs.

    Args:
        folder: Directory containing the files
        extension: File extension (without dot)
        unsanitize_map: Mapping from sanitized_id -> original_id

    Returns:
        Number of files renamed
    """
    if not os.path.isdir(folder):
        print(f"Warning: Folder not found: {folder}")
        return 0

    renamed = 0
    pattern = os.path.join(folder, f"*.{extension}")
    for fpath in glob.glob(pattern):
        basename = os.path.basename(fpath)
        file_id = os.path.splitext(basename)[0]

        if file_id in unsanitize_map:
            original_id = unsanitize_map[file_id]
            new_path = os.path.join(folder, f"{original_id}.{extension}")
            os.rename(fpath, new_path)
            print(f"  Renamed: {basename} -> {original_id}.{extension}")
            renamed += 1

    return renamed


def fix_csv(csv_path: str, unsanitize_map: Dict[str, str]) -> int:
    """
    Fix sanitized IDs in a CSV file's 'id' column.

    Args:
        csv_path: Path to CSV file
        unsanitize_map: Mapping from sanitized_id -> original_id

    Returns:
        Number of IDs fixed
    """
    if not os.path.exists(csv_path):
        print(f"Warning: CSV not found: {csv_path}")
        return 0

    df = pd.read_csv(csv_path)
    if 'id' not in df.columns:
        print(f"Warning: No 'id' column in {csv_path}")
        return 0

    fixed = 0
    original_ids = df['id'].astype(str).tolist()
    new_ids = []
    for oid in original_ids:
        if oid in unsanitize_map:
            new_ids.append(unsanitize_map[oid])
            fixed += 1
        else:
            new_ids.append(oid)

    if fixed > 0:
        df['id'] = new_ids
        df.to_csv(csv_path, index=False)
        print(f"  Fixed {fixed} IDs in {os.path.basename(csv_path)}")

    return fixed


def main():
    parser = argparse.ArgumentParser(
        description='Restore original IDs after external tool sanitization'
    )
    parser.add_argument('--predicted-ids', required=True,
                        help='CSV file with predicted IDs (must have "id" column)')
    parser.add_argument('--rename-files', nargs=2, action='append', default=[],
                        metavar=('FOLDER', 'EXT'),
                        help='Rename files in FOLDER with extension EXT (can be repeated)')
    parser.add_argument('--fix-csv', action='append', default=[],
                        metavar='CSV_PATH',
                        help='Fix IDs in a CSV file (can be repeated)')

    args = parser.parse_args()

    # Load predicted IDs and build mapping
    predicted_ids = load_predicted_ids(args.predicted_ids)
    if not predicted_ids:
        print("No predicted IDs found, nothing to unsanitize")
        return

    unsanitize_map = build_unsanitize_map(predicted_ids)
    if not unsanitize_map:
        print("No IDs require unsanitization (no '+' characters found)")
        return

    print(f"Built unsanitize map for {len(unsanitize_map)} IDs")

    # Rename files
    total_renamed = 0
    for folder, ext in args.rename_files:
        print(f"Renaming .{ext} files in {folder}")
        total_renamed += rename_files(folder, ext, unsanitize_map)

    # Fix CSVs
    total_fixed = 0
    for csv_path in args.fix_csv:
        print(f"Fixing IDs in {csv_path}")
        total_fixed += fix_csv(csv_path, unsanitize_map)

    print(f"\nUnsanitize summary: {total_renamed} files renamed, {total_fixed} CSV IDs fixed")


if __name__ == "__main__":
    main()
