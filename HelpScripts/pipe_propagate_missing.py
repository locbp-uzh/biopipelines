#!/usr/bin/env python3
# Copyright (c) 2026 Gianluca Quargnali @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# This software is freely available for use, modification, and redistribution.
# If you use this software or any derivative of it in published work,
# you must cite the original author and this repository.
"""
Generic helper to propagate missing.csv from upstream tools to downstream tools.

This script reads the missing table from upstream tools and creates a new missing.csv
with updated file paths for the current tool's expected outputs.
"""

import os
import sys
import pandas as pd
from typing import List, Dict, Optional


def find_upstream_missing_csv(input_folders: List[str]) -> Optional[str]:
    """
    Find missing.csv from upstream tool folders.

    Args:
        input_folders: List of potential upstream output folders

    Returns:
        Path to upstream missing.csv, or None if not found
    """
    for folder in input_folders:
        if not folder or not os.path.exists(folder):
            continue

        missing_csv = os.path.join(folder, "missing.csv")
        if os.path.exists(missing_csv):
            return missing_csv

    return None


def load_missing_ids(missing_csv_path: str) -> List[str]:
    """
    Load missing IDs from upstream missing.csv.

    Args:
        missing_csv_path: Path to upstream missing.csv

    Returns:
        List of missing IDs
    """
    if not os.path.exists(missing_csv_path):
        return []

    try:
        df = pd.read_csv(missing_csv_path)
        if 'id' in df.columns:
            return df['id'].astype(str).tolist()
    except Exception as e:
        print(f"Warning: Could not read missing.csv: {e}")

    return []


def propagate_missing_table(
    upstream_folders: List[str],
    output_folder: str,
    file_extensions: Dict[str, str] = None
) -> List[str]:
    """
    Propagate missing table from upstream tools with updated file paths.

    Args:
        upstream_folders: List of potential upstream output folders
        output_folder: Current tool's output folder
        file_extensions: Dictionary mapping column names to file extensions
                        e.g., {"structure": ".pdb", "msa": ".csv"}

    Returns:
        List of missing IDs that were propagated
    """
    # Default file extensions
    if file_extensions is None:
        file_extensions = {
            "structure": ".pdb",
            "msa": ".csv"
        }

    # Find upstream missing.csv
    upstream_missing_csv = find_upstream_missing_csv(upstream_folders)

    if not upstream_missing_csv:
        print("No upstream missing.csv found - creating empty missing.csv")
        output_missing_csv = os.path.join(output_folder, "missing.csv")
        pd.DataFrame(columns=['id'] + list(file_extensions.keys())).to_csv(output_missing_csv, index=False)
        return []

    print(f"Found upstream missing.csv: {upstream_missing_csv}")

    # Load missing IDs
    missing_ids = load_missing_ids(upstream_missing_csv)

    if not missing_ids:
        print("No missing IDs in upstream table - creating empty missing.csv")
        output_missing_csv = os.path.join(output_folder, "missing.csv")
        pd.DataFrame(columns=['id'] + list(file_extensions.keys())).to_csv(output_missing_csv, index=False)
        return []

    print(f"Propagating {len(missing_ids)} missing IDs to current tool")

    # Create missing table with updated paths for current tool
    missing_data = []
    for missing_id in missing_ids:
        row = {"id": missing_id}

        # Add file paths for each column type
        for col_name, extension in file_extensions.items():
            file_path = os.path.join(output_folder, f"{missing_id}{extension}")
            row[col_name] = file_path

        missing_data.append(row)

    # Save missing.csv in current tool's output folder
    missing_df = pd.DataFrame(missing_data)
    output_missing_csv = os.path.join(output_folder, "missing.csv")
    missing_df.to_csv(output_missing_csv, index=False)

    print(f"Created missing.csv: {output_missing_csv}")
    print(f"Propagated {len(missing_ids)} missing IDs")

    return missing_ids


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Propagate missing.csv from upstream tools")
    parser.add_argument("--upstream-folders", nargs="+", required=True,
                       help="Upstream tool output folders")
    parser.add_argument("--output-folder", required=True,
                       help="Current tool's output folder")
    parser.add_argument("--structure-ext", default=".pdb",
                       help="File extension for structure files (default: .pdb)")
    parser.add_argument("--msa-ext", default=".csv",
                       help="File extension for MSA files (default: .csv)")

    args = parser.parse_args()

    file_extensions = {
        "structure": args.structure_ext,
        "msa": args.msa_ext
    }

    missing_ids = propagate_missing_table(
        upstream_folders=args.upstream_folders,
        output_folder=args.output_folder,
        file_extensions=file_extensions
    )

    if missing_ids:
        print(f"\nSuccess: Propagated {len(missing_ids)} missing IDs")
        sys.exit(0)
    else:
        print("\nNo missing IDs to propagate")
        sys.exit(0)
