#!/usr/bin/env python3
# Copyright (c) 2026 Gianluca Quargnali @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Generic helper to propagate missing.csv from upstream tools to downstream tools.

This script copies the missing table from upstream tools as-is to the current
tool's output folder. The missing.csv format is: id, removed_by, cause.
"""

import os
import sys
import pandas as pd
from typing import List, Optional


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
    output_folder: str
) -> List[str]:
    """
    Propagate missing table from upstream tools by copying rows as-is.

    Args:
        upstream_folders: List of potential upstream output folders
        output_folder: Current tool's output folder

    Returns:
        List of missing IDs that were propagated
    """
    output_missing_csv = os.path.join(output_folder, "missing.csv")

    # Find upstream missing.csv
    upstream_missing_csv = find_upstream_missing_csv(upstream_folders)

    if not upstream_missing_csv:
        print("No upstream missing.csv found - creating empty missing.csv")
        pd.DataFrame(columns=['id', 'removed_by', 'cause']).to_csv(output_missing_csv, index=False)
        return []

    print(f"Found upstream missing.csv: {upstream_missing_csv}")

    # Load upstream missing table
    try:
        upstream_df = pd.read_csv(upstream_missing_csv)
    except Exception as e:
        print(f"Warning: Could not read upstream missing.csv: {e}")
        pd.DataFrame(columns=['id', 'removed_by', 'cause']).to_csv(output_missing_csv, index=False)
        return []

    if upstream_df.empty:
        print("No missing IDs in upstream table - creating empty missing.csv")
        pd.DataFrame(columns=['id', 'removed_by', 'cause']).to_csv(output_missing_csv, index=False)
        return []

    print(f"Propagating {len(upstream_df)} missing entries to current tool")

    # Copy upstream rows as-is (preserve original removed_by and cause)
    upstream_df.to_csv(output_missing_csv, index=False)

    print(f"Created missing.csv: {output_missing_csv}")
    print(f"Propagated {len(upstream_df)} missing entries")

    return upstream_df['id'].astype(str).tolist() if 'id' in upstream_df.columns else []


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Propagate missing.csv from upstream tools")
    parser.add_argument("--upstream-folders", nargs="+", required=True,
                       help="Upstream tool output folders")
    parser.add_argument("--output-folder", required=True,
                       help="Current tool's output folder")

    args = parser.parse_args()

    missing_ids = propagate_missing_table(
        upstream_folders=args.upstream_folders,
        output_folder=args.output_folder
    )

    if missing_ids:
        print(f"\nSuccess: Propagated {len(missing_ids)} missing IDs")
        sys.exit(0)
    else:
        print("\nNo missing IDs to propagate")
        sys.exit(0)
