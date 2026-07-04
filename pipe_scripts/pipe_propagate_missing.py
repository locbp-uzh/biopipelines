#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Generic helper to propagate missing.csv from upstream tools to downstream tools.

Merges the upstream missing tables and writes the current tool's own
missing.csv (schema: id, removed_by, kind, cause), in the UPSTREAM (input-axis)
id space. Mapping those ids into the current tool's output id space (for
products like ``prot+lig2``, ``_N`` multipliers, or group keys) is done once,
centrally, by ``pipe_check_completion`` when it excuses missing files — so this
helper stays a uniform raw merge for every tool.
"""

import os
import sys
import pandas as pd
from typing import List, Optional


def find_upstream_missing_csvs(input_folders: List[str]) -> List[str]:
    """
    Find every upstream missing.csv across the given folders.

    Returns the manifests for ALL input axes so a filter on more than one axis
    propagates fully (de-duped per id downstream). Each folder is the directory
    holding ``missing.csv`` (a tool's ``tables/`` dir).
    """
    found = []
    seen = set()
    for folder in input_folders:
        if not folder or not os.path.exists(folder):
            continue
        missing_csv = os.path.join(folder, "missing.csv")
        if os.path.exists(missing_csv) and missing_csv not in seen:
            seen.add(missing_csv)
            found.append(missing_csv)
    return found


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
    missing_csv_path: Optional[str] = None,
    local_missing: Optional[str] = None,
) -> List[str]:
    """Merge upstream missing tables (rows as-is) into this tool's missing.csv.

    Reads every upstream missing.csv across ``upstream_folders`` plus an optional
    ``local_missing`` file this tool produced itself, merges them, and de-dupes
    by id (last wins). Rows stay in the upstream (input-axis) id space; the
    completion check maps them to output ids when excusing files. Always writes
    the destination — even empty — so the completion check finds the path it
    expects.
    """
    columns = ['id', 'removed_by', 'kind', 'cause']
    output_missing_csv = missing_csv_path or os.path.join(output_folder, "missing.csv")
    os.makedirs(os.path.dirname(output_missing_csv), exist_ok=True)

    source_csvs = find_upstream_missing_csvs(upstream_folders)
    if local_missing and os.path.exists(local_missing):
        source_csvs.append(local_missing)
    dfs = []
    for path in source_csvs:
        try:
            df = pd.read_csv(path)
        except Exception as e:
            print(f"Warning: Could not read missing.csv {path}: {e}")
            continue
        if not df.empty:
            dfs.append(df)
            print(f"  Merging {len(df)} missing entries from {path}")

    if dfs:
        merged = pd.concat(dfs, ignore_index=True)
        if 'id' in merged.columns:
            merged = merged.drop_duplicates(subset=['id'], keep='last')
    else:
        merged = pd.DataFrame(columns=columns)

    merged.to_csv(output_missing_csv, index=False)
    print(f"Wrote missing.csv with {len(merged)} entries: {output_missing_csv}")

    return merged['id'].astype(str).tolist() if 'id' in merged.columns else []


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Propagate missing.csv from upstream tools")
    parser.add_argument("--upstream-folders", nargs="+", required=True,
                       help="Upstream tool output folders")
    parser.add_argument("--output-folder", required=True,
                       help="Current tool's output folder")
    parser.add_argument("--missing-csv", default=None,
                       help="Explicit destination for missing.csv")
    parser.add_argument("--local-missing", default=None,
                       help="A missing.csv this tool produced itself, merged in alongside upstream")

    args = parser.parse_args()

    missing_ids = propagate_missing_table(
        upstream_folders=args.upstream_folders,
        output_folder=args.output_folder,
        missing_csv_path=args.missing_csv,
        local_missing=args.local_missing,
    )

    if missing_ids:
        print(f"\nSuccess: Propagated {len(missing_ids)} missing IDs")
        sys.exit(0)
    else:
        print("\nNo missing IDs to propagate")
        sys.exit(0)
