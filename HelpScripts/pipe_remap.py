#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Execution-time script for ReMap tool.

Reads a remap_config.json and:
1. Creates symlinks (or copies on Windows) for stream files with new IDs
2. Rewrites map_table CSVs with remapped IDs in 'id' and all '*.id' provenance columns
3. Rewrites table CSVs with remapped IDs in 'id' column

Usage:
    python pipe_remap.py <config_json>
"""

import sys
import os
import json
import shutil
import platform
import pandas as pd


def create_link_or_copy(source_path, link_path):
    """Create a symlink or copy if symlinks are not supported."""
    # Ensure parent directory exists
    os.makedirs(os.path.dirname(link_path), exist_ok=True)

    # Remove existing target if present
    if os.path.exists(link_path) or os.path.islink(link_path):
        os.remove(link_path)

    if platform.system() == "Windows":
        # On Windows, symlinks may require elevated privileges; fall back to copy
        try:
            os.symlink(os.path.abspath(source_path), link_path)
        except OSError:
            shutil.copy2(source_path, link_path)
    else:
        # Use relative symlinks on Unix for portability
        rel_path = os.path.relpath(source_path, os.path.dirname(link_path))
        os.symlink(rel_path, link_path)


def remap_csv_ids(csv_path, output_path, id_mapping):
    """
    Remap IDs in a CSV file.

    Replaces values in the 'id' column and all '*.id' provenance columns
    according to the id_mapping dict.
    """
    if not os.path.exists(csv_path):
        print(f"  Warning: CSV not found, skipping: {csv_path}")
        return

    df = pd.read_csv(csv_path)

    # Remap 'id' column
    if 'id' in df.columns:
        df['id'] = df['id'].map(lambda x: id_mapping.get(str(x), x))

    # Remap all provenance columns (*.id pattern)
    for col in df.columns:
        if col.endswith('.id') and col != 'id':
            df[col] = df[col].map(lambda x: id_mapping.get(str(x), x))

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    df.to_csv(output_path, index=False)
    print(f"  Wrote remapped CSV: {output_path}")


def main():
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <config_json>")
        sys.exit(1)

    config_path = sys.argv[1]
    with open(config_path, 'r') as f:
        config = json.load(f)

    id_mapping = config["id_mapping"]
    output_folder = config["output_folder"]
    streams = config.get("streams", [])
    tables = config.get("tables", [])
    map_tables = config.get("map_tables", [])

    os.makedirs(output_folder, exist_ok=True)

    print(f"ReMap: {len(id_mapping)} ID mappings")
    for old_id, new_id in id_mapping.items():
        print(f"  {old_id} -> {new_id}")

    # 1. Create symlinks/copies for stream files
    for stream in streams:
        stream_name = stream.get("name", "unknown")
        old_ids = stream.get("ids", [])
        old_files = stream.get("files", [])
        fmt = stream.get("format", "")

        if not old_files:
            print(f"  Stream '{stream_name}': no files (value-based), skipping symlinks")
            continue

        print(f"  Stream '{stream_name}': creating links for {len(old_files)} files")
        for old_id, old_file in zip(old_ids, old_files):
            new_id = id_mapping.get(old_id, old_id)
            ext = os.path.splitext(old_file)[1]
            new_file = os.path.join(output_folder, f"{new_id}{ext}")

            if os.path.exists(old_file):
                create_link_or_copy(old_file, new_file)
                print(f"    {new_id}{ext} -> {old_file}")
            else:
                print(f"    Warning: source file not found: {old_file}")

    # 2. Remap map_table CSVs
    for stream in streams:
        map_table_path = stream.get("map_table", "")
        if not map_table_path:
            continue
        stream_name = stream.get("name", "unknown")
        output_map = os.path.join(output_folder, f"{stream_name}_map.csv")
        print(f"  Remapping map_table: {stream_name}")
        remap_csv_ids(map_table_path, output_map, id_mapping)

    # 3. Remap source table CSVs
    for table in tables:
        table_name = table.get("name", "unknown")
        table_path = table.get("path", "")
        if not table_path:
            continue
        output_table = os.path.join(output_folder, f"{table_name}.csv")
        print(f"  Remapping table: {table_name}")
        remap_csv_ids(table_path, output_table, id_mapping)

    print("ReMap complete.")


if __name__ == "__main__":
    main()
