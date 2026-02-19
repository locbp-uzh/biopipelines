#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Runtime helper script for Load filtering.

Filters CSV tables based on filtered IDs and saves them to Load's output folder.
"""

import os
import sys
import argparse
import json
import pandas as pd
from typing import Dict, List, Any, Set


def filter_csv_file(input_path: str, output_path: str, filtered_ids: Set[str]) -> int:
    """
    Filter a CSV file by IDs and save to output path.

    Args:
        input_path: Original CSV file path
        output_path: Filtered CSV file path
        filtered_ids: Set of IDs to keep

    Returns:
        Number of rows in filtered CSV
    """
    try:
        df = pd.read_csv(input_path)

        if 'id' in df.columns:
            filtered_df = df[df['id'].isin(filtered_ids)]
        else:
            # No 'id' column, copy entire file
            print(f"  Warning: No 'id' column in {os.path.basename(input_path)}, copying as-is")
            filtered_df = df

        # Ensure output directory exists
        os.makedirs(os.path.dirname(output_path), exist_ok=True)

        # Save filtered CSV
        filtered_df.to_csv(output_path, index=False)

        print(f"  ✓ {os.path.basename(output_path)}: {len(filtered_df)}/{len(df)} rows")
        return len(filtered_df)

    except Exception as e:
        print(f"  Error filtering {input_path}: {e}", file=sys.stderr)
        raise


def filter_tables(config_data: Dict[str, Any]) -> int:
    """
    Filter all CSV tables based on configuration.

    Args:
        config_data: Configuration dictionary with filtered_ids and output_structure

    Returns:
        0 if successful, 1 if failed
    """
    filtered_ids = set(config_data['filtered_ids'])
    output_structure = config_data['output_structure']
    output_folder = config_data['output_folder']

    print(f"Filtering tables with {len(filtered_ids)} IDs")
    print(f"Output folder: {output_folder}")

    os.makedirs(output_folder, exist_ok=True)

    # Filter tables
    if 'tables' in output_structure:
        tables = output_structure['tables']
        if isinstance(tables, dict):
            print("\nFiltering tables:")
            for ds_name, ds_info in tables.items():
                # Get original path
                original_path = None
                if isinstance(ds_info, dict) and 'path' in ds_info:
                    original_path = ds_info['path']
                elif isinstance(ds_info, str):
                    original_path = ds_info

                if not original_path or not os.path.exists(original_path):
                    print(f"  Warning: Skipping '{ds_name}' - file not found: {original_path}")
                    continue

                # Filter and save
                filtered_path = os.path.join(output_folder, f"{ds_name}.csv")
                filter_csv_file(original_path, filtered_path, filtered_ids)

    # Filter special CSV lists (compounds, sequences, msas)
    for list_name in ['compounds', 'sequences', 'msas']:
        if list_name in output_structure:
            file_list = output_structure[list_name]
            if isinstance(file_list, list):
                print(f"\nFiltering {list_name}:")
                for file_path in file_list:
                    if isinstance(file_path, str) and file_path.endswith('.csv') and os.path.exists(file_path):
                        filtered_path = os.path.join(output_folder, os.path.basename(file_path))
                        filter_csv_file(file_path, filtered_path, filtered_ids)

    print("\nFiltering complete!")
    return 0


def main():
    parser = argparse.ArgumentParser(description='Filter Load tables')
    parser.add_argument('--config', required=True, help='JSON config file with filter parameters')

    args = parser.parse_args()

    # Load configuration
    if not os.path.exists(args.config):
        print(f"Error: Config file not found: {args.config}", file=sys.stderr)
        sys.exit(1)

    try:
        with open(args.config, 'r') as f:
            config_data = json.load(f)
    except Exception as e:
        print(f"Error loading config: {e}", file=sys.stderr)
        sys.exit(1)

    # Filter tables
    exit_code = filter_tables(config_data)
    sys.exit(exit_code)


if __name__ == "__main__":
    main()
