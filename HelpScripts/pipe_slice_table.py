#!/usr/bin/env python3
"""
Runtime helper script for slicing datasheets to first N rows.

This script takes tool output datasheets and creates new output containing
only the first N rows from each datasheet, along with associated files.
"""

import os
import sys
import argparse
import json
import pandas as pd
import shutil
import glob
from typing import Dict, List, Any, Optional


def slice_datasheet(datasheet_path: str, n_rows: int, output_path: str) -> List[str]:
    """
    Slice a datasheet to first N rows and return the IDs.

    Args:
        datasheet_path: Path to input CSV file
        n_rows: Number of rows to keep from the beginning
        output_path: Path to save sliced CSV

    Returns:
        List of IDs from the sliced datasheet
    """
    print(f"Slicing datasheet: {datasheet_path}")

    if not os.path.exists(datasheet_path):
        print(f"Warning: Datasheet not found: {datasheet_path}")
        return []

    try:
        df = pd.read_csv(datasheet_path)
        print(f"Loaded datasheet: {df.shape}")

        if df.empty:
            print("Warning: Datasheet is empty")
            df.to_csv(output_path, index=False)
            return []

        # Slice to first n_rows
        sliced_df = df.head(n_rows)
        print(f"Sliced to: {sliced_df.shape}")

        # Create output directory if needed
        os.makedirs(os.path.dirname(output_path), exist_ok=True)

        # Save sliced datasheet
        sliced_df.to_csv(output_path, index=False)
        print(f"Saved sliced datasheet: {output_path}")

        # Extract IDs if available
        sliced_ids = []
        if 'id' in sliced_df.columns:
            sliced_ids = sliced_df['id'].astype(str).tolist()
        else:
            # Try other ID columns
            id_columns = [col for col in sliced_df.columns if 'id' in col.lower()]
            if id_columns:
                sliced_ids = sliced_df[id_columns[0]].astype(str).tolist()

        print(f"Extracted {len(sliced_ids)} IDs: {sliced_ids[:5]}{'...' if len(sliced_ids) > 5 else ''}")
        return sliced_ids

    except Exception as e:
        print(f"Error processing datasheet {datasheet_path}: {e}")
        return []


def copy_associated_files(sliced_ids: List[str], input_folder: str, output_folder: str) -> Dict[str, List[str]]:
    """
    Copy structure, compound, and sequence files for sliced IDs.

    Args:
        sliced_ids: List of IDs to copy files for
        input_folder: Input folder containing original files
        output_folder: Output folder to copy files to

    Returns:
        Dictionary mapping file type to list of copied file paths
    """
    copied_files = {"structures": [], "compounds": [], "sequences": []}

    print(f"\nCopying associated files for {len(sliced_ids)} IDs")
    print(f"Input folder: {input_folder}")
    print(f"Output folder: {output_folder}")

    os.makedirs(output_folder, exist_ok=True)

    for i, file_id in enumerate(sliced_ids):
        # Copy structures (PDB/CIF files)
        structure_patterns = [
            os.path.join(input_folder, f"{file_id}.pdb"),
            os.path.join(input_folder, f"*{file_id}*.pdb"),
            os.path.join(input_folder, f"{file_id}.cif"),
            os.path.join(input_folder, f"*{file_id}*.cif"),
            os.path.join(input_folder, "**", f"{file_id}.pdb"),
            os.path.join(input_folder, "**", f"*{file_id}*.pdb"),
        ]

        for pattern in structure_patterns:
            matches = glob.glob(pattern, recursive=True)
            if matches:
                source = matches[0]
                dest = os.path.join(output_folder, f"{file_id}.pdb")
                try:
                    shutil.copy2(source, dest)
                    copied_files["structures"].append(dest)
                    print(f"Copied structure: {os.path.basename(source)} -> {os.path.basename(dest)}")
                    break
                except Exception as e:
                    print(f"Warning: Could not copy structure {source}: {e}")

        # Copy compounds (SDF files)
        compound_patterns = [
            os.path.join(input_folder, f"{file_id}.sdf"),
            os.path.join(input_folder, f"*{file_id}*.sdf"),
            os.path.join(input_folder, "**", f"{file_id}.sdf"),
            os.path.join(input_folder, "**", f"*{file_id}*.sdf"),
        ]

        for pattern in compound_patterns:
            matches = glob.glob(pattern, recursive=True)
            if matches:
                source = matches[0]
                dest = os.path.join(output_folder, f"{file_id}.sdf")
                try:
                    shutil.copy2(source, dest)
                    copied_files["compounds"].append(dest)
                    print(f"Copied compound: {os.path.basename(source)} -> {os.path.basename(dest)}")
                    break
                except Exception as e:
                    print(f"Warning: Could not copy compound {source}: {e}")

        # Copy sequence files (FASTA files)
        sequence_patterns = [
            os.path.join(input_folder, f"{file_id}.fasta"),
            os.path.join(input_folder, f"*{file_id}*.fasta"),
            os.path.join(input_folder, f"{file_id}.fa"),
            os.path.join(input_folder, f"*{file_id}*.fa"),
            os.path.join(input_folder, "**", f"{file_id}.fasta"),
            os.path.join(input_folder, "**", f"*{file_id}*.fasta"),
        ]

        for pattern in sequence_patterns:
            matches = glob.glob(pattern, recursive=True)
            if matches:
                source = matches[0]
                dest = os.path.join(output_folder, f"{file_id}.fasta")
                try:
                    shutil.copy2(source, dest)
                    copied_files["sequences"].append(dest)
                    print(f"Copied sequence: {os.path.basename(source)} -> {os.path.basename(dest)}")
                    break
                except Exception as e:
                    print(f"Warning: Could not copy sequence {source}: {e}")

    # Summary
    for file_type, files in copied_files.items():
        if files:
            print(f"Copied {len(files)} {file_type} files")

    return copied_files


def slice_datasheets(config_data: Dict[str, Any]) -> None:
    """
    Slice all datasheets to first N rows and copy associated files.

    Args:
        config_data: Configuration dictionary with slicing parameters
    """
    input_config = config_data['input_config']
    n_rows = config_data['n_rows']
    output_folder = config_data['output_folder']

    print(f"Slicing datasheets to first {n_rows} rows")
    print(f"Output folder: {output_folder}")

    # Create output directory
    os.makedirs(output_folder, exist_ok=True)

    # Process each datasheet
    all_sliced_ids = []

    datasheets = input_config.get('datasheets', {})
    if not datasheets:
        print("Warning: No datasheets found in input configuration")
        return

    print(f"\nProcessing {len(datasheets)} datasheets:")

    for name, datasheet_info in datasheets.items():
        print(f"\n--- Processing datasheet: {name} ---")

        # Extract path from datasheet_info (handle both dict and object formats)
        if isinstance(datasheet_info, dict):
            datasheet_path = datasheet_info.get('path', '')
        else:
            # Assume it's a path string
            datasheet_path = str(datasheet_info)

        if not datasheet_path:
            print(f"Warning: No path found for datasheet {name}")
            continue

        # Output path uses same filename
        output_filename = os.path.basename(datasheet_path)
        output_path = os.path.join(output_folder, output_filename)

        # Slice the datasheet
        sliced_ids = slice_datasheet(datasheet_path, n_rows, output_path)

        # Keep track of all IDs (use first datasheet's IDs for file copying)
        if not all_sliced_ids and sliced_ids:
            all_sliced_ids = sliced_ids

    # Copy associated files for the sliced IDs
    input_folder = input_config.get('output_folder')
    if input_folder and all_sliced_ids:
        copied_files = copy_associated_files(all_sliced_ids, input_folder, output_folder)
    else:
        print("No input folder or sliced IDs found, skipping file copying")
        copied_files = {"structures": [], "compounds": [], "sequences": []}

    # Summary
    print(f"\n=== Slicing Summary ===")
    print(f"Sliced to first {n_rows} rows from {len(datasheets)} datasheets")
    print(f"Total IDs processed: {len(all_sliced_ids)}")
    for file_type, files in copied_files.items():
        print(f"Copied {len(files)} {file_type} files")
    print(f"Output saved to: {output_folder}")


def main():
    parser = argparse.ArgumentParser(description='Slice datasheets to first N rows')
    parser.add_argument('--config', required=True, help='JSON config file with slicing parameters')

    args = parser.parse_args()

    # Load configuration
    if not os.path.exists(args.config):
        print(f"Error: Config file not found: {args.config}")
        sys.exit(1)

    try:
        with open(args.config, 'r') as f:
            config_data = json.load(f)
    except Exception as e:
        print(f"Error loading config: {e}")
        sys.exit(1)

    # Validate required parameters
    required_params = ['input_config', 'n_rows', 'output_folder']
    for param in required_params:
        if param not in config_data:
            print(f"Error: Missing required parameter: {param}")
            sys.exit(1)

    try:
        slice_datasheets(config_data)
        print("\nDatasheet slicing completed successfully!")

    except Exception as e:
        print(f"Error slicing datasheets: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()