#!/usr/bin/env python3
"""
Helper script for AverageByDatasheet tool.

Computes averages of numeric columns across multiple datasheets.
Each row in the output corresponds to the average from each input datasheet.
Non-numeric values and inconsistent values are discarded.
"""

import pandas as pd
import numpy as np
import argparse
import json
import sys
import os


def load_datasheet(path):
    """Load a datasheet CSV file."""
    try:
        return pd.read_csv(path)
    except Exception as e:
        print(f"Error loading datasheet {path}: {e}")
        return None


def compute_averages(datasheets, datasheet_paths):
    """
    Compute averages across multiple datasheets.

    Args:
        datasheets: List of pandas DataFrames
        datasheet_paths: List of file paths for naming

    Returns:
        pandas DataFrame with averages
    """
    if not datasheets:
        raise ValueError("No valid datasheets to process")

    # Find common numeric columns across all datasheets
    common_columns = None
    numeric_columns = set()

    for df in datasheets:
        if df is None:
            continue

        if common_columns is None:
            common_columns = set(df.columns)
        else:
            common_columns = common_columns.intersection(set(df.columns))

        # Find numeric columns
        for col in df.columns:
            if col in common_columns and pd.api.types.is_numeric_dtype(df[col]):
                numeric_columns.add(col)

    # Keep only numeric columns that exist in all datasheets
    numeric_columns = numeric_columns.intersection(common_columns)

    if not numeric_columns:
        print("Warning: No common numeric columns found across all datasheets")
        return pd.DataFrame()

    print(f"Found {len(numeric_columns)} common numeric columns: {list(numeric_columns)}")

    # Compute averages for each datasheet
    results = []

    for i, (df, path) in enumerate(zip(datasheets, datasheet_paths)):
        if df is None:
            continue

        row = {"datasheet_name": os.path.basename(path)}

        # Compute average for each numeric column
        for col in numeric_columns:
            if col in df.columns:
                # Skip non-numeric values and compute mean
                numeric_values = pd.to_numeric(df[col], errors='coerce').dropna()
                if len(numeric_values) > 0:
                    row[f"avg_{col}"] = numeric_values.mean()
                else:
                    row[f"avg_{col}"] = np.nan
            else:
                row[f"avg_{col}"] = np.nan

        results.append(row)

    return pd.DataFrame(results)


def main():
    parser = argparse.ArgumentParser(description='Compute averages across multiple datasheets')
    parser.add_argument('--config', required=True, help='Configuration file path')
    args = parser.parse_args()

    # Load configuration
    try:
        with open(args.config, 'r') as f:
            config = json.load(f)
    except Exception as e:
        print(f"Error loading configuration: {e}")
        sys.exit(1)

    datasheet_paths = config['datasheets']
    output_path = config['output']

    print(f"Processing {len(datasheet_paths)} datasheets")

    # Load all datasheets
    datasheets = []
    for path in datasheet_paths:
        print(f"Loading: {path}")
        df = load_datasheet(path)
        datasheets.append(df)
        if df is not None:
            print(f"  Shape: {df.shape}")
        else:
            print(f"  Failed to load")

    # Compute averages
    try:
        averages_df = compute_averages(datasheets, datasheet_paths)

        if averages_df.empty:
            print("Warning: No averages computed")
            # Create empty output file
            averages_df = pd.DataFrame(columns=["datasheet_name"])

        # Save results
        averages_df.to_csv(output_path, index=False)
        print(f"Saved averages to: {output_path}")
        print(f"Output shape: {averages_df.shape}")

        if not averages_df.empty:
            print("Sample of computed averages:")
            print(averages_df.head())

    except Exception as e:
        print(f"Error computing averages: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()