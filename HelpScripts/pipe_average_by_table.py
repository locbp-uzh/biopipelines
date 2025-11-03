#!/usr/bin/env python3
"""
Helper script for AverageByTable tool.

Computes averages of numeric columns across multiple tables.
Each row in the output corresponds to the average from each input table.
Non-numeric values and inconsistent values are discarded.
"""

import pandas as pd
import numpy as np
import argparse
import json
import sys
import os


def load_table(path):
    """Load a table CSV file."""
    try:
        return pd.read_csv(path)
    except Exception as e:
        print(f"Error loading table {path}: {e}")
        return None


def compute_averages(tables, table_paths):
    """
    Compute averages across multiple tables.

    Args:
        tables: List of pandas DataFrames
        table_paths: List of file paths for naming

    Returns:
        pandas DataFrame with averages
    """
    if not tables:
        raise ValueError("No valid tables to process")

    # Find common numeric columns across all tables
    common_columns = None
    numeric_columns = set()

    for df in tables:
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

    # Keep only numeric columns that exist in all tables
    numeric_columns = numeric_columns.intersection(common_columns)

    if not numeric_columns:
        print("Warning: No common numeric columns found across all tables")
        return pd.DataFrame()

    print(f"Found {len(numeric_columns)} common numeric columns: {list(numeric_columns)}")

    # Compute averages for each table
    results = []

    for i, (df, path) in enumerate(zip(tables, table_paths)):
        if df is None:
            continue

        row = {"table_name": os.path.basename(path)}

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
    parser = argparse.ArgumentParser(description='Compute averages across multiple tables')
    parser.add_argument('--config', required=True, help='Configuration file path')
    args = parser.parse_args()

    # Load configuration
    try:
        with open(args.config, 'r') as f:
            config = json.load(f)
    except Exception as e:
        print(f"Error loading configuration: {e}")
        sys.exit(1)

    table_paths = config['tables']
    output_path = config['output']

    print(f"Processing {len(table_paths)} tables")

    # Load all tables
    tables = []
    for path in table_paths:
        print(f"Loading: {path}")
        df = load_table(path)
        tables.append(df)
        if df is not None:
            print(f"  Shape: {df.shape}")
        else:
            print(f"  Failed to load")

    # Compute averages
    try:
        averages_df = compute_averages(tables, table_paths)

        if averages_df.empty:
            print("Warning: No averages computed")
            # Create empty output file
            averages_df = pd.DataFrame(columns=["table_name"])

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