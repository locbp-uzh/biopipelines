#!/usr/bin/env python3
"""
Helper script for ExtractMetric tool.

Extracts a specific metric from multiple datasheets and formats the data
for easy import into statistical analysis software like Prism.
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


def extract_metric(datasheets, datasheet_paths, datasheet_names, metric):
    """
    Extract a specific metric from multiple datasheets.

    Args:
        datasheets: List of pandas DataFrames
        datasheet_paths: List of file paths
        datasheet_names: List of column names for output
        metric: Name of metric column to extract

    Returns:
        pandas DataFrame with extracted metric values
    """
    if not datasheets:
        raise ValueError("No valid datasheets to process")

    # Check which datasheets have the metric
    valid_datasheets = []
    valid_names = []

    for df, path, name in zip(datasheets, datasheet_paths, datasheet_names):
        if df is None:
            print(f"Skipping invalid datasheet: {path}")
            continue

        if metric not in df.columns:
            print(f"Warning: Metric '{metric}' not found in {path}")
            print(f"Available columns: {list(df.columns)}")
            continue

        valid_datasheets.append(df)
        valid_names.append(name)

    if not valid_datasheets:
        raise ValueError(f"Metric '{metric}' not found in any datasheet")

    print(f"Found metric '{metric}' in {len(valid_datasheets)} datasheets")

    # Extract metric values from each datasheet
    extracted_data = {}

    for df, name in zip(valid_datasheets, valid_names):
        # Extract metric values, converting to numeric and removing NaN
        metric_values = pd.to_numeric(df[metric], errors='coerce').dropna()
        extracted_data[name] = list(metric_values)

        print(f"  {name}: {len(metric_values)} values")

    # Find the maximum number of values to create DataFrame
    max_length = max(len(values) for values in extracted_data.values()) if extracted_data else 0

    # Pad shorter columns with NaN
    for name in extracted_data:
        while len(extracted_data[name]) < max_length:
            extracted_data[name].append(np.nan)

    # Create DataFrame
    if max_length > 0:
        result_df = pd.DataFrame(extracted_data)
    else:
        # Create empty DataFrame with correct column names
        result_df = pd.DataFrame(columns=valid_names)

    return result_df


def main():
    parser = argparse.ArgumentParser(description='Extract specific metric from multiple datasheets')
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
    metric = config['metric']
    datasheet_names = config['datasheet_names']
    output_path = config['output']

    print(f"Extracting metric '{metric}' from {len(datasheet_paths)} datasheets")

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

    # Extract metric
    try:
        extracted_df = extract_metric(datasheets, datasheet_paths, datasheet_names, metric)

        # Save results
        extracted_df.to_csv(output_path, index=False)
        print(f"Saved extracted metric to: {output_path}")
        print(f"Output shape: {extracted_df.shape}")

        if not extracted_df.empty:
            print(f"Sample of extracted '{metric}' values:")
            print(extracted_df.head())
            print("\nSummary statistics:")
            print(extracted_df.describe())
        else:
            print("No data extracted")

    except Exception as e:
        print(f"Error extracting metric: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()