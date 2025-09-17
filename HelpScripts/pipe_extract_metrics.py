#!/usr/bin/env python3
"""
Helper script for ExtractMetrics tool.

Extracts multiple specific metrics from multiple datasheets and creates
separate CSV files for each metric, formatted for easy import into
statistical analysis software like Prism.
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


def extract_single_metric(datasheets, datasheet_paths, datasheet_names, metric):
    """
    Extract a single metric from multiple datasheets.

    Args:
        datasheets: List of pandas DataFrames
        datasheet_paths: List of file paths
        datasheet_names: List of column names for output
        metric: Name of metric column to extract

    Returns:
        pandas DataFrame with extracted metric values
    """
    # Check which datasheets have the metric
    valid_datasheets = []
    valid_names = []

    for df, path, name in zip(datasheets, datasheet_paths, datasheet_names):
        if df is None:
            print(f"  Skipping invalid datasheet: {path}")
            continue

        if metric not in df.columns:
            print(f"  Warning: Metric '{metric}' not found in {path}")
            continue

        valid_datasheets.append(df)
        valid_names.append(name)

    if not valid_datasheets:
        print(f"  Error: Metric '{metric}' not found in any datasheet")
        return pd.DataFrame(columns=datasheet_names)

    print(f"  Found metric '{metric}' in {len(valid_datasheets)} datasheets")

    # Extract metric values from each datasheet
    extracted_data = {}

    for df, name in zip(valid_datasheets, valid_names):
        # Extract metric values, converting to numeric and removing NaN
        metric_values = pd.to_numeric(df[metric], errors='coerce').dropna()
        extracted_data[name] = list(metric_values)

        print(f"    {name}: {len(metric_values)} values")

    # Find the maximum number of values to create DataFrame
    max_length = max(len(values) for values in extracted_data.values()) if extracted_data else 0

    # Pad shorter columns with NaN
    for name in extracted_data:
        while len(extracted_data[name]) < max_length:
            extracted_data[name].append(np.nan)

    # Add empty columns for datasheets that didn't have this metric
    for name in datasheet_names:
        if name not in extracted_data:
            extracted_data[name] = [np.nan] * max_length

    # Create DataFrame with consistent column order
    if max_length > 0:
        result_df = pd.DataFrame({name: extracted_data[name] for name in datasheet_names})
    else:
        # Create empty DataFrame with correct column names
        result_df = pd.DataFrame(columns=datasheet_names)

    return result_df


def extract_multiple_metrics(datasheets, datasheet_paths, datasheet_names, metrics, output_folder):
    """
    Extract multiple metrics from multiple datasheets.

    Args:
        datasheets: List of pandas DataFrames
        datasheet_paths: List of file paths
        datasheet_names: List of column names for output
        metrics: List of metric names to extract
        output_folder: Directory to save output files

    Returns:
        Dictionary mapping metric names to output file paths
    """
    output_files = {}

    for metric in metrics:
        print(f"Extracting metric: {metric}")

        # Extract this metric
        extracted_df = extract_single_metric(datasheets, datasheet_paths, datasheet_names, metric)

        # Save to separate CSV file
        output_path = os.path.join(output_folder, f"{metric}.csv")
        extracted_df.to_csv(output_path, index=False, na_rep='')

        output_files[metric] = output_path

        print(f"  Saved to: {output_path}")
        print(f"  Output shape: {extracted_df.shape}")

        if not extracted_df.empty:
            print(f"  Sample values:")
            print(f"  {extracted_df.head().to_string()}")
        else:
            print(f"  No data extracted")
        print()

    return output_files


def main():
    parser = argparse.ArgumentParser(description='Extract multiple metrics from multiple datasheets')
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
    metrics = config['metrics']
    datasheet_names = config['datasheet_names']
    output_folder = config['output_folder']

    print(f"Extracting {len(metrics)} metrics from {len(datasheet_paths)} datasheets")
    print(f"Metrics: {', '.join(metrics)}")
    print(f"Output folder: {output_folder}")
    print()

    # Load all datasheets
    datasheets = []
    for path in datasheet_paths:
        print(f"Loading: {path}")
        df = load_datasheet(path)
        datasheets.append(df)
        if df is not None:
            print(f"  Shape: {df.shape}")
            print(f"  Columns: {list(df.columns)}")
        else:
            print(f"  Failed to load")
    print()

    # Extract metrics
    try:
        output_files = extract_multiple_metrics(datasheets, datasheet_paths, datasheet_names, metrics, output_folder)

        print(f"Successfully extracted {len(output_files)} metrics:")
        for metric, filepath in output_files.items():
            print(f"  {metric}: {filepath}")

    except Exception as e:
        print(f"Error extracting metrics: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()