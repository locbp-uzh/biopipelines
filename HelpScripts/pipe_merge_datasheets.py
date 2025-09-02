#!/usr/bin/env python3
"""
Runtime helper script for merging analysis datasheets.

This script merges CSV files from multiple analysis tools into a unified
datasheet, handling metric name collisions with prefixes and adding calculated columns.
"""

import os
import sys
import argparse
import json
import pandas as pd
from typing import Dict, List, Any, Optional


def merge_datasheets(config_data: Dict[str, Any]) -> None:
    """
    Merge multiple analysis CSV files into one unified datasheet.
    
    Args:
        config_data: Configuration dictionary with input files and settings
    """
    input_csvs = config_data['input_csvs']
    merge_key = config_data['merge_key']
    prefixes = config_data.get('prefixes', [])
    calculate = config_data.get('calculate', {})
    output_csv = config_data['output_csv']
    
    print(f"Combining {len(input_csvs)} CSV files on key '{merge_key}'")
    if prefixes:
        print(f"Using prefixes: {prefixes}")
    if calculate:
        print(f"Calculated columns: {list(calculate.keys())}")
    
    # Load all input dataframes
    dataframes = []
    for i, csv_path in enumerate(input_csvs):
        if not os.path.exists(csv_path):
            print(f"Warning: Input file not found: {csv_path}")
            continue
        
        print(f"Loading: {csv_path}")
        try:
            df = pd.read_csv(csv_path)
            print(f"  - Shape: {df.shape}")
            print(f"  - Columns: {list(df.columns)}")
            
            # Check if merge key exists
            if merge_key not in df.columns:
                print(f"Error: Merge key '{merge_key}' not found in {csv_path}")
                print(f"Available columns: {list(df.columns)}")
                continue
            
            # Apply prefix to column names (except merge key and common columns)
            if prefixes and i < len(prefixes) and prefixes[i]:
                prefix = prefixes[i]
                common_cols = {merge_key, 'source_structure', 'id', 'source_id'}
                columns_to_rename = {}
                
                for col in df.columns:
                    if col not in common_cols:
                        new_name = f"{prefix}{col}"
                        columns_to_rename[col] = new_name
                        print(f"  - Renaming column: {col} -> {new_name}")
                
                if columns_to_rename:
                    df = df.rename(columns=columns_to_rename)
            
            dataframes.append(df)
            
        except Exception as e:
            print(f"Error loading {csv_path}: {e}")
            continue
    
    if not dataframes:
        raise ValueError("No valid input dataframes found")
    
    # Start with the first dataframe
    combined_df = dataframes[0]
    print(f"\nStarting with: {combined_df.shape[0]} rows")
    
    # Merge with remaining dataframes
    for i, df in enumerate(dataframes[1:], 1):
        print(f"Merging dataframe {i+1}: {df.shape[0]} rows")
        
        # Check for overlapping columns (excluding merge key)
        overlap = set(combined_df.columns) & set(df.columns) - {merge_key}
        if overlap:
            print(f"  - Warning: Overlapping columns detected: {overlap}")
        
        # Perform outer join to preserve all data
        before_shape = combined_df.shape
        combined_df = pd.merge(combined_df, df, on=merge_key, how='outer')
        after_shape = combined_df.shape
        
        print(f"  - Result: {before_shape} -> {after_shape}")
    
    # Apply calculated columns
    if calculate:
        print(f"\nCalculating derived columns...")
        for col_name, expression in calculate.items():
            try:
                # Use pandas eval for safe expression evaluation
                combined_df[col_name] = combined_df.eval(expression)
                print(f"  - Added column '{col_name}': {expression}")
            except Exception as e:
                print(f"  - Error calculating '{col_name}': {e}")
                # Set to NaN if calculation fails
                combined_df[col_name] = float('nan')
    
    print(f"\nFinal combined dataframe: {combined_df.shape}")
    print(f"Columns: {list(combined_df.columns)}")
    
    # Create output directory
    os.makedirs(os.path.dirname(output_csv), exist_ok=True)
    
    # Save combined dataframe
    combined_df.to_csv(output_csv, index=False)
    print(f"\nCombined datasheet saved: {output_csv}")
    
    # Summary
    print("\nCombination completed successfully!")
    print(f"Input files: {len(input_csvs)}")
    print(f"Final rows: {combined_df.shape[0]}")
    print(f"Final columns: {combined_df.shape[1]}")


def main():
    parser = argparse.ArgumentParser(description='Merge analysis datasheets')
    parser.add_argument('--config', required=True, help='JSON config file with combination parameters')
    
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
    required_params = ['input_csvs', 'merge_key', 'output_csv']
    for param in required_params:
        if param not in config_data:
            print(f"Error: Missing required parameter: {param}")
            sys.exit(1)
    
    try:
        merge_datasheets(config_data)
        
    except Exception as e:
        print(f"Error combining datasheets: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()