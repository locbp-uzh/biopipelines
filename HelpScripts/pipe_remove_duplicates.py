#!/usr/bin/env python3
"""
Runtime helper script for removing duplicate sequences.

This script compares protein sequences in new_pool against reference_pool
and returns only unique sequences from new_pool that don't exist in reference_pool.
"""

import os
import sys
import argparse
import json
import pandas as pd
import shutil
from typing import Dict, List, Any, Optional, Set
from glob import glob


def load_sequences_from_pool(pool_config: Dict[str, Any], compare_column: str) -> pd.DataFrame:
    """
    Load sequences from a pool configuration.
    
    Args:
        pool_config: Pool configuration dictionary
        compare_column: Name of column to compare for duplicates
        
    Returns:
        DataFrame with sequences
    """
    prefix = pool_config['prefix']
    output_folder = pool_config.get('output_folder')
    
    if not output_folder or not os.path.exists(output_folder):
        print(f"Warning: Output folder not found for {prefix}: {output_folder}")
        return pd.DataFrame()
    
    # Try to find sequence CSV file
    csv_path = pool_config.get('sequence_csv')
    
    if not csv_path or not os.path.exists(csv_path):
        # Try alternative sequence file locations
        sequence_files = pool_config.get('sequence_files', [])
        for seq_file in sequence_files:
            if os.path.exists(seq_file):
                csv_path = seq_file
                break
        
        if not csv_path or not os.path.exists(csv_path):
            print(f"Warning: No sequence file found for {prefix}")
            return pd.DataFrame()
    
    try:
        df = pd.read_csv(csv_path)
        print(f"Loaded {prefix} sequences: {df.shape[0]} rows from {csv_path}")
        
        # Check if compare column exists
        if compare_column not in df.columns:
            print(f"Warning: Compare column '{compare_column}' not found in {prefix}")
            print(f"Available columns: {list(df.columns)}")
            return pd.DataFrame()
        
        # Add source information
        df['source_pool'] = prefix
        
        return df
        
    except Exception as e:
        print(f"Error loading sequences from {csv_path}: {e}")
        return pd.DataFrame()


def get_unique_sequences(new_df: pd.DataFrame, reference_df: pd.DataFrame, compare_column: str) -> pd.DataFrame:
    """
    Filter new sequences to only include those not present in reference.
    
    Args:
        new_df: DataFrame with new sequences
        reference_df: DataFrame with reference sequences
        compare_column: Name of column to compare for duplicates
        
    Returns:
        DataFrame with unique sequences from new_df
    """
    if new_df.empty:
        print("No new sequences to process")
        return new_df
    
    if reference_df.empty:
        print("No reference sequences - all new sequences are unique")
        return new_df
    
    print(f"Checking {len(new_df)} new sequences against {len(reference_df)} reference sequences")
    
    # Get reference sequence set
    reference_sequences = set(reference_df[compare_column].dropna())
    print(f"Reference contains {len(reference_sequences)} unique sequences")
    
    # Filter new sequences
    new_sequences = new_df[compare_column].dropna()
    unique_mask = ~new_sequences.isin(reference_sequences)
    unique_df = new_df[unique_mask].copy()
    
    duplicates_count = len(new_df) - len(unique_df)
    print(f"Found {duplicates_count} duplicates, keeping {len(unique_df)} unique sequences")
    
    return unique_df


def save_missing_table(new_df: pd.DataFrame, unique_df: pd.DataFrame, 
                          output_folder: str) -> None:
    """
    Save table of missing/duplicate sequences that were filtered out.
    
    Args:
        new_df: Original DataFrame with all new sequences
        unique_df: DataFrame with unique sequences kept
        output_folder: Output directory
    """
    if 'id' not in new_df.columns:
        print("Warning: No 'id' column found - cannot generate missing table")
        print(f"Available columns in new_df: {list(new_df.columns) if not new_df.empty else 'DataFrame is empty'}")
        return
    
    # Get IDs that were filtered out (missing/duplicates)
    kept_ids = set(unique_df['id']) if not unique_df.empty else set()
    all_ids = set(new_df['id'])
    missing_ids = all_ids - kept_ids
    
    missing_csv = os.path.join(output_folder, "missing.csv")
    os.makedirs(output_folder, exist_ok=True)
    
    if missing_ids:
        # Create rows for missing IDs with their expected file paths (consistent with Filter)
        missing_data = []
        for missing_id in missing_ids:
            structure_path = os.path.join(output_folder, f"{missing_id}.pdb")
            msa_path = os.path.join(output_folder, f"{missing_id}.csv")
            missing_data.append({
                'id': missing_id,
                'structure': structure_path,
                'msa': msa_path
            })
        missing_df = pd.DataFrame(missing_data)
        missing_df.to_csv(missing_csv, index=False)
        print(f"Missing sequences table saved: {missing_csv} ({len(missing_ids)} entries)")
    else:
        # Create empty missing.csv file - pipeline expects it to exist
        empty_df = pd.DataFrame({
            'id': [],
            'structure': [],
            'msa': []
        })
        empty_df.to_csv(missing_csv, index=False)
        print(f"No missing sequences - created empty missing.csv: {missing_csv}")


def remove_duplicates(config_data: Dict[str, Any]) -> None:
    """
    Remove duplicate sequences from pool against history.
    
    Args:
        config_data: Configuration dictionary with pool configurations
    """
    pool_config = config_data['pool_config']
    history_config = config_data['history_config']
    compare = config_data.get('compare', 'sequence')
    
    output_csv = config_data['output_csv']
    output_folder = os.path.dirname(output_csv)
    
    print(f"Removing duplicates by comparing column: {compare}")
    
    # Load sequences from both pools
    print("\nLoading pool sequences...")
    pool_df = load_sequences_from_pool(pool_config, compare)
    
    print("\nLoading history sequences...")
    history_df = load_sequences_from_pool(history_config, compare)
    
    # Filter to get unique sequences
    print("\nFiltering unique sequences...")
    unique_df = get_unique_sequences(pool_df, history_df, compare)
    
    # Save unique sequences (preserve original IDs)
    os.makedirs(os.path.dirname(output_csv), exist_ok=True)
    unique_df.to_csv(output_csv, index=False)
    print(f"\nUnique sequences saved: {output_csv}")
    
    # Generate missing table for filtered sequences
    save_missing_table(pool_df, unique_df, output_folder)
    
    # Verify missing.csv was created
    missing_csv_path = os.path.join(output_folder, "missing.csv")
    if os.path.exists(missing_csv_path):
        print(f"Verified missing.csv exists: {missing_csv_path}")
    else:
        print(f"ERROR: missing.csv not found at expected path: {missing_csv_path}")
    
    print("\nDuplicate removal completed successfully!")
    print(f"Pool sequences: {len(pool_df)}")
    print(f"History sequences: {len(history_df)}")
    print(f"Unique sequences: {len(unique_df)}")
    print(f"Duplicates removed: {len(pool_df) - len(unique_df)}")


def main():
    parser = argparse.ArgumentParser(description='Remove duplicate sequences')
    parser.add_argument('--config', required=True, help='JSON config file with duplicate removal parameters')
    
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
    required_params = ['pool_config', 'history_config', 'output_csv']
    for param in required_params:
        if param not in config_data:
            print(f"Error: Missing required parameter: {param}")
            sys.exit(1)
    
    try:
        remove_duplicates(config_data)
        
    except Exception as e:
        print(f"Error removing duplicates: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()