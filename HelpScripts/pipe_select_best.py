#!/usr/bin/env python3
"""
Runtime helper script for SelectBest tool.

This script applies selection criteria to choose the single best item from
analysis results, with support for multi-objective optimization and tie-breaking.
"""

import os
import sys
import argparse
import json
import pandas as pd
import numpy as np
import shutil
import glob
from typing import Dict, List, Any, Optional


def calculate_composite_score(df: pd.DataFrame, weights: Dict[str, float], function: str) -> pd.Series:
    """
    Calculate composite score from multiple metrics.
    
    Args:
        df: DataFrame with metrics
        weights: Dictionary of metric weights
        function: Composite function type
        
    Returns:
        Series with composite scores
    """
    # Validate that all weighted metrics exist
    missing_metrics = []
    for metric in weights:
        if metric not in df.columns:
            missing_metrics.append(metric)
    
    if missing_metrics:
        raise ValueError(f"Missing metrics for composite score: {missing_metrics}")
    
    # Normalize metrics to [0, 1] range for fair combination
    normalized_metrics = {}
    for metric, weight in weights.items():
        values = df[metric].astype(float)
        
        # Handle missing values
        if values.isnull().any():
            print(f"Warning: Found {values.isnull().sum()} missing values in {metric}")
            values = values.fillna(values.median())
        
        # Normalize to [0, 1] range
        min_val, max_val = values.min(), values.max()
        if max_val > min_val:
            normalized = (values - min_val) / (max_val - min_val)
        else:
            normalized = pd.Series([0.5] * len(values), index=values.index)
        
        normalized_metrics[metric] = normalized * weight
    
    # Apply composite function
    if function == "weighted_sum":
        # Standard weighted sum
        composite = sum(normalized_metrics.values())
    
    elif function == "product":
        # Geometric mean with weights
        composite = np.ones(len(df))
        for metric, weighted_values in normalized_metrics.items():
            # Add small epsilon to avoid zero products
            composite *= (weighted_values + 0.001) ** weights[metric]
    
    elif function == "min":
        # Take minimum of weighted normalized metrics (for robust optimization)
        stacked = pd.concat(normalized_metrics.values(), axis=1)
        composite = stacked.min(axis=1)
    
    elif function == "max":
        # Take maximum of weighted normalized metrics
        stacked = pd.concat(normalized_metrics.values(), axis=1)
        composite = stacked.max(axis=1)
    
    else:
        raise ValueError(f"Unknown composite function: {function}")
    
    return composite



def extract_pool_data_for_id(selected_id: str, pool_folder: str, output_folder: str) -> Dict[str, str]:
    """
    Extract all data types (structures, compounds, sequences) from pool output for selected ID.
    
    Args:
        selected_id: The ID of the selected item
        pool_folder: The pool output folder containing all data types
        output_folder: Where to copy the selected data
        
    Returns:
        Dictionary mapping data type to extracted file path
    """
    extracted_files = {}
    
    # Create output directory
    os.makedirs(output_folder, exist_ok=True)
    
    # 1. Extract structures (PDB/CIF files)
    structure_patterns = [
        os.path.join(pool_folder, f"{selected_id}.pdb"),
        os.path.join(pool_folder, f"*{selected_id}*.pdb"),
        os.path.join(pool_folder, f"{selected_id}.cif"),
        os.path.join(pool_folder, f"*{selected_id}*.cif"),
        os.path.join(pool_folder, "**", f"{selected_id}.pdb"),
        os.path.join(pool_folder, "**", f"*{selected_id}*.pdb"),
    ]
    
    for pattern in structure_patterns:
        matches = glob.glob(pattern, recursive=True)
        if matches:
            source = matches[0]  # Take first match
            dest = os.path.join(output_folder, "best.pdb")
            try:
                shutil.copy2(source, dest)
                extracted_files['structure'] = dest
                print(f"Extracted structure: {source} -> {dest}")
                break
            except Exception as e:
                print(f"Warning: Could not copy structure {source}: {e}")
    
    # 2. Extract compounds (SDF files)
    compound_patterns = [
        os.path.join(pool_folder, f"{selected_id}.sdf"),
        os.path.join(pool_folder, f"*{selected_id}*.sdf"),
        os.path.join(pool_folder, "**", f"{selected_id}.sdf"),
        os.path.join(pool_folder, "**", f"*{selected_id}*.sdf"),
    ]
    
    for pattern in compound_patterns:
        matches = glob.glob(pattern, recursive=True)
        if matches:
            source = matches[0]
            dest = os.path.join(output_folder, "best_compound.sdf")
            try:
                shutil.copy2(source, dest)
                extracted_files['compound'] = dest
                print(f"Extracted compound: {source} -> {dest}")
                break
            except Exception as e:
                print(f"Warning: Could not copy compound {source}: {e}")
    
    # 3. Extract sequences from tables
    sequence_tables = [
        os.path.join(pool_folder, "*sequences*.csv"),
        os.path.join(pool_folder, "*sequence*.csv"),
        os.path.join(pool_folder, "**", "*sequences*.csv"),
    ]
    
    for pattern in sequence_tables:
        matches = glob.glob(pattern, recursive=True)
        if matches:
            source_csv = matches[0]
            try:
                df = pd.read_csv(source_csv)
                if 'id' in df.columns:
                    # Find rows matching the selected ID
                    matching_rows = df[df['id'] == selected_id]
                    if not matching_rows.empty:
                        # Use sequences.csv as the standard filename
                        dest = os.path.join(output_folder, "sequences.csv")
                        matching_rows.to_csv(dest, index=False)
                        extracted_files['sequences'] = dest
                        print(f"Extracted sequences: {len(matching_rows)} rows -> {dest}")
                        break
            except Exception as e:
                print(f"Warning: Could not process sequence table {source_csv}: {e}")
    
    return extracted_files


def select_best_from_arrays(config_data: Dict[str, Any]) -> None:
    """
    Select best item from multiple tables and extract from corresponding pools.
    
    Args:
        config_data: Configuration with pool_folders and table_paths arrays
    """
    pool_folders = config_data['pool_folders']
    table_paths = config_data['table_paths']
    selection_metric = config_data['selection_metric']
    selection_mode = config_data['selection_mode']
    output_csv = config_data['output_csv']
    output_structure = config_data['output_structure']
    
    print(f"Array mode: Selecting best from {len(table_paths)} tables")
    print(f"Selection metric: {selection_metric} ({selection_mode})")
    
    # Load and concatenate all tables
    all_dfs = []
    pool_indices = []  # Track which pool each row comes from
    
    for i, table_path in enumerate(table_paths):
        if not os.path.exists(table_path):
            print(f"Warning: Table not found: {table_path}")
            continue
            
        try:
            df = pd.read_csv(table_path)
            df['_pool_index'] = i  # Add pool index to each row
            all_dfs.append(df)
            pool_indices.extend([i] * len(df))
            print(f"Loaded table {i}: {table_path} ({len(df)} rows)")
        except Exception as e:
            print(f"Error loading {table_path}: {e}")
            continue
    
    if not all_dfs:
        raise ValueError("No valid tables found")
    
    # Concatenate all dataframes
    combined_df = pd.concat(all_dfs, ignore_index=True)
    print(f"Combined dataframe: {combined_df.shape}")
    
    # Handle weights and composite functions
    weights = config_data.get('weights', {})
    tie_breaker = config_data.get('tie_breaker', 'first')
    composite_function = config_data.get('composite_function', 'weighted_sum')
    
    # Convert selection mode
    if selection_mode == "max":
        selection_mode = "maximize"
    elif selection_mode == "min":
        selection_mode = "minimize"
    
    # Determine selection metric
    if weights:
        # Multi-objective selection using composite score
        print(f"Using composite score with weights: {weights}")
        print(f"Composite function: {composite_function}")
        
        composite_scores = calculate_composite_score(combined_df, weights, composite_function)
        combined_df['composite_score'] = composite_scores
        
        # Use composite score as selection metric
        metric_values = composite_scores
        effective_metric = 'composite_score'
    else:
        # Single metric selection
        if selection_metric not in combined_df.columns:
            available_cols = ", ".join(combined_df.columns)
            raise ValueError(f"Selection metric '{selection_metric}' not found. "
                           f"Available columns: {available_cols}")
        
        metric_values = combined_df[selection_metric].astype(float)
        effective_metric = selection_metric
    
    # Handle missing values
    if metric_values.isnull().any():
        print(f"Warning: Found {metric_values.isnull().sum()} missing values in {effective_metric}")
        metric_values = metric_values.fillna(metric_values.median())
    
    print(f"Metric '{effective_metric}' statistics:")
    print(f"  Min: {metric_values.min():.4f}")
    print(f"  Max: {metric_values.max():.4f}")
    print(f"  Mean: {metric_values.mean():.4f}")
    
    # Find best value with tie handling
    if selection_mode == "maximize":
        best_value = metric_values.max()
        best_indices = metric_values[metric_values == best_value].index
    else:  # minimize
        best_value = metric_values.min()
        best_indices = metric_values[metric_values == best_value].index
    
    print(f"Best {effective_metric}: {best_value:.4f}")
    print(f"Number of items with best value: {len(best_indices)}")
    
    # Handle ties
    if len(best_indices) == 1:
        selected_idx = best_indices[0]
        print("No ties - unique best value")
    else:
        print(f"Breaking tie between {len(best_indices)} items using: {tie_breaker}")
        
        if tie_breaker == "first":
            selected_idx = best_indices[0]
            print(f"Selected first item: index {selected_idx}")
        
        elif tie_breaker == "random":
            import random
            selected_idx = random.choice(best_indices.tolist())
            print(f"Randomly selected item: index {selected_idx}")
        
        else:
            # Use another metric as tie breaker
            if tie_breaker not in combined_df.columns:
                print(f"Warning: Tie breaker metric '{tie_breaker}' not found, using first")
                selected_idx = best_indices[0]
            else:
                tie_values = combined_df.loc[best_indices, tie_breaker].astype(float)
                # Maximize tie breaker (could make this configurable)
                best_tie_idx = tie_values.idxmax()
                selected_idx = best_tie_idx
                print(f"Used {tie_breaker} as tie breaker: {tie_values[best_tie_idx]:.4f}")
    
    selected_row = combined_df.iloc[selected_idx]
    pool_index = selected_row['_pool_index']
    
    print(f"Selected best item from pool {pool_index}: {selection_metric} = {selected_row[selection_metric]}")
    
    # Save selected row (without _pool_index column)
    selected_row_clean = selected_row.drop('_pool_index').to_frame().T
    selected_row_clean.to_csv(output_csv, index=False)
    
    # Extract structure from the corresponding pool
    selected_id = selected_row.get('id', 'unknown')
    pool_folder = pool_folders[pool_index]
    output_dir = os.path.dirname(output_structure)
    
    print(f"Extracting structure for ID '{selected_id}' from pool: {pool_folder}")
    extracted_files = extract_pool_data_for_id(selected_id, pool_folder, output_dir)
    
    # Rename structure file to match expected output name
    if 'structure' in extracted_files and os.path.exists(extracted_files['structure']):
        if extracted_files['structure'] != output_structure:
            shutil.move(extracted_files['structure'], output_structure)
            print(f"Renamed structure to: {output_structure}")
    
    print(f"Selection completed successfully!")
    print(f"Selected metric value: {selected_row[selection_metric]}")
    print(f"Output table: {output_csv}")
    print(f"Output structure: {output_structure}")


def select_best_item(config_data: Dict[str, Any]) -> None:
    """
    Select the best item from analysis results using array format.
    
    Args:
        config_data: Configuration dictionary with selection parameters
    """
    # Always use array mode - no legacy support
    select_best_from_arrays(config_data)


def main():
    parser = argparse.ArgumentParser(description='Select best item from analysis results')
    parser.add_argument('--config', required=True, help='JSON config file with selection parameters')
    
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
    
    # Validate required parameters - only support array format
    required_params = ['table_paths', 'pool_folders', 'selection_metric', 'selection_mode', 'output_csv']
    
    for param in required_params:
        if param not in config_data:
            print(f"Error: Missing required parameter: {param}")
            sys.exit(1)
    
    try:
        select_best_item(config_data)
        
    except Exception as e:
        print(f"Error selecting best item: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()