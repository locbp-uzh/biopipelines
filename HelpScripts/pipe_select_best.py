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


def find_matching_structure(row: pd.Series, structures_dir: Optional[str]) -> Optional[str]:
    """
    Find the structure file that matches this datasheet row.
    
    Args:
        row: Row from the datasheet
        structures_dir: Directory containing structure files
        
    Returns:
        Path to matching structure file or None
    """
    if not structures_dir or not os.path.exists(structures_dir):
        return None
    
    # Try different strategies to match structure files
    possible_names = []
    
    # Strategy 1: Use 'id' column if available
    if 'id' in row and pd.notna(row['id']):
        possible_names.append(str(row['id']))
    
    # Strategy 2: Use 'source_structure' column if available
    if 'source_structure' in row and pd.notna(row['source_structure']):
        source = str(row['source_structure'])
        possible_names.append(source)
        # Also try without extension
        if '.' in source:
            possible_names.append(source.rsplit('.', 1)[0])
    
    # Strategy 3: Try structure_id or similar columns
    id_columns = [col for col in row.index if 'id' in col.lower()]
    for col in id_columns:
        if pd.notna(row[col]):
            possible_names.append(str(row[col]))
    
    # Search for matching files
    structure_extensions = ['.pdb', '.cif', '.mmcif']
    
    for name in possible_names:
        for ext in structure_extensions:
            # Try exact match
            candidates = [
                os.path.join(structures_dir, f"{name}{ext}"),
                os.path.join(structures_dir, f"{name}"),
            ]
            
            for candidate in candidates:
                if os.path.exists(candidate):
                    return candidate
            
            # Try glob pattern matching
            pattern = os.path.join(structures_dir, f"*{name}*{ext}")
            matches = glob.glob(pattern)
            if matches:
                return matches[0]  # Return first match
    
    return None


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
    
    # 3. Extract sequences from datasheets
    sequence_datasheets = [
        os.path.join(pool_folder, "*sequences*.csv"),
        os.path.join(pool_folder, "*sequence*.csv"),
        os.path.join(pool_folder, "**", "*sequences*.csv"),
    ]
    
    for pattern in sequence_datasheets:
        matches = glob.glob(pattern, recursive=True)
        if matches:
            source_csv = matches[0]
            try:
                df = pd.read_csv(source_csv)
                if 'id' in df.columns:
                    # Find rows matching the selected ID
                    matching_rows = df[df['id'] == selected_id]
                    if not matching_rows.empty:
                        dest = os.path.join(output_folder, "best_sequences.csv")
                        matching_rows.to_csv(dest, index=False)
                        extracted_files['sequences'] = dest
                        print(f"Extracted sequences: {len(matching_rows)} rows -> {dest}")
                        break
            except Exception as e:
                print(f"Warning: Could not process sequence datasheet {source_csv}: {e}")
    
    return extracted_files


def select_best_from_arrays(config_data: Dict[str, Any]) -> None:
    """
    Select best item from multiple datasheets and extract from corresponding pools.
    
    Args:
        config_data: Configuration with pool_folders and datasheet_paths arrays
    """
    pool_folders = config_data['pool_folders']
    datasheet_paths = config_data['datasheet_paths']
    selection_metric = config_data['selection_metric']
    selection_mode = config_data['selection_mode']
    output_csv = config_data['output_csv']
    output_structure = config_data['output_structure']
    
    print(f"Array mode: Selecting best from {len(datasheet_paths)} datasheets")
    print(f"Selection metric: {selection_metric} ({selection_mode})")
    
    # Load and concatenate all datasheets
    all_dfs = []
    pool_indices = []  # Track which pool each row comes from
    
    for i, datasheet_path in enumerate(datasheet_paths):
        if not os.path.exists(datasheet_path):
            print(f"Warning: Datasheet not found: {datasheet_path}")
            continue
            
        try:
            df = pd.read_csv(datasheet_path)
            df['_pool_index'] = i  # Add pool index to each row
            all_dfs.append(df)
            pool_indices.extend([i] * len(df))
            print(f"Loaded datasheet {i}: {datasheet_path} ({len(df)} rows)")
        except Exception as e:
            print(f"Error loading {datasheet_path}: {e}")
            continue
    
    if not all_dfs:
        raise ValueError("No valid datasheets found")
    
    # Concatenate all dataframes
    combined_df = pd.concat(all_dfs, ignore_index=True)
    print(f"Combined dataframe: {combined_df.shape}")
    
    # Find best row
    if selection_metric not in combined_df.columns:
        raise ValueError(f"Metric '{selection_metric}' not found in combined data")
    
    metric_values = combined_df[selection_metric]
    if selection_mode in ["max", "maximize"]:
        best_idx = metric_values.idxmax()
    else:
        best_idx = metric_values.idxmin()
    
    selected_row = combined_df.iloc[best_idx]
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
    print(f"Output datasheet: {output_csv}")
    print(f"Output structure: {output_structure}")


def select_best_item(config_data: Dict[str, Any]) -> None:
    """
    Select the best item from analysis results.
    
    Args:
        config_data: Configuration dictionary with selection parameters
    """
    # Check for array mode
    if config_data.get('use_array_mode', False):
        select_best_from_arrays(config_data)
        return
    
    # Single mode
    input_csv = config_data['input_csv']
    structures_dir = config_data.get('source_structures_dir')
    selection_metric = config_data['selection_metric']
    selection_mode = config_data['selection_mode']
    
    # Convert mode from "max"/"min" to "maximize"/"minimize"
    if selection_mode == "max":
        selection_mode = "maximize"
    elif selection_mode == "min":
        selection_mode = "minimize"
    weights = config_data.get('weights', {})
    tie_breaker = config_data['tie_breaker']
    composite_function = config_data.get('composite_function', 'weighted_sum')
    output_csv = config_data['output_csv']
    output_structure = config_data['output_structure']
    
    # New pool mode parameters
    use_pool_mode = config_data.get('use_pool_mode', False)
    data_name = config_data.get('data_name')
    pool_output_folder = config_data.get('pool_output_folder')
    
    print(f"Selecting best item using metric: {selection_metric}")
    print(f"Selection mode: {selection_mode}")
    print(f"Input: {input_csv}")
    
    # Load input data
    if not os.path.exists(input_csv):
        raise FileNotFoundError(f"Input CSV not found: {input_csv}")
    
    try:
        df = pd.read_csv(input_csv)
        print(f"Loaded dataframe: {df.shape}")
        print(f"Columns: {list(df.columns)}")
        
    except Exception as e:
        raise ValueError(f"Error loading CSV file: {e}")
    
    if df.empty:
        raise ValueError("Input dataframe is empty - nothing to select")
    
    if len(df) == 1:
        print("Only one item in input - selecting it by default")
        selected_row = df.iloc[0]
    else:
        # Determine selection metric value
        if weights:
            # Multi-objective selection using composite score
            print(f"Using composite score with weights: {weights}")
            print(f"Composite function: {composite_function}")
            
            composite_scores = calculate_composite_score(df, weights, composite_function)
            df['composite_score'] = composite_scores
            
            # Use composite score as selection metric
            metric_values = composite_scores
            effective_metric = 'composite_score'
        else:
            # Single metric selection
            if selection_metric not in df.columns:
                available_cols = ", ".join(df.columns)
                raise ValueError(f"Selection metric '{selection_metric}' not found. "
                               f"Available columns: {available_cols}")
            
            metric_values = df[selection_metric].astype(float)
            effective_metric = selection_metric
        
        # Handle missing values
        if metric_values.isnull().any():
            print(f"Warning: Found {metric_values.isnull().sum()} missing values in {effective_metric}")
            metric_values = metric_values.fillna(metric_values.median())
        
        print(f"Metric '{effective_metric}' statistics:")
        print(f"  Min: {metric_values.min():.4f}")
        print(f"  Max: {metric_values.max():.4f}")
        print(f"  Mean: {metric_values.mean():.4f}")
        
        # Find best value
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
                if tie_breaker not in df.columns:
                    print(f"Warning: Tie breaker metric '{tie_breaker}' not found, using first")
                    selected_idx = best_indices[0]
                else:
                    tie_values = df.loc[best_indices, tie_breaker].astype(float)
                    # Maximize tie breaker (could make this configurable)
                    best_tie_idx = tie_values.idxmax()
                    selected_idx = best_tie_idx
                    print(f"Used {tie_breaker} as tie breaker: {tie_values[best_tie_idx]:.4f}")
        
        selected_row = df.iloc[selected_idx]
    
    # Create output directory
    os.makedirs(os.path.dirname(output_csv), exist_ok=True)
    
    # Save selected item data
    selected_df = pd.DataFrame([selected_row])
    selected_df.to_csv(output_csv, index=False)
    
    print(f"\\nSelected item saved: {output_csv}")
    print("Selected item details:")
    for col, val in selected_row.items():
        print(f"  {col}: {val}")
    
    # Handle data extraction based on mode
    extracted_files = {}
    
    if use_pool_mode and pool_output_folder:
        # Pool mode: extract all data types for selected ID
        selected_id = selected_row.get('id', 'unknown')
        print(f"\\nPool mode: Extracting all data for ID '{selected_id}'")
        
        output_dir = os.path.dirname(output_csv)
        extracted_files = extract_pool_data_for_id(selected_id, pool_output_folder, output_dir)
        
        # Rename structure file to match expected output name
        if 'structure' in extracted_files and output_structure:
            try:
                if os.path.exists(extracted_files['structure']) and extracted_files['structure'] != output_structure:
                    shutil.move(extracted_files['structure'], output_structure)
                    extracted_files['structure'] = output_structure
                    print(f"Renamed structure to: {output_structure}")
            except Exception as e:
                print(f"Warning: Could not rename structure file: {e}")
    
    else:
        # Legacy mode: try to copy structure file only
        structure_copied = False
        if structures_dir and output_structure:
            matching_structure = find_matching_structure(selected_row, structures_dir)
            
            if matching_structure:
                try:
                    shutil.copy2(matching_structure, output_structure)
                    structure_copied = True
                    extracted_files['structure'] = output_structure
                    print(f"\\nCopied structure: {matching_structure} -> {output_structure}")
                except Exception as e:
                    print(f"Warning: Could not copy structure file: {e}")
            else:
                print(f"Warning: Could not find matching structure file for selected item")
    
    # Summary
    print(f"\\nSelection completed successfully!")
    print(f"Selected metric value: {selected_row.get(effective_metric, 'N/A')}")
    print(f"Output datasheet: {output_csv}")
    
    if use_pool_mode:
        print(f"Pool mode - Extracted data types: {list(extracted_files.keys())}")
        for data_type, file_path in extracted_files.items():
            print(f"  {data_type}: {file_path}")
    else:
        if 'structure' in extracted_files:
            print(f"Output structure: {extracted_files['structure']}")
        else:
            print("No structure file available")


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
    
    # Validate required parameters
    required_params = ['input_csv', 'selection_metric', 'selection_mode', 'output_csv']
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