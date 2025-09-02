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


def select_best_item(config_data: Dict[str, Any]) -> None:
    """
    Select the best item from analysis results.
    
    Args:
        config_data: Configuration dictionary with selection parameters
    """
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
    
    # Try to copy the corresponding structure file
    structure_copied = False
    if structures_dir and output_structure:
        matching_structure = find_matching_structure(selected_row, structures_dir)
        
        if matching_structure:
            try:
                shutil.copy2(matching_structure, output_structure)
                structure_copied = True
                print(f"\\nCopied structure: {matching_structure} -> {output_structure}")
            except Exception as e:
                print(f"Warning: Could not copy structure file: {e}")
        else:
            print(f"Warning: Could not find matching structure file for selected item")
    
    # Summary
    print(f"\\nSelection completed successfully!")
    print(f"Selected metric value: {selected_row.get(effective_metric, 'N/A')}")
    print(f"Output datasheet: {output_csv}")
    if structure_copied:
        print(f"Output structure: {output_structure}")
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