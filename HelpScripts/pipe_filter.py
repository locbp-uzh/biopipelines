#!/usr/bin/env python3
"""
Runtime helper script for expression-based filtering.

This script applies pandas query-style expressions to filter CSV files
while preserving all column information.
"""

import os
import sys
import argparse
import json
import pandas as pd
import shutil
import glob
from typing import Dict, List, Any, Optional


def extract_pool_data_for_filtered_ids(filtered_ids: List[str], pool_folder: str, output_folder: str) -> Dict[str, List[str]]:
    """
    Extract structures, compounds, and sequences from pool for filtered IDs.
    
    Args:
        filtered_ids: List of IDs that passed the filter
        pool_folder: Pool output folder containing all data types
        output_folder: Where to copy the filtered data
        
    Returns:
        Dictionary mapping data type to list of extracted file paths
    """
    extracted_files = {"structures": [], "compounds": [], "sequences": []}
    
    # Filter puts structures/compounds directly in output folder (not subdirectories)
    os.makedirs(output_folder, exist_ok=True)
    
    print(f"\nPool mode: Extracting data for {len(filtered_ids)} filtered IDs")
    
    for i, selected_id in enumerate(filtered_ids):
        # Extract structures (PDB/CIF files)
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
                source = matches[0]
                dest = os.path.join(output_folder, f"{selected_id}.pdb")
                try:
                    shutil.copy2(source, dest)
                    extracted_files["structures"].append(dest)
                    print(f"Extracted structure: {source} -> {dest}")
                    break
                except Exception as e:
                    print(f"Warning: Could not copy structure {source}: {e}")
        
        # Extract compounds (SDF files)
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
                dest = os.path.join(output_folder, f"{selected_id}.sdf")
                try:
                    shutil.copy2(source, dest)
                    extracted_files["compounds"].append(dest)
                    print(f"Extracted compound: {source} -> {dest}")
                    break
                except Exception as e:
                    print(f"Warning: Could not copy compound {source}: {e}")
    
    # Extract all tables from pool folder (not just sequences)
    table_patterns = [
        os.path.join(pool_folder, "*.csv"),
        os.path.join(pool_folder, "**", "*.csv")
    ]
    
    processed_files = set()  # Track processed files to avoid duplicates
    
    for pattern in table_patterns:
        matches = glob.glob(pattern, recursive=True)
        for source_csv in matches:
            if source_csv in processed_files:
                continue
            processed_files.add(source_csv)
            
            try:
                df = pd.read_csv(source_csv)
                filename = os.path.basename(source_csv)
                
                # Only process files that have an 'id' column (actual tables)
                # Skip individual MSA files and other non-table CSVs
                if 'id' not in df.columns:
                    continue
                
                if filtered_ids:
                    # Filter rows matching the filtered IDs
                    matching_rows = df[df['id'].isin(filtered_ids)]
                    dest = os.path.join(output_folder, filename)
                    matching_rows.to_csv(dest, index=False)
                    print(f"Filtered table: {filename} ({len(matching_rows)} rows)")
                else:
                    # Create empty table with same columns
                    empty_df = pd.DataFrame(columns=df.columns)
                    dest = os.path.join(output_folder, filename)
                    empty_df.to_csv(dest, index=False)
                    print(f"Created empty table: {filename} with {len(df.columns)} columns")
                    
            except Exception as e:
                print(f"Warning: Could not process table {source_csv}: {e}")
    
    return extracted_files


def apply_filter(config_data: Dict[str, Any]) -> None:
    """
    Apply expression-based filter to a CSV file.
    
    Args:
        config_data: Configuration dictionary with filter parameters
    """
    input_csv = config_data['input_csv']
    expression = config_data['expression']
    max_items = config_data.get('max_items')
    sort_by = config_data.get('sort_by')
    sort_ascending = config_data.get('sort_ascending', True)
    output_csv = config_data['output_csv']
    
    # Pool mode parameters
    use_pool_mode = config_data.get('use_pool_mode', False)
    pool_output_folder = config_data.get('pool_output_folder')
    
    print(f"Applying filter: {expression}")
    print(f"Input: {input_csv}")
    if use_pool_mode:
        print(f"Pool mode: {pool_output_folder}")
    
    # Check input file exists
    if not os.path.exists(input_csv):
        raise FileNotFoundError(f"Input CSV not found: {input_csv}")
    
    # Load input dataframe
    print("Loading input CSV...")
    try:
        df = pd.read_csv(input_csv)
        print(f"Loaded dataframe: {df.shape}")
        print(f"Columns: {list(df.columns)}")
        
    except Exception as e:
        raise ValueError(f"Error loading CSV file: {e}")
    
    if df.empty:
        print("Warning: Input dataframe is empty")
        df.to_csv(output_csv, index=False)
        return
    
    # Apply the filter expression
    print(f"Applying expression: {expression}")
    try:
        # Use pandas query method for safe expression evaluation
        filtered_df = df.query(expression)
        print(f"After filtering: {filtered_df.shape}")
        
    except Exception as e:
        # Try to provide helpful error messages
        if "name '" in str(e) and "' is not defined" in str(e):
            available_cols = ", ".join(df.columns)
            raise ValueError(f"Column not found in expression '{expression}'. "
                           f"Available columns: {available_cols}")
        else:
            raise ValueError(f"Error evaluating expression '{expression}': {e}")
    
    if filtered_df.empty:
        print("Warning: No rows passed the filter")
        filtered_df.to_csv(output_csv, index=False)
        if use_pool_mode and pool_output_folder:
            print("No structures to copy (empty filter result)")
            # Still create empty table files with correct headers
            output_dir = os.path.dirname(output_csv)
            extract_pool_data_for_filtered_ids([], pool_output_folder, output_dir)
            # Create missing_ids.csv for empty filter result
            create_missing_ids_csv(df, filtered_df, pool_output_folder, output_dir)
        return
    
    # Apply sorting if specified
    if sort_by:
        if sort_by not in filtered_df.columns:
            print(f"Warning: Sort column '{sort_by}' not found. Available: {list(filtered_df.columns)}")
        else:
            print(f"Sorting by: {sort_by} ({'ascending' if sort_ascending else 'descending'})")
            filtered_df = filtered_df.sort_values(by=sort_by, ascending=sort_ascending)
    
    # Apply max_items limit if specified
    if max_items and max_items < len(filtered_df):
        print(f"Limiting to top {max_items} items")
        filtered_df = filtered_df.head(max_items)
    
    # Create output directory if needed
    os.makedirs(os.path.dirname(output_csv), exist_ok=True)
    
    # Save filtered results
    filtered_df.to_csv(output_csv, index=False)
    
    # Summary
    original_count = len(df)
    filtered_count = len(filtered_df)
    pass_rate = (filtered_count / original_count) * 100 if original_count > 0 else 0
    
    print(f"\\nFiltering summary:")
    print(f"Original rows: {original_count}")
    print(f"Filtered rows: {filtered_count}")
    print(f"Pass rate: {pass_rate:.1f}%")
    print(f"Output saved: {output_csv}")
    
    # Show column summary
    print(f"\\nOutput columns ({len(filtered_df.columns)}): {', '.join(filtered_df.columns)}")
    
    # Show sample of results if available
    if not filtered_df.empty:
        print("\\nFirst few results:")
        print(filtered_df.head(3).to_string(index=False))
        
        if len(filtered_df) > 3:
            print(f"... and {len(filtered_df) - 3} more rows")
    
    # Handle pool mode: copy structures, compounds, sequences for filtered IDs
    if use_pool_mode and pool_output_folder and not filtered_df.empty:
        # Get filtered IDs for structure copying
        filtered_ids = []
        if 'id' in filtered_df.columns:
            filtered_ids = filtered_df['id'].astype(str).tolist()
        else:
            # Try other ID columns
            id_columns = [col for col in filtered_df.columns if 'id' in col.lower()]
            if id_columns:
                filtered_ids = filtered_df[id_columns[0]].astype(str).tolist()
        
        if filtered_ids:
            output_dir = os.path.dirname(output_csv)
            extracted_files = extract_pool_data_for_filtered_ids(filtered_ids, pool_output_folder, output_dir)
            
            print(f"\\nPool mode summary:")
            print(f"Filtered IDs: {len(filtered_ids)}")
            for data_type, files in extracted_files.items():
                print(f"Extracted {data_type}: {len(files)} files")
        else:
            print("Warning: Could not find ID column for pool mode structure copying")
    
    # Always create missing.csv (shows which IDs were filtered out)
    output_dir = os.path.dirname(output_csv)
    create_missing_ids_csv(df, filtered_df, pool_output_folder or "", output_dir)


def create_missing_ids_csv(original_df: pd.DataFrame, filtered_df: pd.DataFrame, 
                          pool_output_folder: str, output_folder: str) -> None:
    """
    Create missing_ids.csv with IDs that were filtered out and their expected file paths.
    
    Args:
        original_df: Original dataframe before filtering
        filtered_df: Dataframe after filtering
        pool_output_folder: Pool output folder containing original files
        output_folder: Filter output folder
    """
    if 'id' not in original_df.columns:
        print("Warning: No 'id' column found, cannot create missing_ids.csv")
        return
    
    # Get IDs that were filtered out
    all_ids = set(original_df['id'].astype(str))
    passed_ids = set(filtered_df['id'].astype(str)) if not filtered_df.empty else set()
    missing_ids = list(all_ids - passed_ids)
    
    if not missing_ids:
        # Create empty file with headers
        missing_df = pd.DataFrame(columns=['id', 'structure', 'msa'])
    else:
        # Create rows for missing IDs with their expected file paths
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
    
    # Save missing.csv (consistent naming with other tools)
    missing_csv = os.path.join(output_folder, "missing.csv")
    missing_df.to_csv(missing_csv, index=False)
    print(f"Created missing.csv with {len(missing_ids)} filtered out IDs")


def validate_expression_syntax(expression: str, columns: List[str]) -> None:
    """
    Validate that expression syntax is safe and uses valid column names.
    
    Args:
        expression: Filter expression to validate
        columns: List of available column names
    """
    # Check for dangerous patterns
    dangerous_patterns = [
        'import ', '__', 'exec(', 'eval(', 'os.', 'sys.', 'subprocess', 
        'open(', 'file(', 'input(', 'raw_input('
    ]
    
    expr_lower = expression.lower()
    for pattern in dangerous_patterns:
        if pattern in expr_lower:
            raise ValueError(f"Dangerous pattern '{pattern}' not allowed in expression")
    
    # Basic syntax check using pandas query on empty dataframe
    try:
        import pandas as pd
        test_df = pd.DataFrame({col: [] for col in columns})
        test_df.query(expression)
    except Exception as e:
        raise ValueError(f"Invalid expression syntax: {e}")


def main():
    parser = argparse.ArgumentParser(description='Apply expression-based filter to CSV')
    parser.add_argument('--config', required=True, help='JSON config file with filter parameters')
    
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
    required_params = ['input_csv', 'expression', 'output_csv']
    for param in required_params:
        if param not in config_data:
            print(f"Error: Missing required parameter: {param}")
            sys.exit(1)
    
    try:
        apply_filter(config_data)
        print("\nFiltering completed successfully!")
        
    except Exception as e:
        print(f"Error applying filter: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()