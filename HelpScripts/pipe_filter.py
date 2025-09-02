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
from typing import Dict, List, Any, Optional


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
    
    print(f"Applying filter: {expression}")
    print(f"Input: {input_csv}")
    
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
        print("\\nFiltering completed successfully!")
        
    except Exception as e:
        print(f"Error applying filter: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()