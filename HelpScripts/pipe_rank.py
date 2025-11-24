#!/usr/bin/env python3
"""
Runtime helper script for ranking entries by a metric.

This script ranks CSV entries by a specified metric (column or expression),
renames IDs to sequential ranks, and optionally copies structures in ranked order.
"""

import os
import sys
import argparse
import json
import pandas as pd
import shutil
import glob
import re
from typing import Dict, List, Any, Optional


def is_expression(metric: str) -> bool:
    """
    Detect if metric is an expression or a simple column name.

    Args:
        metric: Metric string to check

    Returns:
        True if metric contains operators (is an expression), False if simple column name
    """
    # Check for arithmetic operators or parentheses
    operators = ['+', '-', '*', '/', '(', ')']
    return any(op in metric for op in operators)


def extract_variables_from_expression(metric: str) -> List[str]:
    """
    Extract variable names (column references) from metric expression.

    Args:
        metric: Metric expression

    Returns:
        List of variable names referenced in expression
    """
    # Match identifiers (potential column names)
    potential_vars = re.findall(r'\b[a-zA-Z_][a-zA-Z0-9_]*\b', metric)

    # Filter out Python keywords
    keywords = {'and', 'or', 'not', 'in', 'is', 'if', 'else', 'for', 'while', 'def', 'class'}
    variables = [v for v in potential_vars if v not in keywords]

    return list(set(variables))  # Remove duplicates


def extract_pool_data_for_ranked_ids(ranked_ids: List[str], pool_folder: str,
                                     output_folder: str, prefix: str, num_digits: int) -> Dict[str, List[str]]:
    """
    Extract structures, compounds, and sequences from pool for ranked IDs.
    Copy them in ranked order with sequential naming.

    Args:
        ranked_ids: List of IDs in ranked order
        pool_folder: Pool output folder containing all data types
        output_folder: Where to copy the ranked data
        prefix: Prefix for renamed files (e.g., "rank" -> rank_001.pdb, rank_002.pdb)
        num_digits: Number of digits for zero-padding (determined from pool structure count)

    Returns:
        Dictionary mapping data type to list of extracted file paths
    """
    extracted_files = {"structures": [], "compounds": [], "sequences": []}
    missing_ranked_ids = []  # Track ranked IDs that couldn't be created

    os.makedirs(output_folder, exist_ok=True)

    print(f"\nPool mode: Extracting data for {len(ranked_ids)} ranked IDs")
    print(f"Using {num_digits}-digit zero-padding for ranked IDs")

    for rank_idx, source_id in enumerate(ranked_ids, start=1):
        ranked_id = f"{prefix}_{rank_idx:0{num_digits}d}"

        # Extract structures (PDB/CIF files)
        structure_patterns = [
            os.path.join(pool_folder, f"{source_id}.pdb"),
            os.path.join(pool_folder, f"*{source_id}*.pdb"),
            os.path.join(pool_folder, f"{source_id}.cif"),
            os.path.join(pool_folder, f"*{source_id}*.cif"),
            os.path.join(pool_folder, "**", f"{source_id}.pdb"),
            os.path.join(pool_folder, "**", f"*{source_id}*.pdb"),
        ]

        structure_found = False
        for pattern in structure_patterns:
            matches = glob.glob(pattern, recursive=True)
            if matches:
                source = matches[0]
                ext = os.path.splitext(source)[1]
                dest = os.path.join(output_folder, f"{ranked_id}{ext}")
                try:
                    shutil.copy2(source, dest)
                    extracted_files["structures"].append(dest)
                    print(f"Extracted structure: {source_id} -> {ranked_id}{ext}")
                    structure_found = True
                    break
                except Exception as e:
                    print(f"Warning: Could not copy structure {source}: {e}")

        if not structure_found:
            print(f"Note: Structure not found for {source_id} (may have been filtered out)")
            # Track this ranked ID as missing
            missing_ranked_ids.append(ranked_id)

        # Extract compounds (SDF files)
        compound_patterns = [
            os.path.join(pool_folder, f"{source_id}.sdf"),
            os.path.join(pool_folder, f"*{source_id}*.sdf"),
            os.path.join(pool_folder, "**", f"{source_id}.sdf"),
            os.path.join(pool_folder, "**", f"*{source_id}*.sdf"),
        ]

        compound_found = False
        for pattern in compound_patterns:
            matches = glob.glob(pattern, recursive=True)
            if matches:
                source = matches[0]
                dest = os.path.join(output_folder, f"{ranked_id}.sdf")
                try:
                    shutil.copy2(source, dest)
                    extracted_files["compounds"].append(dest)
                    print(f"Extracted compound: {source_id} -> {ranked_id}.sdf")
                    compound_found = True
                    break
                except Exception as e:
                    print(f"Warning: Could not copy compound {source}: {e}")

        if not compound_found and not structure_found:
            # Only add to missing if BOTH structure and compound are missing
            # (already added above if structure not found)
            pass

    # Extract all tables from pool folder
    table_patterns = [
        os.path.join(pool_folder, "*.csv"),
        os.path.join(pool_folder, "**", "*.csv")
    ]

    processed_files = set()  # Track processed files to avoid duplicates
    upstream_missing_df = None  # Store upstream missing.csv if it exists

    for pattern in table_patterns:
        matches = glob.glob(pattern, recursive=True)
        for source_csv in matches:
            if source_csv in processed_files:
                continue
            processed_files.add(source_csv)

            try:
                df = pd.read_csv(source_csv)
                filename = os.path.basename(source_csv)

                # Special handling for missing.csv - store for later merging
                if filename == "missing.csv":
                    upstream_missing_df = df
                    print(f"Found upstream missing.csv ({len(df)} rows) - will merge with rank missing IDs")
                    continue

                # Only process files that have an 'id' column
                if 'id' not in df.columns:
                    continue

                if ranked_ids:
                    # Filter and reorder rows to match ranked order
                    # Create a mapping from source_id to rank
                    id_to_rank = {source_id: rank_idx for rank_idx, source_id in enumerate(ranked_ids, start=1)}

                    # Filter matching rows
                    matching_rows = df[df['id'].isin(ranked_ids)].copy()

                    if len(matching_rows) == 0:
                        # No matching rows - create empty table with same columns
                        empty_df = pd.DataFrame(columns=df.columns)
                        dest = os.path.join(output_folder, filename)
                        empty_df.to_csv(dest, index=False)
                        print(f"Created empty table: {filename} (no matching IDs in ranked set)")
                        continue

                    # Add rank column for sorting
                    matching_rows['_rank_order'] = matching_rows['id'].map(id_to_rank)

                    # Sort by rank
                    matching_rows = matching_rows.sort_values('_rank_order')

                    # Remove temporary rank column
                    matching_rows = matching_rows.drop('_rank_order', axis=1)

                    # Rename IDs to ranked format with zero-padding (using num_digits passed to this function)
                    matching_rows['id'] = [f"{prefix}_{i:0{num_digits}d}" for i in range(1, len(matching_rows) + 1)]

                    dest = os.path.join(output_folder, filename)
                    matching_rows.to_csv(dest, index=False)
                    print(f"Ranked table: {filename} ({len(matching_rows)} rows)")
                else:
                    # Create empty table with same columns
                    empty_df = pd.DataFrame(columns=df.columns)
                    dest = os.path.join(output_folder, filename)
                    empty_df.to_csv(dest, index=False)
                    print(f"Created empty table: {filename} with {len(df.columns)} columns")

            except Exception as e:
                print(f"Warning: Could not process table {source_csv}: {e}")

    # Create or update missing.csv with ranked IDs that couldn't be created
    if missing_ranked_ids:
        missing_data = []
        for missing_id in missing_ranked_ids:
            structure_path = os.path.join(output_folder, f"{missing_id}.pdb")
            msa_path = os.path.join(output_folder, f"{missing_id}.csv")
            missing_data.append({
                'id': missing_id,
                'structure': structure_path,
                'msa': msa_path
            })
        rank_missing_df = pd.DataFrame(missing_data)

        # Merge with upstream missing.csv if it exists
        if upstream_missing_df is not None and len(upstream_missing_df) > 0:
            # Combine both dataframes
            combined_missing_df = pd.concat([upstream_missing_df, rank_missing_df], ignore_index=True)
        else:
            combined_missing_df = rank_missing_df

        missing_csv_path = os.path.join(output_folder, "missing.csv")
        combined_missing_df.to_csv(missing_csv_path, index=False)
        print(f"\nCreated missing.csv with {len(combined_missing_df)} total missing IDs:")
        print(f"  - {len(upstream_missing_df) if upstream_missing_df is not None else 0} from upstream filter")
        print(f"  - {len(missing_ranked_ids)} ranked IDs that couldn't be created")
    elif upstream_missing_df is not None:
        # Only upstream missing IDs exist - copy as-is
        missing_csv_path = os.path.join(output_folder, "missing.csv")
        upstream_missing_df.to_csv(missing_csv_path, index=False)
        print(f"\nCopied upstream missing.csv ({len(upstream_missing_df)} rows)")

    return extracted_files


def apply_ranking(config_data: Dict[str, Any]) -> None:
    """
    Apply ranking to a CSV file based on a metric.

    Args:
        config_data: Configuration dictionary with ranking parameters
    """
    input_csv = config_data['input_csv']
    metric = config_data['metric']
    ascending = config_data.get('ascending', False)
    prefix = config_data.get('prefix', 'rank')
    top = config_data.get('top')
    output_csv = config_data['output_csv']

    # Pool mode parameters
    use_pool_mode = config_data.get('use_pool_mode', False)
    pool_output_folder = config_data.get('pool_output_folder')

    # Get num_digits from config (determined at pipeline generation time)
    num_digits = config_data.get('num_digits', 1)

    print(f"Ranking by metric: {metric}")
    print(f"Input: {input_csv}")
    print(f"Sort order: {'ascending' if ascending else 'descending'}")
    print(f"Zero-padding: {num_digits} digits")
    if top:
        print(f"Limiting to top {top} entries")
    if use_pool_mode:
        print(f"Pool mode: {pool_output_folder}")

    # Check input file exists
    if not os.path.exists(input_csv):
        raise FileNotFoundError(f"Input CSV not found: {input_csv}")

    # Load input dataframe
    print("\nLoading input CSV...")
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

    # Check if 'id' column exists
    if 'id' not in df.columns:
        raise ValueError("Input CSV must have an 'id' column")

    # Determine if metric is an expression or simple column name
    metric_is_expression = is_expression(metric)

    if metric_is_expression:
        print(f"\nMetric is an expression: {metric}")

        # Extract variables from expression
        variables = extract_variables_from_expression(metric)
        print(f"Variables in expression: {variables}")

        # Validate that all variables exist in dataframe
        missing_vars = [v for v in variables if v not in df.columns]
        if missing_vars:
            available_cols = ", ".join(df.columns)
            raise ValueError(f"Variables not found in dataframe: {missing_vars}. "
                           f"Available columns: {available_cols}")

        # Compute metric using pandas eval
        try:
            df['metric'] = df.eval(metric)
            print(f"Computed metric column with range: [{df['metric'].min():.4f}, {df['metric'].max():.4f}]")
            sort_column = 'metric'
        except Exception as e:
            raise ValueError(f"Error evaluating metric expression '{metric}': {e}")
    else:
        print(f"\nMetric is a column name: {metric}")

        # Check if column exists
        if metric not in df.columns:
            available_cols = ", ".join(df.columns)
            raise ValueError(f"Metric column '{metric}' not found in dataframe. "
                           f"Available columns: {available_cols}")

        sort_column = metric
        print(f"Using existing column with range: [{df[metric].min():.4f}, {df[metric].max():.4f}]")

    # Sort by metric
    print(f"\nSorting by '{sort_column}' ({'ascending' if ascending else 'descending'})")
    sorted_df = df.sort_values(by=sort_column, ascending=ascending).copy()

    # Apply top N limit if specified
    if top and top < len(sorted_df):
        print(f"Limiting to top {top} entries")
        sorted_df = sorted_df.head(top)

    # Store original IDs
    sorted_df['source_id'] = sorted_df['id']

    # Rename IDs to ranked format with zero-padding (num_digits from config)
    sorted_df['id'] = [f"{prefix}_{i:0{num_digits}d}" for i in range(1, len(sorted_df) + 1)]

    # Reorder columns: id, source_id, metric (if computed), variables (if expression), rest
    column_order = ['id', 'source_id']

    if metric_is_expression:
        column_order.append('metric')
        # Add individual variable columns
        variables = extract_variables_from_expression(metric)
        for var in variables:
            if var in sorted_df.columns and var not in column_order:
                column_order.append(var)
    else:
        # Add metric column if not already in order
        if sort_column not in column_order:
            column_order.append(sort_column)

    # Add remaining columns
    for col in sorted_df.columns:
        if col not in column_order:
            column_order.append(col)

    # Reorder dataframe
    sorted_df = sorted_df[column_order]

    # Create output directory if needed
    os.makedirs(os.path.dirname(output_csv), exist_ok=True)

    # Save ranked results
    sorted_df.to_csv(output_csv, index=False)

    # Summary
    total_count = len(df)
    ranked_count = len(sorted_df)

    print(f"\nRanking summary:")
    print(f"Total input rows: {total_count}")
    print(f"Ranked output rows: {ranked_count}")
    print(f"Output saved: {output_csv}")

    # Show column summary
    print(f"\nOutput columns ({len(sorted_df.columns)}): {', '.join(sorted_df.columns)}")

    # Show sample of results
    if not sorted_df.empty:
        print("\nTop ranked entries:")
        display_df = sorted_df.head(5)
        # Show only key columns for readability
        key_cols = ['id', 'source_id']
        if metric_is_expression:
            key_cols.append('metric')
        else:
            key_cols.append(sort_column)
        display_cols = [col for col in key_cols if col in display_df.columns]
        print(display_df[display_cols].to_string(index=False))

        if len(sorted_df) > 5:
            print(f"... and {len(sorted_df) - 5} more rows")

    # Handle pool mode: copy structures, compounds, sequences in ranked order
    if use_pool_mode and pool_output_folder and not sorted_df.empty:
        # Get ranked source IDs for structure copying
        ranked_source_ids = sorted_df['source_id'].astype(str).tolist()

        if ranked_source_ids:
            output_dir = os.path.dirname(output_csv)
            extracted_files = extract_pool_data_for_ranked_ids(
                ranked_source_ids, pool_output_folder, output_dir, prefix, num_digits
            )

            print(f"\nPool mode summary:")
            print(f"Ranked IDs: {len(ranked_source_ids)}")
            for data_type, files in extracted_files.items():
                print(f"Extracted {data_type}: {len(files)} files")
        else:
            print("Warning: No source IDs found for pool mode structure copying")


def main():
    parser = argparse.ArgumentParser(description='Rank entries by metric')
    parser.add_argument('--config', required=True, help='JSON config file with ranking parameters')

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
    required_params = ['input_csv', 'metric', 'output_csv']
    for param in required_params:
        if param not in config_data:
            print(f"Error: Missing required parameter: {param}")
            sys.exit(1)

    try:
        apply_ranking(config_data)
        print("\nRanking completed successfully!")

    except Exception as e:
        print(f"Error applying ranking: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
