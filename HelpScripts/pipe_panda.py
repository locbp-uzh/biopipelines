#!/usr/bin/env python3
"""
Runtime helper script for Panda table transformations.

This script executes pandas-style operations on CSV files including
filtering, sorting, merging, concatenation, and calculated columns.
"""

import os
import sys
import argparse
import json
import glob
import shutil
import pandas as pd
from typing import Dict, List, Any, Optional


def validate_expression(expr: str) -> None:
    """
    Validate that an expression is safe for pandas eval/query.

    Args:
        expr: Expression to validate

    Raises:
        ValueError: If expression contains dangerous patterns
    """
    dangerous_patterns = [
        'import ', '__', 'exec(', 'eval(', 'os.', 'sys.', 'subprocess',
        'open(', 'file(', 'input(', 'raw_input('
    ]

    expr_lower = expr.lower()
    for pattern in dangerous_patterns:
        if pattern in expr_lower:
            raise ValueError(f"Dangerous pattern '{pattern}' not allowed in expression")


def execute_filter(df: pd.DataFrame, params: Dict[str, Any]) -> pd.DataFrame:
    """Apply filter operation using pandas query."""
    expr = params.get("expr", "")
    validate_expression(expr)
    print(f"  Filter: {expr}")
    result = df.query(expr)
    print(f"    -> {len(df)} rows -> {len(result)} rows")
    return result


def execute_head(df: pd.DataFrame, params: Dict[str, Any]) -> pd.DataFrame:
    """Keep first N rows."""
    n = params.get("n", 5)
    print(f"  Head: {n}")
    return df.head(n)


def execute_tail(df: pd.DataFrame, params: Dict[str, Any]) -> pd.DataFrame:
    """Keep last N rows."""
    n = params.get("n", 5)
    print(f"  Tail: {n}")
    return df.tail(n)


def execute_sample(df: pd.DataFrame, params: Dict[str, Any]) -> pd.DataFrame:
    """Randomly sample rows."""
    n = params.get("n")
    frac = params.get("frac")
    random_state = params.get("random_state")

    if n is not None:
        print(f"  Sample: n={n}")
        return df.sample(n=min(n, len(df)), random_state=random_state)
    elif frac is not None:
        print(f"  Sample: frac={frac}")
        return df.sample(frac=frac, random_state=random_state)
    else:
        return df


def execute_drop_duplicates(df: pd.DataFrame, params: Dict[str, Any]) -> pd.DataFrame:
    """Remove duplicate rows."""
    subset = params.get("subset")
    keep = params.get("keep", "first")

    if subset:
        if isinstance(subset, str):
            subset = [subset]
        print(f"  Drop duplicates: subset={subset}, keep={keep}")
    else:
        print(f"  Drop duplicates: all columns, keep={keep}")

    return df.drop_duplicates(subset=subset, keep=keep)


def execute_sort(df: pd.DataFrame, params: Dict[str, Any]) -> pd.DataFrame:
    """Sort rows by column(s)."""
    by = params.get("by")
    ascending = params.get("ascending", True)

    if by:
        if isinstance(by, str):
            by = [by]
        print(f"  Sort: by={by}, ascending={ascending}")
        return df.sort_values(by=by, ascending=ascending)
    return df


def execute_rank(df: pd.DataFrame, params: Dict[str, Any]) -> pd.DataFrame:
    """Assign rank based on a column."""
    by = params.get("by")
    prefix = params.get("prefix", "rank_")
    ascending = params.get("ascending", True)

    if by and by in df.columns:
        rank_col = f"{prefix}{by}"
        print(f"  Rank: by={by}, column={rank_col}")
        df = df.copy()
        df[rank_col] = df[by].rank(ascending=ascending, method='dense').astype(int)
    return df


def execute_select_columns(df: pd.DataFrame, params: Dict[str, Any]) -> pd.DataFrame:
    """Keep only specified columns."""
    cols = params.get("cols", [])
    existing_cols = [c for c in cols if c in df.columns]
    print(f"  Select columns: {existing_cols}")
    return df[existing_cols]


def execute_drop_columns(df: pd.DataFrame, params: Dict[str, Any]) -> pd.DataFrame:
    """Remove specified columns."""
    cols = params.get("cols", [])
    existing_cols = [c for c in cols if c in df.columns]
    print(f"  Drop columns: {existing_cols}")
    return df.drop(columns=existing_cols)


def execute_rename(df: pd.DataFrame, params: Dict[str, Any]) -> pd.DataFrame:
    """Rename columns."""
    mapping = params.get("mapping", {})
    existing_mapping = {k: v for k, v in mapping.items() if k in df.columns}
    print(f"  Rename: {existing_mapping}")
    return df.rename(columns=existing_mapping)


def execute_calculate(df: pd.DataFrame, params: Dict[str, Any]) -> pd.DataFrame:
    """Add calculated columns."""
    exprs = params.get("exprs", {})
    df = df.copy()

    for col_name, expr in exprs.items():
        validate_expression(expr)
        try:
            df[col_name] = df.eval(expr)
            print(f"  Calculate: {col_name} = {expr}")
        except Exception as e:
            print(f"  Warning: Failed to calculate '{col_name}': {e}")
            df[col_name] = float('nan')

    return df


def execute_fillna(df: pd.DataFrame, params: Dict[str, Any]) -> pd.DataFrame:
    """Fill missing values."""
    value = params.get("value")
    column = params.get("column")

    if column and column in df.columns:
        print(f"  Fillna: column={column}, value={value}")
        df = df.copy()
        df[column] = df[column].fillna(value)
    elif value is not None:
        print(f"  Fillna: all columns, value={value}")
        df = df.fillna(value)

    return df


def execute_merge(dataframes: List[pd.DataFrame], params: Dict[str, Any]) -> pd.DataFrame:
    """Merge multiple dataframes horizontally."""
    on = params.get("on", "id")
    how = params.get("how", "outer")
    prefixes = params.get("prefixes", [])
    id_map = params.get("id_map", {})

    print(f"  Merge: on={on}, how={how}")

    # Apply prefixes to each dataframe
    renamed_dfs = []
    for i, df in enumerate(dataframes):
        df = df.copy()

        # Apply ID mapping
        if id_map:
            reverse_map = {}
            for new_id, old_id_list in id_map.items():
                for old_id in old_id_list:
                    reverse_map[old_id] = new_id

            if on in df.columns:
                df[on] = df[on].replace(reverse_map)
                print(f"    Applied ID mapping to dataframe {i}")

        # Apply prefix
        if prefixes and i < len(prefixes) and prefixes[i]:
            prefix = prefixes[i]
            rename_dict = {col: f"{prefix}{col}" for col in df.columns if col != on}
            df = df.rename(columns=rename_dict)
            print(f"    Applied prefix '{prefix}' to dataframe {i}")

        renamed_dfs.append(df)

    # Perform merge
    result = renamed_dfs[0]
    for i, df in enumerate(renamed_dfs[1:], 1):
        # Check for overlapping columns (excluding merge key)
        overlap = set(result.columns) & set(df.columns) - {on}
        if overlap:
            print(f"    Warning: Overlapping columns: {overlap}")
            df = df.drop(columns=list(overlap))

        before_rows = len(result)
        result = pd.merge(result, df, on=on, how=how)
        print(f"    Merged dataframe {i+1}: {before_rows} -> {len(result)} rows")

    return result


def execute_concat(dataframes: List[pd.DataFrame], params: Dict[str, Any]) -> pd.DataFrame:
    """Concatenate multiple dataframes vertically."""
    fill = params.get("fill")
    add_source = params.get("add_source", True)

    print(f"  Concat: fill={fill}, add_source={add_source}")

    if add_source:
        # Add source_table column to each dataframe
        for i, df in enumerate(dataframes):
            dataframes[i] = df.copy()
            dataframes[i]['source_table'] = i

    # Get common columns if fill is None
    if fill is None:
        common_cols = set(dataframes[0].columns)
        for df in dataframes[1:]:
            common_cols &= set(df.columns)
        print(f"    Common columns: {list(common_cols)}")

        # Filter to common columns
        dataframes = [df[list(common_cols)] for df in dataframes]
    else:
        # Keep all columns, fill missing with specified value
        all_cols = set()
        for df in dataframes:
            all_cols.update(df.columns)

        for i, df in enumerate(dataframes):
            missing_cols = all_cols - set(df.columns)
            if missing_cols:
                df = df.copy()
                for col in missing_cols:
                    df[col] = fill
                dataframes[i] = df

    result = pd.concat(dataframes, ignore_index=True)
    print(f"    Concatenated {len(dataframes)} tables: {len(result)} total rows")
    return result


def execute_groupby(df: pd.DataFrame, params: Dict[str, Any]) -> pd.DataFrame:
    """Group by column(s) and aggregate."""
    by = params.get("by")
    agg = params.get("agg", {})

    if not by or not agg:
        return df

    if isinstance(by, str):
        by = [by]

    print(f"  Groupby: by={by}, agg={agg}")
    return df.groupby(by, as_index=False).agg(agg)


def execute_pivot(df: pd.DataFrame, params: Dict[str, Any]) -> pd.DataFrame:
    """Pivot table."""
    index = params.get("index")
    columns = params.get("columns")
    values = params.get("values")
    aggfunc = params.get("aggfunc", "first")

    if not all([index, columns, values]):
        return df

    print(f"  Pivot: index={index}, columns={columns}, values={values}")
    return df.pivot_table(index=index, columns=columns, values=values,
                          aggfunc=aggfunc).reset_index()


def execute_melt(df: pd.DataFrame, params: Dict[str, Any]) -> pd.DataFrame:
    """Unpivot table."""
    id_vars = params.get("id_vars")
    value_vars = params.get("value_vars")
    var_name = params.get("var_name", "variable")
    value_name = params.get("value_name", "value")

    if not id_vars:
        return df

    if isinstance(id_vars, str):
        id_vars = [id_vars]

    print(f"  Melt: id_vars={id_vars}, value_vars={value_vars}")
    return pd.melt(df, id_vars=id_vars, value_vars=value_vars,
                   var_name=var_name, value_name=value_name)


def execute_average_by_source(df: pd.DataFrame, params: Dict[str, Any]) -> pd.DataFrame:
    """
    Average all numeric columns per source table.

    This replaces AverageByTable functionality. After concat with add_source=True,
    this groups by source_table and computes the mean of all numeric columns.
    """
    source_col = params.get("source_col", "source_table")

    if source_col not in df.columns:
        print(f"  Warning: source column '{source_col}' not found. Run concat with add_source=True first.")
        return df

    print(f"  Average by source: grouping by '{source_col}'")

    # Get numeric columns only
    numeric_cols = df.select_dtypes(include=['number']).columns.tolist()

    # Remove the source column from numeric cols if present
    if source_col in numeric_cols:
        numeric_cols.remove(source_col)

    if not numeric_cols:
        print(f"  Warning: No numeric columns found to average")
        return df.groupby(source_col, as_index=False).first()

    print(f"  Averaging columns: {numeric_cols}")

    # Group by source and compute mean of numeric columns
    result = df.groupby(source_col, as_index=False)[numeric_cols].mean()

    # Add a table_name column based on source index
    result['table_name'] = [f"table_{i}" for i in result[source_col]]

    print(f"    -> {len(result)} rows (one per source table)")
    return result


def execute_operation(df_or_dfs, operation: Dict[str, Any], is_multi_table: bool):
    """
    Execute a single operation.

    Args:
        df_or_dfs: Single DataFrame or list of DataFrames
        operation: Operation dictionary with 'type' and 'params'
        is_multi_table: Whether this is a multi-table context

    Returns:
        Resulting DataFrame
    """
    op_type = operation.get("type")
    params = operation.get("params", {})

    # Multi-table operations
    if op_type == "merge":
        if not isinstance(df_or_dfs, list):
            raise ValueError("Merge operation requires multiple dataframes")
        return execute_merge(df_or_dfs, params)

    if op_type == "concat":
        if not isinstance(df_or_dfs, list):
            raise ValueError("Concat operation requires multiple dataframes")
        return execute_concat(df_or_dfs, params)

    # Single-table operations - ensure we have a single DataFrame
    if isinstance(df_or_dfs, list):
        if len(df_or_dfs) == 1:
            df = df_or_dfs[0]
        else:
            raise ValueError(f"Operation '{op_type}' requires a single dataframe. "
                           f"Use 'merge' or 'concat' first to combine tables.")
    else:
        df = df_or_dfs

    operation_funcs = {
        "filter": execute_filter,
        "head": execute_head,
        "tail": execute_tail,
        "sample": execute_sample,
        "drop_duplicates": execute_drop_duplicates,
        "sort": execute_sort,
        "rank": execute_rank,
        "select_columns": execute_select_columns,
        "drop_columns": execute_drop_columns,
        "rename": execute_rename,
        "calculate": execute_calculate,
        "fillna": execute_fillna,
        "groupby": execute_groupby,
        "pivot": execute_pivot,
        "melt": execute_melt,
        "average_by_source": execute_average_by_source,
    }

    func = operation_funcs.get(op_type)
    if func:
        return func(df, params)
    else:
        print(f"  Warning: Unknown operation '{op_type}', skipping")
        return df


def extract_pool_data_for_filtered_ids(filtered_ids: List[str], pool_folder: str,
                                        output_folder: str,
                                        rename_map: Optional[Dict[str, str]] = None) -> Dict[str, List[str]]:
    """
    Extract structures, compounds, and sequences from pool for filtered IDs.

    Args:
        filtered_ids: List of IDs that passed filtering (original IDs for lookup)
        pool_folder: Pool output folder containing all data types
        output_folder: Where to copy the filtered data
        rename_map: Optional mapping from original ID to new ID for renaming output files

    Returns:
        Dictionary mapping data type to list of extracted file paths
    """
    extracted_files = {"structures": [], "compounds": [], "sequences": []}

    os.makedirs(output_folder, exist_ok=True)

    print(f"\nPool mode: Extracting data for {len(filtered_ids)} filtered IDs")
    if rename_map:
        print(f"Renaming files according to rename map")

    for selected_id in filtered_ids:
        # Determine output ID (renamed or original)
        output_id = rename_map.get(str(selected_id), selected_id) if rename_map else selected_id

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
                dest = os.path.join(output_folder, f"{output_id}.pdb")
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
                dest = os.path.join(output_folder, f"{output_id}.sdf")
                try:
                    shutil.copy2(source, dest)
                    extracted_files["compounds"].append(dest)
                    print(f"Extracted compound: {source} -> {dest}")
                    break
                except Exception as e:
                    print(f"Warning: Could not copy compound {source}: {e}")

    # Extract and filter tables from pool folder
    table_patterns = [
        os.path.join(pool_folder, "*.csv"),
        os.path.join(pool_folder, "**", "*.csv")
    ]

    processed_files = set()

    for pattern in table_patterns:
        matches = glob.glob(pattern, recursive=True)
        for source_csv in matches:
            if source_csv in processed_files:
                continue
            processed_files.add(source_csv)

            try:
                df = pd.read_csv(source_csv)
                filename = os.path.basename(source_csv)

                # Only process files that have an 'id' column
                if 'id' not in df.columns:
                    continue

                if filtered_ids:
                    matching_rows = df[df['id'].isin(filtered_ids)].copy()
                    # Apply rename mapping to the id column if provided
                    if rename_map and not matching_rows.empty:
                        matching_rows['original_id'] = matching_rows['id']
                        matching_rows['id'] = matching_rows['id'].astype(str).map(
                            lambda x: rename_map.get(x, x)
                        )
                    dest = os.path.join(output_folder, filename)
                    matching_rows.to_csv(dest, index=False)
                    print(f"Filtered table: {filename} ({len(matching_rows)} rows)")
                else:
                    empty_df = pd.DataFrame(columns=df.columns)
                    dest = os.path.join(output_folder, filename)
                    empty_df.to_csv(dest, index=False)
                    print(f"Created empty table: {filename}")

            except Exception as e:
                print(f"Warning: Could not process table {source_csv}: {e}")

    return extracted_files


def create_missing_csv(original_ids: List[str], filtered_ids: List[str],
                       output_folder: str) -> None:
    """
    Create missing.csv with IDs that were filtered out.

    Args:
        original_ids: All original IDs
        filtered_ids: IDs that passed filtering
        output_folder: Output folder for missing.csv
    """
    all_ids = set(str(i) for i in original_ids)
    passed_ids = set(str(i) for i in filtered_ids)
    missing_ids = list(all_ids - passed_ids)

    if not missing_ids:
        missing_df = pd.DataFrame(columns=['id', 'structure', 'msa'])
    else:
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

    missing_csv = os.path.join(output_folder, "missing.csv")
    missing_df.to_csv(missing_csv, index=False)
    print(f"Created missing.csv with {len(missing_ids)} filtered out IDs")


def run_panda(config_data: Dict[str, Any]) -> None:
    """
    Execute Panda operations on CSV files.

    Args:
        config_data: Configuration dictionary with input files, operations, and settings
    """
    input_csvs = config_data['input_csvs']
    operations = config_data.get('operations', [])
    output_csv = config_data['output_csv']
    use_pool_mode = config_data.get('use_pool_mode', False)
    pool_output_folder = config_data.get('pool_output_folder')
    rename = config_data.get('rename')

    print(f"Panda: Processing {len(input_csvs)} input files")
    print(f"Operations: {len(operations)}")
    if rename:
        print(f"Rename: {rename}_N")

    # Load input dataframes
    dataframes = []
    original_ids = []

    for csv_path in input_csvs:
        if not os.path.exists(csv_path):
            raise FileNotFoundError(f"Input file not found: {csv_path}")

        print(f"Loading: {csv_path}")
        df = pd.read_csv(csv_path)
        print(f"  Shape: {df.shape}")
        print(f"  Columns: {list(df.columns)}")
        dataframes.append(df)

        # Collect original IDs for missing.csv
        if 'id' in df.columns:
            original_ids.extend(df['id'].tolist())

    if not dataframes:
        raise ValueError("No valid input dataframes found")

    # Start with either single df or list of dfs
    is_multi_table = len(dataframes) > 1
    current = dataframes if is_multi_table else dataframes[0]

    # Execute operations sequentially
    print("\nExecuting operations:")
    for i, operation in enumerate(operations):
        print(f"\n[{i+1}/{len(operations)}] {operation['type']}")
        current = execute_operation(current, operation, is_multi_table)

        # After merge/concat, we have a single dataframe
        if operation['type'] in ('merge', 'concat'):
            is_multi_table = False

    # Get final result
    if isinstance(current, list):
        if len(current) == 1:
            result_df = current[0]
        else:
            # If we still have multiple dataframes, concat them
            print("\nAuto-concatenating remaining dataframes")
            result_df = pd.concat(current, ignore_index=True)
    else:
        result_df = current

    # Apply rename if specified - creates new IDs based on row order
    original_to_new_id = {}
    if rename and 'id' in result_df.columns:
        original_ids_ordered = result_df['id'].tolist()
        new_ids = [f"{rename}_{i+1}" for i in range(len(result_df))]
        original_to_new_id = dict(zip([str(x) for x in original_ids_ordered], new_ids))
        result_df = result_df.copy()
        result_df['original_id'] = result_df['id']
        result_df['id'] = new_ids
        print(f"\nApplied rename: {rename}_1 to {rename}_{len(result_df)}")

    # Create output directory
    os.makedirs(os.path.dirname(output_csv), exist_ok=True)

    # Save result
    result_df.to_csv(output_csv, index=False)

    print(f"\nFinal result: {result_df.shape}")
    print(f"Columns: {list(result_df.columns)}")
    print(f"Output saved: {output_csv}")

    # Show sample
    if not result_df.empty:
        print("\nFirst few results:")
        print(result_df.head(3).to_string(index=False))
        if len(result_df) > 3:
            print(f"... and {len(result_df) - 3} more rows")

    # Get filtered IDs for pool mode and missing.csv
    # Use original_id for lookup, new id for output files
    filtered_ids = []
    filtered_ids_for_lookup = []
    if 'id' in result_df.columns:
        filtered_ids = result_df['id'].astype(str).tolist()
        if 'original_id' in result_df.columns:
            filtered_ids_for_lookup = result_df['original_id'].astype(str).tolist()
        else:
            filtered_ids_for_lookup = filtered_ids
    else:
        id_columns = [col for col in result_df.columns if 'id' in col.lower()]
        if id_columns:
            filtered_ids = result_df[id_columns[0]].astype(str).tolist()
            filtered_ids_for_lookup = filtered_ids

    # Handle pool mode
    if use_pool_mode and pool_output_folder and filtered_ids_for_lookup:
        output_dir = os.path.dirname(output_csv)
        extracted = extract_pool_data_for_filtered_ids(
            filtered_ids_for_lookup, pool_output_folder, output_dir,
            rename_map=original_to_new_id if rename else None
        )

        print(f"\nPool mode summary:")
        print(f"Filtered IDs: {len(filtered_ids)}")
        for data_type, files in extracted.items():
            print(f"Extracted {data_type}: {len(files)} files")

    # Create missing.csv
    if original_ids:
        output_dir = os.path.dirname(output_csv)
        create_missing_csv(original_ids, filtered_ids, output_dir)

    print("\nPanda operations completed successfully!")


def main():
    parser = argparse.ArgumentParser(description='Execute Panda table transformations')
    parser.add_argument('--config', required=True, help='JSON config file with operations')

    args = parser.parse_args()

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
    required_params = ['input_csvs', 'output_csv']
    for param in required_params:
        if param not in config_data:
            print(f"Error: Missing required parameter: {param}")
            sys.exit(1)

    try:
        run_panda(config_data)

    except Exception as e:
        print(f"Error executing Panda operations: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
