#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Runtime helper script for Panda table transformations.

This script executes pandas-style operations on CSV files including
filtering, sorting, merging, concatenation, and calculated columns.
"""

import os
import sys
import argparse
import json
import shutil
import pandas as pd
import numpy as np
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
    """Add calculated columns.

    Supports math functions (cos, sin, tan, sqrt, abs, log, exp, radians, degrees,
    arccos, arcsin, arctan, arctan2) via numpy, available directly in expressions.
    """
    exprs = params.get("exprs", {})
    df = df.copy()

    # Expose numpy math functions for use in eval expressions
    math_funcs = {
        'cos': np.cos, 'sin': np.sin, 'tan': np.tan,
        'arccos': np.arccos, 'arcsin': np.arcsin, 'arctan': np.arctan, 'arctan2': np.arctan2,
        'sqrt': np.sqrt, 'abs': np.abs, 'log': np.log, 'exp': np.exp,
        'radians': np.radians, 'degrees': np.degrees,
        'pi': np.pi,
    }

    for col_name, expr in exprs.items():
        validate_expression(expr)
        try:
            df[col_name] = df.eval(expr, engine='python', local_dict=math_funcs)
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
                                        rename_map: Optional[Dict[str, str]] = None,
                                        file_map: Optional[Dict[str, Dict[str, str]]] = None) -> Dict[str, List[str]]:
    """
    Extract data from a single pool for filtered IDs.

    Args:
        filtered_ids: List of IDs that passed filtering (original IDs for lookup)
        pool_folder: Pool output folder containing all data types
        output_folder: Where to copy the filtered data
        rename_map: Optional mapping from original ID to new ID for renaming output files
        file_map: Optional mapping {stream_name: {id: file_path}} for direct file lookup

    Returns:
        Dictionary mapping stream name to list of extracted file paths
    """
    extracted_files = {}

    os.makedirs(output_folder, exist_ok=True)

    print(f"\nPool mode: Extracting data for {len(filtered_ids)} filtered IDs from {pool_folder}")
    if rename_map:
        print(f"Renaming files according to rename map")

    if file_map:
        # Use file map - iterate over all streams dynamically
        print(f"Using file map for direct lookup (streams: {list(file_map.keys())})")

        for stream_name, id_to_file in file_map.items():
            extracted_files[stream_name] = []

            for selected_id in filtered_ids:
                if selected_id not in id_to_file:
                    continue

                source = id_to_file[selected_id]
                output_id = rename_map.get(str(selected_id), selected_id) if rename_map else selected_id
                ext = os.path.splitext(source)[1]
                dest = os.path.join(output_folder, f"{output_id}{ext}")

                shutil.copy2(source, dest)
                extracted_files[stream_name].append(dest)
                print(f"Extracted {stream_name}: {source} -> {dest}")
    else:
        raise ValueError("file_map is required for pool extraction")

    return extracted_files


def extract_from_multiple_pools(result_df: pd.DataFrame, pool_folders: List[str],
                                 output_folder: str,
                                 rename_map: Optional[Dict[str, str]] = None,
                                 file_maps: Optional[List[Dict[str, Dict[str, str]]]] = None) -> Dict[str, List[str]]:
    """
    Extract data from multiple pools based on source_table column.

    When tables are concatenated with add_source=True, each row has a source_table
    column (0, 1, 2...) indicating which table it came from. This function uses that
    to select files from the corresponding pool.

    Args:
        result_df: DataFrame with 'id' and 'source_table' columns
        pool_folders: List of pool folders, indexed by source_table
        output_folder: Where to copy the extracted data
        rename_map: Optional mapping from original ID to new ID
        file_maps: List of file maps, one per pool: [{stream_name: {id: file_path}}]

    Returns:
        Dictionary mapping stream name to list of extracted file paths
    """
    if not file_maps:
        raise ValueError("file_maps is required for multi-pool extraction")

    if 'source_table' not in result_df.columns:
        raise ValueError("source_table column required for multi-pool extraction")

    os.makedirs(output_folder, exist_ok=True)

    print(f"\nMulti-pool mode: Extracting data from {len(pool_folders)} pools")

    # Collect all stream names from all file maps
    all_stream_names = set()
    for fm in file_maps:
        all_stream_names.update(fm.keys())

    extracted_files = {name: [] for name in all_stream_names}

    for _, row in result_df.iterrows():
        # Use original_id for lookup if available (when rename was applied),
        # otherwise use id column
        if 'original_id' in result_df.columns:
            lookup_id = str(row.get('original_id', ''))
        else:
            lookup_id = str(row.get('id', ''))

        # Output ID is always the current 'id' column (may be renamed)
        output_id = str(row.get('id', lookup_id))

        source_idx = int(row.get('source_table', 0))

        if source_idx >= len(file_maps):
            raise ValueError(f"source_table {source_idx} >= num pools {len(file_maps)}")

        file_map = file_maps[source_idx]

        # Extract from all streams in this pool's file map
        for stream_name, id_to_file in file_map.items():
            if lookup_id not in id_to_file:
                continue

            source = id_to_file[lookup_id]
            ext = os.path.splitext(source)[1]
            dest = os.path.join(output_folder, f"{output_id}{ext}")

            shutil.copy2(source, dest)
            extracted_files[stream_name].append(dest)
            print(f"Extracted {stream_name} from pool {source_idx}: {lookup_id} -> {dest}")

    return extracted_files


def create_missing_csv(original_ids: List[str], filtered_ids: List[str],
                       output_folder: str, step_tool_name: str,
                       operations: List[Dict[str, Any]]) -> None:
    """
    Create missing.csv with IDs that were filtered out.

    Args:
        original_ids: All original IDs
        filtered_ids: IDs that passed filtering
        output_folder: Output folder for missing.csv
        step_tool_name: Step and tool name (e.g. "005_Panda")
        operations: List of operation dicts from config
    """
    all_ids = set(str(i) for i in original_ids)
    passed_ids = set(str(i) for i in filtered_ids)
    missing_ids = list(all_ids - passed_ids)

    # Build cause from operations
    filter_exprs = [op['params']['expr'] for op in operations
                    if op.get('type') == 'filter' and 'params' in op and 'expr' in op['params']]
    if filter_exprs:
        cause = "Filtered by: " + ", ".join(filter_exprs)
    else:
        op_summaries = []
        for op in operations:
            op_type = op.get('type', 'unknown')
            if op_type in ('head', 'tail', 'sample') and 'params' in op:
                n = op['params'].get('n', '?')
                op_summaries.append(f"{op_type}({n})")
            else:
                op_summaries.append(op_type)
        cause = "Removed by: " + ", ".join(op_summaries)

    if not missing_ids:
        missing_df = pd.DataFrame(columns=['id', 'removed_by', 'cause'])
    else:
        missing_data = [{'id': mid, 'removed_by': step_tool_name, 'cause': cause}
                        for mid in missing_ids]
        missing_df = pd.DataFrame(missing_data)

    missing_csv = os.path.join(output_folder, "missing.csv")
    missing_df.to_csv(missing_csv, index=False)
    print(f"Created missing.csv with {len(missing_ids)} filtered out IDs")


def filter_and_copy_pool_tables(
    result_df: pd.DataFrame,
    pool_table_maps: List[Dict[str, Dict[str, Any]]],
    pool_folders: List[str],
    output_folder: str,
    id_map_forward: Dict[str, List[str]],
    rename_map: Optional[Dict[str, str]] = None
) -> Dict[str, str]:
    """
    Filter and copy pool tables to output folder, matching rows by ID.

    For each table in the pool(s), this function:
    1. Loads the table
    2. Filters rows to match the IDs in result_df
    3. Applies ID remapping (using id_map_forward for reverse lookup)
    4. Saves the filtered table to output_folder

    Args:
        result_df: DataFrame with filtered IDs (has 'id' column, may have 'source_table')
        pool_table_maps: List of table maps per pool: [{table_name: {"path": str, "columns": list}}]
        pool_folders: List of pool folder paths
        output_folder: Where to save filtered tables
        id_map_forward: Mapping new_id -> [old_ids] for reverse lookup
        rename_map: Optional mapping from lookup_id to output_id

    Returns:
        Dictionary mapping table_name to output path
    """
    if not pool_table_maps:
        return {}

    os.makedirs(output_folder, exist_ok=True)

    # Collect all table names across all pools
    all_table_names = set()
    for tm in pool_table_maps:
        all_table_names.update(tm.keys())

    copied_tables = {}
    has_source_table = 'source_table' in result_df.columns
    has_original_id = 'original_id' in result_df.columns

    print(f"\nFiltering and copying pool tables...")
    print(f"Tables to process: {list(all_table_names)}")

    for table_name in all_table_names:
        # Collect rows from all pools that have this table
        filtered_rows = []

        for pool_idx, table_map in enumerate(pool_table_maps):
            if table_name not in table_map:
                continue

            table_info = table_map[table_name]
            table_path = table_info["path"]

            if not os.path.exists(table_path):
                print(f"  Warning: Table not found: {table_path}")
                continue

            # Load the pool table
            pool_df = pd.read_csv(table_path)

            if 'id' not in pool_df.columns:
                print(f"  Warning: Table {table_name} has no 'id' column, skipping")
                continue

            # For each row in result_df, find matching rows in pool table
            for _, result_row in result_df.iterrows():
                # Determine which pool this row came from
                if has_source_table:
                    row_pool_idx = int(result_row.get('source_table', 0))
                    if row_pool_idx != pool_idx:
                        continue  # This row is from a different pool
                elif len(pool_table_maps) > 1:
                    # Multiple pools but no source_table - can't determine which pool
                    # Try all pools (may get duplicates, but better than nothing)
                    pass

                # Get the ID to look up in the pool table
                if has_original_id:
                    lookup_id = str(result_row.get('original_id', ''))
                else:
                    lookup_id = str(result_row.get('id', ''))

                # Get the output ID (may be renamed)
                output_id = str(result_row.get('id', lookup_id))

                # Try to find matching row in pool table
                # Strategy 1: Direct match
                matching = pool_df[pool_df['id'].astype(str) == lookup_id]

                # Strategy 2: If no match and we have id_map_forward, try reverse lookup
                if matching.empty and id_map_forward:
                    # lookup_id might be the new (merged) ID, find original IDs
                    original_ids = id_map_forward.get(lookup_id, [])
                    for orig_id in original_ids:
                        matching = pool_df[pool_df['id'].astype(str) == orig_id]
                        if not matching.empty:
                            break

                if not matching.empty:
                    # Copy matching row(s) and update ID
                    for _, match_row in matching.iterrows():
                        new_row = match_row.copy()
                        new_row['id'] = output_id
                        filtered_rows.append(new_row)

        # Determine output filename from the pool table path (use basename to match get_output_files)
        output_filename = None
        for tm in pool_table_maps:
            if table_name in tm:
                output_filename = os.path.basename(tm[table_name]["path"])
                break
        if not output_filename:
            raise ValueError(f"Table '{table_name}' not found in any pool_table_maps")

        # Create output table
        if filtered_rows:
            output_df = pd.DataFrame(filtered_rows)
            # Remove duplicate rows (same id might be matched multiple times)
            output_df = output_df.drop_duplicates(subset=['id'], keep='first')
        else:
            # Create empty DataFrame with expected columns
            if pool_table_maps and table_name in pool_table_maps[0]:
                columns = pool_table_maps[0][table_name].get("columns", ["id"])
                if not columns:
                    columns = ["id"]
            else:
                columns = ["id"]
            output_df = pd.DataFrame(columns=columns)

        # Save to output folder (use output_filename from pool table path)
        output_path = os.path.join(output_folder, output_filename)
        output_df.to_csv(output_path, index=False)
        copied_tables[table_name] = output_path
        print(f"  {table_name}: {len(output_df)} rows -> {output_path}")

    return copied_tables


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
    pool_folders = config_data.get('pool_folders', [])
    pool_file_maps = config_data.get('pool_file_maps', [])
    pool_table_maps = config_data.get('pool_table_maps', [])
    id_map_forward = config_data.get('id_map_forward', {})
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
    if use_pool_mode and pool_file_maps:
        output_dir = os.path.dirname(output_csv)

        # Extract files only if we have IDs to extract
        extracted = {}
        if filtered_ids_for_lookup:
            # Use multi-pool extraction if we have source_table column and multiple pools
            if len(pool_file_maps) > 1 and 'source_table' in result_df.columns:
                extracted = extract_from_multiple_pools(
                    result_df, pool_folders, output_dir,
                    rename_map=original_to_new_id if rename else None,
                    file_maps=pool_file_maps
                )
            else:
                # Single pool mode - use first pool
                extracted = extract_pool_data_for_filtered_ids(
                    filtered_ids_for_lookup, pool_folders[0], output_dir,
                    rename_map=original_to_new_id if rename else None,
                    file_map=pool_file_maps[0] if pool_file_maps else None
                )

        print(f"\nPool mode summary:")
        print(f"Filtered IDs: {len(filtered_ids)}")
        print(f"Pool folders: {len(pool_folders)}")
        for stream_name, files in extracted.items():
            print(f"Extracted {stream_name}: {len(files)} files")

        # Filter and copy pool tables (always do this, even with 0 rows - creates empty tables)
        if pool_table_maps:
            copied_tables = filter_and_copy_pool_tables(
                result_df=result_df,
                pool_table_maps=pool_table_maps,
                pool_folders=pool_folders,
                output_folder=output_dir,
                id_map_forward=id_map_forward,
                rename_map=original_to_new_id if rename else None
            )
            print(f"Copied {len(copied_tables)} tables")

    # Create missing.csv
    if original_ids:
        output_dir = os.path.dirname(output_csv)
        step_tool_name = os.path.basename(output_dir)
        create_missing_csv(original_ids, filtered_ids, output_dir,
                           step_tool_name, operations)

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
