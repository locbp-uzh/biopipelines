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
import csv
import json
import shutil
import pandas as pd
import numpy as np
from typing import Dict, List, Any, Optional, Tuple

# Import unified I/O utilities for runtime DataStream expansion
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from biopipelines.biopipelines_io import load_datastream, iterate_files
from biopipelines.id_map_utils import get_mapped_ids, map_table_ids_to_ids, prune_redundant_provenance_columns

# Mirror of biopipelines.panda.Panda.SOURCE — kept inline to avoid a
# config-time import dependency in pipe scripts.
PANDA_SOURCE = "bp_panda_source"


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


def execute_zscore(df: pd.DataFrame, params: Dict[str, Any]) -> pd.DataFrame:
    """Standardize columns to z-scores ((x-mean)/std, population std)."""
    columns = params.get("columns", [])
    by = params.get("by")
    sign = params.get("sign", {}) or {}
    suffix = params.get("suffix", "_z")
    df = df.copy()

    for col in columns:
        if col not in df.columns:
            print(f"  Warning: zscore column '{col}' not found, skipping")
            continue
        s = pd.to_numeric(df[col], errors="coerce")
        if by and by in df.columns:
            # population std within each group (ddof=0)
            grp = s.groupby(df[by])
            z = (s - grp.transform("mean")) / grp.transform("std", ddof=0)
        else:
            std = s.std(ddof=0)
            z = (s - s.mean()) / std if std else s * 0.0
        if sign.get(col, 1) < 0:
            z = -z
        df[f"{col}{suffix}"] = z
        print(f"  Zscore: {col}{suffix}" + (f" (by {by})" if by else "") +
              (" [sign-flipped]" if sign.get(col, 1) < 0 else ""))

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
    """Merge multiple dataframes horizontally.

    Two code paths based on ``on``:
    - ``on`` is not None: exact pandas merge (original behaviour)
    - ``on`` is None: biopipelines ID matching via get_mapped_ids
    """
    on = params.get("on")
    how = params.get("how")
    prefixes = params.get("prefixes", [])
    map_table_paths = params.get("map_table_paths", [])
    grain = params.get("grain", "finest")

    # --- Path A: explicit column merge (on is not None) ---
    if on is not None:
        if how is None:
            how = "outer"
        # Normalise `on` to a per-table list of column names.
        if isinstance(on, list):
            on_per_table = on
            merge_key = on_per_table[0]
        else:
            on_per_table = [on] * len(dataframes)
            merge_key = on

        print(f"  Merge (exact): on={on}, how={how}, merge_key={merge_key}")

        renamed_dfs = []
        for i, df in enumerate(dataframes):
            df = df.copy()

            # Rename per-table key to the canonical merge key
            table_key = on_per_table[i] if i < len(on_per_table) else merge_key
            if table_key != merge_key:
                if table_key in df.columns:
                    if merge_key in df.columns:
                        df = df.drop(columns=[merge_key])
                        print(f"    Dropped original '{merge_key}' column in dataframe {i} (replaced by '{table_key}')")
                    df = df.rename(columns={table_key: merge_key})
                    print(f"    Renamed '{table_key}' -> '{merge_key}' in dataframe {i}")

            # Apply prefix
            if prefixes and i < len(prefixes) and prefixes[i]:
                prefix = prefixes[i]
                rename_dict = {col: f"{prefix}{col}" for col in df.columns if col != merge_key}
                df = df.rename(columns=rename_dict)
                print(f"    Applied prefix '{prefix}' to dataframe {i}")

            renamed_dfs.append(df)

        result = renamed_dfs[0]
        for i, df in enumerate(renamed_dfs[1:], 1):
            overlap = set(result.columns) & set(df.columns) - {merge_key}
            if overlap:
                print(f"    Warning: Overlapping columns: {overlap}")
                df = df.drop(columns=list(overlap))

            before_rows = len(result)
            result = pd.merge(result, df, on=merge_key, how=how)
            print(f"    Merged dataframe {i+1}: {before_rows} -> {len(result)} rows")

        return result

    # --- Path B: ID matching merge (on is None) ---
    merge_key = "id"
    # Derive the join type from grain when not set: both 'finest' and 'coarsest'
    # use a left join so the chosen driver keeps EXACTLY its rows. An outer join
    # would resurrect driver-unmatched rows from the other table with a null id
    # (a coarse design with no fine child, or a dropped fine sibling). Use an
    # explicit how='outer' if you really want those rows.
    if how is None:
        how = "left" if grain in ("finest", "coarsest") else "outer"
    print(f"  Merge (ID matching): how={how}, grain={grain}, merge_key={merge_key}")
    if map_table_paths:
        print(f"    Using {len(map_table_paths)} map_table(s) for provenance")

    # Verify all dataframes have the merge key
    for i, df in enumerate(dataframes):
        if merge_key not in df.columns:
            raise ValueError(f"Dataframe {i} has no '{merge_key}' column. "
                             f"Use on=<column> for explicit column merge.")

    # Pick which table drives the output row set by grain. Grain is a table's
    # depth in the provenance hierarchy (d_1_1 derives from d_1 derives from d),
    # NOT its row count — a coarse table can have more rows than a fine one
    # (many designs, few refolded structures). Rank by the max ID depth, with
    # unique-ID count as a tiebreak and original order as the final tiebreak, so
    # "finest" always lands on the deepest table. Reorder dataframes + their
    # prefixes together so the prefix-to-table correspondence survives.
    order = list(range(len(dataframes)))
    if grain in ("finest", "coarsest"):
        def _depth(df):
            return max((len(map_table_ids_to_ids(i, {"*": "*_<S>"}))
                        for i in df[merge_key].astype(str)), default=0)
        depths = [_depth(df) for df in dataframes]
        counts = [df[merge_key].astype(str).nunique() for df in dataframes]
        sign = 1 if grain == "finest" else -1
        driver = max(range(len(dataframes)),
                     key=lambda i: (sign * depths[i], sign * counts[i], -i))
        order = [driver] + [i for i in range(len(dataframes)) if i != driver]
        if driver != 0:
            print(f"    Grain '{grain}': dataframe {driver} drives rows "
                  f"(depth {depths[driver]}, {counts[driver]} IDs)")
    dataframes = [dataframes[i] for i in order]
    prefixes = ([prefixes[i] for i in order]
                if prefixes and len(prefixes) == len(order) else prefixes)

    # Start with first (driver) dataframe
    result = dataframes[0].copy()
    if prefixes and len(prefixes) > 0 and prefixes[0]:
        prefix = prefixes[0]
        rename_dict = {col: f"{prefix}{col}" for col in result.columns if col != merge_key}
        result = result.rename(columns=rename_dict)
        print(f"    Applied prefix '{prefix}' to dataframe 0")

    # Merge each subsequent dataframe using ID matching
    for i, df_i in enumerate(dataframes[1:], 1):
        df_i = df_i.copy()

        # Get current IDs from the running result and the incoming dataframe
        result_ids = result[merge_key].astype(str).tolist()
        df_i_ids = df_i[merge_key].astype(str).tolist()

        # Use get_mapped_ids to find matches: result_id -> df_i_id
        matched = get_mapped_ids(result_ids, df_i_ids, unique=True,
                                 map_table_paths=map_table_paths)

        # Warn when this fold drops sibling rows: one result row matching many
        # incoming rows means the unique match silently discards the rest. This
        # only happens when the incoming table is finer than the driver (e.g.
        # grain='coarsest'); grain='finest' picks the finest driver so every
        # fold is a coarse->fine broadcast with nothing dropped.
        multi = get_mapped_ids(result_ids, df_i_ids, unique=False,
                               map_table_paths=map_table_paths)
        dropped = sum(max(0, len(v) - 1) for v in multi.values() if v)
        if dropped:
            print(f"    WARNING: dataframe {i} sits at a finer grain than the "
                  f"driver; {dropped} row(s) dropped to keep one match per "
                  f"driver row. Set grain='finest' or aggregate this table "
                  f"first to keep them.")

        n_matched = sum(1 for v in matched.values() if v is not None)
        print(f"    ID matching df {i}: {n_matched}/{len(result_ids)} result IDs matched")

        # Attach each driver row's matched incoming-id as a temp join key, then
        # merge on it. This broadcasts: many driver rows sharing one coarse id
        # all pick up that coarse row's columns (pure rewrite-then-merge would
        # collapse them, since the id->id map is many-to-one).
        join_key = "__merge_join_key__"
        result[join_key] = result[merge_key].astype(str).map(
            lambda x: matched.get(x))
        df_i[join_key] = df_i[merge_key].astype(str)

        # Apply prefix (never to the merge key or the temp join key)
        if prefixes and i < len(prefixes) and prefixes[i]:
            prefix = prefixes[i]
            rename_dict = {col: f"{prefix}{col}" for col in df_i.columns
                           if col not in (merge_key, join_key)}
            df_i = df_i.rename(columns=rename_dict)
            print(f"    Applied prefix '{prefix}' to dataframe {i}")

        # Drop the incoming merge key so the driver's id survives the broadcast;
        # then drop any other overlapping columns.
        df_i = df_i.drop(columns=[merge_key])
        overlap = (set(result.columns) & set(df_i.columns)) - {join_key}
        if overlap:
            print(f"    Warning: Overlapping columns: {overlap}")
            df_i = df_i.drop(columns=list(overlap))

        before_rows = len(result)
        result = pd.merge(result, df_i, on=join_key, how=how)
        result = result.drop(columns=[join_key])
        print(f"    Merged dataframe {i+1}: {before_rows} -> {len(result)} rows")

    return result


def execute_concat(dataframes: List[pd.DataFrame], params: Dict[str, Any]) -> pd.DataFrame:
    """Concatenate multiple dataframes vertically.

    Source tagging is handled implicitly upstream of the operation chain
    (see ``main()`` in this file): when concat appears in the chain over
    multiple inputs, every input dataframe is pre-tagged with
    ``PANDA_SOURCE = i`` before the chain runs, and that column is stripped
    from the final result before writing. So this function just stacks
    rows — the source column is already present.
    """
    fill = params.get("fill")

    print(f"  Concat: fill={fill}")

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

    This replaces AverageByTable functionality. After concat over multiple
    inputs, every row carries an implicit PANDA_SOURCE column; this groups
    by it and computes the mean of all numeric columns.
    """
    source_col = params.get("source_col") or PANDA_SOURCE

    if source_col not in df.columns:
        print(f"  Warning: source column '{source_col}' not found. "
              f"Run concat over multiple inputs first.")
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
        "zscore": execute_zscore,
        "fillna": execute_fillna,
        "groupby": execute_groupby,
        "pivot": execute_pivot,
        "melt": execute_melt,
        "average_by_source": execute_average_by_source,
    }

    # Broadcast per-frame ops over a list of dataframes. Enables the
    # "best per source" idiom (sort -> head -> concat) where the user
    # wants the op applied to each input table independently before
    # they're stacked. Cross-table ops that combine rows across sources
    # (groupby, pivot, melt, average_by_source) require concat first —
    # broadcasting them per-frame would produce surprising results.
    cross_table_ops = {"groupby", "pivot", "melt", "average_by_source"}
    if isinstance(df_or_dfs, list):
        if len(df_or_dfs) == 1:
            df = df_or_dfs[0]
        elif op_type in cross_table_ops:
            raise ValueError(
                f"Operation '{op_type}' aggregates across sources and "
                f"requires a single dataframe. Use 'concat' first."
            )
        else:
            func = operation_funcs.get(op_type)
            if not func:
                print(f"  Warning: Unknown operation '{op_type}', skipping")
                return df_or_dfs
            print(f"  Broadcasting '{op_type}' over {len(df_or_dfs)} dataframes")
            return [func(d, params) for d in df_or_dfs]
    else:
        df = df_or_dfs

    func = operation_funcs.get(op_type)
    if func:
        return func(df, params)
    else:
        print(f"  Warning: Unknown operation '{op_type}', skipping")
        return df


def build_file_map_from_stream_jsons(stream_jsons: Dict[str, str]) -> Dict[str, Dict[str, str]]:
    """Build {stream_name: {id: file_path}} from saved DataStream JSONs.

    Per-id and templated streams populate the inner {id: file_path} dict.
    Shared-file streams (one artifact for the whole stream) are excluded
    here — call ``build_shared_streams_from_stream_jsons`` to discover them.

    Args:
        stream_jsons: {stream_name: json_path} from config

    Returns:
        file_map dict ready for extract_pool_data_for_filtered_ids
    """
    file_map = {}
    for stream_name, json_path in stream_jsons.items():
        ds = load_datastream(json_path)
        if ds.is_shared_file:
            continue
        file_map[stream_name] = {}
        for sid, fpath in iterate_files(ds):
            file_map[stream_name][sid] = fpath
    return file_map


def build_shared_streams_from_stream_jsons(stream_jsons: Dict[str, str]) -> Dict[str, Dict[str, str]]:
    """Build {stream_name: {"src": path, "format": fmt}} for shared-file streams.

    Shared-file streams own a single artifact (e.g. a multi-record FASTA)
    that must be sliced — not copied per id — when only a subset of ids
    survives filtering. This function discovers those streams from the
    saved DataStream JSONs so the extractor can dispatch the format-aware
    slicer once per stream.

    Args:
        stream_jsons: {stream_name: json_path} from config

    Returns:
        Dict mapping stream_name to {"src": absolute source path, "format": stream format}.
    """
    shared = {}
    for stream_name, json_path in stream_jsons.items():
        ds = load_datastream(json_path)
        if ds.is_shared_file:
            shared[stream_name] = {"src": ds.files, "format": ds.format}
    return shared


def extract_pool_data_for_filtered_ids(filtered_ids: List[str], pool_folder: str,
                                        output_folder: str,
                                        rename_map: Optional[Dict[str, str]] = None,
                                        file_map: Optional[Dict[str, Dict[str, str]]] = None,
                                        ignore_missing: bool = False,
                                        stream_folders: Optional[Dict[str, str]] = None,
                                        shared_map: Optional[Dict[str, Dict[str, str]]] = None) -> Dict[str, List[str]]:
    """
    Extract data from a single pool for filtered IDs.

    Args:
        filtered_ids: List of IDs that passed filtering (original IDs for lookup)
        pool_folder: Pool output folder containing all data types
        output_folder: Where to copy the filtered data
        rename_map: Optional mapping from original ID to new ID for renaming output files
        file_map: Optional mapping {stream_name: {id: file_path}} for direct file lookup
        ignore_missing: If True, skip missing files with a warning instead of failing
        shared_map: Optional mapping {stream_name: {"src": path, "format": fmt}}
            for shared-file streams. Each gets sliced once via the
            format-specific function in biopipelines.stream_slicers, then
            written to its stream folder. Filtered IDs after rename_map
            are passed to the slicer (no rename support for shared files —
            the shared artifact's record IDs are the original ones).

    Returns:
        Dictionary mapping stream name to list of extracted file paths
    """
    extracted_files = {}

    os.makedirs(output_folder, exist_ok=True)

    print(f"\nPool mode: Extracting data for {len(filtered_ids)} filtered IDs from {pool_folder}")
    if rename_map:
        print(f"Renaming files according to rename map")

    # Shared-file streams: slice once per stream (not per id). If rename_map
    # is set, the slicer rewrites record IDs so the FASTA/CSV agrees with
    # the downstream stream's ids.
    if shared_map:
        from biopipelines.stream_slicers import get_slicer
        # rename_map keys are *original* IDs (lookup); values are new IDs.
        # The slicer only renames IDs in `filtered_ids` (which are the
        # original IDs at this layer), so pass it through directly.
        slicer_rename = rename_map or None
        for stream_name, info in shared_map.items():
            src = info["src"]
            fmt = info["format"]
            if not os.path.exists(src):
                if ignore_missing:
                    print(f"  Warning: Missing shared file for {stream_name}: {src} (skipped)")
                    extracted_files[stream_name] = []
                    continue
                raise FileNotFoundError(
                    f"Shared file not found for {stream_name}: {src}. "
                    f"Use ignore_missing=True to skip."
                )
            dest_folder = (stream_folders or {}).get(stream_name, output_folder)
            os.makedirs(dest_folder, exist_ok=True)
            dest = os.path.join(dest_folder, os.path.basename(src))
            slicer = get_slicer(fmt)
            slicer(src, dest, filtered_ids, rename_map=slicer_rename)
            extracted_files[stream_name] = [dest]
            print(f"Sliced {stream_name} ({fmt}): {src} -> {dest} (kept {len(filtered_ids)} ids)")

    if file_map:
        # Use file map - iterate over all streams dynamically
        print(f"Using file map for direct lookup (streams: {list(file_map.keys())})")

        for stream_name, id_to_file in file_map.items():
            extracted_files[stream_name] = []

            for selected_id in filtered_ids:
                if selected_id not in id_to_file:
                    # Fallback: use ID matching to find the file
                    match = get_mapped_ids([selected_id], list(id_to_file.keys()), unique=True)
                    matched_id = match.get(selected_id)
                    if not matched_id or matched_id not in id_to_file:
                        continue
                    selected_id_lookup = matched_id
                else:
                    selected_id_lookup = selected_id

                source = id_to_file[selected_id_lookup]

                if not os.path.exists(source):
                    if ignore_missing:
                        print(f"  Warning: Missing file for {stream_name} ID '{selected_id}': {source} (skipped)")
                        continue
                    else:
                        raise FileNotFoundError(
                            f"Pool file not found for {stream_name} ID '{selected_id}': {source}. "
                            f"Use ignore_missing=True to skip missing files."
                        )

                output_id = rename_map.get(str(selected_id), selected_id) if rename_map else selected_id
                ext = os.path.splitext(source)[1]
                dest_folder = (stream_folders or {}).get(stream_name, output_folder)
                dest = os.path.join(dest_folder, f"{output_id}{ext}")

                shutil.copy2(source, dest)
                extracted_files[stream_name].append(dest)
                print(f"Extracted {stream_name}: {source} -> {dest}")
    elif not shared_map:
        raise ValueError("file_map or shared_map is required for pool extraction")

    return extracted_files


def extract_from_multiple_pools(result_df: pd.DataFrame, pool_folders: List[str],
                                 output_folder: str,
                                 rename_map: Optional[Dict[str, str]] = None,
                                 file_maps: Optional[List[Dict[str, Dict[str, str]]]] = None,
                                 ignore_missing: bool = False,
                                 stream_folders: Optional[Dict[str, str]] = None,
                                 shared_maps: Optional[List[Dict[str, Dict[str, str]]]] = None) -> Dict[str, List[str]]:
    """
    Extract data from multiple pools based on the implicit PANDA_SOURCE column.

    When concat runs over multiple inputs, every row carries a PANDA_SOURCE
    column (0, 1, 2...) indicating which table it came from. This function
    uses that to select files from the corresponding pool.

    Args:
        result_df: DataFrame with 'id' and PANDA_SOURCE columns
        pool_folders: List of pool folders, indexed by PANDA_SOURCE
        output_folder: Where to copy the extracted data
        rename_map: Optional mapping from original ID to new ID
        file_maps: List of file maps, one per pool: [{stream_name: {id: file_path}}]
        ignore_missing: If True, skip missing files with a warning instead of failing
        shared_maps: List of shared-file maps, one per pool:
            [{stream_name: {"src": path, "format": fmt}}]. Each shared stream
            is sliced using the IDs present in result_df for the matching
            PANDA_SOURCE.

    Returns:
        Dictionary mapping stream name to list of extracted file paths
    """
    if not file_maps and not shared_maps:
        raise ValueError("file_maps or shared_maps is required for multi-pool extraction")

    if PANDA_SOURCE not in result_df.columns:
        raise ValueError(f"{PANDA_SOURCE} column required for multi-pool extraction")

    os.makedirs(output_folder, exist_ok=True)

    print(f"\nMulti-pool mode: Extracting data from {len(pool_folders)} pools")

    # Collect all stream names from all file maps and shared maps
    all_stream_names = set()
    for fm in (file_maps or []):
        all_stream_names.update(fm.keys())
    for sm in (shared_maps or []):
        all_stream_names.update(sm.keys())

    extracted_files = {name: [] for name in all_stream_names}

    # Shared-file streams across multiple pools: slice each pool into a temp
    # file, then merge them into a single canonical artifact under the
    # stream folder. This keeps the "shared-file means one artifact" rule
    # intact across PANDA_SOURCE values — without merging, pools with the
    # same basename overwrote each other.
    if shared_maps:
        import tempfile
        from biopipelines.stream_slicers import get_slicer, get_merger
        # Group (id, rename_target) by source pool, preserving order.
        ids_by_source: Dict[int, List[str]] = {}
        rename_by_source: Dict[int, Dict[str, str]] = {}
        for _, row in result_df.iterrows():
            lookup_id = str(row.get('original_id', row.get('id', '')))
            out_id = str(row.get('id', lookup_id))
            src_idx = int(row.get(PANDA_SOURCE, 0))
            ids_by_source.setdefault(src_idx, []).append(lookup_id)
            if rename_map and lookup_id != out_id:
                rename_by_source.setdefault(src_idx, {})[lookup_id] = out_id

        # Collect per-stream merge inputs across pools.
        stream_parts: Dict[str, List[str]] = {}
        stream_fmts: Dict[str, str] = {}
        stream_dests: Dict[str, str] = {}
        tmp_root = tempfile.mkdtemp(prefix="panda_shared_")
        for src_idx, sm in enumerate(shared_maps):
            kept = ids_by_source.get(src_idx, [])
            if not kept:
                continue
            rmap = rename_by_source.get(src_idx) or None
            for stream_name, info in sm.items():
                src = info["src"]
                fmt = info["format"]
                if not os.path.exists(src):
                    if ignore_missing:
                        print(f"  Warning: Missing shared file for {stream_name} (pool {src_idx}): {src} (skipped)")
                        continue
                    raise FileNotFoundError(
                        f"Shared file not found for {stream_name} (pool {src_idx}): {src}. "
                        f"Use ignore_missing=True to skip."
                    )
                slicer = get_slicer(fmt)
                part = os.path.join(tmp_root, f"pool{src_idx}_{stream_name}_{os.path.basename(src)}")
                slicer(src, part, kept, rename_map=rmap)
                stream_parts.setdefault(stream_name, []).append(part)
                stream_fmts[stream_name] = fmt
                if stream_name not in stream_dests:
                    dest_folder = (stream_folders or {}).get(stream_name, output_folder)
                    os.makedirs(dest_folder, exist_ok=True)
                    stream_dests[stream_name] = os.path.join(dest_folder, os.path.basename(src))
                print(f"Sliced {stream_name} ({fmt}) from pool {src_idx}: {src} -> {part} (kept {len(kept)} ids)")

        for stream_name, parts in stream_parts.items():
            fmt = stream_fmts[stream_name]
            dest = stream_dests[stream_name]
            merger = get_merger(fmt)
            merger(parts, dest)
            extracted_files[stream_name] = [dest]
            print(f"Merged {stream_name} ({fmt}): {len(parts)} parts -> {dest}")

    if not file_maps:
        return extracted_files

    for _, row in result_df.iterrows():
        # Use original_id for lookup if available (when rename was applied),
        # otherwise use id column
        if 'original_id' in result_df.columns:
            lookup_id = str(row.get('original_id', ''))
        else:
            lookup_id = str(row.get('id', ''))

        # Output ID is always the current 'id' column (may be renamed)
        output_id = str(row.get('id', lookup_id))

        source_idx = int(row.get(PANDA_SOURCE, 0))

        if source_idx >= len(file_maps):
            raise ValueError(f"{PANDA_SOURCE}={source_idx} >= num pools {len(file_maps)}")

        file_map = file_maps[source_idx]

        # Extract from all streams in this pool's file map
        for stream_name, id_to_file in file_map.items():
            if lookup_id not in id_to_file:
                # Fallback: use ID matching to find the file
                match = get_mapped_ids([lookup_id], list(id_to_file.keys()), unique=True)
                matched_id = match.get(lookup_id)
                if not matched_id or matched_id not in id_to_file:
                    continue
                actual_lookup_id = matched_id
            else:
                actual_lookup_id = lookup_id

            source = id_to_file[actual_lookup_id]

            if not os.path.exists(source):
                if ignore_missing:
                    print(f"  Warning: Missing file for {stream_name} ID '{lookup_id}' from pool {source_idx}: {source} (skipped)")
                    continue
                else:
                    raise FileNotFoundError(
                        f"Pool file not found for {stream_name} ID '{lookup_id}' from pool {source_idx}: {source}. "
                        f"Use ignore_missing=True to skip missing files."
                    )

            ext = os.path.splitext(source)[1]
            dest_folder = (stream_folders or {}).get(stream_name, output_folder)
            dest = os.path.join(dest_folder, f"{output_id}{ext}")

            shutil.copy2(source, dest)
            extracted_files[stream_name].append(dest)
            print(f"Extracted {stream_name} from pool {source_idx}: {lookup_id} -> {dest}")

    return extracted_files


def shift_provenance_columns(df: pd.DataFrame, stream_name: str) -> pd.DataFrame:
    """Shift `<stream>.-N.id` provenance columns one generation further back.

    Chained Panda steps each contribute one generation of provenance under
    the same stream alias. Without shifting, every step would overwrite the
    previous `<stream>.id`, erasing the chain ENT -> Panda_1 -> Panda_2.

    The generation tag uses a negative integer (`.-1`, `.-2`, ...) so the
    column name reads as "one step back, two steps back". Negative integers
    cannot collide with Bundle's `<stream>.0.id` / `<stream>.1.id` / ...
    (see combinatorics.py), which use non-negative integers to label
    bundled members of an axis.

    Renaming, from oldest generation outward:

        <stream>.-N.id  -> <stream>.-(N+1).id
        <stream>.-2.id  -> <stream>.-3.id
        <stream>.-1.id  -> <stream>.-2.id
        <stream>.id     -> <stream>.-1.id

    Returns a dataframe with the columns renamed. `<stream>.id` no longer
    appears in the result (it was renamed to `<stream>.-1.id`); the caller
    writes a fresh `<stream>.id` afterwards.
    """
    base = f"{stream_name}.id"
    prefix = f"{stream_name}.-"
    suffix = ".id"

    numbered: List[Tuple[int, str]] = []
    for col in df.columns:
        if not (col.startswith(prefix) and col.endswith(suffix)):
            continue
        middle = col[len(prefix):-len(suffix)]
        if middle.isdigit():
            numbered.append((int(middle), col))

    # Rename from oldest generation downward to avoid collisions.
    numbered.sort(reverse=True)
    rename_map: Dict[str, str] = {}
    for gen, col in numbered:
        rename_map[col] = f"{prefix}{gen + 1}{suffix}"
    if base in df.columns:
        rename_map[base] = f"{prefix}1{suffix}"
    if rename_map:
        df = df.rename(columns=rename_map)
    return df


def order_stream_provenance_columns(df: pd.DataFrame, stream_name: str) -> List[str]:
    """Return same-stream provenance columns in nearest-to-oldest order."""
    axis_col = f"{stream_name}.id"
    gen_cols = [axis_col] if axis_col in df.columns else []
    prev_prefix = f"{stream_name}.-"
    prev_numbered = []
    for col in df.columns:
        if col.startswith(prev_prefix) and col.endswith(".id"):
            mid = col[len(prev_prefix):-len(".id")]
            if mid.isdigit():
                prev_numbered.append((int(mid), col))
    prev_numbered.sort()
    gen_cols.extend(col for _, col in prev_numbered)
    return gen_cols


def load_upstream_provenance_row(
    pool_stream_jsons: List[Dict[str, str]],
    stream_name: str,
    source_idx: int,
    lookup_id: str,
) -> Dict[str, Any]:
    """Load provenance columns for one upstream stream row, if available."""
    if source_idx >= len(pool_stream_jsons):
        return {}
    stream_json = pool_stream_jsons[source_idx].get(stream_name)
    if not stream_json:
        return {}

    try:
        ds = load_datastream(stream_json)
    except Exception:
        return {}

    if not ds.map_table or not os.path.exists(ds.map_table):
        return {}

    try:
        map_df = pd.read_csv(
            ds.map_table,
            keep_default_na=False,
            na_values=['NA', 'N/A', '#N/A'],
        )
    except Exception:
        return {}

    if "id" not in map_df.columns:
        return {}

    ids = map_df["id"].astype(str).tolist()
    matching = map_df[map_df["id"].astype(str) == lookup_id]
    if matching.empty:
        matched = get_mapped_ids(
            [lookup_id],
            ids,
            unique=True,
            map_table_paths=[ds.map_table],
        ).get(lookup_id)
        if matched:
            matching = map_df[map_df["id"].astype(str) == matched]
    if matching.empty:
        return {}

    row = matching.iloc[0].to_dict()
    return {
        col: row[col]
        for col in map_df.columns
        if col != "id" and col.endswith(".id")
    }


def create_missing_csv(original_ids: List[str], filtered_ids: List[str],
                       output_folder: str, step_tool_name: str,
                       operations: List[Dict[str, Any]],
                       rename_map: Optional[Dict[str, str]] = None,
                       removed_by_op: Optional[Dict[str, str]] = None,
                       missing_csv_path: Optional[str] = None) -> None:
    """
    Create missing.csv with IDs that were filtered out or renamed.

    Args:
        original_ids: All original IDs
        filtered_ids: IDs that passed filtering (original IDs, before rename)
        output_folder: Output folder for missing.csv
        step_tool_name: Step and tool name (e.g. "005_Panda")
        operations: List of operation dicts from config
        rename_map: Optional mapping from original ID to new renamed ID
        removed_by_op: Optional per-ID cause mapping from execution tracking
    """
    all_ids = set(str(i) for i in original_ids)
    passed_ids = set(str(i) for i in filtered_ids)
    truly_missing_ids = list(all_ids - passed_ids)

    missing_data = []
    for mid in truly_missing_ids:
        cause = removed_by_op.get(mid, "Unknown") if removed_by_op else "Unknown"
        missing_data.append({'id': mid, 'removed_by': step_tool_name, 'kind': 'filter', 'cause': cause})

    # Add renamed IDs (survived filtering but got a new name)
    if rename_map:
        for orig_id in filtered_ids:
            orig_str = str(orig_id)
            if orig_str in rename_map:
                missing_data.append({
                    'id': orig_str,
                    'removed_by': step_tool_name,
                    'kind': 'filter',
                    'cause': f"Renamed to {rename_map[orig_str]}"
                })

    if missing_data:
        missing_df = pd.DataFrame(missing_data)
    else:
        missing_df = pd.DataFrame(columns=['id', 'removed_by', 'kind', 'cause'])

    missing_csv = missing_csv_path or os.path.join(output_folder, "missing.csv")
    missing_df.to_csv(missing_csv, index=False)
    print(f"Created missing.csv with {len(missing_data)} entries")


def merge_upstream_missing(output_folder: str, upstream_missing_paths: List[str],
                           missing_csv_path: Optional[str] = None) -> None:
    """Merge upstream missing tables into Panda's own missing.csv."""
    missing_csv = missing_csv_path or os.path.join(output_folder, "missing.csv")

    # Load Panda's own missing entries
    if os.path.exists(missing_csv):
        panda_missing = pd.read_csv(missing_csv)
    else:
        panda_missing = pd.DataFrame(columns=['id', 'removed_by', 'kind', 'cause'])

    # Load and merge upstream missing tables
    upstream_dfs = []
    for path in upstream_missing_paths:
        if os.path.exists(path):
            try:
                df = pd.read_csv(path)
                if not df.empty:
                    upstream_dfs.append(df)
                    print(f"  Merging {len(df)} upstream missing entries from {path}")
            except Exception as e:
                print(f"  Warning: Could not read upstream missing {path}: {e}")

    if not upstream_dfs:
        return

    # Concatenate: upstream first, then Panda's own (so Panda's entries win on dedup)
    all_dfs = upstream_dfs + [panda_missing]
    merged = pd.concat(all_dfs, ignore_index=True)
    merged = merged.drop_duplicates(subset=['id'], keep='last')

    merged.to_csv(missing_csv, index=False)
    print(f"  Merged missing.csv: {len(merged)} total entries")


def filter_and_copy_pool_tables(
    result_df: pd.DataFrame,
    pool_table_maps: List[Dict[str, Dict[str, Any]]],
    pool_folders: List[str],
    output_folder: str,
    rename_map: Optional[Dict[str, str]] = None,
    pool_table_targets: Optional[Dict[str, str]] = None,
) -> Dict[str, str]:
    """
    Filter and copy pool tables to output folder, matching rows by ID.

    For each table in the pool(s), this function:
    1. Loads the table
    2. Filters rows to match the IDs in result_df (using get_mapped_ids fallback)
    3. Saves the filtered table to output_folder

    Args:
        result_df: DataFrame with filtered IDs (has 'id' column, may have PANDA_SOURCE)
        pool_table_maps: List of table maps per pool: [{table_name: {"path": str, "columns": list}}]
        pool_folders: List of pool folder paths
        output_folder: Where to save filtered tables
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
    has_source_table = PANDA_SOURCE in result_df.columns
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
            pool_df = pd.read_csv(table_path, keep_default_na=False, na_values=['NA', 'N/A', '#N/A'])

            if 'id' not in pool_df.columns:
                print(f"  Warning: Table {table_name} has no 'id' column, skipping")
                continue

            # For each row in result_df, find matching rows in pool table
            for _, result_row in result_df.iterrows():
                # Determine which pool this row came from
                if has_source_table:
                    row_pool_idx = int(result_row.get(PANDA_SOURCE, 0))
                    if row_pool_idx != pool_idx:
                        continue  # This row is from a different pool
                elif len(pool_table_maps) > 1:
                    # Multiple pools but no source column - can't determine which pool
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

                # Strategy 2: If no direct match, use ID matching fallback
                if matching.empty:
                    pool_table_ids = pool_df['id'].astype(str).tolist()
                    match = get_mapped_ids([lookup_id], pool_table_ids, unique=True)
                    matched_id = match.get(lookup_id)
                    if matched_id:
                        matching = pool_df[pool_df['id'].astype(str) == matched_id]

                if not matching.empty:
                    # Copy matching row(s) and update ID, adding provenance
                    for _, match_row in matching.iterrows():
                        new_row = match_row.copy()
                        original_pool_id = str(match_row['id'])
                        new_row['id'] = output_id
                        # Add provenance column if ID was remapped
                        if output_id != original_pool_id:
                            new_row['pool.id'] = original_pool_id
                        filtered_rows.append(new_row)

        # Use the upstream table's basename as output filename
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

        if pool_table_targets and table_name in pool_table_targets:
            output_path = pool_table_targets[table_name]
        else:
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
    pool_stream_jsons = config_data.get('pool_stream_jsons', [])
    # Build file maps at runtime by expanding DataStream patterns
    pool_file_maps = []
    pool_shared_maps = []
    for stream_jsons in pool_stream_jsons:
        pool_file_maps.append(build_file_map_from_stream_jsons(stream_jsons))
        pool_shared_maps.append(build_shared_streams_from_stream_jsons(stream_jsons))
    pool_table_maps = config_data.get('pool_table_maps', [])
    map_table_paths = config_data.get('map_table_paths', [])
    rename = config_data.get('rename')
    ignore_missing = config_data.get('ignore_missing', False)

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
        df = pd.read_csv(csv_path, keep_default_na=False, na_values=['NA', 'N/A', '#N/A'])
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

    # Implicit source tagging: if the chain contains a concat over multiple
    # inputs, pre-tag every input dataframe with PANDA_SOURCE so downstream
    # ops in the chain (sort, head, groupby, average_by_source, multi-pool
    # extract) can reference it. The column is stripped from the final
    # result before writing — see the strip step after the loop.
    has_concat = any(op.get('type') == 'concat' for op in operations)
    source_tagged = has_concat and is_multi_table
    if source_tagged:
        for i in range(len(current)):
            current[i] = current[i].copy()
            current[i][PANDA_SOURCE] = i

    # Inject map_table_paths into merge operations for ID matching
    for operation in operations:
        if operation.get('type') == 'merge':
            operation.setdefault('params', {})['map_table_paths'] = map_table_paths

    # Execute operations sequentially, tracking which IDs are removed by each
    removed_by_op = {}
    print("\nExecuting operations:")
    for i, operation in enumerate(operations):
        print(f"\n[{i+1}/{len(operations)}] {operation['type']}")

        # Capture IDs before operation. Works for both a single dataframe
        # and a list of dataframes (per-frame broadcast phase) by taking
        # the union of every frame's IDs.
        def _all_ids(d_or_list):
            if isinstance(d_or_list, pd.DataFrame):
                return (set(d_or_list['id'].astype(str).tolist())
                        if 'id' in d_or_list.columns else None)
            if isinstance(d_or_list, list):
                ids = set()
                for d in d_or_list:
                    if isinstance(d, pd.DataFrame) and 'id' in d.columns:
                        ids.update(d['id'].astype(str).tolist())
                return ids if ids else None
            return None

        ids_before = _all_ids(current)

        current = execute_operation(current, operation, is_multi_table)

        ids_after = _all_ids(current)
        if ids_before is not None and ids_after is not None:
            dropped = ids_before - ids_after
            if dropped:
                op_type = operation['type']
                if op_type == 'filter':
                    cause = f"Filtered by: {operation.get('params', {}).get('expr', '?')}"
                elif op_type in ('head', 'tail', 'sample'):
                    n = operation.get('params', {}).get('n', '?')
                    cause = f"Removed by: {op_type}({n})"
                else:
                    cause = f"Removed by: {op_type}"
                for did in dropped:
                    removed_by_op[did] = cause

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

    # Strip the implicit PANDA_SOURCE column from the written CSV — it's a
    # framework-internal column used to drive multi-pool routing and
    # source-grouped operations during the chain. Keep it in result_df so
    # downstream pool extraction can still use it.
    result_to_write = (result_df.drop(columns=[PANDA_SOURCE])
                       if source_tagged and PANDA_SOURCE in result_df.columns
                       else result_df)
    result_to_write.to_csv(output_csv, index=False)

    print(f"\nFinal result: {result_to_write.shape}")
    print(f"Columns: {list(result_to_write.columns)}")
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
    if use_pool_mode:
        output_dir = os.path.dirname(output_csv)

        # Extract files only if we have IDs to extract and file-based streams exist
        stream_folders = config_data.get("stream_folders") or {}
        extracted = {}
        has_per_id = pool_file_maps and any(fm for fm in pool_file_maps)
        has_shared = pool_shared_maps and any(sm for sm in pool_shared_maps)
        if filtered_ids_for_lookup and (has_per_id or has_shared):
            # Use multi-pool extraction if we have the implicit PANDA_SOURCE
            # column (auto-tagged when concat ran over multiple inputs) and
            # multiple pools.
            if len(pool_file_maps) > 1 and PANDA_SOURCE in result_df.columns:
                extracted = extract_from_multiple_pools(
                    result_df, pool_folders, output_dir,
                    rename_map=original_to_new_id if rename else None,
                    file_maps=pool_file_maps,
                    ignore_missing=ignore_missing,
                    stream_folders=stream_folders,
                    shared_maps=pool_shared_maps,
                )
            else:
                # Single pool mode - use first pool
                extracted = extract_pool_data_for_filtered_ids(
                    filtered_ids_for_lookup, pool_folders[0], output_dir,
                    rename_map=original_to_new_id if rename else None,
                    file_map=pool_file_maps[0] if pool_file_maps else None,
                    ignore_missing=ignore_missing,
                    stream_folders=stream_folders,
                    shared_map=pool_shared_maps[0] if pool_shared_maps else None,
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
                rename_map=original_to_new_id if rename else None,
                pool_table_targets=config_data.get("pool_table_targets") or None,
            )
            print(f"Copied {len(copied_tables)} tables")

        # Attach `<stream>.id` provenance for each output stream so
        # downstream tools can link renamed pool IDs back to their upstream
        # parents (e.g. Panda_1 -> a) via the map_table provenance columns.
        #
        # Two cases depending on whether a content table was already
        # written at `map_path` by filter_and_copy_pool_tables above:
        #   (a) Content table exists (Sequence.sequences, Ligand.compounds,
        #       any stream whose "map_table" IS its content). Augment it
        #       in place: insert/override `<stream>.id` right after `id`;
        #       leave every other column (including `pool.id`) untouched.
        #       Do not add a `file` column.
        #   (b) No content table yet (pure per-ID file streams like
        #       structures.pdb). Write a provenance-only CSV with the
        #       canonical `id, file, <stream>.id` schema.
        stream_map_targets = config_data.get('stream_map_targets', [])
        prune_provenance = config_data.get('prune_redundant_provenance', True)
        if stream_map_targets and not result_df.empty and 'id' in result_df.columns:
            has_original = 'original_id' in result_df.columns
            parent_by_id = {}
            for _, r in result_df.iterrows():
                out_id = str(r['id'])
                parent_by_id[out_id] = (
                    str(r['original_id']) if has_original else out_id
                )

            for tgt in stream_map_targets:
                stream_name = tgt.get('stream_name')
                map_path = tgt.get('map_table')
                file_template = tgt.get('file_template', '') or ''
                if not stream_name or not map_path:
                    continue
                axis_col = f"{stream_name}.id"
                os.makedirs(os.path.dirname(map_path), exist_ok=True)

                if os.path.exists(map_path):
                    # (a) Augment the existing content table in place.
                    existing = pd.read_csv(
                        map_path, keep_default_na=False,
                        na_values=['NA', 'N/A', '#N/A'],
                    )
                    if 'id' not in existing.columns:
                        print(
                            f"  Warning: {map_path} has no 'id' column, "
                            f"cannot attach {axis_col} provenance"
                        )
                        continue
                    # Preserve older generations: shift any existing
                    # <stream>.id / <stream>.-N.id one step back before
                    # writing the new immediate-parent column. Without
                    # this, chained Panda steps would overwrite the
                    # previous parent and lose the chain.
                    existing = shift_provenance_columns(existing, stream_name)
                    existing[axis_col] = existing['id'].astype(str).map(
                        parent_by_id
                    ).fillna(existing['id'].astype(str))
                    # Place generation columns immediately after `id` in
                    # generation order — `id, <stream>.id, <stream>.-1.id,
                    # <stream>.-2.id, ...` — then everything else (including
                    # `pool.id`) in its original position.
                    gen_cols = order_stream_provenance_columns(
                        existing, stream_name
                    )
                    other_cols = [c for c in existing.columns if c not in gen_cols]
                    id_pos = other_cols.index('id')
                    new_order = (
                        other_cols[:id_pos + 1]
                        + gen_cols
                        + other_cols[id_pos + 1:]
                    )
                    existing = existing[new_order]
                    if prune_provenance:
                        existing = prune_redundant_provenance_columns(existing)
                    existing.to_csv(map_path, index=False)
                    print(f"Augmented content table with {axis_col}: {map_path}")
                else:
                    # (b) No content table — write a provenance-only map_table.
                    rows_out = []
                    has_source_table = PANDA_SOURCE in result_df.columns
                    for _, r in result_df.iterrows():
                        out_id = str(r['id'])
                        parent = parent_by_id[out_id]
                        lookup_id = (
                            str(r['original_id']) if has_original else out_id
                        )
                        source_idx = int(r.get(PANDA_SOURCE, 0)) if has_source_table else 0
                        file_val = (file_template.replace('<id>', out_id)
                                    if '<id>' in file_template else file_template)
                        row = {
                            'id': out_id,
                            'file': file_val,
                        }
                        row.update(load_upstream_provenance_row(
                            pool_stream_jsons,
                            stream_name,
                            source_idx,
                            lookup_id,
                        ))
                        rows_out.append(row)

                    out_df = pd.DataFrame(rows_out)
                    if not out_df.empty:
                        out_df = shift_provenance_columns(out_df, stream_name)
                        out_df[axis_col] = out_df['id'].astype(str).map(
                            parent_by_id
                        ).fillna(out_df['id'].astype(str))

                        gen_cols = order_stream_provenance_columns(
                            out_df, stream_name
                        )
                        other_cols = [c for c in out_df.columns if c not in gen_cols]
                        id_pos = other_cols.index('id')
                        new_order = (
                            other_cols[:id_pos + 1]
                            + gen_cols
                            + other_cols[id_pos + 1:]
                        )
                        out_df = out_df[new_order]
                        if prune_provenance:
                            out_df = prune_redundant_provenance_columns(out_df)

                    with open(map_path, 'w', newline='') as mf:
                        fieldnames = (
                            list(out_df.columns)
                            if not out_df.empty
                            else ['id', 'file', axis_col]
                        )
                        writer = csv.DictWriter(mf, fieldnames=fieldnames)
                        writer.writeheader()
                        writer.writerows(out_df.to_dict('records'))
                    print(f"Wrote stream map_table: {map_path}")

    # Create missing.csv
    missing_csv_path = config_data.get('missing_csv')
    if original_ids:
        output_dir = os.path.dirname(output_csv)
        step_tool_name = config_data.get('step_tool_name') or os.path.basename(output_dir)
        create_missing_csv(original_ids, filtered_ids_for_lookup, output_dir,
                           step_tool_name, operations,
                           rename_map=original_to_new_id if original_to_new_id else None,
                           removed_by_op=removed_by_op,
                           missing_csv_path=missing_csv_path)

    # Merge upstream missing tables (from pool sources)
    upstream_missing_paths = config_data.get('upstream_missing_paths', [])
    if upstream_missing_paths:
        output_dir = os.path.dirname(output_csv)
        print(f"\nMerging {len(upstream_missing_paths)} upstream missing table(s)")
        merge_upstream_missing(output_dir, upstream_missing_paths,
                               missing_csv_path=missing_csv_path)

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
