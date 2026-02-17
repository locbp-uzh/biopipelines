#!/usr/bin/env python3
# Copyright (c) 2026 Gianluca Quargnali @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
BioPipelines I/O Utilities for SLURM Runtime.

This module provides functions for pipe_*.py scripts to read DataStreams and tables
at SLURM runtime. It handles wildcard resolution, ID mapping, and table lookups.

Functions:
    # DataStream utilities
    load_datastream: Load DataStream from JSON file or dict
    iterate_files: Iterate (id, resolved_file_path) for file-based streams
    iterate_values: Iterate (id, value_dict) for value-based streams
    resolve_file: Resolve single file for an ID (handles wildcards)
    get_value: Get single value from map_table for an ID
    get_all_values: Get all values from map_table for an ID

    # Table reference utilities (for per-structure lookups from external tables)
    load_table: Load a table from path or DATASHEET_REFERENCE string
    lookup_table_value: Look up a value from a table for a given ID
    iterate_table_values: Iterate over (id, value) pairs from a table

Example usage:
    from biopipelines_io import load_datastream, iterate_files, iterate_values

    # Load from JSON or dict
    structures_ds = load_datastream("/path/to/structures.json")

    # Iterate over structure files (wildcards resolved automatically)
    for struct_id, struct_file in iterate_files(structures_ds):
        print(f"Processing {struct_id}: {struct_file}")

    # Iterate over compound values (SMILES from map_table)
    for comp_id, values in iterate_values(compounds_ds, columns=['smiles', 'name']):
        print(f"Compound {comp_id}: {values['smiles']}")

    # Table reference for per-structure data (e.g., fixed positions from DistanceSelector)
    from biopipelines_io import load_table, lookup_table_value, iterate_table_values

    # Load table (supports DATASHEET_REFERENCE:path:column format)
    table, column = load_table("DATASHEET_REFERENCE:/path/to/positions.csv:within")

    # Look up value for a specific structure
    fixed_positions = lookup_table_value(table, "protein_1", column)

    # Or iterate over all structures
    for struct_id, positions in iterate_table_values(table, structure_ids, column):
        print(f"{struct_id}: {positions}")
"""

import glob
import json
import os
from dataclasses import dataclass, field
from typing import Any, Dict, Iterator, List, Optional, Tuple, Union

import pandas as pd


@dataclass
class DataStreamRuntime:
    """
    Runtime representation of a DataStream for SLURM job execution.

    This is a lightweight version of DataStream designed for use in pipe scripts
    at SLURM runtime. It supports lazy loading of map_table data.

    Attributes:
        name: Name of this data stream
        ids: List of unique identifiers for each item
        files: List of file paths (may contain wildcards or be empty)
        map_table: Path to CSV file mapping ids to files/values
        format: Data format ("pdb", "cif", "fasta", "csv", "sdf", "smiles", etc.)
        metadata: Additional tool-specific metadata
        files_contain_wildcards: If True, file paths contain glob patterns
    """
    name: str = ""
    ids: List[str] = field(default_factory=list)
    files: List[str] = field(default_factory=list)
    map_table: str = ""
    format: str = "pdb"
    metadata: Dict[str, Any] = field(default_factory=dict)
    files_contain_wildcards: bool = False

    # Internal cache for map_table data (loaded lazily)
    _map_data: Optional[pd.DataFrame] = field(default=None, repr=False, compare=False)

    def _get_map_data(self) -> pd.DataFrame:
        """
        Lazily load map_table data.

        Returns:
            DataFrame with map_table contents

        Raises:
            FileNotFoundError: If map_table path is set but file doesn't exist
            ValueError: If map_table is not set
        """
        if self._map_data is not None:
            return self._map_data

        if not self.map_table:
            raise ValueError(f"DataStream '{self.name}' has no map_table configured")

        if not os.path.exists(self.map_table):
            raise FileNotFoundError(f"map_table not found: {self.map_table}")

        self._map_data = pd.read_csv(self.map_table)
        return self._map_data


def load_datastream(source: Union[str, Dict[str, Any]]) -> DataStreamRuntime:
    """
    Load a DataStream from a JSON file or dictionary.

    Args:
        source: Either a path to a JSON file or a dictionary with DataStream fields

    Returns:
        DataStreamRuntime instance

    Raises:
        FileNotFoundError: If source is a path that doesn't exist
        ValueError: If required fields (name, ids) are missing
        json.JSONDecodeError: If source is a file with invalid JSON

    Example:
        # From JSON file
        ds = load_datastream("/path/to/structures.json")

        # From dict
        ds = load_datastream({
            "name": "compounds",
            "ids": ["lig_1", "lig_2"],
            "files": [],
            "map_table": "/path/to/compounds_map.csv",
            "format": "smiles"
        })
    """
    if isinstance(source, str):
        # Source is a file path
        if not os.path.exists(source):
            raise FileNotFoundError(f"DataStream JSON file not found: {source}")

        with open(source, 'r') as f:
            data = json.load(f)
    elif isinstance(source, dict):
        data = source
    else:
        raise ValueError(f"source must be a file path (str) or dict, got {type(source)}")

    # Validate required fields
    if 'name' not in data:
        raise ValueError("DataStream missing required field: 'name'")
    if 'ids' not in data:
        raise ValueError("DataStream missing required field: 'ids'")

    # files defaults to empty list
    files = data.get('files', [])

    # Validate files/ids relationship
    if len(files) > 1 and len(files) != len(data['ids']):
        raise ValueError(
            f"Length mismatch: {len(data['ids'])} ids but {len(files)} files. "
            f"Use empty files list or single file for table-based data, "
            f"or provide one file per id."
        )

    return DataStreamRuntime(
        name=data['name'],
        ids=data['ids'],
        files=files,
        map_table=data.get('map_table', ''),
        format=data.get('format', 'pdb'),
        metadata=data.get('metadata', {}),
        files_contain_wildcards=data.get('files_contain_wildcards', False)
    )


def _find_best_match(item_id: str, file_paths: List[str]) -> str:
    """
    Find the best matching file for an ID from a list of file paths.

    Matching priority:
    1. Exact basename match (without extension)
    2. Basename starts with item_id + "_" or "-"
    3. item_id contained in basename
    4. First file (fallback)

    Args:
        item_id: The item identifier to match
        file_paths: List of file paths to search

    Returns:
        Best matching file path

    Raises:
        FileNotFoundError: If file_paths is empty
    """
    if not file_paths:
        raise FileNotFoundError(f"No files found for ID '{item_id}'")

    # Priority 1: Exact basename match
    for fp in file_paths:
        basename = os.path.splitext(os.path.basename(fp))[0]
        if basename == item_id:
            return fp

    # Priority 2: Basename starts with item_id + separator
    for fp in file_paths:
        basename = os.path.splitext(os.path.basename(fp))[0]
        if basename.startswith(f"{item_id}_") or basename.startswith(f"{item_id}-"):
            return fp

    # Priority 3: item_id contained in basename
    for fp in file_paths:
        basename = os.path.splitext(os.path.basename(fp))[0]
        if item_id in basename:
            return fp

    # Priority 4: Fallback to first file
    return file_paths[0]


def iterate_files(ds: DataStreamRuntime) -> Iterator[Tuple[str, str]]:
    """
    Iterate over (id, resolved_file_path) pairs for file-based streams.

    Handles three cases:
    1. files_contain_wildcards=True: Expand globs, match to IDs
    2. len(files) == len(ids): Direct zip iteration
    3. len(files) == 1: Single file for all IDs (bundle case)

    Args:
        ds: DataStreamRuntime instance

    Yields:
        Tuple of (item_id, file_path)

    Raises:
        FileNotFoundError: If wildcard expansion yields no files for an ID
        ValueError: If files list is empty and this is a file-based stream

    Example:
        for struct_id, struct_file in iterate_files(structures_ds):
            process_structure(struct_id, struct_file)
    """
    if not ds.files:
        raise ValueError(f"DataStream '{ds.name}' has no files configured")

    if ds.files_contain_wildcards:
        # Case 1: Wildcards - expand and match
        for i, item_id in enumerate(ds.ids):
            # Get the pattern for this ID (may be per-ID or single pattern)
            if len(ds.files) == len(ds.ids):
                pattern = ds.files[i]
            elif len(ds.files) == 1:
                pattern = ds.files[0]
            else:
                raise ValueError(
                    f"Cannot match {len(ds.files)} patterns to {len(ds.ids)} IDs"
                )

            # Expand glob pattern
            expanded = glob.glob(pattern)
            if not expanded:
                raise FileNotFoundError(
                    f"No files found matching pattern '{pattern}' for ID '{item_id}'"
                )

            # Find best match for this ID
            matched_file = _find_best_match(item_id, expanded)
            yield (item_id, matched_file)

    elif len(ds.files) == len(ds.ids):
        # Case 2: Direct mapping
        for item_id, file_path in zip(ds.ids, ds.files):
            yield (item_id, file_path)

    elif len(ds.files) == 1:
        # Case 3: Single file for all IDs (bundle case)
        single_file = ds.files[0]
        for item_id in ds.ids:
            yield (item_id, single_file)

    else:
        raise ValueError(
            f"Cannot iterate: {len(ds.ids)} ids but {len(ds.files)} files"
        )


def iterate_values(
    ds: DataStreamRuntime,
    columns: Optional[List[str]] = None
) -> Iterator[Tuple[str, Dict[str, Any]]]:
    """
    Iterate over (id, value_dict) pairs for value-based streams.

    Reads values from the map_table CSV and yields dictionaries with
    requested columns for each ID.

    Args:
        ds: DataStreamRuntime instance
        columns: List of column names to include (None = all columns)

    Yields:
        Tuple of (item_id, dict_of_column_values)

    Raises:
        ValueError: If map_table is not configured
        FileNotFoundError: If map_table file doesn't exist
        KeyError: If requested column doesn't exist in map_table

    Example:
        for comp_id, values in iterate_values(compounds_ds, columns=['smiles']):
            smiles = values['smiles']
            process_compound(comp_id, smiles)
    """
    map_data = ds._get_map_data()

    # Validate requested columns exist
    if columns:
        missing_cols = set(columns) - set(map_data.columns)
        if missing_cols:
            raise KeyError(
                f"Columns not found in map_table: {missing_cols}. "
                f"Available columns: {list(map_data.columns)}"
            )

    # Determine which columns to include
    include_cols = columns if columns else list(map_data.columns)

    # Build index for fast lookup
    if 'id' not in map_data.columns:
        raise KeyError("map_table must contain 'id' column")

    id_to_row = {row['id']: row for _, row in map_data.iterrows()}

    for item_id in ds.ids:
        if item_id not in id_to_row:
            raise KeyError(f"ID '{item_id}' not found in map_table")

        row = id_to_row[item_id]
        values = {col: row[col] for col in include_cols}
        yield (item_id, values)


def resolve_file(ds: DataStreamRuntime, item_id: str) -> str:
    """
    Resolve a single file path for an ID.

    Handles wildcard patterns by expanding and matching.

    Args:
        ds: DataStreamRuntime instance
        item_id: The item identifier

    Returns:
        Resolved file path

    Raises:
        ValueError: If ID not found in DataStream or files not configured
        FileNotFoundError: If wildcard expansion yields no files

    Example:
        file_path = resolve_file(structures_ds, "rank0001")
    """
    if not ds.files:
        raise ValueError(f"DataStream '{ds.name}' has no files configured")

    if item_id not in ds.ids:
        raise ValueError(f"ID '{item_id}' not found in DataStream '{ds.name}'")

    idx = ds.ids.index(item_id)

    if ds.files_contain_wildcards:
        # Get pattern for this ID
        if len(ds.files) == len(ds.ids):
            pattern = ds.files[idx]
        elif len(ds.files) == 1:
            pattern = ds.files[0]
        else:
            raise ValueError(
                f"Cannot match {len(ds.files)} patterns to {len(ds.ids)} IDs"
            )

        # Expand and match
        expanded = glob.glob(pattern)
        if not expanded:
            raise FileNotFoundError(
                f"No files found matching pattern '{pattern}' for ID '{item_id}'"
            )

        return _find_best_match(item_id, expanded)

    elif len(ds.files) == len(ds.ids):
        return ds.files[idx]

    elif len(ds.files) == 1:
        return ds.files[0]

    else:
        raise ValueError(
            f"Cannot resolve file: {len(ds.ids)} ids but {len(ds.files)} files"
        )


def get_value(ds: DataStreamRuntime, item_id: str, column: str = "value") -> Any:
    """
    Get a single value from the map_table for an ID.

    Args:
        ds: DataStreamRuntime instance
        item_id: The item identifier
        column: Column name to retrieve (default: "value")

    Returns:
        Value from the specified column

    Raises:
        ValueError: If map_table not configured
        FileNotFoundError: If map_table doesn't exist
        KeyError: If ID or column not found

    Example:
        smiles = get_value(compounds_ds, "ligand_001", column="smiles")
    """
    map_data = ds._get_map_data()

    if column not in map_data.columns:
        raise KeyError(
            f"Column '{column}' not found in map_table. "
            f"Available columns: {list(map_data.columns)}"
        )

    if 'id' not in map_data.columns:
        raise KeyError("map_table must contain 'id' column")

    row = map_data[map_data['id'] == item_id]
    if row.empty:
        raise KeyError(f"ID '{item_id}' not found in map_table")

    return row.iloc[0][column]


def get_all_values(ds: DataStreamRuntime, item_id: str) -> Dict[str, Any]:
    """
    Get all column values from the map_table for an ID.

    Args:
        ds: DataStreamRuntime instance
        item_id: The item identifier

    Returns:
        Dictionary with all column values

    Raises:
        ValueError: If map_table not configured
        FileNotFoundError: If map_table doesn't exist
        KeyError: If ID not found

    Example:
        all_data = get_all_values(compounds_ds, "ligand_001")
        print(f"SMILES: {all_data['smiles']}, Name: {all_data['name']}")
    """
    map_data = ds._get_map_data()

    if 'id' not in map_data.columns:
        raise KeyError("map_table must contain 'id' column")

    row = map_data[map_data['id'] == item_id]
    if row.empty:
        raise KeyError(f"ID '{item_id}' not found in map_table")

    return row.iloc[0].to_dict()


# =============================================================================
# Table Reference Utilities
# =============================================================================
# These functions handle per-structure lookups from external tables, commonly
# used for data like fixed/designed positions from DistanceSelector, or other
# per-structure metadata stored in CSV tables.


# Cache for loaded tables to avoid repeated file reads
_table_cache: Dict[str, pd.DataFrame] = {}


def load_table(
    reference: str,
    id_map: Optional[Dict[str, str]] = None
) -> Tuple[pd.DataFrame, Optional[str]]:
    """
    Load a table from a path or DATASHEET_REFERENCE string.

    Supports two formats:
    1. Direct path: "/path/to/table.csv" - returns (DataFrame, None)
    2. Reference: "DATASHEET_REFERENCE:/path/to/table.csv:column" - returns (DataFrame, column)

    Args:
        reference: Either a file path or DATASHEET_REFERENCE:path:column string
        id_map: Optional ID mapping pattern (currently stored for later use)

    Returns:
        Tuple of (DataFrame, column_name or None)

    Raises:
        FileNotFoundError: If table file doesn't exist
        ValueError: If reference format is invalid

    Example:
        # Direct path
        table, _ = load_table("/path/to/positions.csv")

        # Reference with column
        table, column = load_table("DATASHEET_REFERENCE:/path/to/positions.csv:within")
        value = lookup_table_value(table, "protein_1", column)
    """
    column_name = None

    if reference.startswith("DATASHEET_REFERENCE:"):
        # Parse reference: DATASHEET_REFERENCE:path:column
        # Handle Windows paths with drive letters (e.g., C:\path\to\file.csv)
        # Format: DATASHEET_REFERENCE:<path>:<column>
        # The column is always the LAST colon-separated part
        prefix = "DATASHEET_REFERENCE:"
        remainder = reference[len(prefix):]

        # Find the last colon that's not part of a Windows drive letter
        # Windows drive letters are single letter followed by colon at position 1
        last_colon = remainder.rfind(":")

        # Check if this colon is part of a Windows path (e.g., "C:")
        # A Windows drive letter would be at position 1 (after single letter)
        if last_colon == 1 and remainder[0].isalpha():
            # Only the drive letter colon exists, no column specified
            raise ValueError(
                f"Invalid DATASHEET_REFERENCE format: '{reference}'. "
                f"Expected 'DATASHEET_REFERENCE:path:column'"
            )

        if last_colon == -1:
            raise ValueError(
                f"Invalid DATASHEET_REFERENCE format: '{reference}'. "
                f"Expected 'DATASHEET_REFERENCE:path:column'"
            )

        table_path = remainder[:last_colon]
        column_name = remainder[last_colon + 1:]

        if not table_path or not column_name:
            raise ValueError(
                f"Invalid DATASHEET_REFERENCE format: '{reference}'. "
                f"Expected 'DATASHEET_REFERENCE:path:column'"
            )
    else:
        # Direct path
        table_path = reference

    # Check cache first
    if table_path in _table_cache:
        return _table_cache[table_path], column_name

    # Load table
    if not os.path.exists(table_path):
        raise FileNotFoundError(f"Table not found: {table_path}")

    df = pd.read_csv(table_path)
    _table_cache[table_path] = df

    return df, column_name


def lookup_table_value(
    table: pd.DataFrame,
    item_id: str,
    column: str,
    id_map: Optional[Dict[str, str]] = None,
    id_column: str = "id",
    pdb_column: str = "pdb"
) -> Any:
    """
    Look up a value from a table for a given ID.

    Tries multiple lookup strategies in order:
    1. Match by pdb column (if exists): item_id or item_id.pdb
    2. Match by id column: exact match
    3. Match by id column with ID mapping: strip suffixes according to id_map

    Args:
        table: DataFrame to search
        item_id: The item identifier (structure ID, typically filename without extension)
        column: Column name to retrieve value from
        id_map: Optional ID mapping pattern (e.g., {"*": "*_<N>"}) for suffix stripping
        id_column: Column name for ID lookups (default: "id")
        pdb_column: Column name for PDB filename lookups (default: "pdb")

    Returns:
        Value from the specified column

    Raises:
        KeyError: If ID not found in table or column doesn't exist

    Example:
        table, column = load_table("DATASHEET_REFERENCE:/path/to/positions.csv:within")
        positions = lookup_table_value(table, "protein_1", column)

        # With ID mapping (e.g., "protein_1_2" maps to "protein_1" in table)
        positions = lookup_table_value(table, "protein_1_2", column, id_map={"*": "*_<N>"})
    """
    if column not in table.columns:
        raise KeyError(
            f"Column '{column}' not found in table. "
            f"Available columns: {list(table.columns)}"
        )

    # Track attempted IDs for error reporting
    attempted_ids = []

    # Strategy 1: Try pdb column match (with and without extension)
    if pdb_column in table.columns:
        # Try with .pdb extension
        pdb_name = f"{item_id}.pdb"
        attempted_ids.append(pdb_name)
        matching_rows = table[table[pdb_column] == pdb_name]
        if not matching_rows.empty:
            return matching_rows.iloc[0][column]

        # Try exact match
        attempted_ids.append(item_id)
        matching_rows = table[table[pdb_column] == item_id]
        if not matching_rows.empty:
            return matching_rows.iloc[0][column]

    # Strategy 2: Try id column exact match
    if id_column in table.columns:
        if item_id not in attempted_ids:
            attempted_ids.append(item_id)
        matching_rows = table[table[id_column] == item_id]
        if not matching_rows.empty:
            return matching_rows.iloc[0][column]

        # Strategy 3: Try ID mapping (strip suffixes)
        if id_map:
            try:
                from id_map_utils import map_table_ids_to_ids
                candidate_ids = map_table_ids_to_ids(item_id, id_map)

                for candidate_id in candidate_ids:
                    if candidate_id == item_id:
                        continue  # Already tried
                    attempted_ids.append(candidate_id)
                    matching_rows = table[table[id_column] == candidate_id]
                    if not matching_rows.empty:
                        return matching_rows.iloc[0][column]
            except ImportError:
                # id_map_utils not available, skip this strategy
                pass

    raise KeyError(
        f"ID '{item_id}' not found in table. Tried: {attempted_ids}"
    )


def iterate_table_values(
    table: pd.DataFrame,
    item_ids: List[str],
    column: str,
    id_map: Optional[Dict[str, str]] = None,
    id_column: str = "id",
    pdb_column: str = "pdb"
) -> Iterator[Tuple[str, Any]]:
    """
    Iterate over (id, value) pairs from a table for a list of IDs.

    Args:
        table: DataFrame to search
        item_ids: List of item identifiers to look up
        column: Column name to retrieve values from
        id_map: Optional ID mapping pattern for suffix stripping
        id_column: Column name for ID lookups (default: "id")
        pdb_column: Column name for PDB filename lookups (default: "pdb")

    Yields:
        Tuple of (item_id, value)

    Raises:
        KeyError: If any ID not found in table or column doesn't exist

    Example:
        table, column = load_table("DATASHEET_REFERENCE:/path/to/positions.csv:within")
        for struct_id, positions in iterate_table_values(table, structure_ids, column):
            print(f"{struct_id}: {positions}")
    """
    for item_id in item_ids:
        value = lookup_table_value(
            table, item_id, column, id_map, id_column, pdb_column
        )
        yield (item_id, value)


def clear_table_cache() -> None:
    """
    Clear the table cache.

    Call this if you need to reload tables that may have changed on disk.
    """
    global _table_cache
    _table_cache = {}


# =============================================================================
# Provenance-based ID Resolution
# =============================================================================
# When tools like Panda auto-rename IDs (e.g., Panda_1, Panda_2), the original
# source IDs are preserved in map_table provenance columns (e.g., structures.id).
# These functions trace back through provenance to match renamed IDs to their
# source IDs in external tables.


def resolve_id_by_provenance(
    item_id: str,
    target_ids: set,
    map_table_paths: List[str]
) -> Optional[str]:
    """
    Resolve an ID to a matching target ID using map_table provenance columns.

    When structures pass through Panda with auto-rename (e.g., Panda_1 -> Panda_2),
    their map_table contains provenance columns like 'structures.id' that track
    the original source IDs. This function looks up item_id in those map_tables
    and returns the provenance source ID that exists in target_ids.

    Args:
        item_id: The ID to resolve (e.g., 'Panda_1')
        target_ids: Set of IDs we're trying to match against
        map_table_paths: List of map_table CSV paths to search for provenance

    Returns:
        Matching target ID if found via provenance, None otherwise

    Example:
        # Panda_1's map_table has: id=Panda_1, structures.id=LID_Redesign_001_1
        # Selection table has: id=LID_Redesign_001_1, fixed=10-20+30-40
        matched = resolve_id_by_provenance("Panda_1", {"LID_Redesign_001_1"}, ["/path/to/map.csv"])
        # Returns "LID_Redesign_001_1"
    """
    for map_table_path in map_table_paths:
        if not map_table_path or not os.path.exists(map_table_path):
            continue

        # Use table cache
        if map_table_path in _table_cache:
            df = _table_cache[map_table_path]
        else:
            try:
                df = pd.read_csv(map_table_path)
                _table_cache[map_table_path] = df
            except Exception:
                continue

        if "id" not in df.columns:
            continue

        # Find the row for this item_id
        row = df[df["id"].astype(str) == item_id]
        if row.empty:
            continue

        row = row.iloc[0]

        # Check provenance columns (columns ending in '.id' like 'structures.id')
        provenance_cols = [col for col in df.columns if col.endswith(".id") and col != "id"]

        for prov_col in provenance_cols:
            source_id = str(row[prov_col])
            if source_id in target_ids:
                return source_id

    return None

