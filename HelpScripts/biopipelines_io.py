#!/usr/bin/env python3
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
    from pipe_biopipelines_io import load_datastream, iterate_files, iterate_values

    # Load from JSON or dict
    structures_ds = load_datastream("/path/to/structures.json")

    # Iterate over structure files (wildcards resolved automatically)
    for struct_id, struct_file in iterate_files(structures_ds):
        print(f"Processing {struct_id}: {struct_file}")

    # Iterate over compound values (SMILES from map_table)
    for comp_id, values in iterate_values(compounds_ds, columns=['smiles', 'name']):
        print(f"Compound {comp_id}: {values['smiles']}")

    # Table reference for per-structure data (e.g., fixed positions from DistanceSelector)
    from pipe_biopipelines_io import load_table, lookup_table_value, iterate_table_values

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


if __name__ == "__main__":
    # Run tests when executed directly
    import tempfile
    import shutil

    print("Running pipe_datastream_io tests...\n")

    # Create temporary directory for test files
    test_dir = tempfile.mkdtemp()

    try:
        # Create test files
        test_files_dir = os.path.join(test_dir, "structures")
        os.makedirs(test_files_dir)

        # Create structure files
        for i in range(1, 4):
            with open(os.path.join(test_files_dir, f"protein_{i}.pdb"), 'w') as f:
                f.write(f"ATOM protein_{i}\n")

        # Create wildcard test files (e.g., rank0001_design123.cif)
        wildcard_dir = os.path.join(test_dir, "wildcard_structures")
        os.makedirs(wildcard_dir)
        for i in range(1, 4):
            with open(os.path.join(wildcard_dir, f"rank000{i}_design{i*100}.cif"), 'w') as f:
                f.write(f"ATOM rank000{i}\n")

        # Create map_table CSV
        map_table_path = os.path.join(test_dir, "compounds_map.csv")
        pd.DataFrame({
            'id': ['lig_1', 'lig_2', 'lig_3'],
            'smiles': ['CCO', 'CC(=O)O', 'c1ccccc1'],
            'name': ['ethanol', 'acetic_acid', 'benzene'],
            'value': ['CCO', 'CC(=O)O', 'c1ccccc1']
        }).to_csv(map_table_path, index=False)

        # Test 1: Load from dict
        print("Test 1: load_datastream from dict")
        ds_dict = {
            "name": "structures",
            "ids": ["protein_1", "protein_2", "protein_3"],
            "files": [
                os.path.join(test_files_dir, "protein_1.pdb"),
                os.path.join(test_files_dir, "protein_2.pdb"),
                os.path.join(test_files_dir, "protein_3.pdb"),
            ],
            "format": "pdb"
        }
        ds = load_datastream(ds_dict)
        assert ds.name == "structures"
        assert len(ds.ids) == 3
        print("  PASS: Loaded DataStream from dict\n")

        # Test 2: Load from JSON file
        print("Test 2: load_datastream from JSON file")
        json_path = os.path.join(test_dir, "test_ds.json")
        with open(json_path, 'w') as f:
            json.dump(ds_dict, f)
        ds_from_json = load_datastream(json_path)
        assert ds_from_json.name == "structures"
        print("  PASS: Loaded DataStream from JSON\n")

        # Test 3: iterate_files - direct mapping
        print("Test 3: iterate_files with direct mapping")
        results = list(iterate_files(ds))
        assert len(results) == 3
        assert results[0] == ("protein_1", os.path.join(test_files_dir, "protein_1.pdb"))
        print(f"  PASS: Got {len(results)} file pairs\n")

        # Test 4: iterate_files - single file (bundle case)
        print("Test 4: iterate_files with single file (bundle)")
        bundle_ds = DataStreamRuntime(
            name="bundle",
            ids=["a", "b", "c"],
            files=[os.path.join(test_files_dir, "protein_1.pdb")]
        )
        bundle_results = list(iterate_files(bundle_ds))
        assert len(bundle_results) == 3
        assert all(r[1] == os.path.join(test_files_dir, "protein_1.pdb") for r in bundle_results)
        print("  PASS: Single file returned for all IDs\n")

        # Test 5: iterate_files - with wildcards
        print("Test 5: iterate_files with wildcards")
        wildcard_ds = DataStreamRuntime(
            name="wildcard_structures",
            ids=["rank0001", "rank0002", "rank0003"],
            files=[
                os.path.join(wildcard_dir, "rank0001_*.cif"),
                os.path.join(wildcard_dir, "rank0002_*.cif"),
                os.path.join(wildcard_dir, "rank0003_*.cif"),
            ],
            files_contain_wildcards=True
        )
        wildcard_results = list(iterate_files(wildcard_ds))
        assert len(wildcard_results) == 3
        for item_id, file_path in wildcard_results:
            assert item_id in file_path
            assert os.path.exists(file_path)
        print(f"  PASS: Resolved {len(wildcard_results)} wildcard patterns\n")

        # Test 6: iterate_values
        print("Test 6: iterate_values")
        compounds_ds = DataStreamRuntime(
            name="compounds",
            ids=["lig_1", "lig_2", "lig_3"],
            files=[],
            map_table=map_table_path,
            format="smiles"
        )
        value_results = list(iterate_values(compounds_ds, columns=['smiles', 'name']))
        assert len(value_results) == 3
        assert value_results[0] == ("lig_1", {"smiles": "CCO", "name": "ethanol"})
        print(f"  PASS: Got values for {len(value_results)} compounds\n")

        # Test 7: resolve_file
        print("Test 7: resolve_file")
        resolved = resolve_file(ds, "protein_2")
        assert resolved == os.path.join(test_files_dir, "protein_2.pdb")
        print(f"  PASS: Resolved file: {os.path.basename(resolved)}\n")

        # Test 8: get_value
        print("Test 8: get_value")
        smiles = get_value(compounds_ds, "lig_2", column="smiles")
        assert smiles == "CC(=O)O"
        print(f"  PASS: Got SMILES: {smiles}\n")

        # Test 9: get_all_values
        print("Test 9: get_all_values")
        all_vals = get_all_values(compounds_ds, "lig_3")
        assert all_vals['smiles'] == 'c1ccccc1'
        assert all_vals['name'] == 'benzene'
        print(f"  PASS: Got all values: {all_vals}\n")

        # Test 10: Error cases - missing ID
        print("Test 10: Error handling - missing ID")
        try:
            get_value(compounds_ds, "nonexistent", column="smiles")
            print("  FAIL: Should have raised KeyError")
        except KeyError as e:
            assert "nonexistent" in str(e)
            print(f"  PASS: Raised KeyError as expected\n")

        # Test 11: Error cases - missing column
        print("Test 11: Error handling - missing column")
        try:
            get_value(compounds_ds, "lig_1", column="nonexistent_column")
            print("  FAIL: Should have raised KeyError")
        except KeyError as e:
            assert "nonexistent_column" in str(e)
            print(f"  PASS: Raised KeyError as expected\n")

        # Test 12: Error cases - wildcard with no matches
        print("Test 12: Error handling - wildcard no matches")
        bad_wildcard_ds = DataStreamRuntime(
            name="bad_wildcard",
            ids=["missing"],
            files=[os.path.join(wildcard_dir, "nonexistent_*.xyz")],
            files_contain_wildcards=True
        )
        try:
            list(iterate_files(bad_wildcard_ds))
            print("  FAIL: Should have raised FileNotFoundError")
        except FileNotFoundError as e:
            assert "No files found" in str(e)
            print(f"  PASS: Raised FileNotFoundError as expected\n")

        # =====================================================================
        # Table Reference Tests
        # =====================================================================
        print("=" * 50)
        print("Table Reference Tests\n")

        # Create a positions table (like DistanceSelector output)
        positions_table_path = os.path.join(test_dir, "positions.csv")
        pd.DataFrame({
            'id': ['protein_1', 'protein_2', 'protein_3'],
            'pdb': ['protein_1.pdb', 'protein_2.pdb', 'protein_3.pdb'],
            'within': ['10-20+30-40', '15-25', '5-10'],
            'beyond': ['1-9+21-29+41-100', '1-14+26-100', '1-4+11-100']
        }).to_csv(positions_table_path, index=False)

        # Test 13: load_table - direct path
        print("Test 13: load_table with direct path")
        clear_table_cache()  # Clear cache for clean test
        table, col = load_table(positions_table_path)
        assert col is None
        assert len(table) == 3
        assert 'within' in table.columns
        print("  PASS: Loaded table from direct path\n")

        # Test 14: load_table - DATASHEET_REFERENCE format
        print("Test 14: load_table with DATASHEET_REFERENCE")
        clear_table_cache()
        ref_string = f"DATASHEET_REFERENCE:{positions_table_path}:within"
        table, col = load_table(ref_string)
        assert col == "within"
        assert len(table) == 3
        print(f"  PASS: Parsed reference, column='{col}'\n")

        # Test 15: lookup_table_value - by id
        print("Test 15: lookup_table_value by id")
        value = lookup_table_value(table, "protein_1", "within")
        assert value == "10-20+30-40"
        print(f"  PASS: Got value: {value}\n")

        # Test 16: lookup_table_value - by pdb column
        print("Test 16: lookup_table_value by pdb column")
        # The function should find it via the pdb column when id doesn't match exactly
        value = lookup_table_value(table, "protein_2", "beyond")
        assert value == "1-14+26-100"
        print(f"  PASS: Got value: {value}\n")

        # Test 17: iterate_table_values
        print("Test 17: iterate_table_values")
        results = list(iterate_table_values(
            table, ["protein_1", "protein_3"], "within"
        ))
        assert len(results) == 2
        assert results[0] == ("protein_1", "10-20+30-40")
        assert results[1] == ("protein_3", "5-10")
        print(f"  PASS: Iterated {len(results)} values\n")

        # Test 18: lookup_table_value - missing ID error
        print("Test 18: Error handling - missing ID in table")
        try:
            lookup_table_value(table, "nonexistent_protein", "within")
            print("  FAIL: Should have raised KeyError")
        except KeyError as e:
            assert "nonexistent_protein" in str(e)
            print(f"  PASS: Raised KeyError as expected\n")

        # Test 19: lookup_table_value - missing column error
        print("Test 19: Error handling - missing column")
        try:
            lookup_table_value(table, "protein_1", "nonexistent_column")
            print("  FAIL: Should have raised KeyError")
        except KeyError as e:
            assert "nonexistent_column" in str(e)
            print(f"  PASS: Raised KeyError as expected\n")

        # Test 20: load_table - invalid reference format
        print("Test 20: Error handling - invalid DATASHEET_REFERENCE")
        try:
            load_table("DATASHEET_REFERENCE:only_path_no_column")
            print("  FAIL: Should have raised ValueError")
        except ValueError as e:
            assert "Invalid DATASHEET_REFERENCE" in str(e)
            print(f"  PASS: Raised ValueError as expected\n")

        # Test 21: Table caching works
        print("Test 21: Table caching")
        clear_table_cache()
        table1, _ = load_table(positions_table_path)
        table2, _ = load_table(positions_table_path)
        assert table1 is table2  # Same object due to caching
        print("  PASS: Tables are cached\n")

        print("=" * 50)
        print("All tests passed!")

    finally:
        # Cleanup
        shutil.rmtree(test_dir)
