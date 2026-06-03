#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
BioPipelines I/O Utilities for Execution Time.

This module provides functions for pipe_*.py scripts to read DataStreams and tables
at execution time. It handles wildcard resolution, ID mapping, and table lookups.

Functions:
    # DataStream utilities
    load_datastream: Load DataStream from JSON file or dict
    iterate_files: Iterate (id, resolved_file_path) for file-based streams
    iterate_values: Iterate (id, value_dict) for value-based streams
    resolve_file: Resolve single file for an ID (handles wildcards)
    get_value: Get single value from map_table for an ID
    get_all_values: Get all values from map_table for an ID

    # Table reference utilities (for per-structure lookups from external tables)
    load_table: Load a table from path or TABLE_REFERENCE string
    lookup_table_value: Look up a value from a table for a given ID
    iterate_table_values: Iterate over (id, value) pairs from a table

Example usage:
    from biopipelines.biopipelines_io import load_datastream, iterate_files, iterate_values

    # Load from JSON or dict
    structures_ds = load_datastream("/path/to/structures.json")

    # Iterate over structure files (wildcards resolved automatically)
    for struct_id, struct_file in iterate_files(structures_ds):
        print(f"Processing {struct_id}: {struct_file}")

    # Iterate over compound values (SMILES from map_table)
    for comp_id, values in iterate_values(compounds_ds, columns=['smiles', 'name']):
        print(f"Compound {comp_id}: {values['smiles']}")

    # Table reference for per-structure data (e.g., fixed positions from DistanceSelector)
    from biopipelines.biopipelines_io import load_table, lookup_table_value, iterate_table_values

    # Load table (supports TABLE_REFERENCE:path:column format)
    table, column = load_table("TABLE_REFERENCE:/path/to/positions.csv:within")

    # Look up value for a specific structure
    fixed_positions = lookup_table_value(table, "protein_1", column)

    # Or iterate over all structures
    for struct_id, positions in iterate_table_values(table, structure_ids, column):
        print(f"{struct_id}: {positions}")
"""

import glob
import json
import os
import sys
from typing import Any, Dict, Iterator, List, Optional, Tuple, Union

import pandas as pd

try:
    from . import id_patterns as _id_patterns
except ImportError:
    try:
        import id_patterns as _id_patterns
    except ImportError:
        _id_patterns = None

try:
    from .datastream import DataStream
except ImportError:
    from datastream import DataStream


def load_datastream(source: Union[str, Dict[str, Any]]) -> DataStream:
    """
    Load a DataStream from a JSON file or dictionary.

    Args:
        source: Either a path to a JSON file or a dictionary with DataStream fields

    Returns:
        DataStream instance

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

    # files defaults to empty list. May also be a single str (shared-file form).
    files = data.get('files', [])

    # Validate files/ids relationship (skip for patterns, templates, and
    # shared-file form where one path covers every id).
    is_shared = isinstance(files, str) and bool(files)
    has_patterns = _id_patterns and any(
        _id_patterns.contains_pattern(s) or '[' in s for s in data['ids']
    )
    has_template = (
        isinstance(files, list) and len(files) == 1 and '<id>' in files[0]
    )
    if (not is_shared and not has_patterns and not has_template
            and isinstance(files, list)
            and len(files) > 1 and len(files) != len(data['ids'])):
        raise ValueError(
            f"Length mismatch: {len(data['ids'])} ids but {len(files)} files. "
            f"Use empty files list or single file for table-based data, "
            f"or provide one file per id."
        )

    return DataStream(
        name=data['name'],
        ids=data['ids'],
        files=files,
        map_table=data.get('map_table', ''),
        format=data.get('format', 'pdb'),
        metadata=data.get('metadata', {}),
        _runtime_mode=True
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


def iterate_files(ds: DataStream) -> Iterator[Tuple[str, str]]:
    """
    Iterate over (id, resolved_file_path) pairs for file-based streams.

    Handles three cases:
    1. File paths contain '*' glob: Expand globs, match to IDs
    2. len(files) == len(ids): Direct zip iteration
    3. len(files) == 1: Single file for all IDs (bundle case)

    Missing files are SKIPPED, not fatal: when a wildcard expands to no file for
    an id, that id is skipped with a warning and iteration continues. A declared
    id can legitimately lack a file — e.g. a filtered Pool/Panda still declares
    every original id in its stream but only carries the rows that survived a
    downstream filter. The dropped ids are tracked in the upstream `missing`
    manifest and propagated forward; whether that manifest accounts for every
    absent file (expected) or some are genuinely missing (a real failure) is
    decided by pipe_check_completion, not here. So this iterator never raises on
    a missing file — it just yields the ids whose files are actually present.

    Args:
        ds: DataStream instance

    Yields:
        Tuple of (item_id, file_path) for every id whose file is present.

    Raises:
        ValueError: If files list is empty and this is a file-based stream

    Example:
        for struct_id, struct_file in iterate_files(structures_ds):
            process_structure(struct_id, struct_file)
    """
    # Use expanded IDs and files for iteration
    ids = ds.ids_expanded
    files = ds.files_expanded

    if not files and not ds.files:
        raise ValueError(f"DataStream '{ds.name}' has no files configured")

    def _resolve_glob(item_id: str, pattern: str):
        """Expand a wildcard for one id; skip (warn + return None) on no match."""
        expanded = glob.glob(pattern)
        if not expanded:
            print(f"  Warning: No files found matching pattern '{pattern}' for ID "
                  f"'{item_id}' — skipping (id absent from this stream's files)")
            return None
        return _find_best_match(item_id, expanded)

    # If files_expanded produced results, use them
    if files and len(files) == len(ids):
        # Detect wildcards
        has_wildcards = any('*' in f for f in files)
        if has_wildcards:
            for item_id, file_pattern in zip(ids, files):
                match = _resolve_glob(item_id, file_pattern)
                if match is not None:
                    yield (item_id, match)
        else:
            # Per-id concrete paths (e.g. the '<id>.pdb' template expands to one
            # distinct file per id). Skip an id whose file is absent — a filtered
            # Pool/Panda declares every original id but only carries the surviving
            # files. (Distinct from the shared/bundle branches below, where a single
            # real file legitimately covers every id and must not be per-id-skipped.)
            for item_id, file_path in zip(ids, files):
                if '*' not in file_path and not os.path.exists(file_path):
                    print(f"  Warning: File not found for ID '{item_id}': "
                          f"'{file_path}' — skipping (id absent from this stream's files)")
                    continue
                yield (item_id, file_path)

    elif ds.is_shared_file:
        # Shared-file stream: one artifact path covers every id.
        single = ds.files
        if '*' in single:
            for item_id in ids:
                match = _resolve_glob(item_id, single)
                if match is not None:
                    yield (item_id, match)
        else:
            for item_id in ids:
                yield (item_id, single)

    elif isinstance(ds.files, list) and len(ds.files) == 1:
        # Single file/pattern for all IDs
        single = ds.files[0]
        if '*' in single:
            for item_id in ids:
                match = _resolve_glob(item_id, single)
                if match is not None:
                    yield (item_id, match)
        else:
            for item_id in ids:
                yield (item_id, single)

    else:
        raise ValueError(
            f"Cannot iterate: {len(ids)} ids but {len(files)} files"
        )


def iterate_values(
    ds: DataStream,
    columns: Optional[List[str]] = None
) -> Iterator[Tuple[str, Dict[str, Any]]]:
    """
    Iterate over (id, value_dict) pairs for value-based streams.

    Reads values from the map_table CSV and yields dictionaries with
    requested columns for each ID.

    Args:
        ds: DataStream instance
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

    for item_id in ds.ids_expanded:
        if item_id not in id_to_row:
            raise KeyError(f"ID '{item_id}' not found in map_table")

        row = id_to_row[item_id]
        values = {col: row[col] for col in include_cols}
        yield (item_id, values)


def resolve_file(ds: DataStream, item_id: str) -> str:
    """
    Resolve a single file path for an ID.

    Handles wildcard patterns by expanding and matching.

    Args:
        ds: DataStream instance
        item_id: The item identifier

    Returns:
        Resolved file path

    Raises:
        ValueError: If ID not found in DataStream or files not configured
        FileNotFoundError: If wildcard expansion yields no files

    Example:
        file_path = resolve_file(structures_ds, "rank0001")
    """
    ids = ds.ids_expanded
    files = ds.files_expanded

    if not files and not ds.files:
        raise ValueError(f"DataStream '{ds.name}' has no files configured")

    if item_id not in ids:
        raise ValueError(f"ID '{item_id}' not found in DataStream '{ds.name}'")

    idx = ids.index(item_id)

    if files and len(files) == len(ids):
        file_path_for_id = files[idx]
        if '*' in file_path_for_id:
            expanded = glob.glob(file_path_for_id)
            if not expanded:
                raise FileNotFoundError(
                    f"No files found matching pattern '{file_path_for_id}' for ID '{item_id}'"
                )
            return _find_best_match(item_id, expanded)
        return file_path_for_id

    elif ds.is_shared_file:
        single = ds.files
        if '*' in single:
            expanded = glob.glob(single)
            if not expanded:
                raise FileNotFoundError(
                    f"No files found matching pattern '{single}' for ID '{item_id}'"
                )
            return _find_best_match(item_id, expanded)
        return single

    elif isinstance(ds.files, list) and len(ds.files) == 1:
        single = ds.files[0]
        if '*' in single:
            expanded = glob.glob(single)
            if not expanded:
                raise FileNotFoundError(
                    f"No files found matching pattern '{single}' for ID '{item_id}'"
                )
            return _find_best_match(item_id, expanded)
        return single

    else:
        raise ValueError(
            f"Cannot resolve file: {len(ids)} ids but {len(files)} files"
        )


def get_value(ds: DataStream, item_id: str, column: str = "value") -> Any:
    """
    Get a single value from the map_table for an ID.

    Args:
        ds: DataStream instance
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


def get_all_values(ds: DataStream, item_id: str) -> Dict[str, Any]:
    """
    Get all column values from the map_table for an ID.

    Args:
        ds: DataStream instance
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
) -> Tuple[pd.DataFrame, Optional[str]]:
    """
    Load a table from a path or TABLE_REFERENCE string.

    Supports two formats:
    1. Direct path: "/path/to/table.csv" - returns (DataFrame, None)
    2. Reference: "TABLE_REFERENCE:/path/to/table.csv:column" - returns (DataFrame, column)

    Args:
        reference: Either a file path or TABLE_REFERENCE:path:column string

    Returns:
        Tuple of (DataFrame, column_name or None)

    Raises:
        FileNotFoundError: If table file doesn't exist
        ValueError: If reference format is invalid

    Example:
        # Direct path
        table, _ = load_table("/path/to/positions.csv")

        # Reference with column
        table, column = load_table("TABLE_REFERENCE:/path/to/positions.csv:within")
        value = lookup_table_value(table, "protein_1", column)
    """
    column_name = None

    if reference.startswith("TABLE_REFERENCE:"):
        prefix = "TABLE_REFERENCE:"
        # Parse reference: PREFIX:path:column
        # Handle Windows paths with drive letters (e.g., C:\path\to\file.csv)
        # The column is always the LAST colon-separated part
        remainder = reference[len(prefix):]

        # Find the last colon that's not part of a Windows drive letter
        # Windows drive letters are single letter followed by colon at position 1
        last_colon = remainder.rfind(":")

        # Check if this colon is part of a Windows path (e.g., "C:")
        # A Windows drive letter would be at position 1 (after single letter)
        if last_colon == 1 and remainder[0].isalpha():
            # Only the drive letter colon exists, no column specified
            raise ValueError(
                f"Invalid TABLE_REFERENCE format: '{reference}'. "
                f"Expected 'TABLE_REFERENCE:path:column'"
            )

        if last_colon == -1:
            raise ValueError(
                f"Invalid TABLE_REFERENCE format: '{reference}'. "
                f"Expected 'TABLE_REFERENCE:path:column'"
            )

        table_path = remainder[:last_colon]
        column_name = remainder[last_colon + 1:]

        if not table_path or not column_name:
            raise ValueError(
                f"Invalid TABLE_REFERENCE format: '{reference}'. "
                f"Expected 'TABLE_REFERENCE:path:column'"
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
    id_column: str = "id",
    pdb_column: str = "pdb"
) -> Any:
    """
    Look up a value from a table for a given ID, using the framework's standard
    id matching (``get_mapped_ids``): exact, child-parent (``design+Cy7RR`` ->
    ``design``), and provenance.

    Args:
        table: DataFrame to search
        item_id: The item identifier (structure id).
        column: Column name to retrieve value from
        id_column: Column holding the ids to match against (default: "id";
                   falls back to ``pdb_column`` if no id column exists).

    Returns:
        Value from the specified column

    Raises:
        KeyError: If the id can't be matched, or the column doesn't exist.
    """
    if column not in table.columns:
        raise KeyError(
            f"Column '{column}' not found in table. "
            f"Available columns: {list(table.columns)}"
        )

    key_col = id_column if id_column in table.columns else (
        pdb_column if pdb_column in table.columns else None)
    if key_col is None:
        raise KeyError(
            f"Table has neither '{id_column}' nor '{pdb_column}' column to match "
            f"'{item_id}' against. Columns: {list(table.columns)}"
        )

    # pdb-column ids carry a .pdb extension; strip it so they match design ids.
    table_ids = [str(x) for x in table[key_col].tolist()]
    stripped = [i[:-4] if i.endswith(".pdb") else i for i in table_ids]

    from biopipelines.id_map_utils import get_mapped_ids
    matched = get_mapped_ids([item_id], stripped, unique=True).get(item_id)
    if matched is None:
        raise KeyError(
            f"ID '{item_id}' not found in table (col '{key_col}'). "
            f"Available: {table_ids[:10]}{'...' if len(table_ids) > 10 else ''}"
        )
    # Map the stripped match back to the actual table row.
    row_idx = stripped.index(matched)
    return table.iloc[row_idx][column]


def iterate_table_values(
    table: pd.DataFrame,
    item_ids: List[str],
    column: str,
    id_column: str = "id",
    pdb_column: str = "pdb"
) -> Iterator[Tuple[str, Any]]:
    """
    Iterate over (id, value) pairs from a table for a list of IDs.

    Args:
        table: DataFrame to search
        item_ids: List of item identifiers to look up
        column: Column name to retrieve values from
        id_column: Column name for ID lookups (default: "id")
        pdb_column: Column name for PDB filename lookups (default: "pdb")

    Yields:
        Tuple of (item_id, value)

    Raises:
        KeyError: If any ID not found in table or column doesn't exist

    Example:
        table, column = load_table("TABLE_REFERENCE:/path/to/positions.csv:within")
        for struct_id, positions in iterate_table_values(table, structure_ids, column):
            print(f"{struct_id}: {positions}")
    """
    for item_id in item_ids:
        value = lookup_table_value(table, item_id, column, id_column, pdb_column)
        yield (item_id, value)


def clear_table_cache() -> None:
    """
    Clear the table cache.

    Call this if you need to reload tables that may have changed on disk.
    """
    global _table_cache
    _table_cache = {}


MISSING_COLUMNS = ["id", "removed_by", "kind", "cause"]


def step_id_from_table_path(table_path: str) -> str:
    """The step identifier (``<order>_<Tool>``) a tool stamps into ``removed_by``.

    Derived from a tool-owned table path like ``.../005_XTB/tables/missing.csv``
    -> ``005_XTB`` (the tool's output-folder basename, matching what
    ``pipe_check_completion`` compares against). This lets a same-tool chain
    distinguish an upstream step's failure rows from the current step's own:
    each step writes its OWN step id, so a later step excuses the earlier one
    instead of mistaking it for a local failure.

    Falls back to "" when the path has no recognizable tool folder.
    """
    if not table_path:
        return ""
    # <output_folder>/tables/<name>.csv  ->  basename(<output_folder>)
    tool_folder = os.path.dirname(os.path.dirname(os.path.abspath(table_path)))
    return os.path.basename(tool_folder)


def read_upstream_missing(paths) -> List[Dict[str, Any]]:
    """Merge every upstream ``missing`` manifest into one list of row dicts.

    A consumer may have several id-bearing input axes (e.g. proteins AND
    ligands), each carrying its own ``missing.csv``; passing only the first
    would silently drop the other axis's excused ids. Accepts a single path, a
    list of paths, or None; reads each that exists, concatenates the rows, and
    de-dupes by ``id`` (first kept). Rows stay in the upstream (input-axis) id
    space; mapping them to this tool's output ids is done centrally by
    ``pipe_check_completion`` when excusing files. Used by pipe scripts that
    prepend upstream rows to their own local-failure rows before writing their
    ``missing.csv``.
    """
    if paths is None:
        return []
    if isinstance(paths, str):
        paths = [paths]
    seen = set()
    rows: List[Dict[str, Any]] = []
    for path in paths:
        if not path or not os.path.exists(path):
            continue
        try:
            df = pd.read_csv(path)
        except Exception as e:  # noqa: BLE001 — a malformed upstream file shouldn't abort scoring
            print(f"Warning: could not read upstream missing.csv {path}: {e}", file=sys.stderr)
            continue
        for rec in df.to_dict("records"):
            rid = str(rec.get("id", "")).strip()
            if not rid or rid in seen:
                continue
            seen.add(rid)
            rows.append(rec)
    return rows


# =============================================================================
# Provenance-based ID Resolution
# =============================================================================
# When tools like Panda auto-rename IDs (e.g., Panda_1, Panda_2), the original
# source IDs are preserved in map_table provenance columns (e.g., structures.id).
# These functions trace back through provenance to match renamed IDs to their
# source IDs in external tables.


# Per-process cache of {id -> [provenance_id, ...]} per map_table CSV.
# Mirrors id_map_utils._provenance_lookup_cache; rebuilt lazily once per
# CSV the first time resolve_id_by_provenance hits that path. Keeps
# cross-source lookups O(1) instead of O(rows) per call.
_provenance_index_cache: Dict[str, Dict[str, List[str]]] = {}


def _build_provenance_index(df: pd.DataFrame) -> Dict[str, List[str]]:
    """Build {id -> [provenance values from <stream>.id columns]} from a
    map_table DataFrame. NaN cells are dropped via explicit string
    coercion."""
    index: Dict[str, List[str]] = {}
    if "id" not in df.columns:
        return index
    prov_cols = [c for c in df.columns if c.endswith(".id") and c != "id"]
    if not prov_cols:
        return index
    ids = [str(v) for v in df["id"].tolist()]
    cols = {c: [str(v) for v in df[c].tolist()] for c in prov_cols}
    for i, sid in enumerate(ids):
        provs: List[str] = []
        for c in prov_cols:
            val = cols[c][i]
            if val and val != "nan":
                provs.append(val)
        if provs:
            index[sid] = provs
    return index


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

        # Build / reuse the per-id provenance index for this CSV.
        if map_table_path not in _provenance_index_cache:
            if map_table_path in _table_cache:
                df = _table_cache[map_table_path]
            else:
                try:
                    df = pd.read_csv(map_table_path)
                    _table_cache[map_table_path] = df
                except Exception:
                    continue
            _provenance_index_cache[map_table_path] = _build_provenance_index(df)

        provs = _provenance_index_cache[map_table_path].get(item_id)
        if not provs:
            continue

        for source_id in provs:
            if source_id in target_ids:
                return source_id

    return None


# =============================================================================
# Resolve: Bash snippet generators for runtime resolution
# =============================================================================
# Used at config time by tool classes to generate bash expressions that defer
# file path and table column lookups to execution time.
#
# Usage in tool generate_script methods:
#
#     # Resolve first ID at runtime (safe with lazy streams):
#     f'INPUT_PDB_ID={Resolve.stream_ids(self.pdb_ds_json, index=0)}'
#     f'INPUT_PDB={Resolve.stream_item(self.pdb_ds_json, "$INPUT_PDB_ID")}'
#
#     # In a loop (Resolve.stream_ids without index returns all IDs):
#     f'for sid in {Resolve.stream_ids(self.pdb_ds_json)}; do ...'
#
# The bash function resolve_stream_item is sourced automatically by
# BaseConfig.activate_environment(). No per-tool setup needed.


class Resolve:
    """Static bash snippet generators for runtime resolution."""

    @staticmethod
    def stream_item(ds_json: str, item_id: str, column: Optional[str] = None) -> str:
        """
        Bash expression to resolve one item from a DataStream JSON at runtime.

        For a file-based stream this resolves the item's file path. For a
        value-based stream (files=[]), pass ``column`` to resolve that
        map_table column's value for the id instead.

        Args:
            ds_json: Path to the serialized DataStream JSON file
            item_id: ID of the item to resolve
            column: When set, resolve this map_table column's value (value-based
                stream) rather than a file path.

        Returns:
            Bash subshell expression, e.g. $(resolve_stream_item "..." "...")
        """
        if column is not None:
            return f'$(resolve_stream_item "{ds_json}" "{item_id}" "{column}")'
        return f'$(resolve_stream_item "{ds_json}" "{item_id}")'

    @staticmethod
    def stream_ids(ds_json: str, index: Optional[int] = None) -> str:
        """Bash expression that prints expanded IDs from a DataStream JSON, one per line.

        Handles lazy patterns by matching against map_table at runtime.

        Args:
            ds_json: Path to the serialized DataStream JSON file
            index: If provided, return only the ID at this position (0-based).
        """
        script = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                              "..", "pipe_scripts", "resolve_stream_ids.py")
        if index is not None:
            n = index + 1  # head -n is 1-based
            return f'$(python "{script}" "{ds_json}" | head -n{n} | tail -n1)'
        return f'$(python "{script}" "{ds_json}")'

    @staticmethod
    def table_column(reference, item_id: str, env_name: str = "biopipelines") -> str:
        """
        Bash expression to resolve a table value inline (slow — spawns Python).

        Uses conda/mamba run to execute in the biopipelines environment (on
        Colab / pip mode it falls back to a plain ``python`` call, since no
        ``biopipelines`` env exists there).

        Prefer this only for a one-off, single-id lookup. For a per-id value
        across a loop (e.g. one contig string per input PDB), do NOT call this
        once per id — that spawns Python N times. Instead resolve all ids in a
        single pass with a pipe script that imports ``biopipelines_io``
        (``load_table`` + ``lookup_table_value``), writes a small ``{id: value}``
        JSON, and have the bash loop read each id's value from that JSON. That
        pattern runs under the activated tool env (Colab-safe, no second env)
        and is what the RFdiffusion family uses for per-PDB contigs
        (``pipe_rfdiffusion_contigs.py`` + ``resolve_rfdiffusion_contigs.py``).

        Args:
            reference: TableReference object or TABLE_REFERENCE:path:column string
            item_id: ID of the item to look up
            env_name: Conda environment name (default: "biopipelines")

        Returns:
            Bash subshell expression calling resolve_table_column.py
        """
        from .config_manager import ConfigManager
        cfg = ConfigManager()
        mgr = cfg.get_env_manager()
        script_dir = os.path.dirname(os.path.abspath(__file__))
        script = os.path.join(script_dir, "..", "pipe_scripts", "resolve_table_column.py")
        # On Colab the biopipelines package lives in the base Python (no
        # `biopipelines` conda env exists), so wrapping in `<mgr> run -n
        # biopipelines` would fail. Pip mode is the same case. In both, emit a
        # plain `python` call — the activated tool env already has biopipelines
        # importable via the script's sys.path bootstrap.
        if mgr == "pip" or cfg.get_scheduler() == "colab":
            return f'$(python "{script}" "{reference}" "{item_id}")'
        return f'$({mgr} run -n {env_name} --no-banner python "{script}" "{reference}" "{item_id}")'


class TableReference:
    """Typed reference to a specific column in a CSV table.

    Created by TableInfo attribute access: table.column -> TableReference(path, "column").
    str() produces TABLE_REFERENCE:path:column format consumed by pipe scripts.
    """
    PREFIX = "TABLE_REFERENCE"

    def __init__(self, path: str, column: str):
        self.path = path
        self.column = column

    def __str__(self) -> str:
        return f"{self.PREFIX}:{self.path}:{self.column}"

    def __repr__(self) -> str:
        return f"TableReference({self.path!r}, {self.column!r})"

    def __bool__(self) -> bool:
        return bool(self.path and self.column)

    @classmethod
    def from_string(cls, s: str) -> 'TableReference':
        """Parse a TABLE_REFERENCE:path:column string."""
        if not s.startswith(f"{cls.PREFIX}:"):
            raise ValueError(f"Not a TABLE_REFERENCE string: {s}")
        rest = s[len(cls.PREFIX) + 1:]
        idx = rest.rfind(":")
        if idx <= 0:
            raise ValueError(f"Invalid TABLE_REFERENCE format: {s}")
        return cls(rest[:idx], rest[idx + 1:])

    def resolve(self, item_id: str):
        """Resolve value directly in Python (requires pandas)."""
        table, column = load_table(str(self))
        return lookup_table_value(table, item_id, column)

