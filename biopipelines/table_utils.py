# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Table utilities for accessing tables from tool outputs.

Provides simple functions for navigating tables across different formats.
"""

from typing import Any, Optional, List, Dict


def get_table(source: Any, name: str) -> Any:
    """
    Get a table by name from a source.

    Args:
        source: Dict of tables, ToolOutput, or StandardizedOutput
        name: Table name to retrieve

    Returns:
        TableInfo object

    Raises:
        KeyError: If table not found
    """
    # Handle dict format (most common - direct from get_output_files)
    if isinstance(source, dict):
        if name in source:
            return source[name]
        raise KeyError(f"Table '{name}' not found. Available: {list(source.keys())}")

    # Handle ToolOutput/StandardizedOutput - access their tables attribute
    if hasattr(source, 'tables'):
        return get_table(source.tables, name)

    raise TypeError(f"Cannot get table from type: {type(source).__name__}")


def get_table_path(source: Any, name: str) -> str:
    """
    Get a table's file path by name.

    Args:
        source: Dict of tables, ToolOutput, or StandardizedOutput
        name: Table name to retrieve

    Returns:
        Path string to the table CSV file
    """
    table = get_table(source, name)

    # TableInfo object
    if hasattr(table, 'info'):
        return table.info.path

    # Dict format with 'path' key
    if isinstance(table, dict) and 'path' in table:
        return table['path']

    # Already a path string
    if isinstance(table, str):
        return table

    raise TypeError(f"Cannot extract path from table: {type(table).__name__}")


def list_tables(source: Any) -> List[str]:
    """
    List all available table names from a source.

    Args:
        source: Dict of tables, ToolOutput, or StandardizedOutput

    Returns:
        List of table names
    """
    if isinstance(source, dict):
        return list(source.keys())

    if hasattr(source, 'tables'):
        return list_tables(source.tables)

    raise TypeError(f"Cannot list tables from type: {type(source).__name__}")


def table_exists(source: Any, name: str) -> bool:
    """
    Check if a table exists in the source.

    Args:
        source: Dict of tables, ToolOutput, or StandardizedOutput
        name: Table name to check

    Returns:
        True if table exists
    """
    try:
        get_table(source, name)
        return True
    except (KeyError, TypeError):
        return False
