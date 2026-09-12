# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Table utilities for accessing tables from tool outputs.

These five functions — :func:`get_table`, :func:`get_table_path`,
:func:`get_indexed_table`, :func:`list_tables`, :func:`table_exists` — are the
supported way to read a finished run's tables back into Python. Every one of
them accepts the containers a caller actually holds: a tool's
``StandardizedOutput``, the ``TableContainer`` behind its ``.tables``, a
``ToolOutput``, or the plain ``{name: TableInfo}`` dict that
``Load(...).get_output_files()["tables"]`` returns. Anything else raises a
``TypeError`` naming what was passed and what to pass instead.
"""

from typing import Any, Dict, List

try:
    from .base_config import (
        IndexedTableContainer,
        StandardizedOutput,
        TableContainer,
        TableInfo,
        ToolOutput,
    )
except ImportError:
    import os
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import (
        IndexedTableContainer,
        StandardizedOutput,
        TableContainer,
        TableInfo,
        ToolOutput,
    )


_ACCEPTED = ("a tool output (StandardizedOutput), a ToolOutput, a tables "
             "container (TableContainer), or a plain dict of table name -> TableInfo "
             '(e.g. Load(folder).get_output_files()["tables"])')


def _tables_mapping(source: Any) -> Dict[str, Any]:
    """The ``{name: table}`` mapping behind any accepted container.

    Every entry point routes through here so one container is never accepted by
    one function and rejected by another.
    """
    if isinstance(source, dict):
        return source
    if isinstance(source, TableContainer):
        return source._tables
    if isinstance(source, (StandardizedOutput, ToolOutput)):
        return _tables_mapping(source.tables)
    raise TypeError(
        f"Cannot read tables from {type(source).__name__}. Pass {_ACCEPTED}."
    )


def get_table(source: Any, name: str) -> Any:
    """
    Get a table by name from a source.

    Args:
        source: StandardizedOutput, ToolOutput, TableContainer, or dict of tables
        name: Table name to retrieve

    Returns:
        TableInfo object, or IndexedTableContainer for a per-ID table collection

    Raises:
        KeyError: If the table is not present
        TypeError: If source is not a supported container
    """
    tables = _tables_mapping(source)
    if name in tables:
        return tables[name]
    raise KeyError(f"Table '{name}' not found. Available: {list(tables.keys())}")


def get_table_path(source: Any, name: str) -> str:
    """
    Get a table's file path by name.

    Args:
        source: StandardizedOutput, ToolOutput, TableContainer, or dict of tables
        name: Table name to retrieve

    Returns:
        Path string to the table CSV file

    Raises:
        KeyError: If the table is not present
        TypeError: If source is unsupported, the table is a per-ID collection,
            or no path can be extracted from it
    """
    table = get_table(source, name)

    # A per-ID collection has one path per entry, so there is no single answer.
    if isinstance(table, IndexedTableContainer):
        raise TypeError(
            f"Table '{name}' is an IndexedTableContainer with {len(table)} entries. "
            f"Use get_indexed_table(source, '{name}', entry_id) to get one entry's path."
        )

    if isinstance(table, TableInfo):
        return table.info.path

    if isinstance(table, dict) and 'path' in table:
        return table['path']

    if isinstance(table, str):
        return table

    raise TypeError(
        f"Cannot extract a path from table '{name}' of type {type(table).__name__}"
    )


def list_tables(source: Any) -> List[str]:
    """
    List all available table names from a source.

    Args:
        source: StandardizedOutput, ToolOutput, TableContainer, or dict of tables

    Returns:
        List of table names

    Raises:
        TypeError: If source is not a supported container
    """
    return list(_tables_mapping(source).keys())


def get_indexed_table(source: Any, name: str, entry_id: str) -> Any:
    """
    Get a specific entry from an IndexedTableContainer.

    Args:
        source: StandardizedOutput, ToolOutput, TableContainer, or dict of tables
        name: Name of the indexed table collection
        entry_id: ID of the specific entry

    Returns:
        TableInfo for the given ID

    Raises:
        TypeError: If source is unsupported or the table is not an IndexedTableContainer
        KeyError: If the table or the entry_id is not found
    """
    table = get_table(source, name)
    if isinstance(table, IndexedTableContainer):
        return table[entry_id]
    raise TypeError(
        f"Table '{name}' is a {type(table).__name__}, not an IndexedTableContainer. "
        f"Use get_table_path(source, '{name}') for its single path."
    )


def table_exists(source: Any, name: str) -> bool:
    """
    Check if a table exists in the source.

    An unsupported source raises rather than answering False — a wrong answer
    is worse than an error, since False here reads as "the run has no such table".

    Args:
        source: StandardizedOutput, ToolOutput, TableContainer, or dict of tables
        name: Table name to check

    Returns:
        True if the table is present

    Raises:
        TypeError: If source is not a supported container
    """
    return name in _tables_mapping(source)
