"""
Table Navigator Mixin for unified table access.

Provides elegant navigation across different table formats
(TableContainer, dict, list, legacy formats).
"""

from typing import Union, Optional, List, Any


class TableNavigatorMixin:
    """
    Mixin providing elegant table navigation across all formats.

    Eliminates complex nested if-else chains for accessing tables
    from various sources and formats.
    """

    def get_table(self,
                      source: Any,
                      name: Optional[str] = None,
                      fallback_names: Optional[List[str]] = None) -> Any:
        """
        Universal table getter - works with all formats.

        This method replaces 30+ lines of nested if-else logic for accessing
        tables from different formats with a single function call.

        Args:
            source: ToolOutput, StandardizedOutput, TableContainer, dict, or list
            name: Preferred table name (e.g., 'structures', 'sequences', 'filtered')
            fallback_names: List of fallback names to try if primary name not found
                           Default: ['main', 'structures', 'sequences', 'compounds']

        Returns:
            TableInfo object or path string

        Raises:
            ValueError: If table cannot be found

        Examples:
            >>> # Get 'structures' table, falling back to 'main'
            >>> ds = self.get_table(tool_output, 'structures')

            >>> # Get any available table with custom fallbacks
            >>> ds = self.get_table(
            ...     tool_output.tables,
            ...     name='filtered',
            ...     fallback_names=['merged', 'combined', 'main']
            ... )
        """
        # Default fallback names
        if fallback_names is None:
            fallback_names = ['main', 'structures', 'sequences', 'compounds']

        # Build list of names to try
        names_to_try = []
        if name:
            names_to_try.append(name)
        names_to_try.extend(fallback_names)

        # Remove duplicates while preserving order
        seen = set()
        names_to_try = [n for n in names_to_try if not (n in seen or seen.add(n))]

        # Handle TableContainer (modern format with _tables attribute)
        if hasattr(source, '_tables'):
            for n in names_to_try:
                if n in source._tables:
                    return source._tables[n]
            # If no match found, raise with helpful message
            available = list(source._tables.keys())
            raise ValueError(
                f"Table '{name}' not found in TableContainer. "
                f"Available: {available}"
            )

        # Handle dict format
        if isinstance(source, dict):
            for n in names_to_try:
                if n in source:
                    ds_value = source[n]
                    # Handle nested dict format {'path': '...', 'columns': [...]}
                    if isinstance(ds_value, dict) and 'path' in ds_value:
                        return ds_value
                    # Handle direct path string
                    return ds_value
            # If no match found, raise with helpful message
            available = list(source.keys())
            raise ValueError(
                f"Table '{name}' not found in dict. "
                f"Available: {available}"
            )

        # Handle ToolOutput/StandardizedOutput - recurse on their tables
        if hasattr(source, 'tables'):
            return self.get_table(source.tables, name, fallback_names)

        # Handle legacy list format (first item is always "main")
        if isinstance(source, list) and source:
            if name in ['main', None] or name in fallback_names:
                return source[0]
            raise ValueError(
                f"Legacy list format only supports 'main' table, "
                f"requested: '{name}'"
            )

        # Unsupported type
        raise ValueError(
            f"Cannot access table from type: {type(source)}. "
            f"Expected TableContainer, ToolOutput, StandardizedOutput, dict, or list."
        )

    def get_table_path(self,
                           source: Any,
                           name: Optional[str] = None,
                           fallback_names: Optional[List[str]] = None) -> str:
        """
        Get table path string (convenience method).

        Automatically extracts the path from TableInfo objects.

        Args:
            source: Same as get_table()
            name: Same as get_table()
            fallback_names: Same as get_table()

        Returns:
            Path string to table CSV file

        Example:
            >>> path = self.get_table_path(tool_output, 'structures')
            >>> # Returns: '/path/to/structures.csv'
        """
        table = self.get_table(source, name, fallback_names)

        # Extract path from TableInfo object
        if hasattr(table, 'path'):
            return table.path

        # Handle dict format
        if isinstance(table, dict) and 'path' in table:
            return table['path']

        # Assume it's already a path string
        if isinstance(table, str):
            return table

        raise ValueError(
            f"Cannot extract path from table: {type(table)}"
        )

    def get_all_tables(self, source: Any) -> dict:
        """
        Get all tables as a dictionary.

        Args:
            source: TableContainer, ToolOutput, StandardizedOutput, dict, or list

        Returns:
            Dictionary mapping table names to TableInfo objects or paths

        Example:
            >>> all_ds = self.get_all_tables(tool_output)
            >>> for name, ds in all_ds.items():
            ...     print(f"{name}: {ds.path}")
        """
        # Handle TableContainer
        if hasattr(source, '_tables'):
            return dict(source._tables)

        # Handle dict format
        if isinstance(source, dict):
            return source.copy()

        # Handle ToolOutput/StandardizedOutput
        if hasattr(source, 'tables'):
            return self.get_all_tables(source.tables)

        # Handle legacy list format
        if isinstance(source, list):
            return {'main': source[0]} if source else {}

        raise ValueError(
            f"Cannot get all tables from type: {type(source)}"
        )

    def table_exists(self, source: Any, name: str) -> bool:
        """
        Check if a table exists.

        Args:
            source: Any table source
            name: Table name to check

        Returns:
            True if table exists, False otherwise

        Example:
            >>> if self.table_exists(tool_output, 'structures'):
            ...     ds = self.get_table(tool_output, 'structures')
        """
        try:
            self.get_table(source, name, fallback_names=[])
            return True
        except ValueError:
            return False
