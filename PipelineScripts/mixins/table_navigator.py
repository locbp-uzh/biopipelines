"""
Datasheet Navigator Mixin for unified datasheet access.

Provides elegant navigation across different datasheet formats
(DatasheetContainer, dict, list, legacy formats).
"""

from typing import Union, Optional, List, Any


class DatasheetNavigatorMixin:
    """
    Mixin providing elegant datasheet navigation across all formats.

    Eliminates complex nested if-else chains for accessing datasheets
    from various sources and formats.
    """

    def get_datasheet(self,
                      source: Any,
                      name: Optional[str] = None,
                      fallback_names: Optional[List[str]] = None) -> Any:
        """
        Universal datasheet getter - works with all formats.

        This method replaces 30+ lines of nested if-else logic for accessing
        datasheets from different formats with a single function call.

        Args:
            source: ToolOutput, StandardizedOutput, DatasheetContainer, dict, or list
            name: Preferred datasheet name (e.g., 'structures', 'sequences', 'filtered')
            fallback_names: List of fallback names to try if primary name not found
                           Default: ['main', 'structures', 'sequences', 'compounds']

        Returns:
            DatasheetInfo object or path string

        Raises:
            ValueError: If datasheet cannot be found

        Examples:
            >>> # Get 'structures' datasheet, falling back to 'main'
            >>> ds = self.get_datasheet(tool_output, 'structures')

            >>> # Get any available datasheet with custom fallbacks
            >>> ds = self.get_datasheet(
            ...     tool_output.datasheets,
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

        # Handle DatasheetContainer (modern format with _datasheets attribute)
        if hasattr(source, '_datasheets'):
            for n in names_to_try:
                if n in source._datasheets:
                    return source._datasheets[n]
            # If no match found, raise with helpful message
            available = list(source._datasheets.keys())
            raise ValueError(
                f"Datasheet '{name}' not found in DatasheetContainer. "
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
                f"Datasheet '{name}' not found in dict. "
                f"Available: {available}"
            )

        # Handle ToolOutput/StandardizedOutput - recurse on their datasheets
        if hasattr(source, 'datasheets'):
            return self.get_datasheet(source.datasheets, name, fallback_names)

        # Handle legacy list format (first item is always "main")
        if isinstance(source, list) and source:
            if name in ['main', None] or name in fallback_names:
                return source[0]
            raise ValueError(
                f"Legacy list format only supports 'main' datasheet, "
                f"requested: '{name}'"
            )

        # Unsupported type
        raise ValueError(
            f"Cannot access datasheet from type: {type(source)}. "
            f"Expected DatasheetContainer, ToolOutput, StandardizedOutput, dict, or list."
        )

    def get_datasheet_path(self,
                           source: Any,
                           name: Optional[str] = None,
                           fallback_names: Optional[List[str]] = None) -> str:
        """
        Get datasheet path string (convenience method).

        Automatically extracts the path from DatasheetInfo objects.

        Args:
            source: Same as get_datasheet()
            name: Same as get_datasheet()
            fallback_names: Same as get_datasheet()

        Returns:
            Path string to datasheet CSV file

        Example:
            >>> path = self.get_datasheet_path(tool_output, 'structures')
            >>> # Returns: '/path/to/structures.csv'
        """
        datasheet = self.get_datasheet(source, name, fallback_names)

        # Extract path from DatasheetInfo object
        if hasattr(datasheet, 'path'):
            return datasheet.path

        # Handle dict format
        if isinstance(datasheet, dict) and 'path' in datasheet:
            return datasheet['path']

        # Assume it's already a path string
        if isinstance(datasheet, str):
            return datasheet

        raise ValueError(
            f"Cannot extract path from datasheet: {type(datasheet)}"
        )

    def get_all_datasheets(self, source: Any) -> dict:
        """
        Get all datasheets as a dictionary.

        Args:
            source: DatasheetContainer, ToolOutput, StandardizedOutput, dict, or list

        Returns:
            Dictionary mapping datasheet names to DatasheetInfo objects or paths

        Example:
            >>> all_ds = self.get_all_datasheets(tool_output)
            >>> for name, ds in all_ds.items():
            ...     print(f"{name}: {ds.path}")
        """
        # Handle DatasheetContainer
        if hasattr(source, '_datasheets'):
            return dict(source._datasheets)

        # Handle dict format
        if isinstance(source, dict):
            return source.copy()

        # Handle ToolOutput/StandardizedOutput
        if hasattr(source, 'datasheets'):
            return self.get_all_datasheets(source.datasheets)

        # Handle legacy list format
        if isinstance(source, list):
            return {'main': source[0]} if source else {}

        raise ValueError(
            f"Cannot get all datasheets from type: {type(source)}"
        )

    def datasheet_exists(self, source: Any, name: str) -> bool:
        """
        Check if a datasheet exists.

        Args:
            source: Any datasheet source
            name: Datasheet name to check

        Returns:
            True if datasheet exists, False otherwise

        Example:
            >>> if self.datasheet_exists(tool_output, 'structures'):
            ...     ds = self.get_datasheet(tool_output, 'structures')
        """
        try:
            self.get_datasheet(source, name, fallback_names=[])
            return True
        except ValueError:
            return False
