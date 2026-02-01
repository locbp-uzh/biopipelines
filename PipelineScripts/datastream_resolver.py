"""
DataStream Resolver: Convert any input type into a DataStream.

Provides universal input resolution for ToolOutput, StandardizedOutput,
file paths, lists, and other input formats.
"""

import os
from typing import Any, Optional, List, Dict

from .datastream import DataStream, create_map_table


class DataStreamResolver:
    """
    Resolves any input type into a DataStream.

    This class provides a single entry point for converting various input
    formats into a unified DataStream interface. The caller specifies the
    format - no auto-detection or fallbacks.

    Usage:
        resolver = DataStreamResolver(output_folder="/path/to/tool/output")

        # From ToolOutput - extracts the DataStream directly
        structures = resolver.resolve(tool_output.structures)

        # From file list - caller provides format
        structures = resolver.resolve(["/path/a.pdb", "/path/b.pdb"], format="pdb")

        # From StandardizedOutput attribute
        sequences = resolver.resolve(std_output.sequences)
    """

    def __init__(self, output_folder: str):
        """
        Initialize resolver.

        Args:
            output_folder: Folder where map_table files will be created
        """
        self.output_folder = output_folder

    def resolve(self, input_param: Any, format: str, name: str = "") -> DataStream:
        """
        Universal resolver - converts any input type to DataStream.

        Args:
            input_param: Input of any supported type:
                - DataStream (passed through)
                - List[str] (file paths)
                - str (single file path)
                - Dict with 'ids' and 'files' or 'values' keys
            format: Data format (e.g., "pdb", "fasta", "smiles") - required
            name: Name for the DataStream (e.g., "structures", "sequences")

        Returns:
            DataStream containing the resolved data

        Raises:
            ValueError: If input type is unsupported
        """
        # Already a DataStream - pass through
        if isinstance(input_param, DataStream):
            return input_param

        # None - error, no fallbacks
        if input_param is None:
            raise ValueError("Input cannot be None")

        # Handle list of file paths
        if isinstance(input_param, list):
            return self._resolve_file_list(input_param, format, name)

        # Handle single file path
        if isinstance(input_param, str):
            return self._resolve_string(input_param, format, name)

        # Handle dictionary
        if isinstance(input_param, dict):
            return self._resolve_dict(input_param, format, name)

        raise ValueError(
            f"Unsupported input type: {type(input_param).__name__}. "
            f"Expected DataStream, list, str, or dict."
        )

    def _resolve_file_list(self, files: List[str], format: str, name: str) -> DataStream:
        """Resolve list of file paths to DataStream."""
        if not files:
            raise ValueError("File list cannot be empty")

        # Extract IDs from filenames
        ids = [os.path.splitext(os.path.basename(f))[0] for f in files]

        # Create map_table
        map_table = self._create_map_table(ids, files=files)

        return DataStream(
            name=name,
            ids=ids,
            files=files,
            map_table=map_table,
            format=format
        )

    def _resolve_string(self, input_str: str, format: str, name: str) -> DataStream:
        """Resolve single file path to DataStream."""
        files = [input_str]
        ids = [os.path.splitext(os.path.basename(input_str))[0]]
        map_table = self._create_map_table(ids, files=files)

        return DataStream(
            name=name,
            ids=ids,
            files=files,
            map_table=map_table,
            format=format
        )

    def _resolve_dict(self, input_dict: Dict[str, Any], format: str, name: str) -> DataStream:
        """Resolve dictionary to DataStream."""
        if 'ids' not in input_dict:
            raise ValueError("Dictionary input must contain 'ids' key")

        ids = input_dict['ids']
        files = input_dict.get('files')
        values = input_dict.get('values')

        # Create appropriate map_table
        if values is not None and files is None:
            map_table = self._create_map_table(ids, values=values)
            files_list = []
        elif files is not None:
            map_table = self._create_map_table(ids, files=files)
            files_list = files
        else:
            raise ValueError("Dictionary input must contain either 'files' or 'values' key")

        return DataStream(
            name=name,
            ids=ids,
            files=files_list,
            map_table=map_table,
            format=format,
            metadata=input_dict.get('metadata', {})
        )

    def _create_map_table(
        self,
        ids: List[str],
        files: Optional[List[str]] = None,
        values: Optional[List[str]] = None
    ) -> str:
        """Create map_table CSV file."""
        output_path = os.path.join(self.output_folder, "datastream_map.csv")
        return create_map_table(output_path, ids, files=files, values=values)


def resolve_to_datastream(
    input_param: Any,
    output_folder: str,
    format: str,
    name: str = ""
) -> DataStream:
    """
    Convenience function to resolve any input to DataStream.

    Args:
        input_param: Any supported input type
        output_folder: Where to create map_table
        format: Data format (required)
        name: Name for the DataStream

    Returns:
        DataStream containing the resolved data
    """
    resolver = DataStreamResolver(output_folder=output_folder)
    return resolver.resolve(input_param, format=format, name=name)
