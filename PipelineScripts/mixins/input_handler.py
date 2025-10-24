"""
Input Handler Mixin for standardized input processing across all tools.

Provides universal input resolution for ToolOutput, StandardizedOutput,
strings, lists, and dictionaries.
"""

import os
from typing import List, Optional, Union, Any


class ResolvedInput:
    """Container for resolved input data with consistent interface."""

    def __init__(self, files: List[str], ids: List[str],
                 datasheets: Any, source: Any):
        """
        Initialize resolved input.

        Args:
            files: List of file paths (structures, sequences, or compounds)
            ids: List of IDs corresponding to files
            datasheets: Datasheets associated with input (if any)
            source: Original input source
        """
        self.files = files
        self.ids = ids
        self.datasheets = datasheets
        self.source = source

        # Import here to avoid circular imports
        try:
            from ..base_config import ToolOutput
            self.is_tool_output = isinstance(source, ToolOutput)
        except ImportError:
            self.is_tool_output = False

    def __repr__(self):
        return (f"ResolvedInput(files={len(self.files)}, "
                f"ids={len(self.ids)}, "
                f"is_tool_output={self.is_tool_output})")


class InputHandlerMixin:
    """
    Mixin providing standardized input handling for all data types.

    Eliminates repetitive input type checking and extraction logic
    across all tools.
    """

    @staticmethod
    def resolve_input(input_param: Any,
                      expected_type: str = 'structures') -> ResolvedInput:
        """
        Universal input resolver - handles all input types consistently.

        This method replaces 50+ lines of repetitive if-elif-else chains
        in every tool with a single function call.

        Args:
            input_param: Any input type (ToolOutput, StandardizedOutput, str, list, dict)
            expected_type: What type of data to extract:
                - 'structures': PDB/CIF files
                - 'sequences': FASTA/CSV files
                - 'compounds': SDF/MOL files

        Returns:
            ResolvedInput object with .files, .ids, .datasheets attributes

        Examples:
            >>> resolved = self.resolve_input(self.input_structures, 'structures')
            >>> pdb_files = resolved.files
            >>> structure_ids = resolved.ids
            >>> input_datasheets = resolved.datasheets
        """
        # Import here to avoid circular imports
        try:
            from ..base_config import ToolOutput, StandardizedOutput
        except ImportError:
            ToolOutput = type(None)
            StandardizedOutput = type(None)

        # Determine ID key based on expected type
        id_key = f"{expected_type[:-1]}_ids" if expected_type.endswith('s') else f"{expected_type}_ids"

        # Handle ToolOutput
        if ToolOutput is not type(None) and isinstance(input_param, ToolOutput):
            files = input_param.get_output_files(expected_type) or []
            ids = input_param.get_output_files(id_key) or []

            # Fallback: extract IDs from filenames if not provided
            if not ids and files:
                ids = [os.path.splitext(os.path.basename(f))[0] for f in files]

            return ResolvedInput(
                files=files,
                ids=ids,
                datasheets=input_param.datasheets if hasattr(input_param, 'datasheets') else None,
                source=input_param
            )

        # Handle StandardizedOutput
        if StandardizedOutput is not type(None) and isinstance(input_param, StandardizedOutput):
            files = getattr(input_param, expected_type, [])
            ids = getattr(input_param, id_key, [])

            # Fallback: extract IDs from filenames if not provided
            if not ids and files:
                ids = [os.path.splitext(os.path.basename(f))[0] for f in files]

            return ResolvedInput(
                files=files,
                ids=ids,
                datasheets=input_param.datasheets if hasattr(input_param, 'datasheets') else None,
                source=input_param
            )

        # Handle list of file paths
        if isinstance(input_param, list):
            ids = [os.path.splitext(os.path.basename(f))[0] for f in input_param]
            return ResolvedInput(
                files=input_param,
                ids=ids,
                datasheets=None,
                source=input_param
            )

        # Handle single file path string
        if isinstance(input_param, str):
            # Check if it's a file or directory
            if os.path.isdir(input_param):
                # Directory - list files of expected type
                extensions = {
                    'structures': ['.pdb', '.cif'],
                    'sequences': ['.fasta', '.fa', '.csv'],
                    'compounds': ['.sdf', '.mol', '.mol2']
                }
                exts = extensions.get(expected_type, ['.pdb'])

                files = []
                for ext in exts:
                    files.extend([
                        os.path.join(input_param, f)
                        for f in os.listdir(input_param)
                        if f.endswith(ext)
                    ])

                ids = [os.path.splitext(os.path.basename(f))[0] for f in files]
                return ResolvedInput(files=files, ids=ids, datasheets=None, source=input_param)
            else:
                # Single file
                file_id = os.path.splitext(os.path.basename(input_param))[0]
                return ResolvedInput(
                    files=[input_param],
                    ids=[file_id],
                    datasheets=None,
                    source=input_param
                )

        # Handle dictionary format
        if isinstance(input_param, dict):
            files = input_param.get(expected_type, [])
            ids = input_param.get(id_key, [])

            # Fallback: extract IDs from filenames
            if not ids and files:
                ids = [os.path.splitext(os.path.basename(f))[0] for f in files]

            datasheets = input_param.get('datasheets', None)
            return ResolvedInput(files=files, ids=ids, datasheets=datasheets, source=input_param)

        # Unsupported type
        raise ValueError(
            f"Unsupported input type: {type(input_param)}. "
            f"Expected ToolOutput, StandardizedOutput, str, list, or dict."
        )

    def resolve_multiple_inputs(self, *inputs, expected_type: str = 'structures') -> List[ResolvedInput]:
        """
        Resolve multiple inputs at once.

        Args:
            *inputs: Variable number of input parameters
            expected_type: Type of data to extract

        Returns:
            List of ResolvedInput objects

        Example:
            >>> structures, sequences = self.resolve_multiple_inputs(
            ...     self.input_structures,
            ...     self.input_sequences,
            ...     expected_type='structures'
            ... )
        """
        return [self.resolve_input(inp, expected_type) for inp in inputs]
