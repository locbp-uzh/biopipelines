"""
PyMOL visualization tool with declarative operation-based API.

Creates PyMOL session files (.pse) from structure outputs using a sequence
of operations like Load, Color, Align, etc. Supports ID-based matching
between structures and table columns for per-structure selections.
"""

import os
import json
from typing import Dict, List, Any, Optional, Union, Tuple

try:
    from .base_config import BaseConfig, ToolOutput, StandardizedOutput, TableInfo
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, ToolOutput, StandardizedOutput, TableInfo


class PyMOLOperation:
    """Base class for PyMOL operations."""

    def __init__(self, op_type: str, **kwargs):
        self.op_type = op_type
        self.params = kwargs

    def to_dict(self) -> Dict[str, Any]:
        """Convert operation to dictionary for serialization."""
        return {"op": self.op_type, **self.params}


class PyMOL(BaseConfig):
    """
    PyMOL visualization tool with declarative operation-based API.

    Uses a sequence of operations (Names, Load, Color, Align, etc.) to build
    PyMOL sessions. Supports ID-based matching between structures and table
    columns for per-structure selections and naming.

    Example:
        PyMOL(
            PyMOL.Names(prefix="sensor", basename=fuse.tables.sequences.lengths, suffix="parts"),
            PyMOL.Load(folded),
            PyMOL.Color(folded, selection=fuse.tables.sequences.D1, color="white"),
            PyMOL.Color(folded, selection=fuse.tables.sequences.L1, color="orange"),
            PyMOL.Align("align")
        )
    """

    TOOL_NAME = "PyMOL"
    DEFAULT_ENV = None

    # --- Static methods for creating operations ---

    @staticmethod
    def Names(prefix: str = "", basename=None, suffix: str = "") -> PyMOLOperation:
        """
        Set up ID → PyMOL name mapping for subsequent operations.

        For each ID in the basename table column, creates a mapping:
        ID → "{prefix}_{basename_value}_{suffix}"

        Args:
            prefix: Prefix for PyMOL object names
            basename: Table column reference for the variable part of the name
                      (e.g., fuse.tables.sequences.lengths)
            suffix: Suffix for PyMOL object names

        Returns:
            PyMOLOperation for names mapping
        """
        return PyMOLOperation("names", prefix=prefix, basename=basename, suffix=suffix)

    @staticmethod
    def Load(structures: Union[StandardizedOutput, ToolOutput]) -> PyMOLOperation:
        """
        Load structures into PyMOL with current naming.

        For each structure, uses its ID to look up the PyMOL name from the
        current naming map (set by Names). If no mapping exists, uses the ID
        as the PyMOL object name.

        Args:
            structures: Tool output containing structures to load

        Returns:
            PyMOLOperation for loading structures
        """
        return PyMOLOperation("load", structures=structures)

    @staticmethod
    def Color(structures: Union[StandardizedOutput, ToolOutput],
              selection: Union[Tuple[TableInfo, str], str],
              color: str) -> PyMOLOperation:
        """
        Color a selection on structures.

        For each structure, uses its ID to:
        1. Look up the PyMOL name from the current naming map
        2. Look up the selection value from the table column
        Then applies: color {color}, {pymol_name} and resi {selection}

        Args:
            structures: Tool output containing structures to color
            selection: Table column reference for per-structure selection
                       (e.g., fuse.tables.sequences.D1) or fixed selection string
            color: PyMOL color name (e.g., "white", "marine", "orange")

        Returns:
            PyMOLOperation for coloring
        """
        return PyMOLOperation("color", structures=structures, selection=selection, color=color)

    @staticmethod
    def ColorAF(structures: Union[StandardizedOutput, ToolOutput],
                upper: float = 100) -> PyMOLOperation:
        """
        Color structures by AlphaFold pLDDT (B-factor spectrum).

        Applies spectrum coloring based on B-factor values, which contain
        pLDDT scores in AlphaFold-generated structures.

        Args:
            structures: Tool output containing AlphaFold structures
            upper: Maximum B-factor value for scaling (default: 100 for pLDDT)

        Returns:
            PyMOLOperation for pLDDT coloring
        """
        return PyMOLOperation("coloraf", structures=structures, upper=upper)

    @staticmethod
    def Align(method: str = "align", target: Optional[str] = None) -> PyMOLOperation:
        """
        Align all loaded objects.

        If target is not specified, aligns all objects to the first loaded object.

        Args:
            method: Alignment method - "align", "super", or "cealign"
            target: Target object name for alignment (default: first loaded)

        Returns:
            PyMOLOperation for alignment
        """
        return PyMOLOperation("align", method=method, target=target)

    @staticmethod
    def Show(representation: str = "cartoon",
             selection: Optional[str] = None) -> PyMOLOperation:
        """
        Show a representation for structures.

        Args:
            representation: PyMOL representation (cartoon, sticks, surface, etc.)
            selection: Optional selection to apply to (default: all)

        Returns:
            PyMOLOperation for showing representation
        """
        return PyMOLOperation("show", representation=representation, selection=selection)

    @staticmethod
    def Hide(representation: str = "everything",
             selection: Optional[str] = None) -> PyMOLOperation:
        """
        Hide a representation for structures.

        Args:
            representation: PyMOL representation to hide
            selection: Optional selection to apply to (default: all)

        Returns:
            PyMOLOperation for hiding representation
        """
        return PyMOLOperation("hide", representation=representation, selection=selection)

    @staticmethod
    def Set(setting: str, value: Any, selection: Optional[str] = None) -> PyMOLOperation:
        """
        Set a PyMOL setting.

        Args:
            setting: PyMOL setting name
            value: Setting value
            selection: Optional selection to apply to

        Returns:
            PyMOLOperation for setting
        """
        return PyMOLOperation("set", setting=setting, value=value, selection=selection)

    @staticmethod
    def Save(filename: str = "session.pse") -> PyMOLOperation:
        """
        Save the PyMOL session.

        Args:
            filename: Output filename (relative to output folder)

        Returns:
            PyMOLOperation for saving
        """
        return PyMOLOperation("save", filename=filename)

    # --- Instance methods ---

    def __init__(self, *args, session: str = "session", **kwargs):
        """
        Initialize PyMOL tool with a sequence of operations.

        Args:
            session: Name for the output session file (without .pse extension)
            *args: Sequence of PyMOL operations (Names, Load, Color, etc.)
            **kwargs: Additional configuration parameters

        Example:
            PyMOL(session="my_session",
                  PyMOL.Load(structures),
                  PyMOL.ColorAF(structures))
        """
        self.operations = list(args)
        self.session_name = session

        # Extract structure sources and table references from operations
        self._structure_sources = []
        self._table_references = []
        self._extract_dependencies()

        super().__init__(**kwargs)

    def _extract_dependencies(self):
        """Extract structure sources and table references from operations."""
        for op in self.operations:
            if op.op_type == "load":
                structures = op.params.get("structures")
                if structures is not None:
                    self._structure_sources.append(structures)

            elif op.op_type == "color":
                structures = op.params.get("structures")
                if structures is not None and structures not in self._structure_sources:
                    self._structure_sources.append(structures)

                selection = op.params.get("selection")
                if isinstance(selection, tuple):
                    self._table_references.append(selection)

            elif op.op_type == "coloraf":
                structures = op.params.get("structures")
                if structures is not None and structures not in self._structure_sources:
                    self._structure_sources.append(structures)

            elif op.op_type == "names":
                basename = op.params.get("basename")
                if isinstance(basename, tuple):
                    self._table_references.append(basename)

    def validate_params(self):
        """Validate PyMOL parameters."""
        if not self.operations:
            raise ValueError("At least one operation must be provided")

        # Check that all operations are valid PyMOLOperation objects
        for i, op in enumerate(self.operations):
            if not isinstance(op, PyMOLOperation):
                raise ValueError(f"Operation {i} is not a PyMOLOperation: {type(op)}")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input sources from pipeline context."""
        self.folders = pipeline_folders
        self.input_sources = {}

        # Track structure sources
        for i, source in enumerate(self._structure_sources):
            if hasattr(source, 'output_folder'):
                self.input_sources[f"structures_{i}"] = {
                    'output_folder': source.output_folder,
                    'structures': getattr(source, 'structures', []),
                    'structure_ids': getattr(source, 'structure_ids', [])
                }

    def _serialize_operation(self, op: PyMOLOperation) -> Dict[str, Any]:
        """Serialize an operation to a dictionary for JSON config."""
        result = {"op": op.op_type}

        for key, value in op.params.items():
            if isinstance(value, StandardizedOutput):
                # Serialize StandardizedOutput reference
                result[key] = {
                    "type": "standardized_output",
                    "output_folder": value.output_folder,
                    "structures": value.structures,
                    "structure_ids": value.structure_ids
                }
            elif isinstance(value, tuple) and len(value) == 2:
                # Table column reference: (TableInfo, column_name)
                table_info, column_name = value
                if hasattr(table_info, 'path'):
                    result[key] = {
                        "type": "table_column",
                        "table_path": table_info.path,
                        "column_name": column_name
                    }
                else:
                    result[key] = str(value)
            elif hasattr(value, 'output_folder'):
                # ToolOutput or similar
                result[key] = {
                    "type": "tool_output",
                    "output_folder": value.output_folder,
                    "structures": getattr(value, 'structures', []),
                    "structure_ids": getattr(value, 'structure_ids', [])
                }
            else:
                result[key] = value

        return result

    def generate_script(self, script_path: str) -> str:
        """Generate bash script for PyMOL session creation."""
        # Create config file with all operations
        config = {
            "operations": [self._serialize_operation(op) for op in self.operations],
            "session_name": self.session_name,
            "output_folder": self.output_folder
        }

        config_file = os.path.join(self.output_folder, "pymol_config.json")
        session_file = os.path.join(self.output_folder, f"{self.session_name}.pse")

        # Get HelpScripts path
        help_scripts = self.folders.get("HelpScripts", "HelpScripts")
        pymol_script = os.path.join(help_scripts, "pipe_pymol.py")

        script_content = f"""#!/bin/bash
# PyMOL session creation script
# Generated by BioPipelines

{self.generate_completion_check_header()}

echo "Creating PyMOL session..."
echo "Output: {session_file}"

# Create output directory
mkdir -p "{self.output_folder}"

# Write configuration
cat > "{config_file}" << 'PYMOL_CONFIG_EOF'
{json.dumps(config, indent=2)}
PYMOL_CONFIG_EOF

# Run PyMOL session creation
python "{pymol_script}" --config "{config_file}"

if [ $? -eq 0 ]; then
    echo "PyMOL session created successfully"
else
    echo "ERROR: Failed to create PyMOL session"
    exit 1
fi

{self.generate_completion_check_footer()}
"""

        return script_content

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files."""
        session_file = os.path.join(self.output_folder, f"{self.session_name}.pse")

        return {
            "structures": [],
            "structure_ids": [],
            "compounds": [],
            "compound_ids": [],
            "sequences": [],
            "sequence_ids": [],
            "tables": {},
            "output_folder": self.output_folder,
            "session_file": session_file
        }

    def get_config_display(self) -> List[str]:
        """Get configuration display lines."""
        config_lines = super().get_config_display()

        config_lines.extend([
            f"OPERATIONS: {len(self.operations)}",
            f"SESSION: {self.session_name}.pse"
        ])

        # Show operation summary
        op_types = [op.op_type for op in self.operations]
        config_lines.append(f"SEQUENCE: {' → '.join(op_types)}")

        return config_lines

    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration to dictionary."""
        base_dict = super().to_dict()
        base_dict.update({
            "pymol_params": {
                "session_name": self.session_name,
                "num_operations": len(self.operations),
                "operation_types": [op.op_type for op in self.operations]
            }
        })
        return base_dict

    def __str__(self) -> str:
        """String representation."""
        op_types = [op.op_type for op in self.operations]
        return f"PyMOL({' → '.join(op_types)})"
