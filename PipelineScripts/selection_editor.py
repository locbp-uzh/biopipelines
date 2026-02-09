# Copyright (c) 2026 Gianluca Quargnali @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
SelectionEditor configuration for structure-aware PyMOL selection modification.

Modifies PyMOL-formatted selection strings (e.g., "3-45+58-60") with operations
like expand, shrink, shift, and invert, while respecting actual PDB residue numbering.
"""

import os
import json
from typing import Dict, List, Any, Union

try:
    from .base_config import BaseConfig, StandardizedOutput, TableInfo
    from .file_paths import Path
    from .datastream import DataStream
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo
    from file_paths import Path
    from datastream import DataStream


class SelectionEditor(BaseConfig):
    """
    SelectionEditor configuration for structure-aware selection modification.

    Takes PyMOL-formatted selection strings from a table column and modifies them
    using expand, shrink, shift, or invert operations while validating against actual
    PDB residue numbers.
    """

    TOOL_NAME = "SelectionEditor"

    # Lazy path descriptors
    selections_csv = Path(lambda self: os.path.join(self.output_folder, "selections.csv"))
    config_json = Path(lambda self: os.path.join(self.output_folder, "config.json"))
    structures_ds_json = Path(lambda self: os.path.join(self.output_folder, "structures.json"))
    selection_editor_py = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_selection_editor.py"))

    def __init__(self,
                 structures: Union[DataStream, StandardizedOutput],
                 selection: tuple,
                 expand: int = 0,
                 shrink: int = 0,
                 shift: int = 0,
                 invert: bool = False,
                 **kwargs):
        """
        Initialize SelectionEditor configuration.

        Args:
            structures: Input structures as DataStream or StandardizedOutput
            selection: Table column reference tuple (e.g., tool.tables.structures.designed)
            expand: Number of residues to add on each side of intervals (default: 0)
            shrink: Number of residues to remove from each side of intervals (default: 0)
            shift: Number of residues to shift all intervals (+/-) (default: 0)
            invert: Whether to invert the selection (select complement) (default: False)
            **kwargs: Additional parameters
        """
        # Resolve input to DataStream
        if isinstance(structures, StandardizedOutput):
            self.structures_stream: DataStream = structures.streams.structures
        elif isinstance(structures, DataStream):
            self.structures_stream = structures
        else:
            raise ValueError(f"structures must be DataStream or StandardizedOutput, got {type(structures)}")

        # Store SelectionEditor-specific parameters
        self.selection_ref = selection
        self.expand = expand
        self.shrink = shrink
        self.shift = shift
        self.invert = invert

        # Validate selection reference format
        if not isinstance(selection, tuple) or len(selection) != 2:
            raise ValueError(
                "selection must be a table column reference tuple. "
                "Example: tool.tables.structures.designed"
            )

        # Extract the source table and column from the selection reference
        self.selection_table, self.selection_column = self.selection_ref

        # Initialize base class
        super().__init__(**kwargs)

    def validate_params(self):
        """Validate SelectionEditor-specific parameters."""
        if not self.structures_stream or len(self.structures_stream) == 0:
            raise ValueError("structures parameter is required and must not be empty")

        if not self.selection_ref:
            raise ValueError("selection parameter is required")

        # Validate that at least one operation is specified
        if (self.expand == 0 and self.shrink == 0 and
            self.shift == 0 and not self.invert):
            raise ValueError(
                "At least one operation must be specified: "
                "expand, shrink, shift, or invert"
            )

        # Validate operation values
        if self.expand < 0:
            raise ValueError("expand must be non-negative")
        if self.shrink < 0:
            raise ValueError("shrink must be non-negative")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input sources."""
        self.folders = pipeline_folders

        # Get the selection table path
        if hasattr(self.selection_table, 'path'):
            self.selection_table_path = self.selection_table.path
        else:
            raise ValueError("Invalid selection table reference")

    def get_config_display(self) -> List[str]:
        """Get SelectionEditor configuration display lines."""
        config_lines = super().get_config_display()

        config_lines.extend([
            f"INPUT STRUCTURES: {len(self.structures_stream)} files",
            f"SELECTION COLUMN: {self.selection_column}",
            f"SELECTION SOURCE: {self.selection_table.name}"
        ])

        # Operations
        operations = []
        if self.expand > 0:
            operations.append(f"expand={self.expand}")
        if self.shrink > 0:
            operations.append(f"shrink={self.shrink}")
        if self.shift != 0:
            operations.append(f"shift={self.shift:+d}")
        if self.invert:
            operations.append("invert=True")

        config_lines.append(f"OPERATIONS: {', '.join(operations)}")

        return config_lines

    def generate_script(self, script_path: str) -> str:
        """
        Generate SelectionEditor execution script.

        Args:
            script_path: Path where script should be written

        Returns:
            Script content as string
        """
        # Ensure output folder exists
        os.makedirs(self.output_folder, exist_ok=True)

        # Serialize structures DataStream to JSON for HelpScript to load
        with open(self.structures_ds_json, 'w') as f:
            json.dump(self.structures_stream.to_dict(), f, indent=2)

        # Create config file for helper script
        config_data = {
            "selection_table": self.selection_table_path,
            "selection_column": self.selection_column,
            "structures_json": self.structures_ds_json,
            "expand": self.expand,
            "shrink": self.shrink,
            "shift": self.shift,
            "invert": self.invert,
            "output_csv": self.selections_csv
        }

        with open(self.config_json, 'w') as f:
            json.dump(config_data, f, indent=2)

        # Generate bash script
        script_content = "#!/bin/bash\n"
        script_content += "# SelectionEditor execution script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()

        script_content += f"""echo "Modifying PyMOL selections"
echo "Input structures: {len(self.structures_stream)} files"
echo "Selection column: {self.selection_column}"
echo "Operations: expand={self.expand}, shrink={self.shrink}, shift={self.shift}, invert={self.invert}"

# Run selection modification
python {self.selection_editor_py} \\
    "{self.config_json}"

if [ $? -eq 0 ]; then
    echo "Selection modification completed successfully"
else
    echo "Selection modification failed"
    exit 1
fi

echo "Modified selections saved to: {self.selections_csv}"

"""
        script_content += self.generate_completion_check_footer()

        return script_content

    def get_output_files(self) -> Dict[str, Any]:
        """
        Get expected output files after SelectionEditor execution.

        Returns:
            Dictionary with DataStream objects and tables
        """
        # Create column names
        modified_col = self.selection_column
        original_col = f"original_{self.selection_column}"

        # Organize tables
        tables = {
            "selections": TableInfo(
                name="selections",
                path=self.selections_csv,
                columns=["id", "pdb", modified_col, original_col],
                description=f"Modified PyMOL selections from {self.selection_column}",
                count=len(self.structures_stream)
            )
        }

        return {
            "structures": DataStream.empty("structures", "pdb"),
            "sequences": DataStream.empty("sequences", "fasta"),
            "compounds": DataStream.empty("compounds", "sdf"),
            "tables": tables,
            "output_folder": self.output_folder
        }

    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration including SelectionEditor-specific parameters."""
        base_dict = super().to_dict()
        base_dict.update({
            "selection_editor_params": {
                "selection_column": self.selection_column,
                "selection_table": self.selection_table.name,
                "expand": self.expand,
                "shrink": self.shrink,
                "shift": self.shift,
                "invert": self.invert
            }
        })
        return base_dict
