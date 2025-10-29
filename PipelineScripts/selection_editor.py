"""
SelectionEditor configuration for structure-aware PyMOL selection modification.

Modifies PyMOL-formatted selection strings (e.g., "3-45+58-60") with operations
like expand, shrink, shift, and invert, while respecting actual PDB residue numbering.
"""

import os
from typing import Dict, List, Any, Optional, Union

try:
    from .base_config import BaseConfig, ToolOutput, StandardizedOutput, DatasheetInfo
except ImportError:
    # Fallback for direct execution
    import sys
    import os
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, ToolOutput, StandardizedOutput, DatasheetInfo


class SelectionEditor(BaseConfig):
    """
    SelectionEditor configuration for structure-aware selection modification.

    Takes PyMOL-formatted selection strings from a datasheet column and modifies them
    using expand, shrink, shift, or invert operations while validating against actual
    PDB residue numbers.
    """

    TOOL_NAME = "SelectionEditor"
    DEFAULT_ENV = None  # Loaded from config.yaml

    def __init__(self,
                 structures: Union[ToolOutput, StandardizedOutput, List[str]],
                 selection: tuple,
                 expand: int = 0,
                 shrink: int = 0,
                 shift: int = 0,
                 invert: bool = False,
                 **kwargs):
        """
        Initialize SelectionEditor configuration.

        Args:
            structures: Input structures (ToolOutput/StandardizedOutput from previous tool or list of PDB files)
            selection: Datasheet column reference tuple (e.g., tool.datasheets.structures.designed)
            expand: Number of residues to add on each side of intervals (default: 0)
            shrink: Number of residues to remove from each side of intervals (default: 0)
            shift: Number of residues to shift all intervals (+/-) (default: 0)
            invert: Whether to invert the selection (select complement) (default: False)
            **kwargs: Additional parameters
        """
        # Store SelectionEditor-specific parameters
        self.selection_ref = selection
        self.structures_input = structures
        self.expand = expand
        self.shrink = shrink
        self.shift = shift
        self.invert = invert

        # Track input types
        self.structures_is_tool_output = isinstance(structures, (ToolOutput, StandardizedOutput))

        # Validate selection reference format
        if not isinstance(selection, tuple) or len(selection) != 2:
            raise ValueError(
                "selection must be a datasheet column reference tuple. "
                "Example: tool.datasheets.structures.designed"
            )

        # Initialize base class
        super().__init__(**kwargs)

        # Extract the source tool from the selection reference
        self.selection_datasheet, self.selection_column = self.selection_ref

        # Add dependency on the tool that provides the structures
        if isinstance(structures, (ToolOutput, StandardizedOutput)):
            # Get tool config from StandardizedOutput or ToolOutput
            tool_config = structures.tool if hasattr(structures, 'tool') else (
                structures.config if hasattr(structures, 'config') else None
            )
            if tool_config:
                self.dependencies.append(tool_config)

        # Initialize file paths (will be set in configure_inputs)
        self._initialize_file_paths()

    def validate_params(self):
        """Validate SelectionEditor-specific parameters."""
        if not self.structures_input:
            raise ValueError("structures parameter is required")

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

    def _initialize_file_paths(self):
        """Initialize common file paths used throughout the class."""
        self.selections_csv = None
        self.config_json = None

    def _setup_file_paths(self):
        """Set up all file paths after output_folder is known."""
        # Output files
        self.selections_csv = os.path.join(self.output_folder, "selections.csv")
        self.config_json = os.path.join(self.output_folder, "config.json")

        # Helper script paths
        if hasattr(self, 'folders') and self.folders:
            self.selection_editor_py = os.path.join(
                self.folders["HelpScripts"],
                "pipe_selection_editor.py"
            )
        else:
            self.selection_editor_py = None

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input sources from selection datasheet and structures."""
        self.folders = pipeline_folders
        self._setup_file_paths()

        # Get the selection datasheet path
        if hasattr(self.selection_datasheet, 'path'):
            self.selection_datasheet_path = self.selection_datasheet.path
        else:
            raise ValueError("Invalid selection datasheet reference")

        # Get structures from input parameter
        if self.structures_is_tool_output:
            # structures_input is a StandardizedOutput or ToolOutput
            tool_output = self.structures_input

            # Try to get structures
            source_structures = []
            if isinstance(tool_output, StandardizedOutput):
                source_structures = tool_output.structures
            elif hasattr(tool_output, 'get_output_files'):
                for struct_type in ["structures", "pdbs"]:
                    struct_files = tool_output.get_output_files(struct_type)
                    if struct_files:
                        source_structures = struct_files
                        break

            if not source_structures:
                raise ValueError(f"No structure outputs found from input")

            self.input_sources = {"structures": source_structures}

        elif isinstance(self.structures_input, list):
            self.input_sources = {"structures": self.structures_input}

        else:
            raise ValueError(f"Unsupported structures input type: {type(self.structures_input)}")

    def get_config_display(self) -> List[str]:
        """Get SelectionEditor configuration display lines."""
        config_lines = super().get_config_display()

        # Selection information
        config_lines.append(f"SELECTION COLUMN: {self.selection_column}")
        config_lines.append(f"SELECTION SOURCE: {self.selection_datasheet.name}")

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

        # Structure count
        if hasattr(self, 'input_sources') and "structures" in self.input_sources:
            struct_count = len(self.input_sources["structures"])
            config_lines.append(f"STRUCTURES: {struct_count}")

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

        # Create config file for helper script
        import json
        config_data = {
            "selection_datasheet": self.selection_datasheet_path,
            "selection_column": self.selection_column,
            "structures": self.input_sources["structures"],
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
        script_content += "# Generated by BioPipelines\n\n"
        script_content += self.generate_completion_check_header()

        script_content += f"""echo "Modifying PyMOL selections"
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

    def _predict_structure_ids(self) -> List[str]:
        """
        Predict the structure IDs that will be analyzed.

        Returns:
            List of structure identifiers
        """
        structure_ids = []

        if hasattr(self, 'input_sources') and "structures" in self.input_sources:
            structure_files = self.input_sources["structures"]
            for struct_file in structure_files:
                filename = os.path.basename(struct_file)
                if filename.endswith('.pdb'):
                    structure_ids.append(filename[:-4])
                elif filename.endswith('.cif'):
                    structure_ids.append(filename[:-4])
                else:
                    structure_ids.append(filename)

        return structure_ids

    def get_output_files(self) -> Dict[str, List[str]]:
        """
        Get expected output files after SelectionEditor execution.

        Returns:
            Dictionary mapping output types to file paths with standard keys
        """
        # Ensure file paths are set up
        if not hasattr(self, 'selections_csv') or self.selections_csv is None:
            self._setup_file_paths()

        # Predict structure IDs
        structure_ids = self._predict_structure_ids()

        # Create column names
        modified_col = self.selection_column
        original_col = f"original_{self.selection_column}"

        # Organize datasheets
        datasheets = {
            "selections": DatasheetInfo(
                name="selections",
                path=self.selections_csv,
                columns=["id", "pdb", modified_col, original_col],
                description=f"Modified PyMOL selections from {self.selection_column}",
                count=len(structure_ids)
            )
        }

        return {
            "structures": [],
            "structure_ids": [],
            "compounds": [],
            "compound_ids": [],
            "sequences": [],
            "sequence_ids": [],
            "datasheets": datasheets,
            "output_folder": self.output_folder
        }

    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration including SelectionEditor-specific parameters."""
        base_dict = super().to_dict()
        base_dict.update({
            "selection_editor_params": {
                "selection_column": self.selection_column,
                "selection_datasheet": self.selection_datasheet.name,
                "expand": self.expand,
                "shrink": self.shrink,
                "shift": self.shift,
                "invert": self.invert,
                "structures_provided": self.structures_provided
            }
        })
        return base_dict
