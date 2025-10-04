"""
SliceDatasheet tool for extracting a subset of rows from tool output datasheets.

Takes tool output and creates a new output containing only the first N rows
from each datasheet, useful for sampling or testing subsets.
"""

import os
import pandas as pd
from typing import Dict, List, Any, Optional, Union

try:
    from .base_config import BaseConfig, ToolOutput, StandardizedOutput, DatasheetInfo
except ImportError:
    # Fallback for direct execution
    import sys
    import os
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, ToolOutput, StandardizedOutput, DatasheetInfo


class SliceDatasheet(BaseConfig):
    """
    Tool for extracting first N rows from tool output datasheets.

    Takes any tool output with datasheets and returns a new output containing
    only the first N rows from each datasheet. Also copies associated structure
    files if they exist.
    """

    TOOL_NAME = "SliceDatasheet"
    DEFAULT_ENV = None  # Loaded from config.yaml

    def __init__(self,
                 input: Union[ToolOutput, StandardizedOutput],
                 n_rows: int = 10,
                 **kwargs):
        """
        Initialize SliceDatasheet tool.

        Args:
            input: Tool output containing datasheets to slice
            n_rows: Number of rows to keep from the beginning of each datasheet
            **kwargs: Additional parameters

        Examples:
            # Take first 10 sequences from LigandMPNN output
            first_ten = pipeline.add(SliceDatasheet(
                input=lmpnn,
                n_rows=10
            ))

            # Take first 5 structures for quick testing
            sample = pipeline.add(SliceDatasheet(
                input=boltz,
                n_rows=5
            ))
        """
        if not isinstance(input, (ToolOutput, StandardizedOutput)):
            raise ValueError("input must be a ToolOutput or StandardizedOutput object")

        if n_rows <= 0:
            raise ValueError("n_rows must be positive")

        self.input = input
        self.n_rows = n_rows

        # Initialize base class
        super().__init__(**kwargs)

        # Set up dependencies
        if hasattr(input, 'config'):
            self.dependencies.append(input.config)

    def validate_params(self):
        """Validate SliceDatasheet parameters."""
        if not isinstance(self.input, (ToolOutput, StandardizedOutput)):
            raise ValueError("input must be a ToolOutput or StandardizedOutput object")

        if self.n_rows <= 0:
            raise ValueError("n_rows must be positive")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure inputs from tool output."""
        self.folders = pipeline_folders

        # Extract input configuration
        self.input_config = self._extract_input_config()

        if not self.input_config["datasheets"]:
            raise ValueError(f"Input tool must provide datasheets to slice")

    def _extract_input_config(self) -> Dict[str, Any]:
        """Extract configuration from input tool output."""
        config = {
            "structures": [],
            "structure_ids": [],
            "sequences": [],
            "sequence_ids": [],
            "compounds": [],
            "compound_ids": [],
            "output_folder": None,
            "datasheets": {}
        }

        # Extract basic attributes
        for attr in ["structures", "structure_ids", "sequences", "sequence_ids",
                    "compounds", "compound_ids", "output_folder"]:
            if hasattr(self.input, attr):
                config[attr] = getattr(self.input, attr) or []

        # Extract datasheets
        if hasattr(self.input, 'datasheets'):
            datasheets = self.input.datasheets

            if hasattr(datasheets, '_datasheets'):
                # Standard BioPipelines format
                config["datasheets"] = datasheets._datasheets
            elif isinstance(datasheets, dict):
                # Dict format - convert to DatasheetInfo objects
                for name, info in datasheets.items():
                    if isinstance(info, dict) and 'path' in info:
                        config["datasheets"][name] = DatasheetInfo(
                            name=name,
                            path=info['path'],
                            columns=info.get('columns', []),
                            description=info.get('description', ''),
                            count=info.get('count', 0)
                        )
                    elif hasattr(info, 'path'):
                        config["datasheets"][name] = info
                    else:
                        config["datasheets"][name] = DatasheetInfo(name=name, path=str(info))

        return config

    def get_config_display(self) -> List[str]:
        """Get configuration display lines."""
        config_lines = super().get_config_display()

        config_lines.extend([
            f"N_ROWS: {self.n_rows}",
            f"INPUT DATASHEETS: {len(self.input_config['datasheets'])}"
        ])

        return config_lines

    def generate_script(self, script_path: str) -> str:
        """
        Generate script to slice datasheets.

        Args:
            script_path: Path where script should be written

        Returns:
            Script content as string
        """
        runtime_folder = os.path.dirname(script_path)
        output_folder = self.output_folder
        os.makedirs(output_folder, exist_ok=True)

        # Create config file for slicing
        config_file = os.path.join(output_folder, "slice_config.json")
        config_data = {
            "input_config": {
                # Convert DatasheetInfo objects to dicts for JSON serialization
                "datasheets": {
                    name: info.to_dict() if hasattr(info, 'to_dict') else str(info)
                    for name, info in self.input_config["datasheets"].items()
                },
                "structures": self.input_config["structures"],
                "structure_ids": self.input_config["structure_ids"],
                "sequences": self.input_config["sequences"],
                "sequence_ids": self.input_config["sequence_ids"],
                "compounds": self.input_config["compounds"],
                "compound_ids": self.input_config["compound_ids"],
                "output_folder": self.input_config["output_folder"]
            },
            "n_rows": self.n_rows,
            "output_folder": output_folder
        }

        import json
        with open(config_file, 'w') as f:
            json.dump(config_data, f, indent=2)

        # Generate script content
        script_content = f"""#!/bin/bash
# SliceDatasheet execution script
# Generated by BioPipelines pipeline system

{self.generate_completion_check_header()}

echo "Slicing datasheets to first {self.n_rows} rows"
echo "Input folder: {self.input_config['output_folder']}"
echo "Output folder: {output_folder}"

# Run Python slicing script
python "{os.path.join(self.folders['HelpScripts'], 'pipe_slice_datasheet.py')}" \\
  --config "{config_file}"

if [ $? -eq 0 ]; then
    echo "Successfully sliced datasheets"
    echo "Results written to: {output_folder}"
else
    echo "Error: Failed to slice datasheets"
    exit 1
fi

{self.generate_completion_check_footer()}
"""

        return script_content

    def get_output_files(self) -> Dict[str, List[str]]:
        """
        Get expected output files after slicing.

        Returns:
            Dictionary with output file paths
        """
        output_config = {
            "structures": [],
            "structure_ids": [],
            "sequences": [],
            "sequence_ids": [],
            "compounds": [],
            "compound_ids": [],
            "datasheets": {},
            "output_folder": self.output_folder
        }

        # Copy structure files (first n_rows only)
        if self.input_config["structures"]:
            for i, struct_path in enumerate(self.input_config["structures"][:self.n_rows]):
                filename = os.path.basename(struct_path)
                output_config["structures"].append(os.path.join(self.output_folder, filename))

        # Copy sequence files (first n_rows only)
        if self.input_config["sequences"]:
            for i, seq_path in enumerate(self.input_config["sequences"][:self.n_rows]):
                filename = os.path.basename(seq_path)
                output_config["sequences"].append(os.path.join(self.output_folder, filename))

        # Copy compound files (first n_rows only)
        if self.input_config["compounds"]:
            for i, comp_path in enumerate(self.input_config["compounds"][:self.n_rows]):
                filename = os.path.basename(comp_path)
                output_config["compounds"].append(os.path.join(self.output_folder, filename))

        # Slice ID arrays
        output_config["structure_ids"] = self.input_config["structure_ids"][:self.n_rows]
        output_config["sequence_ids"] = self.input_config["sequence_ids"][:self.n_rows]
        output_config["compound_ids"] = self.input_config["compound_ids"][:self.n_rows]

        # Process datasheets
        for name, info in self.input_config["datasheets"].items():
            output_filename = os.path.basename(info.path if hasattr(info, 'path') else str(info))
            output_path = os.path.join(self.output_folder, output_filename)

            # Create new DatasheetInfo with sliced description
            output_config["datasheets"][name] = DatasheetInfo(
                name=name,
                path=output_path,
                columns=info.columns if hasattr(info, 'columns') else [],
                description=f"first {self.n_rows} rows from {info.description if hasattr(info, 'description') else name}",
                count=self.n_rows
            )

        return output_config

    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration."""
        base_dict = super().to_dict()
        base_dict.update({
            "tool_params": {
                "n_rows": self.n_rows
            }
        })
        return base_dict