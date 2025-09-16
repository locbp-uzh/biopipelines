"""
AverageByDatasheet tool for calculating averages across multiple datasheets.

Takes a list of datasheets and computes the average of numeric columns across all datasheets.
Each row in the output corresponds to the average from each input datasheet.
Non-numeric values and values that are not consistent across all datasheets are discarded.
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


class AverageByDatasheet(BaseConfig):
    """
    Pipeline tool for calculating averages across multiple datasheets.

    Takes a list of datasheets (typically from different cycles) and computes
    the average of numeric columns. Each row in the output corresponds to the
    average values from one input datasheet.

    Commonly used for:
    - Analyzing cycle progression in iterative design
    - Computing average metrics across different conditions
    - Summarizing results from multiple experimental runs
    """

    # Tool identification
    TOOL_NAME = "AverageByDatasheet"
    DEFAULT_ENV = "ProteinEnv"
    COMPATIBLE_ENVS = ["ProteinEnv", "MutationEnv"]
    DEFAULT_RESOURCES = {"gpu": "T4", "memory": "4GB", "time": "1:00:00"}

    def __init__(self,
                 datasheets: List[Union[str, Dict, ToolOutput, StandardizedOutput]],
                 **kwargs):
        """
        Initialize AverageByDatasheet tool.

        Args:
            datasheets: List of datasheets to average across
            **kwargs: Additional parameters

        Examples:
            # Average metrics across cycles
            averages = pipeline.add(AverageByDatasheet(
                datasheets=[cycle0.output.datasheets.merged,
                           cycle1.output.datasheets.merged,
                           cycle2.output.datasheets.merged]
            ))
        """
        self.datasheets_input = datasheets

        # Initialize base class
        super().__init__(**kwargs)

        # Set up dependencies
        for ds in datasheets:
            if hasattr(ds, 'config'):
                self.dependencies.append(ds.config)

    def validate_params(self):
        """Validate AverageByDatasheet parameters."""
        if not self.datasheets_input:
            raise ValueError("datasheets parameter is required and cannot be empty")

        if not isinstance(self.datasheets_input, list):
            raise ValueError("datasheets must be a list")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input datasheets from previous tools."""
        self.folders = pipeline_folders

        # Extract datasheet paths
        self.datasheet_paths = []
        for i, ds in enumerate(self.datasheets_input):
            path = self._extract_datasheet_path(ds, f"datasheet_{i}")
            self.datasheet_paths.append(path)

    def _extract_datasheet_path(self, input_obj: Union[str, Dict, ToolOutput, StandardizedOutput], input_type: str) -> str:
        """Extract datasheet path from various input formats."""
        if isinstance(input_obj, str):
            # Direct file path
            return input_obj
        elif isinstance(input_obj, dict):
            # Dictionary with path
            if 'path' in input_obj:
                return input_obj['path']
            else:
                raise ValueError(f"Dictionary input for {input_type} must have 'path' key")
        elif hasattr(input_obj, 'path'):
            # DatasheetInfo object
            return input_obj.path
        else:
            raise ValueError(f"Could not extract path from {input_type}: {type(input_obj)}")

    def get_config_display(self) -> List[str]:
        """Get configuration display lines."""
        config_lines = super().get_config_display()

        config_lines.extend([
            f"NUM_DATASHEETS: {len(self.datasheets_input)}",
            f"DATASHEET_TYPES: {[type(ds).__name__ for ds in self.datasheets_input[:3]]}{'...' if len(self.datasheets_input) > 3 else ''}"
        ])

        return config_lines

    def generate_script(self, script_path: str) -> str:
        """
        Generate script to compute averages across datasheets.

        Args:
            script_path: Path where script should be written

        Returns:
            Script content as string
        """
        runtime_folder = os.path.dirname(script_path)
        output_folder = self.output_folder
        os.makedirs(output_folder, exist_ok=True)

        # Output file
        averages_csv = os.path.join(output_folder, "averages.csv")

        # Create config file for the helper script
        config_file = os.path.join(output_folder, "average_by_datasheet_config.json")
        config_data = {
            "datasheets": self.datasheet_paths,
            "output": averages_csv
        }

        import json
        with open(config_file, 'w') as f:
            json.dump(config_data, f, indent=2)

        # Generate script content
        script_content = f"""#!/bin/bash
# AverageByDatasheet execution script
# Generated by BioPipelines pipeline system

{self.generate_completion_check_header()}

echo "Computing averages across datasheets"
echo "Number of datasheets: {len(self.datasheet_paths)}"
echo "Output: {averages_csv}"

# Run Python helper script
python "{os.path.join(self.folders['HelpScripts'], 'pipe_average_by_datasheet.py')}" \\
  --config "{config_file}"

if [ $? -eq 0 ]; then
    echo "Successfully computed averages"
    echo "Output file: {averages_csv}"
else
    echo "Error: Failed to compute averages"
    exit 1
fi

{self.generate_completion_check_footer()}
"""

        return script_content

    def get_output_files(self) -> Dict[str, List[str]]:
        """
        Get expected output files after averaging.

        Returns:
            Dictionary with output file paths and datasheet information
        """
        averages_csv = os.path.join(self.output_folder, "averages.csv")

        # Define datasheet that will be created
        datasheets = {
            "averages": DatasheetInfo(
                name="averages",
                path=averages_csv,
                columns=["datasheet_name"] + ["avg_*"],  # Will be determined at runtime
                description="Average values for each numeric column across all datasheets",
                count=len(self.datasheets_input)  # One row per input datasheet
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
        """Serialize configuration."""
        base_dict = super().to_dict()
        base_dict.update({
            "tool_params": {
                "num_datasheets": len(self.datasheets_input)
            }
        })
        return base_dict