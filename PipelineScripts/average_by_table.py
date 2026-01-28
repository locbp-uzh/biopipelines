"""
AverageByTable tool for calculating averages across multiple tables.

Takes a list of tables and computes the average of numeric columns across all tables.
Each row in the output corresponds to the average from each input table.
Non-numeric values and values that are not consistent across all tables are discarded.
"""

import os
from typing import Dict, List, Any, Optional, Union

try:
    from .base_config import BaseConfig, ToolOutput, StandardizedOutput, TableInfo
except ImportError:
    # Fallback for direct execution
    import sys
    import os
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, ToolOutput, StandardizedOutput, TableInfo


class AverageByTable(BaseConfig):
    """
    Pipeline tool for calculating averages across multiple tables.

    Takes a list of tables (typically from different cycles) and computes
    the average of numeric columns. Each row in the output corresponds to the
    average values from one input table.

    Commonly used for:
    - Analyzing cycle progression in iterative design
    - Computing average metrics across different conditions
    - Summarizing results from multiple experimental runs
    """

    # Tool identification
    TOOL_NAME = "AverageByTable"
    

    def __init__(self,
                 tables: List[Union[str, Dict, ToolOutput, StandardizedOutput]],
                 **kwargs):
        """
        Initialize AverageByTable tool.

        Args:
            tables: List of tables to average across
            **kwargs: Additional parameters

        Examples:
            # Average metrics across cycles
            averages = pipeline.add(AverageByTable(
                tables=[cycle0.tables.merged,
                           cycle1.tables.merged,
                           cycle2.tables.merged]
            ))
        """
        self.tables_input = tables

        # Initialize base class
        super().__init__(**kwargs)

        # Set up dependencies
        for ds in tables:
            if hasattr(ds, 'config'):
                self.dependencies.append(ds.config)

    def validate_params(self):
        """Validate AverageByTable parameters."""
        if not self.tables_input:
            raise ValueError("tables parameter is required and cannot be empty")

        if not isinstance(self.tables_input, list):
            raise ValueError("tables must be a list")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input tables from previous tools."""
        self.folders = pipeline_folders

        # Extract table paths
        self.table_paths = []
        for i, ds in enumerate(self.tables_input):
            path = self._extract_table_path(ds, f"table_{i}")
            self.table_paths.append(path)

    def _extract_table_path(self, input_obj: Union[str, Dict, ToolOutput, StandardizedOutput], input_type: str) -> str:
        """Extract table path from various input formats."""
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
            # TableInfo object
            return input_obj.path
        else:
            raise ValueError(f"Could not extract path from {input_type}: {type(input_obj)}")

    def get_config_display(self) -> List[str]:
        """Get configuration display lines."""
        config_lines = super().get_config_display()

        config_lines.extend([
            f"NUM_DATASHEETS: {len(self.tables_input)}",
            f"DATASHEET_TYPES: {[type(ds).__name__ for ds in self.tables_input[:3]]}{'...' if len(self.tables_input) > 3 else ''}"
        ])

        return config_lines

    def generate_script(self, script_path: str) -> str:
        """Generate script to compute averages across tables."""
        os.makedirs(self.output_folder, exist_ok=True)

        script_content = "#!/bin/bash\n"
        script_content += "# AverageByTable execution script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()
        script_content += self.generate_script_run_average()
        script_content += self.generate_completion_check_footer()

        return script_content

    def generate_script_run_average(self) -> str:
        """Generate the averaging part of the script."""
        averages_csv = os.path.join(self.output_folder, "averages.csv")

        config_file = os.path.join(self.output_folder, "average_by_table_config.json")
        config_data = {
            "tables": self.table_paths,
            "output": averages_csv
        }

        import json
        with open(config_file, 'w') as f:
            json.dump(config_data, f, indent=2)

        return f"""echo "Computing averages across tables"
echo "Number of tables: {len(self.table_paths)}"
echo "Output: {averages_csv}"

python "{os.path.join(self.folders['HelpScripts'], 'pipe_average_by_table.py')}" \\
  --config "{config_file}"

if [ $? -eq 0 ]; then
    echo "Successfully computed averages"
    echo "Output file: {averages_csv}"
else
    echo "Error: Failed to compute averages"
    exit 1
fi

"""

    def get_output_files(self) -> Dict[str, List[str]]:
        """
        Get expected output files after averaging.

        Returns:
            Dictionary with output file paths and table information
        """
        averages_csv = os.path.join(self.output_folder, "averages.csv")

        # Define table that will be created
        tables = {
            "averages": TableInfo(
                name="averages",
                path=averages_csv,
                columns=["table_name"] + ["avg_*"],  # Will be determined at runtime
                description="Average values for each numeric column across all tables",
                count=len(self.tables_input)  # One row per input table
            )
        }

        return {
            "structures": [],
            "structure_ids": [],
            "compounds": [],
            "compound_ids": [],
            "sequences": [],
            "sequence_ids": [],
            "tables": tables,
            "output_folder": self.output_folder
        }

    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration."""
        base_dict = super().to_dict()
        base_dict.update({
            "tool_params": {
                "num_tables": len(self.tables_input)
            }
        })
        return base_dict