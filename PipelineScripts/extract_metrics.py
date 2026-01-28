"""
ExtractMetrics tool for extracting multiple metrics from multiple tables.

Takes a list of tables and multiple metric names, then extracts all values
for each metric from each table. Creates separate CSV files for each metric,
perfect for copying into Prism for column graph analysis.
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


class ExtractMetrics(BaseConfig):
    """
    Pipeline tool for extracting multiple metrics from multiple tables.

    Takes a list of tables and extracts all values for specified metrics.
    Creates separate CSV files for each metric, with output format optimized
    for statistical analysis tools like Prism.

    Commonly used for:
    - Extracting multiple metrics across cycles for trend analysis
    - Preparing data for statistical plotting with separate files per metric
    - Comparing metric distributions between conditions
    - Generating input files for Prism column graphs
    """

    # Tool identification
    TOOL_NAME = "ExtractMetrics"
    

    def __init__(self,
                 tables: List[Union[str, Dict, ToolOutput, StandardizedOutput]],
                 metrics: List[str],
                 table_names: Optional[List[str]] = None,
                 **kwargs):
        """
        Initialize ExtractMetrics tool.

        Args:
            tables: List of tables to extract metrics from
            metrics: List of metric column names to extract
            table_names: Optional custom names for columns (defaults to table indices)
            **kwargs: Additional parameters

        Examples:
            # Extract multiple metrics across cycles
            metrics_extract = pipeline.add(ExtractMetrics(
                tables=[cycle0.tables.merged,
                           cycle1.tables.merged,
                           cycle2.tables.merged],
                metrics=["affinity_delta", "affinity_delta_R", "affinity_delta_S"],
                table_names=["Cycle0", "Cycle1", "Cycle2"]
            ))
        """
        self.tables_input = tables
        self.metrics = metrics
        self.table_names = table_names

        # Initialize base class
        super().__init__(**kwargs)

        # Set up dependencies
        for ds in tables:
            if hasattr(ds, 'config'):
                self.dependencies.append(ds.config)

    def validate_params(self):
        """Validate ExtractMetrics parameters."""
        if not self.tables_input:
            raise ValueError("tables parameter is required and cannot be empty")

        if not isinstance(self.tables_input, list):
            raise ValueError("tables must be a list")

        if not self.metrics:
            raise ValueError("metrics parameter is required and cannot be empty")

        if not isinstance(self.metrics, list):
            raise ValueError("metrics must be a list")

        if self.table_names and len(self.table_names) != len(self.tables_input):
            raise ValueError("table_names length must match tables length")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input tables from previous tools."""
        self.folders = pipeline_folders

        # Extract table paths
        self.table_paths = []
        for i, ds in enumerate(self.tables_input):
            path = self._extract_table_path(ds, f"table_{i}")
            self.table_paths.append(path)

        # Set default column names if not provided
        if not self.table_names:
            self.table_names = [f"Table_{i}" for i in range(len(self.tables_input))]

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
            f"METRICS: {', '.join(self.metrics[:3])}{'...' if len(self.metrics) > 3 else ''} ({len(self.metrics)} total)",
            f"NUM_DATASHEETS: {len(self.tables_input)}",
            f"COLUMN_NAMES: {', '.join(self.table_names[:3])}{'...' if len(self.table_names) > 3 else ''}"
        ])

        return config_lines

    def generate_script(self, script_path: str) -> str:
        """
        Generate script to extract metrics from tables.

        Args:
            script_path: Path where script should be written

        Returns:
            Script content as string
        """
        runtime_folder = os.path.dirname(script_path)
        output_folder = self.output_folder
        os.makedirs(output_folder, exist_ok=True)

        # Create config file for the helper script
        config_file = os.path.join(output_folder, "extract_metrics_config.json")
        config_data = {
            "tables": self.table_paths,
            "metrics": self.metrics,
            "table_names": self.table_names,
            "output_folder": output_folder
        }

        import json
        with open(config_file, 'w') as f:
            json.dump(config_data, f, indent=2)

        # Generate script content
        script_content = f"""#!/bin/bash
# ExtractMetrics execution script
# Generated by BioPipelines pipeline system

{self.generate_completion_check_header()}

echo "Extracting {len(self.metrics)} metrics from {len(self.table_paths)} tables"
echo "Metrics: {', '.join(self.metrics)}"
echo "Output folder: {output_folder}"

# Run Python helper script
python "{os.path.join(self.folders['HelpScripts'], 'pipe_extract_metrics.py')}" \\
  --config "{config_file}"

if [ $? -eq 0 ]; then
    echo "Successfully extracted all metrics"
    echo "Output files created in: {output_folder}"
    echo "Ready for copy-paste into Prism or other analysis software"
else
    echo "Error: Failed to extract metrics"
    exit 1
fi

{self.generate_completion_check_footer()}
"""

        return script_content

    def get_output_files(self) -> Dict[str, List[str]]:
        """
        Get expected output files after metric extraction.

        Returns:
            Dictionary with output file paths and table information
        """
        # Create separate CSV file for each metric
        tables = {}

        for metric in self.metrics:
            csv_path = os.path.join(self.output_folder, f"{metric}.csv")

            tables[metric] = TableInfo(
                name=metric,
                path=csv_path,
                columns=self.table_names,  # Column names from table names
                description=f"Extracted '{metric}' values from all tables, formatted for statistical analysis",
                count=None  # Will be determined at runtime based on data
            )

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
                "metrics": self.metrics,
                "num_tables": len(self.tables_input),
                "table_names": self.table_names
            }
        })
        return base_dict