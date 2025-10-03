"""
ExtractMetrics tool for extracting multiple metrics from multiple datasheets.

Takes a list of datasheets and multiple metric names, then extracts all values
for each metric from each datasheet. Creates separate CSV files for each metric,
perfect for copying into Prism for column graph analysis.
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


class ExtractMetrics(BaseConfig):
    """
    Pipeline tool for extracting multiple metrics from multiple datasheets.

    Takes a list of datasheets and extracts all values for specified metrics.
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
    DEFAULT_ENV = "ProteinEnv"

    def __init__(self,
                 datasheets: List[Union[str, Dict, ToolOutput, StandardizedOutput]],
                 metrics: List[str],
                 datasheet_names: Optional[List[str]] = None,
                 **kwargs):
        """
        Initialize ExtractMetrics tool.

        Args:
            datasheets: List of datasheets to extract metrics from
            metrics: List of metric column names to extract
            datasheet_names: Optional custom names for columns (defaults to datasheet indices)
            **kwargs: Additional parameters

        Examples:
            # Extract multiple metrics across cycles
            metrics_extract = pipeline.add(ExtractMetrics(
                datasheets=[cycle0.datasheets.merged,
                           cycle1.datasheets.merged,
                           cycle2.datasheets.merged],
                metrics=["affinity_delta", "affinity_delta_R", "affinity_delta_S"],
                datasheet_names=["Cycle0", "Cycle1", "Cycle2"]
            ))
        """
        self.datasheets_input = datasheets
        self.metrics = metrics
        self.datasheet_names = datasheet_names

        # Initialize base class
        super().__init__(**kwargs)

        # Set up dependencies
        for ds in datasheets:
            if hasattr(ds, 'config'):
                self.dependencies.append(ds.config)

    def validate_params(self):
        """Validate ExtractMetrics parameters."""
        if not self.datasheets_input:
            raise ValueError("datasheets parameter is required and cannot be empty")

        if not isinstance(self.datasheets_input, list):
            raise ValueError("datasheets must be a list")

        if not self.metrics:
            raise ValueError("metrics parameter is required and cannot be empty")

        if not isinstance(self.metrics, list):
            raise ValueError("metrics must be a list")

        if self.datasheet_names and len(self.datasheet_names) != len(self.datasheets_input):
            raise ValueError("datasheet_names length must match datasheets length")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input datasheets from previous tools."""
        self.folders = pipeline_folders

        # Extract datasheet paths
        self.datasheet_paths = []
        for i, ds in enumerate(self.datasheets_input):
            path = self._extract_datasheet_path(ds, f"datasheet_{i}")
            self.datasheet_paths.append(path)

        # Set default column names if not provided
        if not self.datasheet_names:
            self.datasheet_names = [f"Datasheet_{i}" for i in range(len(self.datasheets_input))]

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
            f"METRICS: {', '.join(self.metrics[:3])}{'...' if len(self.metrics) > 3 else ''} ({len(self.metrics)} total)",
            f"NUM_DATASHEETS: {len(self.datasheets_input)}",
            f"COLUMN_NAMES: {', '.join(self.datasheet_names[:3])}{'...' if len(self.datasheet_names) > 3 else ''}"
        ])

        return config_lines

    def generate_script(self, script_path: str) -> str:
        """
        Generate script to extract metrics from datasheets.

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
            "datasheets": self.datasheet_paths,
            "metrics": self.metrics,
            "datasheet_names": self.datasheet_names,
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

echo "Extracting {len(self.metrics)} metrics from {len(self.datasheet_paths)} datasheets"
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
            Dictionary with output file paths and datasheet information
        """
        # Create separate CSV file for each metric
        datasheets = {}

        for metric in self.metrics:
            csv_path = os.path.join(self.output_folder, f"{metric}.csv")

            datasheets[metric] = DatasheetInfo(
                name=metric,
                path=csv_path,
                columns=self.datasheet_names,  # Column names from datasheet names
                description=f"Extracted '{metric}' values from all datasheets, formatted for statistical analysis",
                count=None  # Will be determined at runtime based on data
            )

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
                "metrics": self.metrics,
                "num_datasheets": len(self.datasheets_input),
                "datasheet_names": self.datasheet_names
            }
        })
        return base_dict