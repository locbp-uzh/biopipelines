# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
ExtractMetrics tool for extracting multiple metrics from multiple tables.

Takes a list of tables and multiple metric names, then extracts all values
for each metric from each table. Creates separate CSV files for each metric,
perfect for copying into Prism for column graph analysis.
"""

import os
from typing import Dict, List, Any, Optional, Union

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

    TOOL_NAME = "ExtractMetrics"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        return """echo "=== ExtractMetrics ==="
echo "Uses biopipelines environment (no additional installation needed)."
echo "=== ExtractMetrics ready ==="
"""

    # Lazy path descriptors
    config_file = Path(lambda self: os.path.join(self.output_folder, "extract_metrics_config.json"))
    extract_py = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_extract_metrics.py"))

    def __init__(self,
                 tables: List[Union[str, TableInfo]],
                 metrics: List[str],
                 table_names: Optional[List[str]] = None,
                 **kwargs):
        """
        Initialize ExtractMetrics tool.

        Args:
            tables: List of tables to extract metrics from (TableInfo or path strings)
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

        super().__init__(**kwargs)

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

    def _extract_table_path(self, input_obj: Union[str, TableInfo], input_type: str) -> str:
        """Extract table path from various input formats."""
        if isinstance(input_obj, str):
            return input_obj
        elif isinstance(input_obj, TableInfo):
            return input_obj.info.path
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
        """Generate script to extract metrics from tables."""
        os.makedirs(self.output_folder, exist_ok=True)

        script_content = "#!/bin/bash\n"
        script_content += "# ExtractMetrics execution script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()
        script_content += self.generate_script_run_extraction()
        script_content += self.generate_completion_check_footer()

        return script_content

    def generate_script_run_extraction(self) -> str:
        """Generate the metrics extraction part of the script."""
        import json

        config_data = {
            "tables": self.table_paths,
            "metrics": self.metrics,
            "table_names": self.table_names,
            "output_folder": self.output_folder
        }

        with open(self.config_file, 'w') as f:
            json.dump(config_data, f, indent=2)

        return f"""echo "Extracting {len(self.metrics)} metrics from {len(self.table_paths)} tables"
echo "Metrics: {', '.join(self.metrics)}"
echo "Output folder: {self.output_folder}"

python "{self.extract_py}" --config "{self.config_file}"

"""

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files after metric extraction."""
        tables = {}

        for metric in self.metrics:
            csv_path = os.path.join(self.output_folder, f"{metric}.csv")

            tables[metric] = TableInfo(
                name=metric,
                path=csv_path,
                columns=self.table_names,
                description=f"Extracted '{metric}' values from all tables",
                count=None
            )

        return {
            "structures": DataStream.empty("structures", "pdb"),
            "sequences": DataStream.empty("sequences", "fasta"),
            "compounds": DataStream.empty("compounds", "sdf"),
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