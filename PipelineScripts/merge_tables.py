"""
MergeTables tool for merging multiple analysis CSV files.

Merges analysis results from multiple tools into a unified table for filtering.
Handles metric name collisions with prefixes and maintains all information.
"""

import os
import json
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


class MergeTables(BaseConfig):
    """
    Tool for combining multiple analysis tables into a unified CSV.

    Merges CSV files from multiple analysis tools on a common key column,
    handling metric name collisions with prefixes and preserving all data.
    """

    TOOL_NAME = "MergeTables"

    # Lazy path descriptors
    merged_csv = Path(lambda self: os.path.join(self.output_folder, "merged.csv"))
    config_file = Path(lambda self: os.path.join(self.output_folder, "merge_config.json"))
    helper_script = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_merge_tables.py"))

    def __init__(self,
                 tables: List[Union[StandardizedOutput, TableInfo, str]],
                 prefixes: Optional[List[str]] = None,
                 key: str = "id",
                 calculate: Optional[Dict[str, str]] = None,
                 id_map: Optional[Dict[str, List[str]]] = None,
                 **kwargs):
        """
        Initialize MergeTables tool.

        Args:
            tables: List of tables (file paths, StandardizedOutput.tables entries, TableInfo objects)
            prefixes: Optional list of prefixes for each table (must match tables length)
            key: Column name to merge on (default: "id")
            calculate: Optional dict of {new_column: expression} for calculated columns
            id_map: Optional dict mapping new_id -> [old_id1, old_id2, ...] to consolidate different IDs
            **kwargs: Additional parameters

        Examples:
            # Simple combination
            merged = pipeline.add(MergeTables(
                tables=[confidence_results.tables.analysis,
                           distance_results.tables.analysis],
                key="structure_id"
            ))

            # With prefixes and explicit table specification
            merged = pipeline.add(MergeTables(
                tables=[boltz_apo.tables.affinity, boltz_holo.tables.affinity],
                prefixes=["apo_", "holo_"],
                calculate={"affinity_diff": "holo_affinity - apo_affinity"}
            ))

            # With ID mapping to consolidate different IDs into common entities
            merged = pipeline.add(MergeTables(
                tables=[open_results.tables.affinity, close_results.tables.affinity],
                prefixes=["open_", "close_"],
                id_map={"original": ["HT_Cy7_C_R", "HT_Cy7_C_RR"]},
                calculate={"affinity_delta": "open_affinity_pred_value - close_affinity_pred_value"}
            ))
        """
        self.tables_input = tables
        self.prefixes = prefixes or []
        self.merge_key = key
        self.calculate = calculate or {}
        self.id_map = id_map or {}

        # Validate inputs
        if not self.tables_input:
            raise ValueError("At least one table must be provided")

        if self.prefixes and len(self.prefixes) != len(self.tables_input):
            raise ValueError("Number of prefixes must match number of tables")

        # Initialize base class
        super().__init__(**kwargs)

    def validate_params(self):
        """Validate MergeTables parameters."""
        if not self.merge_key:
            raise ValueError("merge_key cannot be empty")

        if not self.tables_input:
            raise ValueError("At least one table is required")

        if self.prefixes and len(self.prefixes) != len(self.tables_input):
            raise ValueError("Number of prefixes must match number of tables")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input tables from analysis tools."""
        self.folders = pipeline_folders
        self.input_csv_paths = []

        for i, table in enumerate(self.tables_input):
            csv_path = self._get_csv_path_from_table(table)
            if csv_path is None:
                raise ValueError(f"Could not determine CSV path from table: {table}")
            self.input_csv_paths.append(csv_path)

    def _get_csv_path_from_table(self, table: Any) -> Optional[str]:
        """
        Get predicted CSV file path from table reference.

        Args:
            table: Table reference (path, TableInfo, or StandardizedOutput)

        Returns:
            Predicted path to CSV file or None if cannot determine
        """
        # Handle tuple from column reference (table_info, column_name)
        if isinstance(table, tuple) and len(table) == 2:
            # Extract the TableInfo object from the tuple
            table_info, column_name = table
            # Use the TableInfo path
            if hasattr(table_info, 'path'):
                return table_info.path

        # Direct path string
        if isinstance(table, str):
            return table

        # TableInfo object
        if isinstance(table, TableInfo):
            return table.path

        # StandardizedOutput object
        if hasattr(table, 'tables'):
            tables = table.tables

            if hasattr(tables, '_tables'):
                if len(tables._tables) == 1:
                    # Single table - use it automatically
                    ds_info = next(iter(tables._tables.values()))
                    if hasattr(ds_info, 'path'):
                        return ds_info.path
                else:
                    # Multiple tables - ambiguous
                    table_names = list(tables._tables.keys())
                    raise ValueError(f"Ambiguous: StandardizedOutput has multiple tables {table_names}. Use explicit specification like tables=[tool.tables.affinity, ...]")

        # Predict based on output folder
        if hasattr(table, 'output_folder'):
            output_folder = table.output_folder
            return os.path.join(output_folder, 'results.csv')

        return None

    def get_config_display(self) -> List[str]:
        """Get configuration display lines."""
        config_lines = super().get_config_display()

        config_lines.extend([
            f"INPUTS: {len(self.tables_input)} tables",
            f"MERGE KEY: {self.merge_key}",
        ])

        if self.prefixes:
            prefixes_str = ", ".join(f"{i}: '{prefix}'" for i, prefix in enumerate(self.prefixes))
            config_lines.append(f"PREFIXES: {prefixes_str}")

        if self.calculate:
            calc_str = ", ".join(f"{k}={v}" for k, v in self.calculate.items())
            config_lines.append(f"CALCULATED: {calc_str}")

        if self.id_map:
            id_map_str = ", ".join(f"{new_id}={old_ids}" for new_id, old_ids in self.id_map.items())
            config_lines.append(f"ID_MAP: {id_map_str}")

        return config_lines

    def generate_script(self, script_path: str) -> str:
        """Generate MergeTables execution script."""
        os.makedirs(self.output_folder, exist_ok=True)

        # Create config file for the combination
        config_data = {
            "input_csvs": self.input_csv_paths,
            "merge_key": self.merge_key,
            "prefixes": self.prefixes,
            "calculate": self.calculate,
            "id_map": self.id_map,
            "output_csv": self.merged_csv
        }

        with open(self.config_file, 'w') as f:
            json.dump(config_data, f, indent=2)

        script_content = "#!/bin/bash\n"
        script_content += "# MergeTables execution script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()

        script_content += f"""echo "Combining analysis tables"
echo "Input tables: {len(self.input_csv_paths)}"
echo "Merge key: {self.merge_key}"
echo "Output: {self.merged_csv}"

python "{self.helper_script}" \\
  --config "{self.config_file}"

if [ $? -eq 0 ]; then
    echo "Successfully merged {len(self.input_csv_paths)} tables"
    echo "Output written to: {self.merged_csv}"
else
    echo "Error: Failed to merge tables"
    exit 1
fi

"""
        script_content += self.generate_completion_check_footer()

        return script_content

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files after combining tables."""
        # Only predict columns that are guaranteed: merge key and calculated columns
        expected_columns = [self.merge_key]

        # Add calculated columns
        if self.calculate:
            expected_columns.extend(self.calculate.keys())

        tables = {
            "merged": TableInfo(
                name="merged",
                path=self.merged_csv,
                columns=expected_columns,
                description="Merged analysis results from multiple tools",
                count="variable"
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
        """Serialize configuration."""
        base_dict = super().to_dict()
        base_dict.update({
            "tool_params": {
                "num_tables": len(self.tables_input),
                "merge_key": self.merge_key,
                "prefixes": self.prefixes,
                "calculate": self.calculate,
                "id_map": self.id_map
            }
        })
        return base_dict
