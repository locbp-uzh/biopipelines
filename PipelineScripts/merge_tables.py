"""
MergeTables tool for merging multiple analysis CSV files.

Merges analysis results from multiple tools into a unified table for filtering.
Handles metric name collisions with prefixes and maintains all information.
"""

import os
import pandas as pd
from typing import Dict, List, Any, Optional, Union

try:
    from .base_config import BaseConfig, ToolOutput, StandardizedOutput, TableInfo
except ImportError:
    # Fallback for direct execution
    import sys
    import os
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, ToolOutput, StandardizedOutput, TableInfo


class MergeTables(BaseConfig):
    """
    Tool for combining multiple analysis tables into a unified CSV.
    
    Merges CSV files from multiple analysis tools on a common key column,
    handling metric name collisions with prefixes and preserving all data.
    """
    
    TOOL_NAME = "MergeTables"
    DEFAULT_ENV = None  # Loaded from config.yaml
    
    def __init__(self,
                 tables: List[Any],
                 prefixes: Optional[List[str]] = None,
                 key: str = "id",
                 calculate: Optional[Dict[str, str]] = None,
                 id_map: Optional[Dict[str, List[str]]] = None,
                 **kwargs):
        """
        Initialize MergeTables tool.
        
        Args:
            tables: List of tables (file paths, ToolOutput.tables entries, etc.)
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
        self.tables = tables
        self.prefixes = prefixes or []
        self.merge_key = key
        self.calculate = calculate or {}
        self.id_map = id_map or {}
        
        # Validate inputs
        if not self.tables:
            raise ValueError("At least one table must be provided")
        
        if self.prefixes and len(self.prefixes) != len(self.tables):
            raise ValueError("Number of prefixes must match number of tables")
        
        
        # Initialize base class
        super().__init__(**kwargs)
        
        # Set up dependencies
        for table in self.tables:
            if hasattr(table, 'config'):
                self.dependencies.append(table.config)
    
    def validate_params(self):
        """Validate MergeTables parameters."""
        if not self.merge_key:
            raise ValueError("merge_key cannot be empty")
        
        if not self.tables:
            raise ValueError("At least one table is required")
        
        if self.prefixes and len(self.prefixes) != len(self.tables):
            raise ValueError("Number of prefixes must match number of tables")
    
    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input tables from analysis tools."""
        self.folders = pipeline_folders
        self.input_csv_paths = []
        
        for i, table in enumerate(self.tables):
            csv_path = self._get_csv_path_from_table(table)
            if csv_path is None:
                raise ValueError(f"Could not determine CSV path from table: {table}")
            self.input_csv_paths.append(csv_path)
    
    def _get_csv_path_from_table(self, table: Any) -> Optional[str]:
        """
        Get predicted CSV file path from table reference.
        
        Args:
            table: Table reference (path, dict with 'path', or ToolOutput)
            
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
        
        # Dict with path key
        if isinstance(table, dict) and 'path' in table:
            return table['path']

        # TableInfo object (from LoadOutput)
        if hasattr(table, 'path') and hasattr(table, 'name') and hasattr(table, 'columns'):
            # This is likely a TableInfo object - use its path directly
            return table.path

        # ToolOutput object
        if hasattr(table, 'tables'):
            tables = table.tables
            if isinstance(tables, dict):
                if len(tables) == 1:
                    # Single table - use it automatically
                    ds_info = next(iter(tables.values()))
                else:
                    # Multiple tables - ambiguous
                    table_names = list(tables.keys())
                    raise ValueError(f"Ambiguous: ToolOutput has multiple tables {table_names}. Use explicit specification like tables=[tool.tables.affinity, ...]")
                
                if isinstance(ds_info, dict) and 'path' in ds_info:
                    return ds_info['path']
                else:
                    return str(ds_info)
        
        # Predict based on output folder (simplified - no assumptions about names)
        if hasattr(table, 'output_folder'):
            output_folder = table.output_folder
            # Just use a generic prediction - no assumptions about specific tool outputs
            return os.path.join(output_folder, 'results.csv')
        
        return None
    
    def get_config_display(self) -> List[str]:
        """Get configuration display lines."""
        config_lines = super().get_config_display()
        
        config_lines.extend([
            f"INPUTS: {len(self.tables)} tables",
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
        """
        Generate script to merge analysis tables.
        
        Args:
            script_path: Path where script should be written
            
        Returns:
            Script content as string
        """
        runtime_folder = os.path.dirname(script_path)
        output_folder = self.output_folder
        os.makedirs(output_folder, exist_ok=True)
        
        # Output CSV path
        merged_csv = os.path.join(output_folder, "merged.csv")
        
        # Create config file for the combination
        config_file = os.path.join(output_folder, "merge_config.json")
        config_data = {
            "input_csvs": self.input_csv_paths,
            "merge_key": self.merge_key,
            "prefixes": self.prefixes,
            "calculate": self.calculate,
            "id_map": self.id_map,
            "output_csv": merged_csv
        }
        
        import json
        with open(config_file, 'w') as f:
            json.dump(config_data, f, indent=2)
        
        # Generate script content
        script_content = f"""#!/bin/bash
# MergeTables execution script
# Generated by BioPipelines pipeline system

{self.generate_completion_check_header()}

echo "Combining analysis tables"
echo "Input tables: {len(self.input_csv_paths)}"
echo "Merge key: {self.merge_key}"
echo "Output: {merged_csv}"

# Run Python combination script
python "{os.path.join(self.folders['HelpScripts'], 'pipe_merge_tables.py')}" \\
  --config "{config_file}"

if [ $? -eq 0 ]; then
    echo "Successfully merged {len(self.input_csv_paths)} tables"
    echo "Output written to: {merged_csv}"
else
    echo "Error: Failed to merge tables"
    exit 1
fi

{self.generate_completion_check_footer()}
"""
        
        return script_content
    
    def get_output_files(self) -> Dict[str, List[str]]:
        """
        Get expected output files after combining tables.
        
        Returns:
            Dictionary with output file paths
        """
        merged_csv = os.path.join(self.output_folder, "merged.csv")
        
        # Determine expected columns by combining input tables info
        expected_columns = [self.merge_key, "source_structure"]  # Common columns
        
        # Predict prefixed metric columns based on known tool outputs
        if self.prefixes:
            # Predict common metric names with prefixes
            common_metrics = [
                "affinity_pred_value", "affinity_confidence",  # Boltz2 affinity
                "pLDDT", "confidence",  # General confidence metrics
                "distance", "energy"     # Analysis metrics
            ]
            
            for i, prefix in enumerate(self.prefixes):
                if prefix:  # Only add prefixed columns for non-empty prefixes
                    for metric in common_metrics:
                        expected_columns.append(f"{prefix}{metric}")
                else:
                    # For empty prefix, add the original metric names
                    expected_columns.extend(common_metrics)
        
        # Add calculated columns
        if self.calculate:
            expected_columns.extend(self.calculate.keys())
        
        tables = {
            "merged": TableInfo(
                name="merged",
                path=merged_csv,
                columns=expected_columns,
                description="merged analysis results from multiple tools",
                count="variable"  # Depends on input data
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
                "num_tables": len(self.tables),
                "merge_key": self.merge_key,
                "prefixes": self.prefixes,
                "calculate": self.calculate,
                "id_map": self.id_map
            }
        })
        return base_dict