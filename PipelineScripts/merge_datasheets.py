"""
MergeDatasheets tool for merging multiple analysis CSV files.

Merges analysis results from multiple tools into a unified datasheet for filtering.
Handles metric name collisions with prefixes and maintains all information.
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


class MergeDatasheets(BaseConfig):
    """
    Tool for combining multiple analysis datasheets into a unified CSV.
    
    Merges CSV files from multiple analysis tools on a common key column,
    handling metric name collisions with prefixes and preserving all data.
    """
    
    TOOL_NAME = "MergeDatasheets"
    DEFAULT_ENV = "ProteinEnv"
    
    def __init__(self,
                 datasheets: List[Any],
                 prefixes: Optional[List[str]] = None,
                 key: str = "id",
                 calculate: Optional[Dict[str, str]] = None,
                 id_map: Optional[Dict[str, List[str]]] = None,
                 **kwargs):
        """
        Initialize MergeDatasheets tool.
        
        Args:
            datasheets: List of datasheets (file paths, ToolOutput.datasheets entries, etc.)
            prefixes: Optional list of prefixes for each datasheet (must match datasheets length)
            key: Column name to merge on (default: "id")
            calculate: Optional dict of {new_column: expression} for calculated columns
            id_map: Optional dict mapping new_id -> [old_id1, old_id2, ...] to consolidate different IDs
            **kwargs: Additional parameters
            
        Examples:
            # Simple combination
            merged = pipeline.add(MergeDatasheets(
                datasheets=[confidence_results.datasheets.analysis,
                           distance_results.datasheets.analysis],
                key="structure_id"
            ))
            
            # With prefixes and explicit datasheet specification
            merged = pipeline.add(MergeDatasheets(
                datasheets=[boltz_apo.datasheets.affinity, boltz_holo.datasheets.affinity],
                prefixes=["apo_", "holo_"],
                calculate={"affinity_diff": "holo_affinity - apo_affinity"}
            ))
            
            # With ID mapping to consolidate different IDs into common entities
            merged = pipeline.add(MergeDatasheets(
                datasheets=[open_results.datasheets.affinity, close_results.datasheets.affinity],
                prefixes=["open_", "close_"],
                id_map={"original": ["HT_Cy7_C_R", "HT_Cy7_C_RR"]},
                calculate={"affinity_delta": "open_affinity_pred_value - close_affinity_pred_value"}
            ))
        """
        self.datasheets = datasheets
        self.prefixes = prefixes or []
        self.merge_key = key
        self.calculate = calculate or {}
        self.id_map = id_map or {}
        
        # Validate inputs
        if not self.datasheets:
            raise ValueError("At least one datasheet must be provided")
        
        if self.prefixes and len(self.prefixes) != len(self.datasheets):
            raise ValueError("Number of prefixes must match number of datasheets")
        
        
        # Initialize base class
        super().__init__(**kwargs)
        
        # Set up dependencies
        for datasheet in self.datasheets:
            if hasattr(datasheet, 'config'):
                self.dependencies.append(datasheet.config)
    
    def validate_params(self):
        """Validate MergeDatasheets parameters."""
        if not self.merge_key:
            raise ValueError("merge_key cannot be empty")
        
        if not self.datasheets:
            raise ValueError("At least one datasheet is required")
        
        if self.prefixes and len(self.prefixes) != len(self.datasheets):
            raise ValueError("Number of prefixes must match number of datasheets")
    
    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input datasheets from analysis tools."""
        self.folders = pipeline_folders
        self.input_csv_paths = []
        
        for i, datasheet in enumerate(self.datasheets):
            csv_path = self._get_csv_path_from_datasheet(datasheet)
            if csv_path is None:
                raise ValueError(f"Could not determine CSV path from datasheet: {datasheet}")
            self.input_csv_paths.append(csv_path)
    
    def _get_csv_path_from_datasheet(self, datasheet: Any) -> Optional[str]:
        """
        Get predicted CSV file path from datasheet reference.
        
        Args:
            datasheet: Datasheet reference (path, dict with 'path', or ToolOutput)
            
        Returns:
            Predicted path to CSV file or None if cannot determine
        """
        # Handle tuple from column reference (datasheet_info, column_name)
        if isinstance(datasheet, tuple) and len(datasheet) == 2:
            # Extract the DatasheetInfo object from the tuple
            datasheet_info, column_name = datasheet
            # Use the DatasheetInfo path
            if hasattr(datasheet_info, 'path'):
                return datasheet_info.path

        # Direct path string
        if isinstance(datasheet, str):
            return datasheet
        
        # Dict with path key
        if isinstance(datasheet, dict) and 'path' in datasheet:
            return datasheet['path']

        # DatasheetInfo object (from LoadOutput)
        if hasattr(datasheet, 'path') and hasattr(datasheet, 'name') and hasattr(datasheet, 'columns'):
            # This is likely a DatasheetInfo object - use its path directly
            return datasheet.path

        # ToolOutput object
        if hasattr(datasheet, 'datasheets'):
            datasheets = datasheet.datasheets
            if isinstance(datasheets, dict):
                if len(datasheets) == 1:
                    # Single datasheet - use it automatically
                    ds_info = next(iter(datasheets.values()))
                else:
                    # Multiple datasheets - ambiguous
                    datasheet_names = list(datasheets.keys())
                    raise ValueError(f"Ambiguous: ToolOutput has multiple datasheets {datasheet_names}. Use explicit specification like datasheets=[tool.datasheets.affinity, ...]")
                
                if isinstance(ds_info, dict) and 'path' in ds_info:
                    return ds_info['path']
                else:
                    return str(ds_info)
        
        # Predict based on output folder (simplified - no assumptions about names)
        if hasattr(datasheet, 'output_folder'):
            output_folder = datasheet.output_folder
            # Just use a generic prediction - no assumptions about specific tool outputs
            return os.path.join(output_folder, 'results.csv')
        
        return None
    
    def get_config_display(self) -> List[str]:
        """Get configuration display lines."""
        config_lines = super().get_config_display()
        
        config_lines.extend([
            f"INPUTS: {len(self.datasheets)} datasheets",
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
        Generate script to merge analysis datasheets.
        
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
# MergeDatasheets execution script
# Generated by BioPipelines pipeline system

{self.generate_completion_check_header()}

echo "Combining analysis datasheets"
echo "Input datasheets: {len(self.input_csv_paths)}"
echo "Merge key: {self.merge_key}"
echo "Output: {merged_csv}"

# Run Python combination script
python "{os.path.join(self.folders['HelpScripts'], 'pipe_merge_datasheets.py')}" \\
  --config "{config_file}"

if [ $? -eq 0 ]; then
    echo "Successfully merged {len(self.input_csv_paths)} datasheets"
    echo "Output written to: {merged_csv}"
else
    echo "Error: Failed to merge datasheets"
    exit 1
fi

{self.generate_completion_check_footer()}
"""
        
        return script_content
    
    def get_output_files(self) -> Dict[str, List[str]]:
        """
        Get expected output files after combining datasheets.
        
        Returns:
            Dictionary with output file paths
        """
        merged_csv = os.path.join(self.output_folder, "merged.csv")
        
        # Determine expected columns by combining input datasheets info
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
        
        datasheets = {
            "merged": DatasheetInfo(
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
            "datasheets": datasheets,
            "output_folder": self.output_folder
        }
    
    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration."""
        base_dict = super().to_dict()
        base_dict.update({
            "tool_params": {
                "num_datasheets": len(self.datasheets),
                "merge_key": self.merge_key,
                "prefixes": self.prefixes,
                "calculate": self.calculate,
                "id_map": self.id_map
            }
        })
        return base_dict