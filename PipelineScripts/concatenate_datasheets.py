"""
ConcatenateDatasheets tool for combining multiple tool outputs.

Concatenates tool outputs from multiple sources, preserving all structures,
compounds, sequences, and datasheets for cyclic pipeline operations.
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


class ConcatenateDatasheets(BaseConfig):
    """
    Tool for concatenating multiple tool outputs into a unified result.
    
    Combines structures, sequences, compounds and datasheets from multiple tool outputs,
    useful for iterative optimization cycles where previous results need to be carried forward.
    """
    
    TOOL_NAME = "ConcatenateDatasheets"
    DEFAULT_ENV = "ProteinEnv"
    DEFAULT_RESOURCES = {"gpu": "T4", "memory": "4GB", "time": "1:00:00"}
    
    def __init__(self,
                 datasheets: List[Any],
                 fill: Optional[str] = None,
                 **kwargs):
        """
        Initialize ConcatenateDatasheets tool.
        
        Args:
            datasheets: List of datasheets to concatenate vertically (like SQL UNION)
            fill: How to handle missing columns - None removes non-common columns, 
                  string value fills missing columns with that value
            **kwargs: Additional parameters
            
        Examples:
            # Simple concatenation removing non-common columns
            combined = pipeline.add(ConcatenateDatasheets(
                datasheets=[cycle0_sequences.output.datasheets.sequences,
                           cycle1_sequences.output.datasheets.sequences]
            ))
            
            # Fill missing columns with empty string
            all_results = pipeline.add(ConcatenateDatasheets(
                datasheets=[tool1.output.datasheets.analysis,
                           tool2.output.datasheets.analysis],
                fill=""
            ))
        """
        self.datasheets = datasheets
        self.fill = fill
        
        # Validate inputs
        if not self.datasheets:
            raise ValueError("At least one datasheet must be provided")
        
        # Initialize base class
        super().__init__(**kwargs)
        
        # Set up dependencies - datasheets can be tool outputs or direct datasheet references
        for datasheet in self.datasheets:
            if hasattr(datasheet, 'config'):
                self.dependencies.append(datasheet.config)
    
    def validate_params(self):
        """Validate ConcatenateDatasheets parameters."""
        if not self.datasheets:
            raise ValueError("At least one datasheet is required")
        
        # fill parameter can be any string or None - no validation needed
    
    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure inputs from datasheets."""
        self.folders = pipeline_folders
        
        # Extract information from each datasheet
        self.datasheet_configs = []
        
        for i, datasheet in enumerate(self.datasheets):
            datasheet_config = self._extract_datasheet_config(datasheet, f"datasheet_{i}")
            self.datasheet_configs.append(datasheet_config)
    
    def _extract_datasheet_config(self, datasheet: Any, prefix: str) -> Dict[str, Any]:
        """Extract configuration from a datasheet reference."""
        config = {
            "prefix": prefix,
            "datasheet_path": None
        }
        
        # Handle different datasheet input formats
        if isinstance(datasheet, str):
            # Direct file path
            config["datasheet_path"] = datasheet
        elif hasattr(datasheet, 'path'):
            # DatasheetInfo object
            config["datasheet_path"] = datasheet.path
        elif hasattr(datasheet, '_datasheets'):
            # Tool output datasheets collection - get first one
            for name, info in datasheet._datasheets.items():
                if hasattr(info, 'path'):
                    config["datasheet_path"] = info.path
                else:
                    config["datasheet_path"] = str(info)
                break
        else:
            # Try to convert to string as fallback
            config["datasheet_path"] = str(datasheet)
        
        return config
    
    def get_config_display(self) -> List[str]:
        """Get configuration display lines."""
        config_lines = super().get_config_display()
        
        config_lines.extend([
            f"DATASHEETS: {len(self.datasheets)} inputs",
            f"FILL: {self.fill if self.fill is not None else 'remove non-common columns'}",
        ])
        
        return config_lines
    
    def generate_script(self, script_path: str) -> str:
        """
        Generate script to concatenate tool outputs.
        
        Args:
            script_path: Path where script should be written
            
        Returns:
            Script content as string
        """
        runtime_folder = os.path.dirname(script_path)
        output_folder = self.output_folder
        os.makedirs(output_folder, exist_ok=True)
        
        # Output paths - only what we actually need
        concatenated_csv = os.path.join(output_folder, "concatenated.csv")
        
        # Create config file for concatenation
        config_file = os.path.join(output_folder, "concatenate_config.json")
        config_data = {
            "datasheet_configs": self.datasheet_configs,
            "fill": self.fill,
            "output_csv": concatenated_csv,
            "output_folder": output_folder
        }
        
        import json
        with open(config_file, 'w') as f:
            json.dump(config_data, f, indent=2)
        
        # Generate script content
        script_content = f"""#!/bin/bash
# ConcatenateDatasheets execution script
# Generated by BioPipelines pipeline system

{self.generate_completion_check_header()}

echo "Concatenating datasheets"
echo "Input datasheets: {len(self.datasheets)}"
echo "Fill strategy: {self.fill if self.fill is not None else 'remove non-common columns'}"
echo "Output: {output_folder}"

# Run Python concatenation script
python "{os.path.join(self.folders['HelpScripts'], 'pipe_concatenate_datasheets.py')}" \\
  --config "{config_file}"

if [ $? -eq 0 ]; then
    echo "Successfully concatenated {len(self.datasheets)} datasheets"
    echo "Output written to: {output_folder}"
else
    echo "Error: Failed to concatenate datasheets"
    exit 1
fi

{self.generate_completion_check_footer()}
"""
        
        return script_content
    
    def get_output_files(self) -> Dict[str, List[str]]:
        """
        Get expected output files after concatenation.
        
        Returns:
            Dictionary with output file paths
        """
        # ConcatenateDatasheets is a datasheet-only operation like MergeDatasheets
        # It combines datasheets from multiple tools but doesn't handle structures/compounds
        
        concatenated_csv = os.path.join(self.output_folder, "concatenated.csv")
        
        # Try to determine output columns from input datasheets
        expected_columns = ["id", "source_datasheet"]
        if self.datasheets and hasattr(self.datasheets[0], 'columns'):
            # Try to get columns from first datasheet if available
            try:
                expected_columns = list(self.datasheets[0].columns) + ["source_datasheet"]
            except:
                # Fallback to default columns
                expected_columns = ["id", "source_datasheet"]
        
        datasheets = {
            "concatenated": DatasheetInfo(
                name="concatenated", 
                path=concatenated_csv,
                columns=expected_columns,
                description="concatenated results from multiple tool outputs",
                count="variable"
            )
        }
        
        # ConcatenateDatasheets operates only on datasheets - return empty lists like MergeDatasheets
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
                "fill": self.fill
            }
        })
        return base_dict