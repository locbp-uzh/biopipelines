"""
ConcatenateTables tool for combining multiple tool outputs.

Concatenates tool outputs from multiple sources, preserving all structures,
compounds, sequences, and tables for cyclic pipeline operations.
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


class ConcatenateTables(BaseConfig):
    """
    Tool for concatenating multiple tool outputs into a unified result.
    
    Combines structures, sequences, compounds and tables from multiple tool outputs,
    useful for iterative optimization cycles where previous results need to be carried forward.
    """
    
    TOOL_NAME = "ConcatenateTables"
    
    
    def __init__(self,
                 tables: List[Any],
                 fill: Optional[str] = None,
                 **kwargs):
        """
        Initialize ConcatenateTables tool.
        
        Args:
            tables: List of tables to concatenate vertically (like SQL UNION)
            fill: How to handle missing columns - None removes non-common columns, 
                  string value fills missing columns with that value
            **kwargs: Additional parameters
            
        Examples:
            # Simple concatenation removing non-common columns
            combined = pipeline.add(ConcatenateTables(
                tables=[cycle0_sequences.tables.sequences,
                           cycle1_sequences.tables.sequences]
            ))
            
            # Fill missing columns with empty string
            all_results = pipeline.add(ConcatenateTables(
                tables=[tool1.tables.analysis,
                           tool2.tables.analysis],
                fill=""
            ))
        """
        self.tables = tables
        self.fill = fill
        
        # Validate inputs
        if not self.tables:
            raise ValueError("At least one table must be provided")
        
        # Initialize base class
        super().__init__(**kwargs)
        
        # Set up dependencies - tables can be tool outputs or direct table references
        for table in self.tables:
            if hasattr(table, 'config'):
                self.dependencies.append(table.config)
    
    def validate_params(self):
        """Validate ConcatenateTables parameters."""
        if not self.tables:
            raise ValueError("At least one table is required")
        
        # fill parameter can be any string or None - no validation needed
    
    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure inputs from tables."""
        self.folders = pipeline_folders
        
        # Extract information from each table
        self.table_configs = []
        
        for i, table in enumerate(self.tables):
            table_config = self._extract_table_config(table, f"table_{i}")
            self.table_configs.append(table_config)
    
    def _extract_table_config(self, table: Any, prefix: str) -> Dict[str, Any]:
        """Extract configuration from a table reference."""
        config = {
            "prefix": prefix,
            "table_path": None
        }
        
        # Handle different table input formats
        if isinstance(table, str):
            # Direct file path
            config["table_path"] = table
        elif hasattr(table, 'path'):
            # TableInfo object
            config["table_path"] = table.path
        elif hasattr(table, '_tables'):
            # Tool output tables collection - get first one
            for name, info in table._tables.items():
                if hasattr(info, 'path'):
                    config["table_path"] = info.path
                else:
                    config["table_path"] = str(info)
                break
        else:
            # Try to convert to string as fallback
            config["table_path"] = str(table)
        
        return config
    
    def get_config_display(self) -> List[str]:
        """Get configuration display lines."""
        config_lines = super().get_config_display()
        
        config_lines.extend([
            f"DATASHEETS: {len(self.tables)} inputs",
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
            "table_configs": self.table_configs,
            "fill": self.fill,
            "output_csv": concatenated_csv,
            "output_folder": output_folder
        }
        
        import json
        with open(config_file, 'w') as f:
            json.dump(config_data, f, indent=2)
        
        # Generate script content
        script_content = f"""#!/bin/bash
# ConcatenateTables execution script
# Generated by BioPipelines pipeline system

{self.generate_completion_check_header()}

echo "Concatenating tables"
echo "Input tables: {len(self.tables)}"
echo "Fill strategy: {self.fill if self.fill is not None else 'remove non-common columns'}"
echo "Output: {output_folder}"

# Run Python concatenation script
python "{os.path.join(self.folders['HelpScripts'], 'pipe_concatenate_tables.py')}" \\
  --config "{config_file}"

if [ $? -eq 0 ]; then
    echo "Successfully concatenated {len(self.tables)} tables"
    echo "Output written to: {output_folder}"
else
    echo "Error: Failed to concatenate tables"
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
        # ConcatenateTables is a table-only operation like MergeTables
        # It combines tables from multiple tools but doesn't handle structures/compounds
        
        concatenated_csv = os.path.join(self.output_folder, "concatenated.csv")
        
        # Try to determine output columns from input tables
        expected_columns = ["id", "source_table"]
        if self.tables and hasattr(self.tables[0], 'columns'):
            # Try to get columns from first table if available
            try:
                expected_columns = list(self.tables[0].columns) + ["source_table"]
            except:
                # Fallback to default columns
                expected_columns = ["id", "source_table"]
        
        tables = {
            "concatenated": TableInfo(
                name="concatenated", 
                path=concatenated_csv,
                columns=expected_columns,
                description="concatenated results from multiple tool outputs",
                count="variable"
            )
        }
        
        # ConcatenateTables operates only on tables - return empty lists like MergeTables
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
                "fill": self.fill
            }
        })
        return base_dict