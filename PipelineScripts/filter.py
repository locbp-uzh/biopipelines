"""
Filter tool for expression-based filtering of unified analysis datasheets.

Takes unified datasheets from CombineDatasheets and applies pandas query-style
expressions to filter rows while preserving all column information.
"""

import os
import pandas as pd
from typing import Dict, List, Any, Optional, Union

try:
    from .base_config import BaseConfig, ToolOutput, StandardizedOutput
except ImportError:
    # Fallback for direct execution
    import sys
    import os
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, ToolOutput, StandardizedOutput


class FilterResult:
    """
    Result from applying expression-based filtering.
    """
    
    def __init__(self, 
                 input_csv: str,
                 output_csv: str,
                 expression: str,
                 total_input: int,
                 kept_count: int,
                 filtered_count: int):
        """
        Initialize filter result.
        
        Args:
            input_csv: Path to input CSV file
            output_csv: Path to output filtered CSV file  
            expression: Filter expression that was applied
            total_input: Total number of input rows
            kept_count: Number of rows that passed the filter
            filtered_count: Number of rows that were filtered out
        """
        self.input_csv = input_csv
        self.output_csv = output_csv
        self.expression = expression
        self.total_input = total_input
        self.kept_count = kept_count
        self.filtered_count = filtered_count
        self.pass_rate = kept_count / total_input if total_input > 0 else 0.0
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for serialization."""
        return {
            "input_csv": self.input_csv,
            "output_csv": self.output_csv,
            "expression": self.expression,
            "total_input": self.total_input,
            "kept_count": self.kept_count,
            "filtered_count": self.filtered_count,
            "pass_rate": self.pass_rate
        }
    
    def summary(self) -> str:
        """Get human-readable summary."""
        return (f"Filter ({self.expression}): {self.kept_count}/{self.total_input} items kept "
                f"({self.pass_rate:.1%} pass rate)")


class Filter(BaseConfig):
    """
    Expression-based filter tool for unified analysis datasheets.
    
    Takes a unified datasheet (typically from CombineDatasheets) and applies
    pandas query-style expressions to filter rows while preserving all columns.
    """
    
    # Tool identification
    TOOL_NAME = "Filter"
    DEFAULT_ENV = "ProteinEnv"
    COMPATIBLE_ENVS = ["ProteinEnv", "Boltz2Env", "ligandmpnn_env"]
    DEFAULT_RESOURCES = {"gpu": "T4", "memory": "4GB", "time": "1:00:00"}
    
    def __init__(self,
                 input: Union[ToolOutput, StandardizedOutput],
                 expression: str,
                 max_items: Optional[int] = None,
                 sort_by: Optional[str] = None,
                 sort_ascending: bool = True,
                 **kwargs):
        """
        Initialize Filter tool.
        
        Args:
            input: Input from previous tool (typically CombineDatasheets)
            expression: Pandas query-style expression (e.g., "pLDDT>80 and distance < 5.0")
            max_items: Maximum number of items to keep after filtering
            sort_by: Column name to sort by before applying max_items limit
            sort_ascending: Sort order (True for ascending, False for descending)
            **kwargs: Additional parameters
            
        Examples:
            # Simple filtering
            filtered = pipeline.add(Filter(
                input=combined_analysis,
                expression="pLDDT > 80"
            ))
            
            # Complex filtering with multiple conditions
            filtered = pipeline.add(Filter(
                input=combined_analysis,
                expression="pLDDT > 80 and distance < 5.0 and confidence > 0.9"
            ))
            
            # Filtering with item limit and sorting
            filtered = pipeline.add(Filter(
                input=combined_analysis,
                expression="pLDDT > 70",
                max_items=10,
                sort_by="pLDDT",
                sort_ascending=False  # Best pLDDT first
            ))
        """
        self.filter_input = input
        self.expression = expression
        self.max_items = max_items
        self.sort_by = sort_by
        self.sort_ascending = sort_ascending
        
        # Validate expression
        self._validate_expression()
        
        # Initialize base class
        super().__init__(**kwargs)
        
        # Set up dependency
        if hasattr(input, 'config'):
            self.dependencies.append(input.config)
    
    def _validate_expression(self):
        """Validate that the expression is safe for pandas query."""
        if not self.expression.strip():
            raise ValueError("Filter expression cannot be empty")
        
        # Basic safety check - ensure only safe characters and operations
        import re
        allowed_pattern = r'^[a-zA-Z_][a-zA-Z0-9_\s\.\+\-\*\/\(\)<>=!&|and\sor\snot\s\d]+$'
        
        if not re.match(allowed_pattern, self.expression):
            raise ValueError(f"Invalid characters in expression: {self.expression}")
        
        # Check for dangerous keywords
        dangerous_keywords = ['import', 'exec', 'eval', '__', 'os.', 'sys.', 'subprocess']
        expr_lower = self.expression.lower()
        for keyword in dangerous_keywords:
            if keyword in expr_lower:
                raise ValueError(f"Dangerous keyword '{keyword}' not allowed in expression")
    
    def validate_params(self):
        """Validate Filter parameters."""
        if not isinstance(self.filter_input, (ToolOutput, StandardizedOutput)):
            raise ValueError("Input must be a ToolOutput or StandardizedOutput object")
        
        if self.max_items is not None and self.max_items <= 0:
            raise ValueError("max_items must be positive")
    
    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input datasheet from previous tool."""
        self.folders = pipeline_folders
        
        # Get the input CSV file from the previous tool
        self.input_csv_path = None
        
        if hasattr(self.filter_input, 'datasheets'):
            datasheets = self.filter_input.datasheets
            if isinstance(datasheets, dict):
                # Find the main datasheet (look for 'combined' or first available)
                if 'combined' in datasheets:
                    ds_info = datasheets['combined']
                else:
                    ds_info = next(iter(datasheets.values()))
                
                if isinstance(ds_info, dict) and 'path' in ds_info:
                    self.input_csv_path = ds_info['path']
                else:
                    self.input_csv_path = str(ds_info)
        
        # Fallback: predict path based on output folder
        if not self.input_csv_path and hasattr(self.filter_input, 'output_folder'):
            output_folder = self.filter_input.output_folder
            # Predict common CSV names that tools would generate (don't check existence)
            common_names = [
                'combined_analysis.csv',    # MergeDatasheets output
                'analysis_results.csv',     # General analysis output
                'filtered_results.csv',     # Previous filter output
                'results.csv'               # Generic results
            ]
            
            # Use first predicted path (don't check existence)
            self.input_csv_path = os.path.join(output_folder, common_names[0])
        
        if not self.input_csv_path:
            raise ValueError(f"Could not predict input CSV path from previous tool: {self.filter_input}")
    
    def get_config_display(self) -> List[str]:
        """Get configuration display lines."""
        config_lines = super().get_config_display()
        
        config_lines.extend([
            f"EXPRESSION: {self.expression}",
            f"MAX ITEMS: {self.max_items if self.max_items else 'unlimited'}",
        ])
        
        if self.sort_by:
            order = "ascending" if self.sort_ascending else "descending"
            config_lines.append(f"SORT BY: {self.sort_by} ({order})")
        
        return config_lines
    
    def generate_script(self, script_path: str) -> str:
        """
        Generate filter execution script.
        
        Args:
            script_path: Path where script should be written
            
        Returns:
            Script content as string
        """
        runtime_folder = os.path.dirname(script_path)
        output_folder = self.output_folder
        os.makedirs(output_folder, exist_ok=True)
        
        # Output CSV path
        filtered_csv = os.path.join(output_folder, "filtered_results.csv")
        
        # Create config file for the filter
        config_file = os.path.join(runtime_folder, "filter_config.json")
        config_data = {
            "input_csv": self.input_csv_path,
            "expression": self.expression,
            "max_items": self.max_items,
            "sort_by": self.sort_by,
            "sort_ascending": self.sort_ascending,
            "output_csv": filtered_csv
        }
        
        import json
        with open(config_file, 'w') as f:
            json.dump(config_data, f, indent=2)
        
        # Generate script content
        script_content = f"""#!/bin/bash
# Filter execution script
# Generated by BioPipelines pipeline system

{self.generate_completion_check_header()}

echo "Applying filter expression"
echo "Input: {self.input_csv_path}"
echo "Expression: {self.expression}"
echo "Output: {filtered_csv}"

# Run Python filtering script
python "{os.path.join(self.folders['HelpScripts'], 'pipe_filter.py')}" \\
  --config "{config_file}"

if [ $? -eq 0 ]; then
    echo "Filtering completed successfully"
    echo "Results written to: {filtered_csv}"
else
    echo "Error: Filtering failed"
    exit 1
fi

{self.generate_completion_check_footer()}
"""
        
        return script_content
    
    def get_output_files(self) -> Dict[str, List[str]]:
        """
        Get expected output files after filtering.
        
        Returns:
            Dictionary with output file paths
        """
        filtered_csv = os.path.join(self.output_folder, "filtered_results.csv")
        
        datasheets = {
            "filtered": {
                "path": filtered_csv,
                "columns": "preserved_from_input",
                "description": f"Filtered results using expression: {self.expression}",
                "count": "variable"
            }
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
                "expression": self.expression,
                "max_items": self.max_items,
                "sort_by": self.sort_by,
                "sort_ascending": self.sort_ascending
            }
        })
        return base_dict