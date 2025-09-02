"""
SelectBest tool for choosing the single best item from analysis results.

Takes analysis datasheets and applies selection criteria to pick exactly one
best structure/sequence/compound for use in iterative optimization cycles.
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


class SelectBest(BaseConfig):
    """
    Tool for selecting the single best item from analysis results.
    
    Takes analysis datasheets and applies selection criteria to return exactly
    one best structure/sequence/compound with its associated data.
    """
    
    TOOL_NAME = "SelectBest"
    DEFAULT_ENV = "ProteinEnv"
    COMPATIBLE_ENVS = ["ProteinEnv", "Boltz2Env", "ligandmpnn_env"]
    DEFAULT_RESOURCES = {"gpu": "T4", "memory": "4GB", "time": "1:00:00"}
    
    def __init__(self,
                 input: Union[ToolOutput, StandardizedOutput],
                 metric: str,
                 mode: str = "max",
                 weights: Optional[Dict[str, float]] = None,
                 tie_breaker: str = "first",
                 composite_function: str = "weighted_sum",
                 name: str = "best",
                 **kwargs):
        """
        Initialize SelectBest tool.
        
        Args:
            input: Input from previous tool (analysis results or filtered data)
            metric: Primary metric to optimize for selection
            mode: "max" or "min" the metric
            weights: Dict of {metric_name: weight} for multi-metric selection
            tie_breaker: How to break ties ("first", "random", or metric name)
            composite_function: How to combine metrics ("weighted_sum", "product", "min", "max")
            name: Name for the output structure file (default: "best")
            **kwargs: Additional parameters
            
        Examples:
            # Simple: best binding affinity
            best = pipeline.add(SelectBest(
                input=analysis_results,
                metric="binding_affinity",
                mode="max"
            ))
            
            # Multi-objective with weights
            best = pipeline.add(SelectBest(
                input=combined_results,
                metric="composite_score",
                weights={"binding_affinity": 0.6, "pLDDT": 0.4},
                mode="max"
            ))
            
            # Custom tie breaking
            best = pipeline.add(SelectBest(
                input=filtered_data,
                metric="energy",
                mode="min",
                tie_breaker="pLDDT"  # Use pLDDT to break ties
            ))
        """
        self.selection_input = input
        self.metric = metric
        self.mode = mode
        self.weights = weights or {}
        self.tie_breaker = tie_breaker
        self.composite_function = composite_function
        self.output_name = name
        
        # Validate parameters
        self._validate_selection_params()
        
        # Initialize base class
        super().__init__(**kwargs)
        
        # Set up dependency
        if hasattr(input, 'config'):
            self.dependencies.append(input.config)
    
    def _validate_selection_params(self):
        """Validate selection parameters."""
        if self.mode not in ["max", "min"]:
            raise ValueError("mode must be 'max' or 'min'")
        
        if self.tie_breaker not in ["first", "random"] and not isinstance(self.tie_breaker, str):
            raise ValueError("tie_breaker must be 'first', 'random', or a metric name")
        
        if self.composite_function not in ["weighted_sum", "product", "min", "max"]:
            raise ValueError("composite_function must be one of: weighted_sum, product, min, max")
        
        if not self.metric:
            raise ValueError("metric cannot be empty")
    
    def validate_params(self):
        """Validate SelectBest parameters."""
        if not isinstance(self.selection_input, (ToolOutput, StandardizedOutput)):
            raise ValueError("Input must be a ToolOutput or StandardizedOutput object")
    
    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input datasheet from previous tool."""
        self.folders = pipeline_folders
        
        # Get the input CSV file from the previous tool
        self.input_csv_path = None
        
        if hasattr(self.selection_input, 'datasheets'):
            datasheets = self.selection_input.datasheets
            if isinstance(datasheets, dict):
                # Find the main datasheet (look for 'combined', 'filtered', or first available)
                priority_names = ['combined', 'filtered', 'analysis', 'structures']
                
                for name in priority_names:
                    if name in datasheets:
                        ds_info = datasheets[name]
                        break
                else:
                    # Take first available
                    ds_info = next(iter(datasheets.values()))
                
                if isinstance(ds_info, dict) and 'path' in ds_info:
                    self.input_csv_path = ds_info['path']
                else:
                    self.input_csv_path = str(ds_info)
        
        # Fallback: predict path based on output folder
        if not self.input_csv_path and hasattr(self.selection_input, 'output_folder'):
            output_folder = self.selection_input.output_folder
            # Predict common CSV names that tools would generate (don't check existence)
            common_names = [
                'filtered_results.csv',     # Filter output
                'combined_analysis.csv',    # MergeDatasheets output  
                'analysis_results.csv',     # Analysis output
                'results.csv'               # Generic results
            ]
            
            # Use first predicted path (don't check existence)
            self.input_csv_path = os.path.join(output_folder, common_names[0])
        
        if not self.input_csv_path:
            raise ValueError(f"Could not predict input CSV path from previous tool: {self.selection_input}")
        
        # Also predict source structures directory
        self.source_structures_dir = None
        if hasattr(self.selection_input, 'output_folder'):
            # Predict common structure directory locations (don't check existence)
            potential_dirs = [
                self.selection_input.output_folder,
                os.path.join(self.selection_input.output_folder, 'structures'),
                os.path.join(self.selection_input.output_folder, 'pdbs')
            ]
            
            # Use first predicted directory (don't verify contents)
            self.source_structures_dir = potential_dirs[0]
    
    def get_config_display(self) -> List[str]:
        """Get configuration display lines."""
        config_lines = super().get_config_display()
        
        config_lines.extend([
            f"SELECTION METRIC: {self.metric}",
            f"SELECTION MODE: {self.mode}",
            f"TIE BREAKER: {self.tie_breaker}",
        ])
        
        if self.weights:
            weights_str = ", ".join(f"{k}:{v}" for k, v in self.weights.items())
            config_lines.extend([
                f"COMPOSITE FUNCTION: {self.composite_function}",
                f"WEIGHTS: {weights_str}"
            ])
        
        return config_lines
    
    def generate_script(self, script_path: str) -> str:
        """
        Generate SelectBest execution script.
        
        Args:
            script_path: Path where script should be written
            
        Returns:
            Script content as string
        """
        runtime_folder = os.path.dirname(script_path)
        output_folder = self.output_folder
        os.makedirs(output_folder, exist_ok=True)
        
        # Output files
        selected_csv = os.path.join(output_folder, "selected_best.csv")
        selected_structure = os.path.join(output_folder, f"{self.output_name}.pdb")
        
        # Create config file for selection
        config_file = os.path.join(runtime_folder, "select_best_config.json")
        config_data = {
            "input_csv": self.input_csv_path,
            "source_structures_dir": self.source_structures_dir,
            "selection_metric": self.metric,
            "selection_mode": self.mode,
            "weights": self.weights,
            "tie_breaker": self.tie_breaker,
            "composite_function": self.composite_function,
            "output_csv": selected_csv,
            "output_structure": selected_structure
        }
        
        import json
        with open(config_file, 'w') as f:
            json.dump(config_data, f, indent=2)
        
        # Generate script content
        script_content = f"""#!/bin/bash
# SelectBest execution script
# Generated by BioPipelines pipeline system

{self.generate_completion_check_header()}

echo "Selecting best item from analysis results"
echo "Input: {self.input_csv_path}"
echo "Selection metric: {self.metric} ({self.mode})"
echo "Output: {selected_structure}"

# Run Python selection script
python "{os.path.join(self.folders['HelpScripts'], 'pipe_select_best.py')}" \\
  --config "{config_file}"

if [ $? -eq 0 ]; then
    echo "Selection completed successfully"
    echo "Best structure: {selected_structure}"
    echo "Best data: {selected_csv}"
else
    echo "Error: Selection failed"
    exit 1
fi

{self.generate_completion_check_footer()}
"""
        
        return script_content
    
    def get_output_files(self) -> Dict[str, List[str]]:
        """
        Get expected output files after selection.
        
        Returns:
            Dictionary with output file paths (single structure + datasheet)
        """
        selected_structure = os.path.join(self.output_folder, f"{self.output_name}.pdb")
        selected_csv = os.path.join(self.output_folder, "selected_best.csv")
        
        # Structure ID should match the structure name without extension
        structure_id = self.output_name
        
        datasheets = {
            "selected": {
                "path": selected_csv,
                "columns": "preserved_from_input",
                "description": f"Single best item selected using {self.metric}",
                "count": 1
            }
        }
        
        return {
            "structures": [selected_structure],
            "structure_ids": [structure_id],
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
                "metric": self.metric,
                "mode": self.mode,
                "weights": self.weights,
                "tie_breaker": self.tie_breaker,
                "composite_function": self.composite_function,
                "output_name": self.output_name
            }
        })
        return base_dict