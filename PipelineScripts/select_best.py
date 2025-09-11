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
                 input: Union[ToolOutput, StandardizedOutput] = None,
                 pool: Union[ToolOutput, StandardizedOutput, List[Union[ToolOutput, StandardizedOutput]]] = None,
                 data: Union[str, ToolOutput, StandardizedOutput] = None,
                 datasheets: Union[List[Union[ToolOutput, StandardizedOutput]], List[str]] = None,
                 metric: str = None,
                 mode: str = "max",
                 weights: Optional[Dict[str, float]] = None,
                 tie_breaker: str = "first",
                 composite_function: str = "weighted_sum",
                 name: str = "best",
                 **kwargs):
        """
        Initialize SelectBest tool.
        
        Args:
            input: [LEGACY] Input from previous tool (analysis results or filtered data)
            pool: Single or list of tool outputs to select structures from
            data: [LEGACY] Which datasheet to evaluate for selection (string name) or tool output to evaluate
            datasheets: List of datasheets/tool outputs to concatenate for evaluation (new multi-cycle approach)
            metric: Primary metric to optimize for selection
            mode: "max" or "min" the metric
            weights: Dict of {metric_name: weight} for multi-metric selection
            tie_breaker: How to break ties ("first", "random", or metric name)
            composite_function: How to combine metrics ("weighted_sum", "product", "min", "max")
            name: Name for the output structure file (default: "best")
            **kwargs: Additional parameters
            
        Examples:
            # Multi-cycle approach: pools + datasheets
            best = pipeline.add(SelectBest(
                pool=[boltz1.output, boltz2.output],     # List of structure pools to search
                datasheets=[analysis1.output.datasheets.merged, analysis2.output.datasheets.merged],  # List of analysis datasheets
                data="compounds",            # Evaluate compounds datasheet
                metric="binding_affinity",
                mode="max"
            ))
            
            # Legacy approach: keyword input
            best = pipeline.add(SelectBest(
                input=analysis_results,
                metric="binding_affinity",
                mode="max"
            ))
            
            # Multi-objective with weights
            best = pipeline.add(SelectBest(
                pool=combined_results,
                data="structures", 
                metric="composite_score",
                weights={"binding_affinity": 0.6, "pLDDT": 0.4},
                mode="max"
            ))
        """
        # Determine which approach is being used
        if pool is not None and datasheets is not None:
            # New multi-cycle approach: pool + datasheets
            self.use_multi_cycle_mode = True
            self.use_pool_mode = True
            
            # Handle pool - can be single or list
            if isinstance(pool, list):
                self.pool_outputs = pool
            else:
                self.pool_outputs = [pool]
            
            # Handle datasheets - list of datasheets/tool outputs
            self.datasheets = datasheets
            self.data_name = None
            self.data_output = None
            self.selection_input = None
            
            # Set up dependencies
            dependencies = []
            for p in self.pool_outputs:
                if hasattr(p, 'config'):
                    dependencies.append(p.config)
            for ds in self.datasheets:
                if hasattr(ds, 'config'):
                    dependencies.append(ds.config)
            dependency_source = self.pool_outputs[0]
            
        elif pool is not None and data is not None:
            # Legacy pool + data approach
            self.use_multi_cycle_mode = False
            self.use_pool_mode = True
            self.pool_outputs = [pool] if not isinstance(pool, list) else pool
            self.datasheets = None
            
            # Handle data parameter - can be string (datasheet name) or tool output
            if isinstance(data, str):
                self.data_name = data
                self.data_output = None
            else:
                # data is a tool output - use it for evaluation
                self.data_name = None
                self.data_output = data
            
            self.selection_input = None
            
            # Set up dependencies
            dependencies = []
            for p in self.pool_outputs:
                if hasattr(p, 'config'):
                    dependencies.append(p.config)
            if self.data_output and hasattr(data, 'config'):
                dependencies.append(data.config)
            dependency_source = self.pool_outputs[0]
            
        elif input is not None:
            # Legacy input approach (including positional)
            self.use_multi_cycle_mode = False
            self.use_pool_mode = False
            self.pool_outputs = None
            self.datasheets = None
            self.data_name = None
            self.data_output = None
            self.selection_input = input
            dependencies = []
            dependency_source = input
        else:
            raise ValueError("Must specify either 'input' (legacy), 'pool'+'data' (single cycle), or 'pool'+'datasheets' (multi-cycle)")
        
        if metric is None:
            raise ValueError("metric parameter is required")
            
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
        
        # Set up dependencies
        if self.use_pool_mode:
            # Pool mode: add pool(s) and data/datasheet dependencies
            if self.use_multi_cycle_mode:
                # Multi-cycle mode: add all pools and datasheets
                for pool_output in self.pool_outputs:
                    if hasattr(pool_output, 'config'):
                        self.dependencies.append(pool_output.config)
                for datasheet in self.datasheets:
                    if hasattr(datasheet, 'config'):
                        self.dependencies.append(datasheet.config)
            else:
                # Single cycle mode: add pool and data
                for pool_output in self.pool_outputs:
                    if hasattr(pool_output, 'config'):
                        self.dependencies.append(pool_output.config)
                if self.data_output and hasattr(self.data_output, 'config'):
                    self.dependencies.append(self.data_output.config)
        else:
            # Legacy mode
            if hasattr(dependency_source, 'config'):
                self.dependencies.append(dependency_source.config)
    
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
        if self.use_pool_mode:
            # New pool + data mode
            if not isinstance(self.pool_output, (ToolOutput, StandardizedOutput)):
                raise ValueError("pool must be a ToolOutput or StandardizedOutput object")
            
            # Validate data parameter
            if self.data_name:
                # String datasheet name
                if not isinstance(self.data_name, str) or not self.data_name:
                    raise ValueError("data must be a non-empty string specifying the datasheet name")
            elif self.data_output:
                # Tool output for evaluation
                if not isinstance(self.data_output, (ToolOutput, StandardizedOutput)):
                    raise ValueError("data must be a ToolOutput or StandardizedOutput object")
            else:
                raise ValueError("data parameter is required in pool mode")
        else:
            # Legacy input mode
            if not isinstance(self.selection_input, (ToolOutput, StandardizedOutput)):
                raise ValueError("input must be a ToolOutput or StandardizedOutput object")
    
    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input datasheet from previous tool."""
        self.folders = pipeline_folders
        
        if self.use_pool_mode:
            # New pool + data mode: extract specific datasheet for evaluation
            self._configure_pool_mode()
        else:
            # Legacy mode: use existing logic
            self._configure_legacy_mode()
    
    def _configure_pool_mode(self):
        """Configure inputs for pool + data selection mode."""
        # Get the specified datasheet for evaluation
        self.input_csv_path = None
        
        # Check if data parameter is a tool output (cross-tool evaluation)
        if self.data_output:
            # Data comes from a separate tool output
            if hasattr(self.data_output, 'datasheets'):
                datasheets = self.data_output.datasheets
                
                # For tool outputs, use the main/combined datasheet
                if hasattr(datasheets, 'combined'):
                    ds_info = datasheets.combined
                elif hasattr(datasheets, 'filtered'):
                    ds_info = datasheets.filtered
                elif hasattr(datasheets, '_datasheets'):
                    # Get first available datasheet
                    ds_info = next(iter(datasheets._datasheets.values()))
                elif isinstance(datasheets, dict):
                    ds_info = next(iter(datasheets.values()))
                else:
                    ds_info = datasheets
                
                # Extract path from datasheet info
                if hasattr(ds_info, 'path'):
                    self.input_csv_path = ds_info.path
                elif isinstance(ds_info, dict) and 'path' in ds_info:
                    self.input_csv_path = ds_info['path']
                elif isinstance(ds_info, str):
                    self.input_csv_path = ds_info
                else:
                    self.input_csv_path = str(ds_info)
            
            if not self.input_csv_path:
                raise ValueError(f"Could not predict input CSV path from data_output: {self.data_output}. "
                               f"The input tool must provide proper datasheet path predictions.")
                
        else:
            # Data comes from a named datasheet within the pool output
            if hasattr(self.pool_output, 'datasheets'):
                datasheets = self.pool_output.datasheets
                
                # Look for the specified datasheet name
                if hasattr(datasheets, self.data_name):
                    ds_info = getattr(datasheets, self.data_name)
                    if hasattr(ds_info, 'path'):
                        self.input_csv_path = ds_info.path
                    elif isinstance(ds_info, str):
                        self.input_csv_path = ds_info
                    else:
                        self.input_csv_path = str(ds_info)
                elif isinstance(datasheets, dict) and self.data_name in datasheets:
                    ds_info = datasheets[self.data_name]
                    if isinstance(ds_info, dict) and 'path' in ds_info:
                        self.input_csv_path = ds_info['path']
                    else:
                        self.input_csv_path = str(ds_info)
            
            # Fallback: predict path based on pool output folder and data name
            if not self.input_csv_path and hasattr(self.pool_output, 'output_folder'):
                output_folder = self.pool_output.output_folder
                predicted_name = f"{self.data_name}_results.csv"
                self.input_csv_path = os.path.join(output_folder, predicted_name)
        
        # Validation
        if not self.input_csv_path:
            if self.data_output:
                raise ValueError(f"Could not predict input CSV path from data tool output: {self.data_output}")
            else:
                available_datasheets = []
                if hasattr(self.pool_output, 'datasheets'):
                    if hasattr(self.pool_output.datasheets, '_datasheets'):
                        available_datasheets = list(self.pool_output.datasheets._datasheets.keys())
                    elif isinstance(self.pool_output.datasheets, dict):
                        available_datasheets = list(self.pool_output.datasheets.keys())
                
                raise ValueError(f"Could not find datasheet '{self.data_name}' in pool output. "
                               f"Available datasheets: {available_datasheets}")
        
        # Store the full pool output for later data extraction
        self.pool_source = self.pool_output
        
        # Predict source data directories from pool output
        self.source_structures_dir = None
        if hasattr(self.pool_output, 'output_folder'):
            potential_dirs = [
                self.pool_output.output_folder,
                os.path.join(self.pool_output.output_folder, 'structures'),
                os.path.join(self.pool_output.output_folder, 'pdbs')
            ]
            self.source_structures_dir = potential_dirs[0]
    
    def _configure_legacy_mode(self):
        """Configure inputs for legacy input mode."""
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
        
        if not self.input_csv_path:
            raise ValueError(f"Could not predict input CSV path from previous tool: {self.selection_input}. "
                           f"The input tool must provide proper datasheet path predictions.")
        
        # Store for legacy compatibility
        self.pool_source = self.selection_input
        
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
    
    def _get_input_columns(self) -> List[str]:
        """Get column names from the input datasheet."""
        if self.use_pool_mode:
            # Pool mode: get columns from data source
            if self.data_output:
                # Data comes from a separate tool output
                if hasattr(self.data_output, 'datasheets') and hasattr(self.data_output.datasheets, '_datasheets'):
                    # Look through datasheets to find one with column info
                    for name, info in self.data_output.datasheets._datasheets.items():
                        if hasattr(info, 'columns') and info.columns:
                            return info.columns
            else:
                # Data comes from a named datasheet within the pool output
                if hasattr(self.pool_output, 'datasheets') and hasattr(self.pool_output.datasheets, '_datasheets'):
                    if self.data_name in self.pool_output.datasheets._datasheets:
                        info = self.pool_output.datasheets._datasheets[self.data_name]
                        if hasattr(info, 'columns') and info.columns:
                            return info.columns
        else:
            # Legacy mode: get columns from selection input
            if hasattr(self.selection_input, 'datasheets') and hasattr(self.selection_input.datasheets, '_datasheets'):
                # Look through the datasheets to find one with column info
                for name, info in self.selection_input.datasheets._datasheets.items():
                    if hasattr(info, 'columns') and info.columns:
                        return info.columns
        
        # Return empty list if columns cannot be determined
        return []
    
    def get_config_display(self) -> List[str]:
        """Get configuration display lines."""
        config_lines = super().get_config_display()
        
        if self.use_pool_mode:
            if self.data_output:
                config_lines.extend([
                    f"MODE: Cross-tool Pool + Data Selection",
                    f"POOL SOURCE: Tool Output",
                    f"DATA SOURCE: Tool Output",
                ])
            else:
                config_lines.extend([
                    f"MODE: Pool + Data Selection",
                    f"DATA SOURCE: {self.data_name}",
                ])
        else:
            config_lines.append("MODE: Legacy Input Selection")
        
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
        config_file = os.path.join(output_folder, "select_best_config.json")
        config_data = {
            "input_csv": self.input_csv_path,
            "source_structures_dir": self.source_structures_dir,
            "selection_metric": self.metric,
            "selection_mode": self.mode,
            "weights": self.weights,
            "tie_breaker": self.tie_breaker,
            "composite_function": self.composite_function,
            "output_csv": selected_csv,
            "output_structure": selected_structure,
            # New pool mode configuration
            "use_pool_mode": self.use_pool_mode,
            "data_name": self.data_name if self.use_pool_mode else None,
            "pool_output_folder": getattr(self.pool_source, 'output_folder', None) if hasattr(self, 'pool_source') else None
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
            Dictionary with output file paths for selected item
        """
        if self.use_pool_mode:
            return self._get_pool_mode_outputs()
        else:
            return self._get_legacy_mode_outputs()
    
    def _get_pool_mode_outputs(self) -> Dict[str, List[str]]:
        """Get outputs for pool mode - includes all data types for selected ID."""
        selected_csv = os.path.join(self.output_folder, "selected_best.csv")
        selected_structure = os.path.join(self.output_folder, f"{self.output_name}.pdb")
        selected_compound = os.path.join(self.output_folder, f"{self.output_name}_compound.sdf")
        selected_sequence = os.path.join(self.output_folder, f"{self.output_name}_sequences.csv")
        
        # Use output name as the selected ID
        selected_id = self.output_name
        
        # Get actual column names from evaluation datasheet
        input_columns = self._get_input_columns()
        
        # Define datasheets that will be created
        datasheets = {
            "selected": {
                "path": selected_csv,
                "columns": input_columns,  # Use actual columns
                "description": f"Best item selected from {self.data_name} using {self.metric}",
                "count": 1
            }
        }
        
        return {
            "structures": [selected_structure],
            "structure_ids": [selected_id],
            "compounds": [selected_compound], 
            "compound_ids": [selected_id],
            "sequences": [selected_sequence],
            "sequence_ids": [selected_id],
            "datasheets": datasheets,
            "output_folder": self.output_folder
        }
    
    def _get_legacy_mode_outputs(self) -> Dict[str, List[str]]:
        """Get outputs for legacy mode - single structure + datasheet only."""
        selected_structure = os.path.join(self.output_folder, f"{self.output_name}.pdb")
        selected_csv = os.path.join(self.output_folder, "selected_best.csv")
        
        # Structure ID should match the structure name without extension
        structure_id = self.output_name
        
        # Get actual column names from input
        input_columns = self._get_input_columns()
        
        datasheets = {
            "selected": {
                "path": selected_csv,
                "columns": input_columns,  # Use actual columns
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
                "use_pool_mode": self.use_pool_mode,
                "data_name": self.data_name if self.use_pool_mode else None,
                "metric": self.metric,
                "mode": self.mode,
                "weights": self.weights,
                "tie_breaker": self.tie_breaker,
                "composite_function": self.composite_function,
                "output_name": self.output_name
            }
        })
        return base_dict