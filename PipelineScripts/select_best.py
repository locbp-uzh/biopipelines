"""
SelectBest tool for choosing the single best item from analysis results.

Takes analysis datasheets and applies selection criteria to pick exactly one
best structure/sequence/compound for use in iterative optimization cycles.
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


class SelectBest(BaseConfig):
    """
    Tool for selecting the single best item from analysis results.
    
    Takes analysis datasheets and applies selection criteria to return exactly
    one best structure/sequence/compound with its associated data.
    """
    
    TOOL_NAME = "SelectBest"
    DEFAULT_ENV = None  # Loaded from config.yaml
    
    def __init__(self,
                 pool: Union[ToolOutput, StandardizedOutput, List[Union[ToolOutput, StandardizedOutput]]],
                 datasheets: Union[List[Union[ToolOutput, StandardizedOutput, DatasheetInfo, str]], List[str]],
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
            pool: Single or list of tool outputs to select structures from
            datasheets: Single or list of datasheets to evaluate for selection
            metric: Primary metric to optimize for selection
            mode: "max" or "min" the metric
            weights: Dict of {metric_name: weight} for multi-metric selection
            tie_breaker: How to break ties ("first", "random", or metric name)
            composite_function: How to combine metrics ("weighted_sum", "product", "min", "max")
            name: Name for the output structure file (default: "best")
            **kwargs: Additional parameters
            
        Examples:
            # Compare across multiple pools (e.g., engineered vs original)
            best = pipeline.add(SelectBest(
                pool=[original, engineered],
                datasheets=[original_analysis.datasheets.merged, engineered_analysis.datasheets.merged],
                metric="binding_affinity",
                mode="max"
            ))
            
            # Single pool selection
            best = pipeline.add(SelectBest(
                pool=boltz,
                datasheets=analysis.datasheets.merged,
                metric="pLDDT",
                mode="max"
            ))
            
            # Multi-objective with weights
            best = pipeline.add(SelectBest(
                pool=combined_results,
                datasheets=analysis.datasheets.merged,
                metric="composite_score",
                weights={"binding_affinity": 0.6, "pLDDT": 0.4},
                mode="max"
            ))
        """
        # Handle pool - can be single or list, always convert to list
        if isinstance(pool, list):
            self.pool_outputs = pool
        else:
            self.pool_outputs = [pool]
        
        # Handle datasheets - always expect list, convert if single
        if isinstance(datasheets, list):
            self.datasheets = datasheets
        else:
            self.datasheets = [datasheets]
        
        # Set up dependencies
        dependencies = []
        for p in self.pool_outputs:
            if hasattr(p, 'config'):
                dependencies.append(p.config)
        for ds in self.datasheets:
            if hasattr(ds, 'config'):
                dependencies.append(ds.config)
            
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
        
        # Set up dependencies (already handled above in each case)
        for dep in dependencies:
            self.dependencies.append(dep)
    
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
        # Validate pool outputs
        if not self.pool_outputs:
            raise ValueError("pool must be provided")
        
        for i, pool_output in enumerate(self.pool_outputs):
            if not isinstance(pool_output, (ToolOutput, StandardizedOutput)):
                raise ValueError(f"pool[{i}] is {type(pool_output)}, must be a ToolOutput or StandardizedOutput object. Value: {pool_output}")
        
        # Validate datasheets
        if not self.datasheets:
            raise ValueError("datasheets must be provided")
        
        if len(self.datasheets) != len(self.pool_outputs):
            raise ValueError(f"Number of datasheets ({len(self.datasheets)}) must match number of pools ({len(self.pool_outputs)})")
        
        for i, datasheet in enumerate(self.datasheets):
            if not isinstance(datasheet, (ToolOutput, StandardizedOutput, DatasheetInfo, str)):
                raise ValueError(f"datasheet[{i}] is {type(datasheet)}, must be a ToolOutput, StandardizedOutput, DatasheetInfo, or str object. Value: {datasheet}")
    
    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input datasheet from previous tool."""
        self.folders = pipeline_folders
        self._configure_pool_mode()
    
    def _configure_pool_mode(self):
        """Configure inputs for pool + datasheets selection."""
        # Extract pool folders and datasheet paths
        self.pool_folders = []
        self.datasheet_paths = []
        
        for pool in self.pool_outputs:
            if hasattr(pool, 'output_folder'):
                self.pool_folders.append(pool.output_folder)
            else:
                raise ValueError(f"Pool {pool} must have output_folder")
        
        for datasheet in self.datasheets:
            # Handle different datasheet input types - check DatasheetInfo first
            if hasattr(datasheet, 'path') and hasattr(datasheet, 'name') and hasattr(datasheet, 'columns'):
                # DatasheetInfo object - use path directly
                self.datasheet_paths.append(datasheet.path)
            elif isinstance(datasheet, str):
                # String path - use directly
                self.datasheet_paths.append(datasheet)
            elif hasattr(datasheet, 'datasheets'):
                # ToolOutput object - extract datasheet path
                ds_obj = datasheet.datasheets
                datasheet_path = None

                if hasattr(ds_obj, 'merged'):
                    datasheet_path = ds_obj.merged.path
                elif hasattr(ds_obj, 'combined'):
                    datasheet_path = ds_obj.combined.path
                elif hasattr(ds_obj, 'filtered'):
                    datasheet_path = ds_obj.filtered.path
                elif hasattr(ds_obj, '_datasheets'):
                    # Get first available datasheet
                    first_ds = next(iter(ds_obj._datasheets.values()))
                    datasheet_path = first_ds.path
                else:
                    raise ValueError(f"Could not find datasheet in {datasheet}")

                self.datasheet_paths.append(datasheet_path)
            else:
                raise ValueError(f"Unsupported datasheet type: {type(datasheet)}")
    
    def get_config_display(self) -> List[str]:
        """Get configuration display lines."""
        config_lines = super().get_config_display()
        
        config_lines.extend([
            f"POOLS: {len(self.pool_outputs)}",
            f"DATASHEETS: {len(self.datasheets)}",
            f"METRIC: {self.metric} ({self.mode})"
        ])
        
        return config_lines
    
    def generate_script(self, script_path: str) -> str:
        """
        Generate script to select the best item from analysis results.
        
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
            "selection_metric": self.metric,
            "selection_mode": self.mode,
            "weights": self.weights,
            "tie_breaker": self.tie_breaker,
            "composite_function": self.composite_function,
            "output_csv": selected_csv,
            "output_structure": selected_structure,
            "output_sequences": os.path.join(output_folder, "sequences.csv"),
            "pool_folders": self.pool_folders,
            "datasheet_paths": self.datasheet_paths
        }
        
        import json
        with open(config_file, 'w') as f:
            json.dump(config_data, f, indent=2)
        
        # Generate script content
        script_content = f"""#!/bin/bash
# SelectBest execution script
# Generated by BioPipelines pipeline system

{self.generate_completion_check_header()}

echo "Selecting best item from {len(self.datasheet_paths)} datasheets"
echo "Selection metric: {self.metric} ({self.mode})"
echo "Output: {selected_structure}"

# Run Python selection script
python "{os.path.join(self.folders['HelpScripts'], 'pipe_select_best.py')}" \\
  --config "{config_file}"

if [ $? -eq 0 ]; then
    echo "Successfully selected best item"
    echo "Selected structure: {selected_structure}"
    echo "Selected data: {selected_csv}"
else
    echo "Error: Failed to select best item"
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
        selected_csv = os.path.join(self.output_folder, "selected_best.csv")
        selected_structure = os.path.join(self.output_folder, f"{self.output_name}.pdb")
        selected_compound = os.path.join(self.output_folder, f"{self.output_name}_compound.sdf")
        selected_sequence = os.path.join(self.output_folder, "sequences.csv")
        
        # Use output name as the selected ID
        selected_id = self.output_name
        
        # Define datasheets that will be created
        datasheets = {
            "selected": {
                "path": selected_csv,
                "columns": ["id", self.metric],  # Basic columns
                "description": f"Best item selected using {self.metric}",
                "count": 1
            },
            "sequences": {
                "path": selected_sequence,
                "columns": ["id", "sequence"],  # Standard sequence format
                "description": "Selected best sequence",
                "count": 1
            }
        }
        
        # Only predict files we're certain will be created
        # Structures are always created (extracted from pools)
        # Compounds and sequences are only created if they exist in source pools
        return {
            "structures": [selected_structure],
            "structure_ids": [selected_id],
            "compounds": [],  # Only created if available in source pools
            "compound_ids": [],
            "sequences": [selected_sequence],  # Include sequences CSV  
            "sequence_ids": [selected_id],
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
