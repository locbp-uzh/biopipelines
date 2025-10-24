"""
SelectBest tool for choosing the single best item from analysis results (REFACTORED).

Takes analysis datasheets and applies selection criteria to pick exactly one
best structure/sequence/compound for use in iterative optimization cycles.

REFACTORED VERSION using DatasheetNavigatorMixin for cleaner datasheet access.
"""

import os
import pandas as pd
from typing import Dict, List, Any, Optional, Union

try:
    from .base_config import BaseConfig, ToolOutput, StandardizedOutput, DatasheetInfo
    from .mixins import DatasheetNavigatorMixin, FilePathDescriptor
except ImportError:
    import sys
    import os
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, ToolOutput, StandardizedOutput, DatasheetInfo
    from mixins import DatasheetNavigatorMixin, FilePathDescriptor


class SelectBest(DatasheetNavigatorMixin, BaseConfig):
    """
    Tool for selecting the single best item from analysis results.

    REFACTORED to use DatasheetNavigatorMixin for elegant datasheet access.
    """

    TOOL_NAME = "SelectBest"
    DEFAULT_ENV = None

    # ============================================================================
    # AUTOMATIC FILE PATH MANAGEMENT
    # ============================================================================
    selected_csv = FilePathDescriptor("selected_best.csv")
    config_file = FilePathDescriptor("select_best_config.json")
    select_script = FilePathDescriptor("pipe_select_best.py", folder_key="HelpScripts")

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
        """Initialize SelectBest tool."""
        # Handle pool - can be single or list
        if isinstance(pool, list):
            self.pool_outputs = pool
        else:
            self.pool_outputs = [pool]

        # Handle datasheets - always expect list
        if isinstance(datasheets, list):
            self.datasheets = datasheets
        else:
            self.datasheets = [datasheets]

        # Store parameters
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
        for p in self.pool_outputs:
            if hasattr(p, 'config'):
                self.dependencies.append(p.config)
        for ds in self.datasheets:
            if hasattr(ds, 'config'):
                self.dependencies.append(ds.config)

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
        if not self.pool_outputs:
            raise ValueError("pool must be provided")

        for i, pool_output in enumerate(self.pool_outputs):
            if not isinstance(pool_output, (ToolOutput, StandardizedOutput)):
                raise ValueError(
                    f"pool[{i}] is {type(pool_output)}, must be ToolOutput or StandardizedOutput"
                )

        if not self.datasheets:
            raise ValueError("datasheets must be provided")

        if len(self.datasheets) != len(self.pool_outputs):
            raise ValueError(
                f"Number of datasheets ({len(self.datasheets)}) must match "
                f"number of pools ({len(self.pool_outputs)})"
            )

        for i, datasheet in enumerate(self.datasheets):
            if not isinstance(datasheet, (ToolOutput, StandardizedOutput, DatasheetInfo, str)):
                raise ValueError(
                    f"datasheet[{i}] is {type(datasheet)}, must be ToolOutput, "
                    f"StandardizedOutput, DatasheetInfo, or str"
                )

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """
        Configure input datasheet from previous tool.

        REFACTORED: Uses DatasheetNavigatorMixin for elegant datasheet access.
        OLD: 45+ lines of nested if-else
        NEW: 20 lines with get_datasheet_path() mixin!
        """
        self.folders = pipeline_folders
        self._configure_pool_mode()

    def _configure_pool_mode(self):
        """
        Configure inputs for pool + datasheets selection.

        REFACTORED: Uses get_datasheet_path() for cleaner datasheet extraction.
        """
        # Extract pool folders
        self.pool_folders = []
        for pool in self.pool_outputs:
            if hasattr(pool, 'output_folder'):
                self.pool_folders.append(pool.output_folder)
            else:
                raise ValueError(f"Pool {pool} must have output_folder")

        # ========================================================================
        # ELEGANT DATASHEET PATH EXTRACTION using mixin
        # ========================================================================
        self.datasheet_paths = []
        for datasheet in self.datasheets:
            # Direct string path
            if isinstance(datasheet, str):
                self.datasheet_paths.append(datasheet)
                continue

            # DatasheetInfo object
            if hasattr(datasheet, 'path') and hasattr(datasheet, 'name') and hasattr(datasheet, 'columns'):
                self.datasheet_paths.append(datasheet.path)
                continue

            # ToolOutput or StandardizedOutput - use mixin!
            if hasattr(datasheet, 'datasheets'):
                try:
                    # Try to get common datasheet names
                    path = self.get_datasheet_path(
                        datasheet,
                        name=None,  # No preference
                        fallback_names=['merged', 'combined', 'filtered', 'structures', 'main']
                    )
                    self.datasheet_paths.append(path)
                except ValueError as e:
                    raise ValueError(f"Could not find datasheet in {datasheet}: {e}")
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
        """Generate script to select the best item from analysis results."""
        output_folder = self.output_folder
        os.makedirs(output_folder, exist_ok=True)

        # Output files
        selected_structure = os.path.join(output_folder, f"{self.output_name}.pdb")

        # Create config file for selection
        config_data = {
            "selection_metric": self.metric,
            "selection_mode": self.mode,
            "weights": self.weights,
            "tie_breaker": self.tie_breaker,
            "composite_function": self.composite_function,
            "output_csv": self.selected_csv,  # Auto-managed path!
            "output_structure": selected_structure,
            "output_sequences": os.path.join(output_folder, "sequences.csv"),
            "pool_folders": self.pool_folders,
            "datasheet_paths": self.datasheet_paths
        }

        import json
        with open(self.config_file, 'w') as f:  # Auto-managed path!
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
python "{self.select_script}" \\
  --config "{self.config_file}"

if [ $? -eq 0 ]; then
    echo "Successfully selected best item"
    echo "Selected structure: {selected_structure}"
    echo "Selected data: {self.selected_csv}"
else
    echo "Error: Failed to select best item"
    exit 1
fi

{self.generate_completion_check_footer()}
"""

        return script_content

    def get_output_files(self) -> Dict[str, List[str]]:
        """Get expected output files after selection."""
        selected_structure = os.path.join(self.output_folder, f"{self.output_name}.pdb")
        selected_sequence = os.path.join(self.output_folder, "sequences.csv")

        # Use output name as the selected ID
        selected_id = self.output_name

        # Define datasheets
        datasheets = {
            "selected": DatasheetInfo(
                name="selected",
                path=self.selected_csv,  # Auto-managed path!
                columns=["id", self.metric],
                description=f"Best item selected using {self.metric}",
                count=1
            ),
            "sequences": DatasheetInfo(
                name="sequences",
                path=selected_sequence,
                columns=["id", "sequence"],
                description="Selected best sequence",
                count=1
            )
        }

        return {
            "structures": [selected_structure],
            "structure_ids": [selected_id],
            "compounds": [],
            "compound_ids": [],
            "sequences": [selected_sequence],
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
