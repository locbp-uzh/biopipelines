"""
SelectBest tool for choosing the single best item from analysis results.

Takes analysis tables and applies selection criteria to pick exactly one
best structure/sequence/compound for use in iterative optimization cycles.
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


class SelectBest(BaseConfig):
    """
    Tool for selecting the single best item from analysis results.

    Takes analysis tables and applies selection criteria to return exactly
    one best structure/sequence/compound with its associated data.
    """

    TOOL_NAME = "SelectBest"

    # Lazy path descriptors
    selected_csv = Path(lambda self: os.path.join(self.output_folder, "selected_best.csv"))
    selected_structure = Path(lambda self: os.path.join(self.output_folder, f"{self.output_name}.pdb"))
    selected_compound = Path(lambda self: os.path.join(self.output_folder, f"{self.output_name}_compound.sdf"))
    sequences_csv = Path(lambda self: os.path.join(self.output_folder, "sequences.csv"))
    config_file = Path(lambda self: os.path.join(self.output_folder, "select_best_config.json"))
    helper_script = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_select_best.py"))

    def __init__(self,
                 pool: Union[StandardizedOutput, List[StandardizedOutput]],
                 tables: Union[List[Union[StandardizedOutput, TableInfo, str]], StandardizedOutput, TableInfo, str],
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
            pool: Single or list of StandardizedOutput to select structures from
            tables: Single or list of tables to evaluate for selection
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
                tables=[original_analysis.tables.merged, engineered_analysis.tables.merged],
                metric="binding_affinity",
                mode="max"
            ))

            # Single pool selection
            best = pipeline.add(SelectBest(
                pool=boltz,
                tables=analysis.tables.merged,
                metric="pLDDT",
                mode="max"
            ))

            # Multi-objective with weights
            best = pipeline.add(SelectBest(
                pool=combined_results,
                tables=analysis.tables.merged,
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

        # Handle tables - always expect list, convert if single
        if isinstance(tables, list):
            self.tables_input = tables
        else:
            self.tables_input = [tables]

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
            if not isinstance(pool_output, StandardizedOutput):
                raise ValueError(f"pool[{i}] is {type(pool_output)}, must be a StandardizedOutput object. Value: {pool_output}")

        # Validate tables
        if not self.tables_input:
            raise ValueError("tables must be provided")

        if len(self.tables_input) != len(self.pool_outputs):
            raise ValueError(f"Number of tables ({len(self.tables_input)}) must match number of pools ({len(self.pool_outputs)})")

        for i, table in enumerate(self.tables_input):
            if not isinstance(table, (StandardizedOutput, TableInfo, str)):
                raise ValueError(f"table[{i}] is {type(table)}, must be a StandardizedOutput, TableInfo, or str object. Value: {table}")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input table from previous tool."""
        self.folders = pipeline_folders
        self._configure_pool_mode()

    def _configure_pool_mode(self):
        """Configure inputs for pool + tables selection."""
        # Extract pool folders and table paths
        self.pool_folders = []
        self.table_paths = []

        for pool in self.pool_outputs:
            if hasattr(pool, 'output_folder'):
                self.pool_folders.append(pool.output_folder)
            else:
                raise ValueError(f"Pool {pool} must have output_folder")

        for table in self.tables_input:
            # Handle different table input types - check TableInfo first
            if hasattr(table, 'path') and hasattr(table, 'name') and hasattr(table, 'columns'):
                # TableInfo object - use path directly
                self.table_paths.append(table.path)
            elif isinstance(table, str):
                # String path - use directly
                self.table_paths.append(table)
            elif hasattr(table, 'tables'):
                # StandardizedOutput object - extract table path
                ds_obj = table.tables
                table_path = None

                if hasattr(ds_obj, 'merged'):
                    table_path = ds_obj.merged.path
                elif hasattr(ds_obj, 'combined'):
                    table_path = ds_obj.combined.path
                elif hasattr(ds_obj, 'filtered'):
                    table_path = ds_obj.filtered.path
                elif hasattr(ds_obj, '_tables'):
                    # Get first available table
                    first_ds = next(iter(ds_obj._tables.values()))
                    table_path = first_ds.path
                else:
                    raise ValueError(f"Could not find table in {table}")

                self.table_paths.append(table_path)
            else:
                raise ValueError(f"Unsupported table type: {type(table)}")

    def get_config_display(self) -> List[str]:
        """Get configuration display lines."""
        config_lines = super().get_config_display()

        config_lines.extend([
            f"POOLS: {len(self.pool_outputs)}",
            f"DATASHEETS: {len(self.tables_input)}",
            f"METRIC: {self.metric} ({self.mode})"
        ])

        return config_lines

    def generate_script(self, script_path: str) -> str:
        """Generate SelectBest execution script."""
        os.makedirs(self.output_folder, exist_ok=True)

        # Create config file for selection
        config_data = {
            "selection_metric": self.metric,
            "selection_mode": self.mode,
            "weights": self.weights,
            "tie_breaker": self.tie_breaker,
            "composite_function": self.composite_function,
            "output_csv": self.selected_csv,
            "output_structure": self.selected_structure,
            "output_sequences": self.sequences_csv,
            "pool_folders": self.pool_folders,
            "table_paths": self.table_paths
        }

        with open(self.config_file, 'w') as f:
            json.dump(config_data, f, indent=2)

        script_content = "#!/bin/bash\n"
        script_content += "# SelectBest execution script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()

        script_content += f"""echo "Selecting best item from {len(self.table_paths)} tables"
echo "Selection metric: {self.metric} ({self.mode})"
echo "Output: {self.selected_structure}"

python "{self.helper_script}" \\
  --config "{self.config_file}"

if [ $? -eq 0 ]; then
    echo "Successfully selected best item"
    echo "Selected structure: {self.selected_structure}"
    echo "Selected data: {self.selected_csv}"
else
    echo "Error: Failed to select best item"
    exit 1
fi

"""
        script_content += self.generate_completion_check_footer()

        return script_content

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files after selection."""
        # Use output name as the selected ID
        selected_id = self.output_name

        # Define tables that will be created
        tables = {
            "selected": TableInfo(
                name="selected",
                path=self.selected_csv,
                columns=["id", self.metric],
                description=f"Best item selected using {self.metric}",
                count=1
            ),
            "sequences": TableInfo(
                name="sequences",
                path=self.sequences_csv,
                columns=["id", "sequence"],
                description="Selected best sequence",
                count=1
            )
        }

        # Create DataStreams for output
        structures = DataStream(
            name="structures",
            ids=[selected_id],
            files=[self.selected_structure],
            format="pdb"
        )

        sequences = DataStream(
            name="sequences",
            ids=[selected_id],
            files=[self.sequences_csv],
            map_table=self.sequences_csv,
            format="csv"
        )

        # Compounds only created if available in source pools
        compounds = DataStream.empty("compounds", "sdf")

        return {
            "structures": structures,
            "sequences": sequences,
            "compounds": compounds,
            "tables": tables,
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
