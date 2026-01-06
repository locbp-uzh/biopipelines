"""
Rank tool for sorting and ranking entries based on a metric.

Takes tables and ranks entries by a specified metric (column or computed expression),
optionally renames IDs to sequential ranks (rank_1, rank_2, ...), and can limit to top N.
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


class Rank(BaseConfig):
    """
    Rank tool for sorting entries by a metric and optionally renaming IDs.

    Ranks entries in a table based on a metric (existing column or computed expression),
    creates sequential ranked IDs (e.g., rank_1, rank_2), and optionally copies
    structures/compounds in ranked order.
    """

    # Tool identification
    TOOL_NAME = "Rank"
    DEFAULT_ENV = None  # Loaded from config.yaml

    def __init__(self,
                 data: Union[ToolOutput, StandardizedOutput, TableInfo] = None,
                 pool: Union[ToolOutput, StandardizedOutput] = None,
                 metric: str = None,
                 ascending: bool = False,
                 prefix: str = "rank",
                 top: Optional[int] = None,
                 **kwargs):
        """
        Initialize Rank tool.

        Args:
            data: Table input to rank (required)
            pool: Structure pool for copying ranked structures (optional)
            metric: Column name or expression for ranking (e.g., "pLDDT" or "0.8*pLDDT+0.2*affinity") (required)
            ascending: Sort order - False for descending (highest first), True for ascending
            prefix: Prefix for renamed IDs (default: "rank", produces rank_1, rank_2, ...)
            top: Limit to top N entries after ranking (optional)
            **kwargs: Additional parameters

        Examples:
            # Rank by single column
            ranked = pipeline.add(Rank(
                data=analysis,
                metric="pLDDT",
                ascending=False  # Higher pLDDT is better
            ))

            # Rank by computed metric with renamed IDs
            ranked = pipeline.add(Rank(
                data=analysis,
                metric="0.8*pLDDT + 0.2*affinity",
                prefix="model",
                top=10  # Keep only top 10
            ))

            # Rank with pool mode (copy structures in ranked order)
            ranked = pipeline.add(Rank(
                pool=boltz_results,
                data=analysis,
                metric="pLDDT",
                top=5
            ))
        """
        # Validate required parameters
        if data is None:
            raise ValueError("'data' parameter is required")
        if metric is None:
            raise ValueError("'metric' parameter is required")

        # Determine mode and set up inputs
        if pool is not None:
            # Pool mode: rank data and copy structures from pool
            self.use_pool_mode = True
            self.pool_output = pool
            self.data_input = data
        else:
            # Data mode: rank data only
            self.use_pool_mode = False
            self.pool_output = None
            self.data_input = data

        self.metric = metric
        self.ascending = ascending
        self.prefix = prefix
        self.top = top

        # Validate metric
        self._validate_metric()

        # Initialize base class
        super().__init__(**kwargs)

        # Set up dependencies
        dependencies = []
        if hasattr(data, 'config'):
            dependencies.append(data.config)
        if pool and hasattr(pool, 'config'):
            dependencies.append(pool.config)
        self.dependencies.extend(dependencies)

    def _validate_metric(self):
        """Validate that the metric is safe for evaluation."""
        if not self.metric.strip():
            raise ValueError("Metric cannot be empty")

        # Basic safety check - ensure only safe characters and operations
        import re
        allowed_pattern = r'^[a-zA-Z_][a-zA-Z0-9_\s\.\+\-\*\/\(\)<>=!&|and\sor\snot\s\d]+$'

        if not re.match(allowed_pattern, self.metric):
            raise ValueError(f"Invalid characters in metric: {self.metric}")

        # Check for dangerous keywords
        dangerous_keywords = ['import', 'exec', 'eval', '__', 'os.', 'sys.', 'subprocess']
        metric_lower = self.metric.lower()
        for keyword in dangerous_keywords:
            if keyword in metric_lower:
                raise ValueError(f"Dangerous keyword '{keyword}' not allowed in metric")

    def validate_params(self):
        """Validate Rank parameters."""
        if not isinstance(self.data_input, (ToolOutput, StandardizedOutput, TableInfo)):
            raise ValueError("data must be a ToolOutput, StandardizedOutput, or TableInfo object")

        if self.pool_output and not isinstance(self.pool_output, (ToolOutput, StandardizedOutput)):
            raise ValueError("pool must be a ToolOutput or StandardizedOutput object")

        if self.top is not None and self.top <= 0:
            raise ValueError("top must be positive")

        if not self.prefix:
            raise ValueError("prefix cannot be empty")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input table and pool from previous tools."""
        self.folders = pipeline_folders

        # Configure data input (required)
        self.input_csv_path = None
        self.input_table_name = "ranked"  # Default name
        self.input_table_info = None
        self.input_row_count = None  # Will be determined at runtime

        # Check if data_input is a TableInfo object directly
        if isinstance(self.data_input, TableInfo):
            self.input_csv_path = self.data_input.path
            self.input_table_name = self.data_input.name
            self.input_table_info = self.data_input
            if hasattr(self.data_input, 'count') and isinstance(self.data_input.count, int):
                self.input_row_count = self.data_input.count
        elif hasattr(self.data_input, 'tables'):
            tables = self.data_input.tables

            # Handle TableContainer objects
            if hasattr(tables, '_tables'):
                # Get the first available TableInfo object and its name
                first_name, ds_info = next(iter(tables._tables.items()))
                self.input_csv_path = ds_info.path
                self.input_table_name = first_name
                self.input_table_info = ds_info
                # Try to get count if available
                if hasattr(ds_info, 'count') and isinstance(ds_info.count, int):
                    self.input_row_count = ds_info.count
            elif isinstance(tables, dict):
                # Handle raw dict (legacy format)
                first_name, ds_info = next(iter(tables.items()))
                self.input_table_name = first_name

                if isinstance(ds_info, dict) and 'path' in ds_info:
                    self.input_csv_path = ds_info['path']
                elif hasattr(ds_info, 'path'):
                    # Handle TableInfo objects
                    self.input_csv_path = ds_info.path
                    self.input_table_info = ds_info
                else:
                    self.input_csv_path = str(ds_info)

        if not self.input_csv_path:
            raise ValueError(f"Could not predict input CSV path from data tool: {self.data_input}. "
                           f"The data tool must provide proper table path predictions.")

        # Configure pool input (optional)
        if self.use_pool_mode:
            self.pool_folder = None
            if hasattr(self.pool_output, 'output_folder'):
                self.pool_folder = self.pool_output.output_folder

            if not self.pool_folder:
                raise ValueError(f"Could not predict pool folder from pool tool: {self.pool_output}. "
                               f"The pool tool must provide output_folder.")

    def _get_input_columns(self) -> List[str]:
        """Get column names from the input table."""
        # If data_input is a TableInfo directly, get columns from it
        if isinstance(self.data_input, TableInfo):
            if hasattr(self.data_input, 'columns') and self.data_input.columns:
                return self.data_input.columns
        # Try to get columns from data tool's table info
        elif hasattr(self.data_input, 'tables') and hasattr(self.data_input.tables, '_tables'):
            # Look through the tables to find one with column info
            for name, info in self.data_input.tables._tables.items():
                if hasattr(info, 'columns') and info.columns:
                    return info.columns

        # No fallbacks - columns will be determined at runtime
        return []

    def _detect_metric_variables(self) -> List[str]:
        """
        Detect variable names used in the metric expression.
        Returns list of column names referenced in the expression.
        """
        import re
        # Extract potential column names (identifiers)
        # Match word characters that could be column names
        potential_vars = re.findall(r'\b[a-zA-Z_][a-zA-Z0-9_]*\b', self.metric)

        # Filter out Python keywords and operators
        keywords = {'and', 'or', 'not', 'in', 'is', 'if', 'else', 'for', 'while', 'def', 'class'}
        variables = [v for v in potential_vars if v not in keywords]

        return list(set(variables))  # Remove duplicates

    def get_config_display(self) -> List[str]:
        """Get configuration display lines."""
        config_lines = super().get_config_display()

        order = "ascending (low to high)" if self.ascending else "descending (high to low)"

        config_lines.extend([
            f"MODE: {'Pool + Data' if self.use_pool_mode else 'Data only'}",
            f"METRIC: {self.metric}",
            f"SORT ORDER: {order}",
            f"PREFIX: {self.prefix}",
            f"TOP N: {self.top if self.top else 'all'}",
        ])

        return config_lines

    def generate_script(self, script_path: str) -> str:
        """
        Generate rank execution script.

        Args:
            script_path: Path where script should be written

        Returns:
            Script content as string
        """
        runtime_folder = os.path.dirname(script_path)
        output_folder = self.output_folder
        os.makedirs(output_folder, exist_ok=True)

        # Output CSV path
        ranked_csv_name = "ranked.csv"
        ranked_csv = os.path.join(output_folder, ranked_csv_name)

        # Determine num_digits from pool structure_ids count (only structures can be ranked)
        predicted_count = 0
        if self.use_pool_mode:
            if hasattr(self.pool_output, 'structure_ids') and self.pool_output.structure_ids:
                predicted_count = len(self.pool_output.structure_ids)
            elif hasattr(self.pool_output, 'structures') and self.pool_output.structures:
                predicted_count = len(self.pool_output.structures)

        # Apply top limit if specified
        if self.top and predicted_count > 0:
            predicted_count = min(predicted_count, self.top)

        # Calculate num_digits
        num_digits = len(str(predicted_count)) if predicted_count > 0 else 1

        # Create config file for the rank operation
        config_file = os.path.join(output_folder, "rank_config.json")
        config_data = {
            "input_csv": self.input_csv_path,
            "metric": self.metric,
            "ascending": self.ascending,
            "prefix": self.prefix,
            "top": self.top,
            "output_csv": ranked_csv,
            "use_pool_mode": self.use_pool_mode,
            "pool_output_folder": self.pool_folder if self.use_pool_mode else None,
            "num_digits": num_digits  # Pass the calculated num_digits
        }

        import json
        with open(config_file, 'w') as f:
            json.dump(config_data, f, indent=2)

        # Generate script content
        script_content = f"""#!/bin/bash
# Rank execution script
# Generated by BioPipelines pipeline system

{self.generate_completion_check_header()}

echo "Ranking entries by metric"
echo "Input: {self.input_csv_path}"
echo "Metric: {self.metric}"
echo "Output: {ranked_csv}"

# Run Python ranking script
python "{os.path.join(self.folders['HelpScripts'], 'pipe_rank.py')}" \\
  --config "{config_file}"

if [ $? -eq 0 ]; then
    echo "Results written to: {ranked_csv}"
else
    echo "Error: Ranking failed"
    exit 1
fi

{self.generate_completion_check_footer()}
"""

        return script_content

    def _calculate_num_digits(self, count: int) -> int:
        """Calculate number of digits needed for zero-padding based on count."""
        if count == 0:
            return 1
        return len(str(count))

    def get_output_files(self) -> Dict[str, List[str]]:
        """
        Get expected output files after ranking.

        Returns:
            Dictionary with output file paths
        """
        ranked_csv = os.path.join(self.output_folder, "ranked.csv")

        # Predict output columns
        input_columns = self._get_input_columns()

        # Add new columns that will be created
        output_columns = ["id", "source_id"]  # Renamed id and original id

        # Detect if metric is an expression or simple column name
        metric_vars = self._detect_metric_variables()

        # Add metric column (if computed expression)
        if len(metric_vars) > 1 or any(op in self.metric for op in ['+', '-', '*', '/', '(', ')']):
            # It's an expression - add computed metric column
            output_columns.append("metric")
            # Add individual variable columns
            output_columns.extend(metric_vars)
        else:
            # It's a simple column reference - will be preserved
            output_columns.append(self.metric)

        # Add remaining original columns (excluding 'id')
        for col in input_columns:
            if col not in output_columns and col != 'id':
                output_columns.append(col)

        tables = {
            "ranked": TableInfo(
                name="ranked",
                path=ranked_csv,
                columns=output_columns if output_columns else ["id", "source_id", "metric"],
                description=f"Ranked results by metric: {self.metric}",
                count="variable"
            )
        }

        if self.use_pool_mode:
            # Pool mode: predict copying structures in ranked order
            pool_output_dict = self.pool_output._data.copy() if hasattr(self.pool_output, '_data') else {}

            # Update paths to point to our output folder
            updated_structures = []
            updated_compounds = []
            updated_sequences = []

            # Use the SAME num_digits calculation as in generate_script()
            # This ensures prediction matches runtime output
            predicted_count = 0
            if hasattr(self.pool_output, 'structure_ids') and self.pool_output.structure_ids:
                predicted_count = len(self.pool_output.structure_ids)
            elif hasattr(self.pool_output, 'structures') and self.pool_output.structures:
                predicted_count = len(self.pool_output.structures)

            # Apply top limit if specified
            if self.top and predicted_count > 0:
                predicted_count = min(predicted_count, self.top)

            # Calculate num_digits (same logic as in generate_script)
            num_digits = len(str(predicted_count)) if predicted_count > 0 else 1

            # Copy structure file paths with ranked names
            if hasattr(self.pool_output, 'structures') and self.pool_output.structures:
                for i, struct_path in enumerate(self.pool_output.structures):
                    if self.top and i >= self.top:
                        break
                    ext = os.path.splitext(struct_path)[1]
                    ranked_name = f"{self.prefix}_{i+1:0{num_digits}d}{ext}"
                    updated_structures.append(os.path.join(self.output_folder, ranked_name))

            # Copy compound file paths with ranked names
            if hasattr(self.pool_output, 'compounds') and self.pool_output.compounds:
                for i, comp_path in enumerate(self.pool_output.compounds):
                    if self.top and i >= self.top:
                        break
                    ext = os.path.splitext(comp_path)[1]
                    ranked_name = f"{self.prefix}_{i+1:0{num_digits}d}{ext}"
                    updated_compounds.append(os.path.join(self.output_folder, ranked_name))

            # Copy sequence file paths with ranked names
            if hasattr(self.pool_output, 'sequences') and self.pool_output.sequences:
                for i, seq_path in enumerate(self.pool_output.sequences):
                    if self.top and i >= self.top:
                        break
                    ext = os.path.splitext(seq_path)[1]
                    ranked_name = f"{self.prefix}_{i+1:0{num_digits}d}{ext}"
                    updated_sequences.append(os.path.join(self.output_folder, ranked_name))

            # Combine pool tables with ranked table
            combined_tables = tables.copy()
            if hasattr(self.pool_output, 'tables') and hasattr(self.pool_output.tables, '_tables'):
                for name, info in self.pool_output.tables._tables.items():
                    # Skip the 'ranked' table to avoid duplication (we create our own)
                    if name == 'ranked':
                        continue
                    filename = os.path.basename(info.path)
                    combined_tables[name] = TableInfo(
                        name=name,
                        path=os.path.join(self.output_folder, filename),
                        columns=info.columns,
                        description=info.description,
                        count=info.count
                    )

            # Generate ranked IDs with zero-padding
            structure_ids = [f"{self.prefix}_{i+1:0{num_digits}d}" for i in range(len(updated_structures))]
            compound_ids = [f"{self.prefix}_{i+1:0{num_digits}d}" for i in range(len(updated_compounds))]
            sequence_ids = [f"{self.prefix}_{i+1:0{num_digits}d}" for i in range(len(updated_sequences))]

            return {
                "structures": updated_structures,
                "structure_ids": structure_ids,
                "compounds": updated_compounds,
                "compound_ids": compound_ids,
                "sequences": updated_sequences,
                "sequence_ids": sequence_ids,
                "tables": combined_tables,
                "output_folder": self.output_folder
            }
        else:
            # Data mode: only ranked CSV
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
                "metric": self.metric,
                "ascending": self.ascending,
                "prefix": self.prefix,
                "top": self.top
            }
        })
        return base_dict
