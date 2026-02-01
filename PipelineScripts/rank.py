"""
Rank tool for sorting and ranking entries based on a metric.

Takes tables and ranks entries by a specified metric (column or computed expression),
optionally renames IDs to sequential ranks (rank_1, rank_2, ...), and can limit to top N.
"""

import os
import json
import re
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


class Rank(BaseConfig):
    """
    Rank tool for sorting entries by a metric and optionally renaming IDs.

    Ranks entries in a table based on a metric (existing column or computed expression),
    creates sequential ranked IDs (e.g., rank_1, rank_2), and optionally copies
    structures/compounds in ranked order.
    """

    TOOL_NAME = "Rank"

    # Lazy path descriptors
    ranked_csv = Path(lambda self: os.path.join(self.output_folder, "ranked.csv"))
    metrics_csv = Path(lambda self: os.path.join(self.output_folder, "metrics.csv"))
    config_file = Path(lambda self: os.path.join(self.output_folder, "rank_config.json"))
    helper_script = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_rank.py"))

    def __init__(self,
                 data: Union[StandardizedOutput, TableInfo],
                 metric: str,
                 pool: Union[StandardizedOutput, None] = None,
                 ascending: bool = False,
                 prefix: str = "rank",
                 top: Optional[int] = None,
                 **kwargs):
        """
        Initialize Rank tool.

        Args:
            data: Table input to rank (required) - StandardizedOutput or TableInfo
            metric: Column name or expression for ranking (e.g., "pLDDT" or "0.8*pLDDT+0.2*affinity") (required)
            pool: Structure pool for copying ranked structures (optional)
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
        if data is None:
            raise ValueError("'data' parameter is required")

        # Handle tuple notation: data=tool.tables.table_name.column_name
        # This returns (TableInfo, column_name) tuple
        if isinstance(data, tuple) and len(data) == 2:
            table_info, column_name = data
            if hasattr(table_info, 'path'):
                # Extract metric from tuple if not explicitly provided
                if metric is None:
                    metric = column_name
                # Use TableInfo as data
                data = table_info

        if metric is None:
            raise ValueError("'metric' parameter is required")

        # Determine mode and set up inputs
        self.use_pool_mode = pool is not None
        self.pool_output = pool
        self.data_input = data

        self.metric = metric
        self.ascending = ascending
        self.prefix = prefix
        self.top = top

        # Validate metric
        self._validate_metric()

        # Initialize base class
        super().__init__(**kwargs)

    def _validate_metric(self):
        """Validate that the metric is safe for evaluation."""
        if not self.metric.strip():
            raise ValueError("Metric cannot be empty")

        # Basic safety check - ensure only safe characters and operations
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
        if not isinstance(self.data_input, (StandardizedOutput, TableInfo)):
            raise ValueError("data must be a StandardizedOutput or TableInfo object")

        if self.pool_output and not isinstance(self.pool_output, StandardizedOutput):
            raise ValueError("pool must be a StandardizedOutput object")

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
            raise ValueError(f"Could not determine input CSV path from data: {self.data_input}")

        # Configure pool input (optional)
        if self.use_pool_mode:
            self.pool_folder = None
            if hasattr(self.pool_output, 'output_folder'):
                self.pool_folder = self.pool_output.output_folder

            if not self.pool_folder:
                raise ValueError(f"Could not determine pool folder from pool: {self.pool_output}")

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
        """Generate rank execution script."""
        os.makedirs(self.output_folder, exist_ok=True)

        # Determine num_digits from pool structure_ids count
        predicted_count = 0
        if self.use_pool_mode:
            if self.pool_output.structures and len(self.pool_output.structures) > 0:
                predicted_count = len(self.pool_output.structures)

        if self.top and predicted_count > 0:
            predicted_count = min(predicted_count, self.top)

        num_digits = len(str(predicted_count)) if predicted_count > 0 else 1

        config_data = {
            "input_csv": self.input_csv_path,
            "metric": self.metric,
            "ascending": self.ascending,
            "prefix": self.prefix,
            "top": self.top,
            "output_csv": self.ranked_csv,
            "use_pool_mode": self.use_pool_mode,
            "pool_output_folder": self.pool_folder if self.use_pool_mode else None,
            "num_digits": num_digits
        }

        with open(self.config_file, 'w') as f:
            json.dump(config_data, f, indent=2)

        script_content = "#!/bin/bash\n"
        script_content += "# Rank execution script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()

        script_content += f"""echo "Ranking entries by metric"
echo "Input: {self.input_csv_path}"
echo "Metric: {self.metric}"
echo "Output: {self.ranked_csv}"

python "{self.helper_script}" \\
  --config "{self.config_file}"

if [ $? -eq 0 ]; then
    echo "Results written to: {self.ranked_csv}"
else
    echo "Error: Ranking failed"
    exit 1
fi

"""
        script_content += self.generate_completion_check_footer()

        return script_content

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files after ranking."""
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
            metric_col_name = "metric"
        else:
            # It's a simple column reference - will be preserved
            output_columns.append(self.metric)
            metric_col_name = self.metric

        # Add remaining original columns (excluding 'id')
        for col in input_columns:
            if col not in output_columns and col != 'id':
                output_columns.append(col)

        if not output_columns:
            raise ValueError("Failed to determine output columns - this indicates a logic error")

        tables = {
            "ranked": TableInfo(
                name="ranked",
                path=self.ranked_csv,
                columns=output_columns,
                description=f"Ranked results by metric: {self.metric}",
                count="variable"
            ),
            "metrics": TableInfo(
                name="metrics",
                path=self.metrics_csv,
                columns=["id", "source_id", metric_col_name],
                description=f"Summary table with id, source_id, and {metric_col_name}",
                count="variable"
            )
        }

        if self.use_pool_mode:
            # Pool mode: predict copying structures in ranked order
            updated_structures = []
            updated_compounds = []
            updated_sequences = []

            # Use the SAME num_digits calculation as in generate_script()
            predicted_count = 0
            if self.pool_output.structures and len(self.pool_output.structures) > 0:
                predicted_count = len(self.pool_output.structures)

            # Apply top limit if specified
            if self.top and predicted_count > 0:
                predicted_count = min(predicted_count, self.top)

            # Calculate num_digits (same logic as in generate_script)
            num_digits = len(str(predicted_count)) if predicted_count > 0 else 1

            # Copy structure file paths with ranked names
            if self.pool_output.structures and len(self.pool_output.structures) > 0:
                for i, struct_path in enumerate(self.pool_output.structures.files):
                    if self.top and i >= self.top:
                        break
                    ext = os.path.splitext(struct_path)[1]
                    ranked_name = f"{self.prefix}_{i+1:0{num_digits}d}{ext}"
                    updated_structures.append(os.path.join(self.output_folder, ranked_name))

            # Copy compound file paths with ranked names
            if self.pool_output.compounds and len(self.pool_output.compounds) > 0:
                for i, comp_path in enumerate(self.pool_output.compounds.files):
                    if self.top and i >= self.top:
                        break
                    ext = os.path.splitext(comp_path)[1]
                    ranked_name = f"{self.prefix}_{i+1:0{num_digits}d}{ext}"
                    updated_compounds.append(os.path.join(self.output_folder, ranked_name))

            # Copy sequence file paths with ranked names
            if self.pool_output.sequences and len(self.pool_output.sequences) > 0:
                for i, seq_path in enumerate(self.pool_output.sequences.files):
                    if self.top and i >= self.top:
                        break
                    ext = os.path.splitext(seq_path)[1]
                    ranked_name = f"{self.prefix}_{i+1:0{num_digits}d}{ext}"
                    updated_sequences.append(os.path.join(self.output_folder, ranked_name))

            # Combine pool tables with ranked table
            if hasattr(self.pool_output, 'tables') and hasattr(self.pool_output.tables, '_tables'):
                for name, info in self.pool_output.tables._tables.items():
                    # Skip the 'ranked' table to avoid duplication (we create our own)
                    if name == 'ranked':
                        continue
                    filename = os.path.basename(info.path)
                    tables[name] = TableInfo(
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

            structures = DataStream(
                name="structures",
                ids=structure_ids,
                files=updated_structures,
                format="pdb"
            )
            compounds = DataStream(
                name="compounds",
                ids=compound_ids,
                files=updated_compounds,
                format="sdf"
            )
            sequences = DataStream(
                name="sequences",
                ids=sequence_ids,
                files=updated_sequences,
                format="fasta"
            )

            return {
                "structures": structures,
                "compounds": compounds,
                "sequences": sequences,
                "tables": tables,
                "output_folder": self.output_folder
            }
        else:
            # Data mode: only ranked CSV
            return {
                "structures": DataStream.empty("structures", "pdb"),
                "sequences": DataStream.empty("sequences", "fasta"),
                "compounds": DataStream.empty("compounds", "sdf"),
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
