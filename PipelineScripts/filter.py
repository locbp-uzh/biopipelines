"""
Filter tool for expression-based filtering of unified analysis tables.

Takes unified tables from CombineTables and applies pandas query-style
expressions to filter rows while preserving all column information.
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
    Expression-based filter tool for unified analysis tables.

    Takes a unified table (typically from CombineTables) and applies
    pandas query-style expressions to filter rows while preserving all columns.
    """

    TOOL_NAME = "Filter"

    # Lazy path descriptors
    config_file = Path(lambda self: os.path.join(self.output_folder, "filter_config.json"))
    missing_csv = Path(lambda self: os.path.join(self.output_folder, "missing.csv"))
    helper_script = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_filter.py"))

    def __init__(self,
                 data: Union[StandardizedOutput, TableInfo],
                 expression: str,
                 pool: Union[StandardizedOutput, None] = None,
                 max_items: Optional[int] = None,
                 sort_by: Optional[str] = None,
                 sort_ascending: bool = True,
                 **kwargs):
        """
        Initialize Filter tool.

        Args:
            data: Table input to filter (required) - StandardizedOutput or TableInfo
            expression: Pandas query-style expression (e.g., "pLDDT>80 and distance < 5.0") (required)
            pool: Tool output having same ids as the data table (optional). Structures passing the filter will be copied
            max_items: Maximum number of items to keep after filtering
            sort_by: Column name to sort by before applying max_items limit
            sort_ascending: Sort order (True for ascending, False for descending)
            **kwargs: Additional parameters

        Examples:
            # Data mode - filter table only
            filtered = pipeline.add(Filter(
                data=combined_analysis,
                expression="pLDDT > 80"
            ))

            # Pool mode - filter table and copy structures
            filtered = pipeline.add(Filter(
                pool=boltz_results,
                data=combined_analysis,
                expression="pLDDT > 80 and distance < 5.0"
            ))
        """
        if data is None:
            raise ValueError("'data' parameter is required")
        if expression is None:
            raise ValueError("'expression' parameter is required")

        # Determine mode and set up inputs
        self.use_pool_mode = pool is not None
        self.pool_output = pool
        self.data_input = data

        self.expression = expression
        self.max_items = max_items
        self.sort_by = sort_by
        self.sort_ascending = sort_ascending

        # Validate expression
        self._validate_expression()

        # Initialize base class
        super().__init__(**kwargs)

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
        if not isinstance(self.data_input, (StandardizedOutput, TableInfo)):
            raise ValueError("data must be a StandardizedOutput or TableInfo object")

        if self.pool_output and not isinstance(self.pool_output, StandardizedOutput):
            raise ValueError("pool must be a StandardizedOutput object")

        if self.max_items is not None and self.max_items <= 0:
            raise ValueError("max_items must be positive")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input table and pool from previous tools."""
        self.folders = pipeline_folders

        # Configure data input (required)
        self.input_csv_path = None
        self.input_table_name = "filtered"  # Default name
        self.input_table_info = None

        # Check if data_input is a TableInfo object directly
        if isinstance(self.data_input, TableInfo):
            self.input_csv_path = self.data_input.path
            self.input_table_name = self.data_input.name
            self.input_table_info = self.data_input
        elif hasattr(self.data_input, 'tables'):
            tables = self.data_input.tables

            # Handle TableContainer objects
            if hasattr(tables, '_tables'):
                # Get the first available TableInfo object and its name
                first_name, ds_info = next(iter(tables._tables.items()))
                self.input_csv_path = ds_info.path
                self.input_table_name = first_name
                self.input_table_info = ds_info
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

    def _get_filtered_csv_path(self) -> str:
        """Get the path for the filtered CSV output."""
        output_table_name = getattr(self, 'input_table_name', 'filtered')
        return os.path.join(self.output_folder, f"{output_table_name}.csv")

    def get_config_display(self) -> List[str]:
        """Get configuration display lines."""
        config_lines = super().get_config_display()

        config_lines.extend([
            f"MODE: {'Pool + Data' if self.use_pool_mode else 'Data only'}",
            f"EXPRESSION: {self.expression}",
            f"MAX ITEMS: {self.max_items if self.max_items else 'unlimited'}",
        ])

        if self.sort_by:
            order = "ascending" if self.sort_ascending else "descending"
            config_lines.append(f"SORT BY: {self.sort_by} ({order})")

        return config_lines

    def generate_script(self, script_path: str) -> str:
        """Generate filter execution script."""
        os.makedirs(self.output_folder, exist_ok=True)

        filtered_csv = self._get_filtered_csv_path()

        config_data = {
            "input_csv": self.input_csv_path,
            "expression": self.expression,
            "max_items": self.max_items,
            "sort_by": self.sort_by,
            "sort_ascending": self.sort_ascending,
            "output_csv": filtered_csv,
            "use_pool_mode": self.use_pool_mode,
            "pool_output_folder": self.pool_folder if self.use_pool_mode else None
        }

        with open(self.config_file, 'w') as f:
            json.dump(config_data, f, indent=2)

        script_content = "#!/bin/bash\n"
        script_content += "# Filter execution script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()

        script_content += f"""echo "Applying filter expression"
echo "Input: {self.input_csv_path}"
echo "Expression: {self.expression}"
echo "Output: {filtered_csv}"

python "{self.helper_script}" \\
  --config "{self.config_file}"

if [ $? -eq 0 ]; then
    echo "Results written to: {filtered_csv}"
else
    echo "Error: Filtering failed"
    exit 1
fi

"""
        script_content += self.generate_completion_check_footer()

        return script_content

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files after filtering."""
        filtered_csv = self._get_filtered_csv_path()

        # Get actual column names from input
        input_columns = self._get_input_columns()

        # Use the same table name as the input to preserve structure
        output_table_name = getattr(self, 'input_table_name', 'filtered')

        tables = {
            output_table_name: TableInfo(
                name=output_table_name,
                path=filtered_csv,
                columns=input_columns,
                description=f"Filtered results using expression: {self.expression}",
                count="variable"
            ),
            "missing": TableInfo(
                name="missing",
                path=self.missing_csv,
                columns=["id", "structure", "msa"],
                description="IDs that were filtered out and their expected file paths",
                count="variable"
            )
        }

        if self.use_pool_mode:
            # Pool mode: structures, compounds, sequences from pool
            updated_structures = []
            updated_compounds = []
            updated_sequences = []

            # Copy structure file paths with original names
            if self.pool_output.structures and len(self.pool_output.structures) > 0:
                for struct_path in self.pool_output.structures.files:
                    filename = os.path.basename(struct_path)
                    updated_structures.append(os.path.join(self.output_folder, filename))

            # Copy compound file paths with original names
            if self.pool_output.compounds and len(self.pool_output.compounds) > 0:
                for comp_path in self.pool_output.compounds.files:
                    filename = os.path.basename(comp_path)
                    updated_compounds.append(os.path.join(self.output_folder, filename))

            # Copy sequence file paths with original names
            if self.pool_output.sequences and len(self.pool_output.sequences) > 0:
                for seq_path in self.pool_output.sequences.files:
                    filename = os.path.basename(seq_path)
                    updated_sequences.append(os.path.join(self.output_folder, filename))

            # Combine pool tables with filtered table
            if hasattr(self.pool_output, 'tables') and hasattr(self.pool_output.tables, '_tables'):
                for name, info in self.pool_output.tables._tables.items():
                    filename = os.path.basename(info.path)
                    tables[name] = TableInfo(
                        name=name,
                        path=os.path.join(self.output_folder, filename),
                        columns=info.columns,
                        description=info.description,
                        count=info.count
                    )

            structures = DataStream(
                name="structures",
                ids=self.pool_output.structures.ids if self.pool_output.structures else [],
                files=updated_structures,
                format="pdb"
            )
            compounds = DataStream(
                name="compounds",
                ids=self.pool_output.compounds.ids if self.pool_output.compounds else [],
                files=updated_compounds,
                format="sdf"
            )
            sequences = DataStream(
                name="sequences",
                ids=self.pool_output.sequences.ids if self.pool_output.sequences else [],
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
            # Data mode: only filtered CSV
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
                "expression": self.expression,
                "max_items": self.max_items,
                "sort_by": self.sort_by,
                "sort_ascending": self.sort_ascending
            }
        })
        return base_dict
