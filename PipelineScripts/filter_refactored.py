"""
Filter tool for expression-based filtering of unified analysis datasheets (REFACTORED).

Takes unified datasheets from CombineDatasheets and applies pandas query-style
expressions to filter rows while preserving all column information.

REFACTORED VERSION using new mixin architecture for cleaner, more maintainable code.
"""

import os
import pandas as pd
from typing import Dict, List, Any, Optional, Union

try:
    from .base_config import BaseConfig, ToolOutput, StandardizedOutput, DatasheetInfo
    from .mixins import InputHandlerMixin, DatasheetNavigatorMixin, FilePathDescriptor
except ImportError:
    import sys
    import os
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, ToolOutput, StandardizedOutput, DatasheetInfo
    from mixins import InputHandlerMixin, DatasheetNavigatorMixin, FilePathDescriptor


class FilterResult:
    """Result from applying expression-based filtering."""

    def __init__(self, input_csv: str, output_csv: str, expression: str,
                 total_input: int, kept_count: int, filtered_count: int):
        self.input_csv = input_csv
        self.output_csv = output_csv
        self.expression = expression
        self.total_input = total_input
        self.kept_count = kept_count
        self.filtered_count = filtered_count
        self.pass_rate = kept_count / total_input if total_input > 0 else 0.0

    def to_dict(self) -> Dict[str, Any]:
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
        return (f"Filter ({self.expression}): {self.kept_count}/{self.total_input} items kept "
                f"({self.pass_rate:.1%} pass rate)")


class Filter(DatasheetNavigatorMixin, BaseConfig):
    """
    Expression-based filter tool for unified analysis datasheets.

    REFACTORED to use DatasheetNavigatorMixin for elegant datasheet access.
    """

    TOOL_NAME = "Filter"
    DEFAULT_ENV = None

    # ============================================================================
    # AUTOMATIC FILE PATH MANAGEMENT
    # ============================================================================
    config_file = FilePathDescriptor("filter_config.json")
    filter_script = FilePathDescriptor("pipe_filter.py", folder_key="HelpScripts")

    def __init__(self,
                 data: Union[ToolOutput, StandardizedOutput] = None,
                 pool: Union[ToolOutput, StandardizedOutput] = None,
                 expression: str = None,
                 max_items: Optional[int] = None,
                 sort_by: Optional[str] = None,
                 sort_ascending: bool = True,
                 **kwargs):
        """Initialize Filter tool."""
        # Validate required parameters
        if data is None:
            raise ValueError("'data' parameter is required")
        if expression is None:
            raise ValueError("'expression' parameter is required")

        # Store parameters
        self.data_input = data
        self.pool_output = pool
        self.expression = expression
        self.max_items = max_items
        self.sort_by = sort_by
        self.sort_ascending = sort_ascending

        # Determine mode
        self.use_pool_mode = (pool is not None)

        # Validate expression
        self._validate_expression()

        # Initialize base class
        super().__init__(**kwargs)

        # Set up dependencies
        dependencies = []
        if hasattr(data, 'config'):
            dependencies.append(data.config)
        if pool and hasattr(pool, 'config'):
            dependencies.append(pool.config)
        self.dependencies.extend(dependencies)

    def _validate_expression(self):
        """Validate that the expression is safe for pandas query."""
        if not self.expression.strip():
            raise ValueError("Filter expression cannot be empty")

        import re
        allowed_pattern = r'^[a-zA-Z_][a-zA-Z0-9_\s\.\+\-\*\/\(\)<>=!&|and\sor\snot\s\d]+$'

        if not re.match(allowed_pattern, self.expression):
            raise ValueError(f"Invalid characters in expression: {self.expression}")

        dangerous_keywords = ['import', 'exec', 'eval', '__', 'os.', 'sys.', 'subprocess']
        expr_lower = self.expression.lower()
        for keyword in dangerous_keywords:
            if keyword in expr_lower:
                raise ValueError(f"Dangerous keyword '{keyword}' not allowed in expression")

    def validate_params(self):
        """Validate Filter parameters."""
        if not isinstance(self.data_input, (ToolOutput, StandardizedOutput)):
            raise ValueError("data must be a ToolOutput or StandardizedOutput object")

        if self.pool_output and not isinstance(self.pool_output, (ToolOutput, StandardizedOutput)):
            raise ValueError("pool must be a ToolOutput or StandardizedOutput object")

        if self.max_items is not None and self.max_items <= 0:
            raise ValueError("max_items must be positive")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """
        Configure input datasheet and pool from previous tools.

        REFACTORED: Uses DatasheetNavigatorMixin for elegant datasheet access.
        OLD: 40+ lines of nested if-else
        NEW: 8 lines with get_datasheet() mixin!
        """
        self.folders = pipeline_folders

        # ========================================================================
        # ELEGANT DATASHEET ACCESS using mixin
        # ========================================================================
        try:
            # Get datasheet info (handles all formats automatically)
            self.input_datasheet_info = self.get_datasheet(
                self.data_input,
                name='filtered',
                fallback_names=['structures', 'sequences', 'main']
            )

            # Get path (works with DatasheetInfo, dict, or string)
            self.input_csv_path = self.get_datasheet_path(
                self.data_input,
                name='filtered',
                fallback_names=['structures', 'sequences', 'main']
            )

            # Get datasheet name
            if hasattr(self.input_datasheet_info, 'name'):
                self.input_datasheet_name = self.input_datasheet_info.name
            else:
                self.input_datasheet_name = 'filtered'

        except ValueError as e:
            raise ValueError(f"Could not find datasheet in data input: {e}")

        # Configure pool if in pool mode
        if self.use_pool_mode:
            if hasattr(self.pool_output, 'output_folder'):
                self.pool_folder = self.pool_output.output_folder
            else:
                raise ValueError("Pool must provide output_folder")

    def _get_input_columns(self) -> List[str]:
        """Get column names from the input datasheet."""
        if hasattr(self.input_datasheet_info, 'columns') and self.input_datasheet_info.columns:
            return self.input_datasheet_info.columns
        return []

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
        output_folder = self.output_folder
        os.makedirs(output_folder, exist_ok=True)

        # Output CSV path - use same name as input datasheet
        output_datasheet_name = getattr(self, 'input_datasheet_name', 'filtered')
        filtered_csv_name = f"{output_datasheet_name}.csv"
        filtered_csv = os.path.join(output_folder, filtered_csv_name)

        # Create config file for the filter
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

        import json
        with open(self.config_file, 'w') as f:  # Auto-managed path!
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
python "{self.filter_script}" \\
  --config "{self.config_file}"

if [ $? -eq 0 ]; then
    echo "Results written to: {filtered_csv}"
else
    echo "Error: Filtering failed"
    exit 1
fi

{self.generate_completion_check_footer()}
"""

        return script_content

    def get_output_files(self) -> Dict[str, List[str]]:
        """Get expected output files after filtering."""
        # Get input columns
        input_columns = self._get_input_columns()

        # Output CSV paths
        output_datasheet_name = getattr(self, 'input_datasheet_name', 'filtered')
        filtered_csv_name = f"{output_datasheet_name}.csv"
        filtered_csv = os.path.join(self.output_folder, filtered_csv_name)
        missing_csv = os.path.join(self.output_folder, "missing.csv")

        # Create datasheets
        datasheets = {
            output_datasheet_name: DatasheetInfo(
                name=output_datasheet_name,
                path=filtered_csv,
                columns=input_columns,
                description=f"Filtered results using expression: {self.expression}",
                count="variable"
            ),
            "missing": DatasheetInfo(
                name="missing",
                path=missing_csv,
                columns=["id", "structure", "msa"],
                description="IDs that were filtered out",
                count="variable"
            )
        }

        if self.use_pool_mode:
            # Pool mode: predict copying ALL pool files
            updated_structures = []
            updated_compounds = []
            updated_sequences = []

            # Copy structure file paths
            if hasattr(self.pool_output, 'structures') and self.pool_output.structures:
                for struct_path in self.pool_output.structures:
                    filename = os.path.basename(struct_path)
                    updated_structures.append(os.path.join(self.output_folder, filename))

            # Copy compound file paths
            if hasattr(self.pool_output, 'compounds') and self.pool_output.compounds:
                for comp_path in self.pool_output.compounds:
                    filename = os.path.basename(comp_path)
                    updated_compounds.append(os.path.join(self.output_folder, filename))

            # Copy sequence file paths
            if hasattr(self.pool_output, 'sequences') and self.pool_output.sequences:
                for seq_path in self.pool_output.sequences:
                    filename = os.path.basename(seq_path)
                    updated_sequences.append(os.path.join(self.output_folder, filename))

            # Combine datasheets
            combined_datasheets = datasheets.copy()
            if hasattr(self.pool_output, 'datasheets') and hasattr(self.pool_output.datasheets, '_datasheets'):
                for name, info in self.pool_output.datasheets._datasheets.items():
                    filename = os.path.basename(info.path)
                    combined_datasheets[name] = DatasheetInfo(
                        name=name,
                        path=os.path.join(self.output_folder, filename),
                        columns=info.columns,
                        description=info.description,
                        count=info.count
                    )

            return {
                "structures": updated_structures,
                "structure_ids": self.pool_output.structure_ids if hasattr(self.pool_output, 'structure_ids') else [],
                "compounds": updated_compounds,
                "compound_ids": self.pool_output.compound_ids if hasattr(self.pool_output, 'compound_ids') else [],
                "sequences": updated_sequences,
                "sequence_ids": self.pool_output.sequence_ids if hasattr(self.pool_output, 'sequence_ids') else [],
                "datasheets": combined_datasheets,
                "output_folder": self.output_folder
            }
        else:
            # Data mode: only filtered CSV
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
