# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Panda tool for unified pandas-style table transformations.

Provides a single tool that accepts TableInfo objects directly and supports
pandas-style transformations including filtering, sorting, merging, concatenation,
and calculated columns.
"""

import os
import json
from dataclasses import dataclass, field
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


@dataclass
class Operation:
    """
    Represents a single pandas operation to be applied to a dataframe.

    Operations are created via static methods on the Panda class:
        Panda.filter("pLDDT > 80")
        Panda.sort("affinity", ascending=False)
        Panda.merge(prefixes=["apo_", "holo_"])
    """

    type: str
    params: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        """Convert operation to dictionary for JSON serialization."""
        return {"type": self.type, "params": self.params}

    def __repr__(self) -> str:
        params_str = ", ".join(f"{k}={v!r}" for k, v in self.params.items())
        return f"Operation({self.type}, {params_str})"


class Panda(BaseConfig):
    """
    Unified tool for pandas-style table transformations.

    Supports single table operations (filter, sort, head, tail, etc.) and
    multi-table operations (merge, concat) with a consistent operation-based API.

    Examples:
        # Single table operations
        result = Panda(
            tables=boltz.tables.confidence,
            operations=[
                Panda.filter("pLDDT > 80"),
                Panda.sort("affinity", ascending=False),
                Panda.head(10)
            ]
        )

        # Multi-table merge (replaces MergeTables)
        merged = Panda(
            tables=[apo.tables.affinity, holo.tables.affinity],
            operations=[
                Panda.merge(prefixes=["apo_", "holo_"]),
                Panda.calculate({"delta": "holo_affinity - apo_affinity"})
            ]
        )

        # Calculate with math functions (cos, sin, sqrt, log, exp, radians, degrees, pi, ...)
        fret = Panda(
            tables=[distances.tables.result, angles.tables.angles],
            operations=[
                Panda.merge(),
                Panda.calculate({
                    "kappa2": "cos(orientation) ** 2",
                    "R0_eff": "49.0 * (kappa2 / 0.6667) ** (1.0 / 6.0)",
                    "efficiency": "1 / (1 + (distance / R0_eff) ** 6)"
                })
            ]
        )

        # Multi-table concat (replaces ConcatenateTables)
        combined = Panda(
            tables=[cycle0.tables.sequences, cycle1.tables.sequences],
            operations=[Panda.concat(fill="")]
        )
    """

    TOOL_NAME = "Panda"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        return """echo "=== Panda ==="
echo "Uses biopipelines environment (no additional installation needed)."
echo "=== Panda ready ==="
"""

    # Path descriptors
    output_csv = Path(lambda self: os.path.join(self.output_folder, f"{self._get_output_name()}.csv"))
    config_file = Path(lambda self: os.path.join(self.output_folder, "panda_config.json"))
    panda_py = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_panda.py"))
    missing_csv = Path(lambda self: os.path.join(self.output_folder, "missing.csv"))
    sequences_csv = Path(lambda self: os.path.join(self.output_folder, "sequences.csv"))

    # ========== Static methods for creating operations ==========

    @staticmethod
    def filter(expr: str) -> Operation:
        """
        Filter rows using a pandas query expression.

        Args:
            expr: Pandas query expression (e.g., "pLDDT > 80 and distance < 5.0")

        Returns:
            Operation to filter rows

        Example:
            Panda.filter("pLDDT > 80 and contacts >= 3")
        """
        return Operation(type="filter", params={"expr": expr})

    @staticmethod
    def head(n: int) -> Operation:
        """
        Keep only the first N rows.

        Args:
            n: Number of rows to keep

        Returns:
            Operation to keep first N rows

        Example:
            Panda.head(10)
        """
        return Operation(type="head", params={"n": n})

    @staticmethod
    def tail(n: int) -> Operation:
        """
        Keep only the last N rows.

        Args:
            n: Number of rows to keep

        Returns:
            Operation to keep last N rows

        Example:
            Panda.tail(5)
        """
        return Operation(type="tail", params={"n": n})

    @staticmethod
    def sample(n: Optional[int] = None, frac: Optional[float] = None,
               random_state: Optional[int] = None) -> Operation:
        """
        Randomly sample rows.

        Args:
            n: Number of rows to sample (mutually exclusive with frac)
            frac: Fraction of rows to sample (mutually exclusive with n)
            random_state: Random seed for reproducibility

        Returns:
            Operation to sample rows

        Example:
            Panda.sample(n=100)
            Panda.sample(frac=0.1, random_state=42)
        """
        return Operation(type="sample", params={"n": n, "frac": frac, "random_state": random_state})

    @staticmethod
    def drop_duplicates(subset: Optional[Union[str, List[str]]] = None,
                        keep: str = "first") -> Operation:
        """
        Remove duplicate rows.

        Args:
            subset: Column(s) to consider for duplicates (None = all columns)
            keep: Which duplicate to keep ("first", "last", or False for none)

        Returns:
            Operation to remove duplicates

        Example:
            Panda.drop_duplicates(subset="sequence", keep="first")
        """
        return Operation(type="drop_duplicates", params={"subset": subset, "keep": keep})

    @staticmethod
    def sort(by: Union[str, List[str]], ascending: Union[bool, List[bool]] = True) -> Operation:
        """
        Sort rows by column(s).

        Args:
            by: Column name(s) to sort by
            ascending: Sort order (True for ascending)

        Returns:
            Operation to sort rows

        Example:
            Panda.sort("affinity", ascending=False)
            Panda.sort(["group", "score"], ascending=[True, False])
        """
        return Operation(type="sort", params={"by": by, "ascending": ascending})

    @staticmethod
    def rank(by: str, prefix: str = "rank_", ascending: bool = True) -> Operation:
        """
        Assign rank IDs based on a column.

        Args:
            by: Column name to rank by
            prefix: Prefix for the rank column name
            ascending: Rank order (True for ascending)

        Returns:
            Operation to add rank column

        Example:
            Panda.rank(by="affinity", ascending=False)
        """
        return Operation(type="rank", params={"by": by, "prefix": prefix, "ascending": ascending})

    @staticmethod
    def select_columns(cols: List[str]) -> Operation:
        """
        Keep only specified columns.

        Args:
            cols: List of column names to keep

        Returns:
            Operation to select columns

        Example:
            Panda.select_columns(["id", "sequence", "score"])
        """
        return Operation(type="select_columns", params={"cols": cols})

    @staticmethod
    def drop_columns(cols: List[str]) -> Operation:
        """
        Remove specified columns.

        Args:
            cols: List of column names to drop

        Returns:
            Operation to drop columns

        Example:
            Panda.drop_columns(["temp_col", "debug_info"])
        """
        return Operation(type="drop_columns", params={"cols": cols})

    @staticmethod
    def rename(mapping: Dict[str, str]) -> Operation:
        """
        Rename columns.

        Args:
            mapping: Dictionary mapping old names to new names

        Returns:
            Operation to rename columns

        Example:
            Panda.rename({"old_name": "new_name", "score": "affinity"})
        """
        return Operation(type="rename", params={"mapping": mapping})

    @staticmethod
    def calculate(exprs: Dict[str, str]) -> Operation:
        """
        Add calculated columns using pandas eval expressions.

        Supports standard arithmetic (+, -, *, /, **, %) and math functions:
        cos, sin, tan, arccos, arcsin, arctan, arctan2, sqrt, abs, log, exp,
        radians, degrees, pi.

        Expressions can reference columns defined earlier in the same calculate call.

        Args:
            exprs: Dictionary mapping new column names to expressions

        Returns:
            Operation to add calculated columns

        Example:
            Panda.calculate({
                "delta": "holo_affinity - apo_affinity",
                "normalized": "score / max_score",
                "angle_rad": "radians(angle_deg)",
                "kappa2": "cos(angle_rad) ** 2",
                "efficiency": "1 / (1 + (distance / R0_eff) ** 6)"
            })
        """
        return Operation(type="calculate", params={"exprs": exprs})

    @staticmethod
    def fillna(value: Any = None, column: Optional[str] = None) -> Operation:
        """
        Fill missing values.

        Args:
            value: Value to fill NA/NaN with
            column: Specific column to fill (None = all columns)

        Returns:
            Operation to fill missing values

        Example:
            Panda.fillna(0)
            Panda.fillna("unknown", column="category")
        """
        return Operation(type="fillna", params={"value": value, "column": column})

    @staticmethod
    def merge(on: Optional[Union[str, List[str]]] = None, how: str = "outer",
              prefixes: Optional[List[str]] = None) -> Operation:
        """
        Merge multiple tables horizontally (like SQL JOIN).

        This operation combines columns from multiple tables based on a common key.
        Requires multi-table input (tables parameter instead of table).

        Args:
            on: Column name(s) to merge on. Can be:
                - None (default): uses biopipelines ID matching (get_mapped_ids)
                  which handles exact, child/parent, sibling, and provenance matching
                - A single string: exact pandas merge on that column (e.g., "id")
                - A list of strings (one per table): each table is joined on its
                  own column, which is renamed to the first entry before merging
                  (e.g., ["id", "id", "structures.id"])
            how: Join type ("inner", "outer", "left", "right")
            prefixes: List of prefixes for each table's columns (prevents name collisions)

        Returns:
            Operation to merge tables

        Example:
            Panda.merge(prefixes=["apo_", "holo_"])
            Panda.merge(on=["id", "id", "structures.id"])
        """
        return Operation(type="merge", params={
            "on": on, "how": how, "prefixes": prefixes
        })

    @staticmethod
    def concat(fill: Optional[str] = "", add_source: bool = True) -> Operation:
        """
        Concatenate multiple tables vertically (like SQL UNION).

        This operation stacks rows from multiple tables.
        Requires multi-table input (tables parameter instead of table).

        Args:
            fill: Value to fill missing columns (None = remove non-common columns)
            add_source: Add a "source_table" column indicating origin

        Returns:
            Operation to concatenate tables

        Example:
            Panda.concat(fill="")
            Panda.concat(fill=0, add_source=False)
        """
        return Operation(type="concat", params={"fill": fill, "add_source": add_source})

    @staticmethod
    def groupby(by: Union[str, List[str]], agg: Dict[str, str]) -> Operation:
        """
        Group by column(s) and aggregate.

        Args:
            by: Column(s) to group by
            agg: Dictionary mapping column names to aggregation functions
                 (e.g., "sum", "mean", "min", "max", "count", "first", "last")

        Returns:
            Operation to group and aggregate

        Example:
            Panda.groupby("category", {"score": "mean", "count": "count"})
        """
        return Operation(type="groupby", params={"by": by, "agg": agg})

    @staticmethod
    def pivot(index: str, columns: str, values: str,
              aggfunc: str = "first") -> Operation:
        """
        Pivot table (wide format).

        Args:
            index: Column to use as index
            columns: Column to use for new column headers
            values: Column to use for values
            aggfunc: Aggregation function if duplicates exist

        Returns:
            Operation to pivot table

        Example:
            Panda.pivot(index="sample_id", columns="metric", values="value")
        """
        return Operation(type="pivot", params={
            "index": index, "columns": columns, "values": values, "aggfunc": aggfunc
        })

    @staticmethod
    def melt(id_vars: Union[str, List[str]],
             value_vars: Optional[List[str]] = None,
             var_name: str = "variable",
             value_name: str = "value") -> Operation:
        """
        Unpivot table (long format).

        Args:
            id_vars: Column(s) to keep as identifiers
            value_vars: Column(s) to unpivot (None = all non-id columns)
            var_name: Name for the variable column
            value_name: Name for the value column

        Returns:
            Operation to melt table

        Example:
            Panda.melt(id_vars="id", value_vars=["score1", "score2"])
        """
        return Operation(type="melt", params={
            "id_vars": id_vars, "value_vars": value_vars,
            "var_name": var_name, "value_name": value_name
        })

    @staticmethod
    def average_by_source(source_col: str = "source_table") -> Operation:
        """
        Average all numeric columns per source table (replaces AverageByTable).

        This operation is typically used after concat to compute the mean of all
        numeric columns for each input table. Each source table becomes one row
        in the output.

        Args:
            source_col: Column name identifying the source table (default: "source_table")

        Returns:
            Operation to average by source

        Example:
            # Replaces AverageByTable
            Panda(
                tables=[cycle0.tables.merged, cycle1.tables.merged, cycle2.tables.merged],
                operations=[
                    Panda.concat(add_source=True),
                    Panda.average_by_source()
                ]
            )
        """
        return Operation(type="average_by_source", params={"source_col": source_col})

    # ========== Tool implementation ==========

    def __init__(self,
                 tables: Union[TableInfo, StandardizedOutput, str,
                               List[Union[TableInfo, StandardizedOutput, str]], None] = None,
                 operations: Optional[List[Operation]] = None,
                 pool: Optional[Union[StandardizedOutput, List[StandardizedOutput]]] = None,
                 rename: Optional[str] = None,
                 ignore_missing: bool = True,
                 **kwargs):
        """
        Initialize Panda tool.

        Args:
            tables: Table input(s). Can be a single table or a list of tables.
                    Single table: Panda(tables=boltz.tables.confidence, ...)
                    Multiple tables: Panda(tables=[apo.tables.affinity, holo.tables.affinity], ...)
            operations: List of operations to apply sequentially
            pool: Tool output(s) for pool mode - structures matching filtered IDs will be copied.
                  Can be a single pool or a list of pools matching `tables` for multi-pool selection.
                  When using multiple pools with concat(add_source=True), files are copied from
                  the pool corresponding to each row's source_table index.
            rename: If provided, output IDs will be renamed to {rename}_1, {rename}_2, etc.
                    Useful after sorting to get ranked IDs (e.g., rename="best" -> best_1, best_2, ...)
            ignore_missing: If True (default), skip missing pool files with a warning instead
                    of failing. Set to False to raise an error on any missing file.
            **kwargs: Additional parameters

        Output:
            Streams: inherits all streams from pool input (if pool mode)
            Tables:
                result: columns derived from input + applied operations (dynamic)
                missing: id | removed_by | cause
        """
        if tables is None:
            raise ValueError("Must specify 'tables' (single table or list of tables)")

        # Normalize: wrap single table in a list
        if isinstance(tables, list):
            self.tables_input = tables
        else:
            self.tables_input = [tables]

        # Keep table_input for internal backwards compat (single table = first element)
        self.table_input = self.tables_input[0] if len(self.tables_input) == 1 else None

        self.operations = operations or []
        self.rename = rename
        self.ignore_missing = ignore_missing

        # Handle pool - can be single or list
        if pool is None:
            self.pool_outputs = []
            self.use_pool_mode = False
        elif isinstance(pool, list):
            self.pool_outputs = pool
            self.use_pool_mode = True
        else:
            self.pool_outputs = [pool]
            self.use_pool_mode = True

        # Validate operations
        self._validate_operations()

        super().__init__(**kwargs)

    def _validate_operations(self):
        """Validate that operations are compatible with input type."""
        if not self.operations:
            return

        has_merge = any(op.type == "merge" for op in self.operations)
        has_concat = any(op.type == "concat" for op in self.operations)

        if (has_merge or has_concat) and len(self.tables_input) < 2:
            raise ValueError(
                f"Operations 'merge' and 'concat' require multiple tables. "
                f"Use 'tables=[table1, table2, ...]'"
            )

    def validate_params(self):
        """Validate Panda parameters."""
        # Validate operation types
        valid_types = {
            "filter", "head", "tail", "sample", "drop_duplicates",
            "sort", "rank", "select_columns", "drop_columns", "rename",
            "calculate", "fillna", "merge", "concat", "groupby", "pivot", "melt",
            "average_by_source"
        }

        for op in self.operations:
            if op.type not in valid_types:
                raise ValueError(f"Unknown operation type: {op.type}")

        # Validate filter expressions for safety
        for op in self.operations:
            if op.type == "filter":
                self._validate_expression(op.params.get("expr", ""))
            elif op.type == "calculate":
                for expr in op.params.get("exprs", {}).values():
                    self._validate_expression(expr)

    def _validate_expression(self, expr: str):
        """Validate that an expression is safe for pandas eval/query."""
        if not expr or not expr.strip():
            raise ValueError("Expression cannot be empty")

        # Check for dangerous keywords
        dangerous_keywords = ['import', 'exec', 'eval', '__', 'os.', 'sys.', 'subprocess', 'open(']
        expr_lower = expr.lower()
        for keyword in dangerous_keywords:
            if keyword in expr_lower:
                raise ValueError(f"Dangerous keyword '{keyword}' not allowed in expression")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input tables from previous tools."""
        self.folders = pipeline_folders

        self.input_csv_paths = [self._resolve_table_path(t) for t in self.tables_input]

        # Collect map_table paths from pool outputs for ID matching at runtime
        self.map_table_paths = []
        if self.use_pool_mode:
            for pool in self.pool_outputs:
                if hasattr(pool, 'streams') and pool.streams:
                    for stream_name in pool.streams.keys():
                        stream = pool.streams.get(stream_name)
                        if stream and stream.map_table:
                            self.map_table_paths.append(stream.map_table)

        # Auto-rename for pool mode:
        # - sort/head/tail/sample present: rename (output IDs are unpredictable at configuration time)
        # - Otherwise: no rename, predict all input IDs (deduplicated across pools),
        #   missing ones from filter/etc. are tracked in missing.csv
        if self.use_pool_mode and not self.rename:
            needs_rename = any(op.type in ("sort", "head", "tail", "sample") for op in self.operations)

            if needs_rename:
                # Derive rename prefix from step folder name (e.g., "010_Panda_Cycle1" -> "Panda_Cycle1")
                folder_name = os.path.basename(self.output_folder)
                parts = folder_name.split("_", 1)
                if len(parts) > 1 and parts[0].isdigit():
                    self.rename = parts[1]  # e.g., "Panda_Cycle1"
                else:
                    self.rename = folder_name

        # Configure pool mode - collect all pool folders and file mappings
        if self.use_pool_mode:
            self.pool_folders = []
            self.pool_stream_jsons = []  # List of {stream_name: json_path} dicts per pool
            self.pool_table_maps = []  # List of {table_name: {"path": str, "columns": list}} for each pool
            for pool_idx, pool in enumerate(self.pool_outputs):
                if hasattr(pool, 'output_folder'):
                    self.pool_folders.append(pool.output_folder)
                else:
                    raise ValueError("Each pool must have an output_folder attribute")

                # Save each pool stream's DataStream to JSON for runtime expansion
                stream_jsons = {}
                if hasattr(pool, 'streams'):
                    for stream_name in pool.streams.keys():
                        stream = pool.streams.get(stream_name)
                        if stream and len(stream) > 0 and len(stream.files) > 0:
                            json_path = os.path.join(
                                self.output_folder,
                                f"pool_{pool_idx}_{stream_name}_ds.json"
                            )
                            stream.save_json(json_path)
                            stream_jsons[stream_name] = json_path
                self.pool_stream_jsons.append(stream_jsons)

                # Build table map from pool tables (for filtering/copying at execution time)
                table_map = {}
                if hasattr(pool, 'tables') and hasattr(pool.tables, '_tables'):
                    for table_name, table_info in pool.tables._tables.items():
                        if table_name not in ('result', 'missing'):  # Skip Panda's own tables
                            table_map[table_name] = {
                                "path": table_info.info.path,
                                "columns": list(table_info.info.columns) if table_info.info.columns else []
                            }
                self.pool_table_maps.append(table_map)

    def _resolve_table_path(self, table_input: Any) -> str:
        """
        Resolve table input to a CSV file path.

        Args:
            table_input: TableInfo, StandardizedOutput, or string path

        Returns:
            Path to CSV file
        """
        # TableInfo object (must check before the duck-type TableReference check
        # below, because TableInfo.__getattr__ returns TableReference for any
        # attribute access including .path and .column)
        if isinstance(table_input, TableInfo):
            return table_input.info.path

        # Handle TableReference from column access (table.column)
        if hasattr(table_input, 'path') and hasattr(table_input, 'column'):
            return table_input.path

        # Direct path string
        if isinstance(table_input, str):
            return table_input

        # StandardizedOutput object
        if hasattr(table_input, 'tables'):
            tables = table_input.tables
            if hasattr(tables, '_tables'):
                if len(tables._tables) == 1:
                    ds_info = next(iter(tables._tables.values()))
                    if hasattr(ds_info, 'info'):
                        return ds_info.info.path
                else:
                    table_names = list(tables._tables.keys())
                    raise ValueError(
                        f"Ambiguous: StandardizedOutput has multiple tables {table_names}. "
                        f"Did you mean <tool>.tables.{table_names[0]}?"
                    )

        # Try output folder prediction
        if hasattr(table_input, 'output_folder'):
            return os.path.join(table_input.output_folder, 'results.csv')

        raise ValueError(f"Cannot resolve table path from: {type(table_input)}")

    def _get_output_name(self) -> str:
        """Auto-generate output name from operations."""
        if not self.operations:
            return "result"

        op_names = [op.type for op in self.operations[:3]]
        return "_".join(op_names)

    def _get_expected_columns(self) -> List[str]:
        """Predict output columns based on input and operations."""
        # Start with input columns if available
        columns = []

        if self.tables_input:
            first_table = self.tables_input[0]
            if isinstance(first_table, TableInfo) and first_table.info.columns:
                columns = list(first_table.info.columns)

        # Apply operation effects
        for op in self.operations:
            if op.type == "calculate":
                # Add calculated column names
                columns.extend(op.params.get("exprs", {}).keys())
            elif op.type == "select_columns":
                columns = op.params.get("cols", columns)
            elif op.type == "drop_columns":
                cols_to_drop = set(op.params.get("cols", []))
                columns = [c for c in columns if c not in cols_to_drop]
            elif op.type == "rename":
                mapping = op.params.get("mapping", {})
                columns = [mapping.get(c, c) for c in columns]
            elif op.type == "rank":
                prefix = op.params.get("prefix", "rank_")
                by = op.params.get("by", "")
                columns.append(f"{prefix}{by}")
            elif op.type == "concat":
                if op.params.get("add_source", True):
                    if "source_table" not in columns:
                        columns.append("source_table")

        return columns if columns else ["id"]

    def get_config_display(self) -> List[str]:
        """Get configuration display lines."""
        config_lines = super().get_config_display()

        n_tables = len(self.tables_input)
        config_lines.append(f"INPUT{'S' if n_tables > 1 else ''}: {n_tables} table{'s' if n_tables > 1 else ''}")

        config_lines.append(f"OPERATIONS: {len(self.operations)}")
        for i, op in enumerate(self.operations):
            params_str = ", ".join(f"{k}={v}" for k, v in op.params.items() if v is not None)
            config_lines.append(f"  {i+1}. {op.type}({params_str})")

        if self.use_pool_mode:
            config_lines.append(f"POOL MODE: enabled")

        if self.rename:
            config_lines.append(f"RENAME: {self.rename}_N")

        return config_lines

    def generate_script(self, script_path: str) -> str:
        """Generate Panda execution script."""
        os.makedirs(self.output_folder, exist_ok=True)

        # Create config file
        step_tool_name = os.path.basename(self.output_folder)
        config_data = {
            "input_csvs": self.input_csv_paths,
            "operations": [op.to_dict() for op in self.operations],
            "output_csv": self.output_csv,
            "use_pool_mode": self.use_pool_mode,
            "pool_folders": getattr(self, 'pool_folders', []),
            "pool_stream_jsons": getattr(self, 'pool_stream_jsons', []),
            "pool_table_maps": getattr(self, 'pool_table_maps', []),
            "map_table_paths": getattr(self, 'map_table_paths', []),
            "rename": self.rename,
            "ignore_missing": self.ignore_missing,
            "step_tool_name": step_tool_name
        }

        with open(self.config_file, 'w') as f:
            json.dump(config_data, f, indent=2)

        script_content = "#!/bin/bash\n"
        script_content += "# Panda execution script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()

        ops_summary = ", ".join(op.type for op in self.operations)
        script_content += f"""echo "Panda table transformations"
echo "Input tables: {len(self.input_csv_paths)}"
echo "Operations: {ops_summary}"
echo "Output: {self.output_csv}"

python "{self.panda_py}" \\
  --config "{self.config_file}"

if [ $? -eq 0 ]; then
    echo "Successfully applied {len(self.operations)} operations"
    echo "Output written to: {self.output_csv}"
else
    echo "Error: Panda operations failed"
fi

"""
        script_content += self.generate_completion_check_footer()

        return script_content

    def _get_predicted_output_count(self) -> Optional[int]:
        """
        Predict the number of output rows based on operations.

        Returns:
            Predicted count if determinable (e.g., head/tail/sample with n), None otherwise.
        """
        for op in self.operations:
            if op.type == "head":
                return op.params.get("n")
            elif op.type == "tail":
                return op.params.get("n")
            elif op.type == "sample":
                return op.params.get("n")  # Only if n is specified, not frac
        return None

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files after transformation."""
        expected_columns = self._get_expected_columns()

        tables = {
            "result": TableInfo(
                name="result",
                path=self.output_csv,
                columns=expected_columns,
                description=f"Transformed table with operations: {', '.join(op.type for op in self.operations)}",
                count="variable"
            ),
            "missing": TableInfo(
                name="missing",
                path=self.missing_csv,
                columns=["id", "removed_by", "cause"],
                description="IDs removed by this tool with removal reason",
                count="variable"
            )
        }

        if self.use_pool_mode and self.pool_outputs:
            # Pool mode: dynamically collect all streams from all pools
            # Collect stream info: {stream_name: {"ids": [], "files": [], "format": str, "map_table": str}}
            stream_data = {}

            for pool in self.pool_outputs:
                if hasattr(pool, 'streams'):
                    for stream_name in pool.streams.keys():
                        stream = pool.streams.get(stream_name)
                        if not stream or len(stream) == 0:
                            continue
                        if stream_name not in stream_data:
                            stream_data[stream_name] = {
                                "ids": [],
                                "files": [],
                                "format": stream.format,
                                "map_table": stream.map_table or ""
                            }
                        stream_data[stream_name]["ids"].extend(list(stream.ids))
                        # Only collect files for file-based streams
                        if stream.is_file_based():
                            stream_data[stream_name]["files"].extend(stream.files)

            # Deduplicate IDs across pools (first occurrence wins)
            for sdata in stream_data.values():
                seen = set()
                deduped_ids = []
                deduped_files = []
                is_template = (len(sdata["files"]) == 1 and '<id>' in sdata["files"][0])
                if is_template:
                    # Template files: dedup IDs only, keep template as-is
                    for sid in sdata["ids"]:
                        if sid not in seen:
                            seen.add(sid)
                            deduped_ids.append(sid)
                    deduped_files = sdata["files"]
                else:
                    pairs = zip(sdata["ids"], sdata["files"]) if sdata["files"] else ((sid, None) for sid in sdata["ids"])
                    for sid, sf in pairs:
                        if sid not in seen:
                            seen.add(sid)
                            deduped_ids.append(sid)
                            if sf is not None:
                                deduped_files.append(sf)
                sdata["ids"] = deduped_ids
                sdata["files"] = deduped_files

            # Check if head/tail/sample limits the output count
            predicted_count = self._get_predicted_output_count()

            # Build output streams
            output_streams = {}
            for stream_name, data in stream_data.items():
                pool_count = len(data["ids"])
                # Use predicted count if available, otherwise use pool count
                num_items = min(predicted_count, pool_count) if predicted_count is not None else pool_count

                # Value-based streams (no files): propagate as-is with original map_table
                if not data["files"]:
                    output_streams[stream_name] = DataStream(
                        name=stream_name,
                        ids=data["ids"][:num_items],
                        files=[],
                        map_table=data.get("map_table") or "",
                        format=data["format"]
                    )
                    continue

                if self.rename:
                    # Predict renamed IDs and file paths based on predicted output count
                    new_ids = [f"{self.rename}_{i+1}" for i in range(num_items)]
                    new_files = []
                    ext = os.path.splitext(data["files"][0])[1]
                    for new_id in new_ids:
                        new_files.append(os.path.join(self.output_folder, f"{new_id}{ext}"))
                else:
                    new_ids = data["ids"][:num_items]
                    new_files = []
                    ext = os.path.splitext(data["files"][0])[1]
                    for new_id in new_ids:
                        new_files.append(os.path.join(self.output_folder, f"{new_id}{ext}"))

                # Preserve pool's map_table, pointing to the copy in output_folder
                pool_map_table = data.get("map_table", "")
                if pool_map_table:
                    map_table = os.path.join(self.output_folder, os.path.basename(pool_map_table))
                else:
                    map_table = None

                output_streams[stream_name] = DataStream(
                    name=stream_name,
                    ids=new_ids,
                    files=new_files,
                    map_table=map_table,
                    format=data["format"]
                )


            # Include tables from first pool (they typically have same schema)
            first_pool = self.pool_outputs[0]
            if hasattr(first_pool, 'tables') and hasattr(first_pool.tables, '_tables'):
                for name, info in first_pool.tables._tables.items():
                    if name not in tables:  # Don't overwrite our result table
                        filename = os.path.basename(info.info.path)
                        tables[name] = TableInfo(
                            name=name,
                            path=os.path.join(self.output_folder, filename),
                            columns=info.info.columns,
                            description=info.info.description,
                            count=info.info.count
                        )

            output_streams["tables"] = tables
            output_streams["output_folder"] = self.output_folder
            return output_streams
        else:
            return {
                "tables": tables,
                "output_folder": self.output_folder
            }

    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration."""
        base_dict = super().to_dict()
        base_dict.update({
            "tool_params": {
                "num_inputs": len(self.input_csv_paths) if hasattr(self, 'input_csv_paths') else 1,
                "operations": [op.to_dict() for op in self.operations],
                "pool_mode": self.use_pool_mode
            }
        })
        return base_dict
