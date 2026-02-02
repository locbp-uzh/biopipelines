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
        Panda.merge(on="id", prefixes=["apo_", "holo_"])
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
            table=boltz.tables.confidence,
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
                Panda.merge(on="id", prefixes=["apo_", "holo_"]),
                Panda.calculate({"delta": "holo_affinity - apo_affinity"})
            ]
        )

        # Multi-table concat (replaces ConcatenateTables)
        combined = Panda(
            tables=[cycle0.tables.sequences, cycle1.tables.sequences],
            operations=[Panda.concat(fill="")]
        )
    """

    TOOL_NAME = "Panda"

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

        Args:
            exprs: Dictionary mapping new column names to expressions

        Returns:
            Operation to add calculated columns

        Example:
            Panda.calculate({
                "delta": "holo_affinity - apo_affinity",
                "normalized": "score / max_score"
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
    def merge(on: str = "id", how: str = "outer",
              prefixes: Optional[List[str]] = None,
              id_map: Optional[Dict[str, List[str]]] = None) -> Operation:
        """
        Merge multiple tables horizontally (like SQL JOIN).

        This operation combines columns from multiple tables based on a common key.
        Requires multi-table input (tables parameter instead of table).

        Args:
            on: Column name to merge on (must exist in all tables)
            how: Join type ("inner", "outer", "left", "right")
            prefixes: List of prefixes for each table's columns (prevents name collisions)
            id_map: Dictionary mapping new_id -> [old_id1, old_id2, ...] to consolidate IDs

        Returns:
            Operation to merge tables

        Example:
            Panda.merge(on="id", prefixes=["apo_", "holo_"])
            Panda.merge(on="structure_id", how="inner", id_map={"common": ["id_a", "id_b"]})
        """
        return Operation(type="merge", params={
            "on": on, "how": how, "prefixes": prefixes, "id_map": id_map
        })

    @staticmethod
    def concat(fill: Optional[str] = None, add_source: bool = True) -> Operation:
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
                 table: Union[TableInfo, StandardizedOutput, str, None] = None,
                 tables: Optional[List[Union[TableInfo, StandardizedOutput, str]]] = None,
                 operations: Optional[List[Operation]] = None,
                 pool: Optional[Union[StandardizedOutput, List[StandardizedOutput]]] = None,
                 rename: Optional[str] = None,
                 **kwargs):
        """
        Initialize Panda tool.

        Args:
            table: Single table input (mutually exclusive with tables)
            tables: Multiple table inputs for merge/concat operations
            operations: List of operations to apply sequentially
            pool: Tool output(s) for pool mode - structures matching filtered IDs will be copied.
                  Can be a single pool or a list of pools matching `tables` for multi-pool selection.
                  When using multiple pools with concat(add_source=True), files are copied from
                  the pool corresponding to each row's source_table index.
            rename: If provided, output IDs will be renamed to {rename}_1, {rename}_2, etc.
                    Useful after sorting to get ranked IDs (e.g., rename="best" -> best_1, best_2, ...)
            **kwargs: Additional parameters

        Examples:
            # Single table with operations
            result = Panda(
                table=boltz.tables.confidence,
                operations=[Panda.filter("pLDDT > 80"), Panda.sort("score")]
            )

            # Multi-table merge
            merged = Panda(
                tables=[apo.tables.affinity, holo.tables.affinity],
                operations=[Panda.merge(on="id", prefixes=["apo_", "holo_"])]
            )

            # With pool mode - copy structures matching filtered IDs
            filtered = Panda(
                table=combined.tables.merged,
                operations=[Panda.filter("delta > 0")],
                pool=boltz_output
            )

            # Sort and rename to get ranked output
            ranked = Panda(
                table=boltz.tables.confidence,
                operations=[Panda.sort("confidence_score", ascending=False)],
                rename="best",  # Output will have IDs: best_1, best_2, ...
                pool=boltz
            )

            # Multi-pool selection (replaces SelectBest) - select best from multiple cycles
            best = Panda(
                tables=[cycle1.tables.result, cycle2.tables.result, cycle3.tables.result],
                operations=[
                    Panda.concat(add_source=True),
                    Panda.sort("metric", ascending=True),
                    Panda.head(1)
                ],
                pool=[cycle1_output, cycle2_output, cycle3_output],  # Pools match tables
                rename="best"
            )
        """
        # Validate mutual exclusivity
        if table is not None and tables is not None:
            raise ValueError("Cannot specify both 'table' and 'tables'. Use 'table' for single table, 'tables' for multiple.")

        if table is None and tables is None:
            raise ValueError("Must specify either 'table' (single) or 'tables' (multiple)")

        self.table_input = table
        self.tables_input = tables or []
        self.operations = operations or []
        self.rename = rename

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

        if (has_merge or has_concat) and self.table_input is not None:
            raise ValueError(
                f"Operations 'merge' and 'concat' require multiple tables. "
                f"Use 'tables=[...]' instead of 'table=...'"
            )

        if not (has_merge or has_concat) and self.tables_input:
            # If multiple tables but no merge/concat, use first table only
            # (user might want to do other operations after merge/concat)
            pass

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

        if self.table_input is not None:
            # Single table mode
            self.input_csv_paths = [self._resolve_table_path(self.table_input)]
        else:
            # Multi-table mode
            self.input_csv_paths = [self._resolve_table_path(t) for t in self.tables_input]

        # Configure pool mode - collect all pool folders
        if self.use_pool_mode:
            self.pool_folders = []
            for pool in self.pool_outputs:
                if hasattr(pool, 'output_folder'):
                    self.pool_folders.append(pool.output_folder)
                else:
                    raise ValueError("Each pool must have an output_folder attribute")

    def _resolve_table_path(self, table_input: Any) -> str:
        """
        Resolve table input to a CSV file path.

        Args:
            table_input: TableInfo, StandardizedOutput, or string path

        Returns:
            Path to CSV file
        """
        # Handle tuple from column reference (table_info, column_name)
        if isinstance(table_input, tuple) and len(table_input) == 2:
            table_info, _ = table_input
            if hasattr(table_info, 'path'):
                return table_info.path

        # Direct path string
        if isinstance(table_input, str):
            return table_input

        # TableInfo object
        if isinstance(table_input, TableInfo):
            return table_input.path

        # StandardizedOutput object
        if hasattr(table_input, 'tables'):
            tables = table_input.tables
            if hasattr(tables, '_tables'):
                if len(tables._tables) == 1:
                    ds_info = next(iter(tables._tables.values()))
                    if hasattr(ds_info, 'path'):
                        return ds_info.path
                else:
                    table_names = list(tables._tables.keys())
                    raise ValueError(
                        f"Ambiguous: StandardizedOutput has multiple tables {table_names}. "
                        f"Use explicit specification like table=tool.tables.affinity"
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

        if self.table_input is not None:
            if isinstance(self.table_input, TableInfo) and self.table_input.columns:
                columns = list(self.table_input.columns)
        elif self.tables_input:
            # For merge, we'd have columns from all tables
            # For concat, we'd have common columns
            first_table = self.tables_input[0]
            if isinstance(first_table, TableInfo) and first_table.columns:
                columns = list(first_table.columns)

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

        if self.table_input is not None:
            config_lines.append(f"INPUT: 1 table")
        else:
            config_lines.append(f"INPUTS: {len(self.tables_input)} tables")

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
        config_data = {
            "input_csvs": self.input_csv_paths,
            "operations": [op.to_dict() for op in self.operations],
            "output_csv": self.output_csv,
            "use_pool_mode": self.use_pool_mode,
            "pool_folders": getattr(self, 'pool_folders', []),
            "rename": self.rename
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
    exit 1
fi

"""
        script_content += self.generate_completion_check_footer()

        return script_content

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
                columns=["id", "structure", "msa"],
                description="IDs that were filtered out and their expected file paths",
                count="variable"
            )
        }

        if self.use_pool_mode and self.pool_outputs:
            # Pool mode: predict structures, compounds, sequences from pool(s)
            # Collect all IDs from all pools
            all_struct_ids = []
            all_comp_ids = []
            all_seq_ids = []

            for pool in self.pool_outputs:
                if pool.structures:
                    all_struct_ids.extend(pool.structures.ids)
                if pool.compounds:
                    all_comp_ids.extend(pool.compounds.ids)
                if pool.sequences:
                    all_seq_ids.extend(pool.sequences.ids)

            pool_struct_count = len(all_struct_ids)
            pool_comp_count = len(all_comp_ids)
            pool_seq_count = len(all_seq_ids)

            if self.rename:
                # Predict renamed IDs: {rename}_1, {rename}_2, ...
                struct_ids = [f"{self.rename}_{i+1}" for i in range(pool_struct_count)]
                comp_ids = [f"{self.rename}_{i+1}" for i in range(pool_comp_count)]
                seq_ids = [f"{self.rename}_{i+1}" for i in range(pool_seq_count)]

                # Predict renamed file paths
                updated_structures = [os.path.join(self.output_folder, f"{sid}.pdb") for sid in struct_ids]
                updated_compounds = [os.path.join(self.output_folder, f"{cid}.sdf") for cid in comp_ids]
                updated_sequences = [os.path.join(self.output_folder, f"{seqid}.fasta") for seqid in seq_ids]
            else:
                # Keep original IDs
                struct_ids = all_struct_ids
                comp_ids = all_comp_ids
                seq_ids = all_seq_ids

                # Keep original filenames
                updated_structures = []
                updated_compounds = []
                updated_sequences = []

                for pool in self.pool_outputs:
                    if pool.structures and len(pool.structures) > 0:
                        for struct_path in pool.structures.files:
                            filename = os.path.basename(struct_path)
                            updated_structures.append(os.path.join(self.output_folder, filename))

                    if pool.compounds and len(pool.compounds) > 0:
                        for comp_path in pool.compounds.files:
                            filename = os.path.basename(comp_path)
                            updated_compounds.append(os.path.join(self.output_folder, filename))

                    if pool.sequences and len(pool.sequences) > 0:
                        for seq_path in pool.sequences.files:
                            filename = os.path.basename(seq_path)
                            updated_sequences.append(os.path.join(self.output_folder, filename))

            # Include tables from first pool (they typically have same schema)
            first_pool = self.pool_outputs[0]
            if hasattr(first_pool, 'tables') and hasattr(first_pool.tables, '_tables'):
                for name, info in first_pool.tables._tables.items():
                    if name not in tables:  # Don't overwrite our result table
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
                ids=struct_ids,
                files=updated_structures,
                format="pdb"
            )
            compounds = DataStream(
                name="compounds",
                ids=comp_ids,
                files=updated_compounds,
                format="sdf"
            )
            # For sequences, use the result CSV as map_table (contains id and sequence columns)
            # This allows Boltz2 and other tools to read sequences from the filtered result
            sequences = DataStream(
                name="sequences",
                ids=seq_ids,
                files=updated_sequences,
                map_table=self.output_csv,  # Result CSV serves as sequence map
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
                "num_inputs": len(self.input_csv_paths) if hasattr(self, 'input_csv_paths') else 1,
                "operations": [op.to_dict() for op in self.operations],
                "pool_mode": self.use_pool_mode
            }
        })
        return base_dict
