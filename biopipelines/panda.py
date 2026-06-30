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
    TOOL_VERSION = "1.1"

    # Internal column name used to track which input table a row came from
    # when concat runs over multiple inputs. Auto-added before the operation
    # chain starts and stripped from the final result by pipe_panda. Tools
    # that need to reference it during the chain (groupby, etc.) should use
    # this constant.
    SOURCE = "bp_panda_source"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        return """echo "=== Panda ==="
echo "Uses biopipelines environment (no additional installation needed)."
touch "$INSTALL_SUCCESS"
echo "=== Panda ready ==="
"""

    # Path descriptors — result table lives in tables/ (standalone output);
    # config JSON in configuration/. In pool mode the per-stream files and
    # their map_tables land inside their respective stream folders.
    output_csv = Path(lambda self: self.table_path(self._get_output_name()))
    config_file = Path(lambda self: self.configuration_path("panda_config.json"))
    panda_py = Path(lambda self: self.pipe_script_path("pipe_panda.py"))
    missing_csv = Path(lambda self: self.table_path("missing"))
    sequences_csv = Path(lambda self: self.stream_path("sequences", "sequences.csv"))

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
    def zscore(columns: Union[str, List[str]], by: Optional[str] = None,
               sign: Optional[Dict[str, int]] = None, suffix: str = "_z") -> Operation:
        """
        Standardize column(s) to z-scores: (x - mean) / std (population std,
        ddof=0). Writes one new column per input, named ``<column><suffix>``.

        Use to put metrics on a common scale before combining them (e.g. a
        weighted sum via calculate), since raw metrics live on different scales.

        Args:
            columns: Column name or list of names to standardize.
            by: Optional grouping column — standardize within each group (e.g.
                z within each parent design) instead of globally.
            sign: Optional {column: -1} to flip a "lower is better" metric so
                higher z = better (e.g. aggregation score, affinity_pred_value).
            suffix: Suffix for the output columns (default "_z").

        Returns:
            Operation to add z-score column(s).

        Example:
            # combine on a common scale: higher = better, both equally weighted
            Panda.zscore(["plddt", "aggrescan_avg"], sign={"aggrescan_avg": -1})
            Panda.calculate({"score": "plddt_z + aggrescan_avg_z"})
        """
        cols = [columns] if isinstance(columns, str) else list(columns)
        return Operation(type="zscore", params={
            "columns": cols, "by": by, "sign": sign or {}, "suffix": suffix})

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
    def merge(on: Optional[Union[str, List[str]]] = None, how: Optional[str] = None,
              prefixes: Optional[List[str]] = None, grain: str = "finest") -> Operation:
        """
        Merge multiple tables horizontally (like SQL JOIN).

        This operation combines columns from multiple tables based on a common key.
        Requires multi-table input (tables parameter instead of table).

        Args:
            on: Column name(s) to merge on. Usually leave as None — biopipelines
                ID matching pairs the rows for you; you rarely need to pass this.
                - None (default): uses biopipelines ID matching (get_mapped_ids)
                  which handles exact, child/parent, sibling, and provenance matching
                - A single string: exact pandas merge on that column (e.g., "id")
                - A list of strings (one per table): each table is joined on its
                  own column, which is renamed to the first entry before merging
                  (e.g., ["id", "id", "structures.id"])
            grain: Which input decides the output row level when tables sit at
                different fan-out levels (e.g. one design has many sequences).
                Only applies to ID matching (on=None).
                - "finest" (default): the table at the deepest provenance ID level
                  drives the rows (ID count breaks ties); coarser tables broadcast
                  their values onto each fine-grain row (one design's mean_plddt
                  repeats across its sequences). Order-independent and lossless — the
                  natural choice for plotting.
                - "coarsest": the table at the shallowest provenance ID level drives
                  the rows (ID count breaks ties); finer tables collapse to one match
                  per coarse row. Aggregate the finer table first if you don't want
                  siblings dropped — a WARNING is emitted whenever rows are dropped.
            how: Join type ("inner", "outer", "left", "right"). Decides what happens
                to non-overlapping keys once the grain is fixed; rarely needs setting.
                Default None derives it from grain: a "left" join for both finest and
                coarsest so the driver keeps exactly its rows (unmatched rows from the
                other table are not resurrected with a null id). Pass how="outer"
                explicitly to keep them. For on=<column> it falls back to "outer".
            prefixes: List of prefixes for each table's columns (prevents name collisions)

        Returns:
            Operation to merge tables

        Example:
            Panda.merge(prefixes=["apo_", "holo_"])
            Panda.merge(on=["id", "id", "structures.id"])
        """
        return Operation(type="merge", params={
            "on": on, "how": how, "prefixes": prefixes, "grain": grain
        })

    @staticmethod
    def concat(fill: Optional[str] = "") -> Operation:
        """
        Concatenate multiple tables vertically (like SQL UNION).

        This operation stacks rows from multiple tables.
        Requires multi-table input (tables parameter instead of table).

        Source tracking is implicit: when there are >1 input tables and the
        operation chain contains a concat, every row is tagged with its
        origin table index in the internal ``Panda.SOURCE`` column for the
        duration of the chain. The column is stripped from the final result
        before it is written. To reference it inside the chain (e.g. for a
        groupby), use ``Panda.SOURCE`` as the column name.

        Args:
            fill: Value to fill missing columns (None = remove non-common columns)

        Returns:
            Operation to concatenate tables

        Example:
            Panda.concat()
            Panda.concat(fill=0)
            # Best per source: groupby Panda.SOURCE then take the head.
            Panda(tables=[a, b, c], operations=[
                Panda.concat(),
                Panda.sort("score", ascending=False),
                Panda.groupby(Panda.SOURCE, {"id": "first", "score": "first"}),
            ])
        """
        return Operation(type="concat", params={"fill": fill})

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
    def average_by_source(source_col: Optional[str] = None) -> Operation:
        """
        Average all numeric columns per source table (replaces AverageByTable).

        This operation is typically used after concat to compute the mean of all
        numeric columns for each input table. Each source table becomes one row
        in the output.

        Args:
            source_col: Column name identifying the source table. Defaults to
                       Panda.SOURCE — the implicit column auto-added when
                       concat runs over multiple inputs.

        Returns:
            Operation to average by source

        Example:
            # Replaces AverageByTable
            Panda(
                tables=[cycle0.tables.merged, cycle1.tables.merged, cycle2.tables.merged],
                operations=[
                    Panda.concat(),
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
                 prune_redundant_provenance: bool = True,
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
                  When using multiple pools with concat in the operation chain,
                  files are copied from the pool corresponding to each row's
                  origin (tracked implicitly via Panda.SOURCE).
            rename: If provided, output IDs will be renamed to {rename}_1, {rename}_2, etc.
                    Useful after sorting to get ranked IDs (e.g., rename="best" -> best_1, best_2, ...)
            ignore_missing: If True (default), skip missing pool files with a warning instead
                    of failing. Set to False to raise an error on any missing file.
            prune_redundant_provenance: If True (default), drop plain ``<axis>.id``
                    provenance columns from emitted stream map_tables when every
                    non-empty value is recoverable from the row id by id-semantics
                    alone (the same matcher downstream joins use). Lineage the
                    matcher cannot reconstruct (Panda-rename ``.-N.id`` columns,
                    ``pool.id`` / ``original.id``) is always kept. Set False to
                    retain every provenance column verbatim.
            **kwargs: Additional parameters

        Output:
            Streams: inherits all streams from pool input (if pool mode)
            Tables:
                result: columns derived from input + applied operations (dynamic)
                missing: id | removed_by | kind | cause
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
        self.prune_redundant_provenance = prune_redundant_provenance

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
            "calculate", "zscore", "fillna", "merge", "concat", "groupby", "pivot", "melt",
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
                # Derive rename prefix from step folder name, keeping the
                # execution-order int as a leading namespace so two unnamed Panda
                # steps don't collide on the same Panda_<N> ids (which suffix/
                # exact matching would cross-link before provenance). Strip the
                # zero-padding: "010_Panda_Cycle1" -> "10_Panda_Cycle1".
                folder_name = os.path.basename(self.output_folder)
                parts = folder_name.split("_", 1)
                if len(parts) > 1 and parts[0].isdigit():
                    self.rename = f"{int(parts[0])}_{parts[1]}"  # e.g., "10_Panda_Cycle1"
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
                        if stream and len(stream) > 0 and (stream.is_shared_file or len(stream.files) > 0):
                            json_path = self.configuration_path(
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

            # Collect upstream missing table paths for propagation
            self.upstream_missing_paths = []
            for pool in self.pool_outputs:
                missing_path = self._get_upstream_missing_table_path(pool)
                if missing_path:
                    self.upstream_missing_paths.append(missing_path)

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
            # Note: concat injects Panda.SOURCE into the working frame for
            # the rest of the chain, but pipe_panda strips it from the
            # final result, so it's not part of the output column schema.

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
        # configuration/, tables/, stream sub-dirs all auto-created by the pipeline.

        # Create config file
        step_tool_name = os.path.basename(self.output_folder)

        # In pool mode, compute per-stream map_table targets so pipe_panda
        # can emit `<stream>_map.csv` with `<stream>.id` (upstream ID)
        # provenance columns. Also collect per-stream destination folders
        # so extracted files land in the correct stream folder rather than
        # flat under output_folder.
        stream_map_targets = []
        stream_folders_map = {}
        pool_table_targets = {}
        if self.use_pool_mode:
            try:
                out = self.get_output_files()
                from .datastream import DataStream as _DS
                for name, spec in out.items():
                    if isinstance(spec, _DS):
                        # Always record the stream folder so extracted files
                        # land in their proper home, even for streams that
                        # don't carry a dedicated map_table.
                        stream_folders_map[name] = self.stream_folder(name)
                        if spec.map_table:
                            if spec.is_shared_file:
                                # Shared artifact: pass the path as-is, no
                                # <id> template substitution at runtime.
                                file_template = spec.files
                            elif isinstance(spec.files, list) and spec.files:
                                file_template = spec.files[0]
                            else:
                                file_template = ""
                            stream_map_targets.append({
                                "stream_name": name,
                                "map_table": spec.map_table,
                                "file_template": file_template,
                            })
                # Propagated pool tables now live under our tables/; map
                # source basenames to their new canonical destinations so
                # pipe_panda copies them into the right place.
                # TableInfo's ``__getattr__`` masks ``.path`` with a
                # TableReference for column access, so go through ``.info.path``.
                from .base_config import TableInfo as _TI
                for table_name, info in out.get("tables", {}).items():
                    if isinstance(info, _TI):
                        pool_table_targets[table_name] = info.info.path
            except Exception:
                stream_map_targets = []
                stream_folders_map = {}
                pool_table_targets = {}

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
            "step_tool_name": step_tool_name,
            "upstream_missing_paths": getattr(self, 'upstream_missing_paths', []),
            "stream_map_targets": stream_map_targets,
            "stream_folders": stream_folders_map,
            "pool_table_targets": pool_table_targets,
            "missing_csv": self.missing_csv,
            "prune_redundant_provenance": self.prune_redundant_provenance,
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

        Returns the minimum n across all head/tail/sample operations,
        since each one can only reduce the count further.

        Returns:
            Predicted count if determinable (e.g., head/tail/sample with n), None otherwise.
        """
        counts = []
        for op in self.operations:
            if op.type in ("head", "tail", "sample"):
                n = op.params.get("n")
                if n is not None:
                    counts.append(n)
        return min(counts) if counts else None

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files after transformation."""
        expected_columns = self._get_expected_columns()

        tables = {
            "result": TableInfo(
                name="result",
                path=self.output_csv,
                columns=expected_columns,
                description=f"Transformed table with operations: {', '.join(op.type for op in self.operations)}"
            ),
            "missing": TableInfo(
                name="missing",
                path=self.missing_csv,
                columns=["id", "removed_by", "kind", "cause"],
                description="IDs removed by this tool with removal reason"
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
                                # Track each pool's format separately so we can
                                # merge them into a "a|b" union format below
                                # when pools disagree (e.g. one produces .pdb,
                                # another .cif). Mirrors pdb.py/rcsb.py.
                                "formats": [],
                                "map_table": stream.map_table or "",
                                # Shared-file streams (e.g. multi-record FASTA):
                                # remember the source path; the output stream's
                                # files stays a single str pointing at the
                                # sliced copy under our stream folder.
                                "shared_src": None,
                            }
                        stream_data[stream_name]["ids"].extend(list(stream.ids))
                        for f in (stream.format or "").split("|"):
                            f = f.strip()
                            if f and f not in stream_data[stream_name]["formats"]:
                                stream_data[stream_name]["formats"].append(f)
                        # Only collect files for file-based streams. A
                        # value-based stream has files=[] (content lives in its
                        # map_table), so emptiness is the signal — not the
                        # format string (a per-id MSA can be format="csv" yet
                        # file-based).
                        if stream.files:
                            if stream.is_shared_file:
                                # First pool's basename wins for the predicted
                                # dest. At runtime pipe_panda slices each pool
                                # then merges the parts into this single dest
                                # via stream_slicers.get_merger(fmt).
                                if stream_data[stream_name]["shared_src"] is None:
                                    stream_data[stream_name]["shared_src"] = stream.files
                            else:
                                stream_data[stream_name]["files"].extend(stream.files)

            # Collapse the per-pool formats into a single "a|b" string —
            # downstream code reads sdata["format"] as before.
            for sdata in stream_data.values():
                fmts = sdata.pop("formats", [])
                sdata["format"] = "|".join(fmts) if fmts else ""

            # Deduplicate IDs across pools (first occurrence wins)
            for sdata in stream_data.values():
                seen = set()
                deduped_ids = []
                deduped_files = []
                # Shared-file streams: there's no per-id files list to dedup,
                # only ids — the slicer will produce one shared output file.
                if sdata["shared_src"] is not None:
                    for sid in sdata["ids"]:
                        if sid not in seen:
                            seen.add(sid)
                            deduped_ids.append(sid)
                    sdata["ids"] = deduped_ids
                    # Leave sdata["files"] as [] — shared output is built below
                    # from shared_src and emitted as a str.
                    continue
                # A stream is template-shaped if every collected files entry
                # carries an <id> placeholder. With one pool that's a single
                # element; with N pools contributing different extensions
                # (e.g. <id>.pdb + <id>.cif) it's N elements that all
                # template-expand. Either way: dedup IDs only, collapse the
                # files list to a single canonical template (with .* when
                # extensions disagree, mirroring pdb.py:716-718).
                files_list = sdata["files"]
                is_template = bool(files_list) and all('<id>' in f for f in files_list)
                if is_template:
                    exts = {os.path.splitext(f)[1] for f in files_list}
                    if len(exts) == 1:
                        canonical_template = files_list[0]
                    else:
                        # Mixed extensions across pools — use wildcard.
                        # All templates share the same dirname/<id> stem,
                        # so take dirname(files_list[0]) for the prefix.
                        template_dir = os.path.dirname(files_list[0])
                        canonical_template = os.path.join(template_dir, "<id>.*")
                    sdata["files"] = [canonical_template]
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

            # Build output streams — per-stream files + map_tables land in
            # their own stream folder; standalone tables (carried forward
            # from the first pool) go under tables/.
            output_streams = {}
            for stream_name, data in stream_data.items():
                stream_dir = self.stream_folder(stream_name)
                # The content-bearing map_table convention: if the pool's
                # map_table is named like a content file (e.g. sequences.csv),
                # keep its filename; otherwise emit <stream>_map.csv.
                pool_map_table = data.get("map_table", "")
                if pool_map_table:
                    map_basename = os.path.basename(pool_map_table)
                    map_table = os.path.join(stream_dir, map_basename)
                else:
                    map_table = None

                # Shared-file streams: emit a single str path pointing at
                # the sliced artifact under our stream folder. The slicer
                # runs at execution time in pipe_panda — here we just
                # predict the output path.
                if data.get("shared_src"):
                    if self.rename and predicted_count is not None:
                        new_ids = [f"{self.rename}_{i+1}" for i in range(predicted_count)]
                    elif self.rename:
                        new_ids = [f"{self.rename}_[<N>]"]
                    else:
                        new_ids = data["ids"]
                    shared_basename = os.path.basename(data["shared_src"])
                    shared_dest = os.path.join(stream_dir, shared_basename)
                    output_streams[stream_name] = DataStream(
                        name=stream_name,
                        ids=new_ids,
                        files=shared_dest,   # str — shared form preserved
                        map_table=map_table or "",
                        format=data["format"]
                    )
                    continue

                # Value-based streams (no files): propagate as-is with original map_table
                if not data["files"]:
                    if self.rename and predicted_count is not None:
                        new_ids = [f"{self.rename}_{i+1}" for i in range(predicted_count)]
                    elif self.rename:
                        new_ids = [f"{self.rename}_[<N>]"]
                    else:
                        new_ids = data["ids"]
                    output_streams[stream_name] = DataStream(
                        name=stream_name,
                        ids=new_ids,
                        files=[],
                        map_table=map_table or "",
                        format=data["format"]
                    )
                    continue

                ext = os.path.splitext(data["files"][0])[1]

                if self.rename and predicted_count is not None:
                    # Rename with known count: generate concrete IDs and file paths
                    new_ids = [f"{self.rename}_{i+1}" for i in range(predicted_count)]
                    new_files = [os.path.join(stream_dir, f"{nid}{ext}") for nid in new_ids]
                elif self.rename:
                    # Rename with unknown count: lazy pattern
                    new_ids = [f"{self.rename}_[<N>]"]
                    new_files = [os.path.join(stream_dir, f"<id>{ext}")]
                else:
                    # No rename: copy pool IDs as-is, use template for files
                    new_ids = data["ids"]
                    new_files = [os.path.join(stream_dir, f"<id>{ext}")]

                output_streams[stream_name] = DataStream(
                    name=stream_name,
                    ids=new_ids,
                    files=new_files,
                    map_table=map_table,
                    format=data["format"]
                )


            # Include tables from first pool (they typically have same schema).
            # When a table's name matches a declared stream whose map_table
            # IS the content table (Sequence.sequences, Ligand.compounds, …),
            # route the propagated copy to the stream folder so the pool
            # copy lands at the same path as the stream's map_table — one
            # file, one source of truth. Otherwise, standalone TableInfos
            # live under tables/.
            first_pool = self.pool_outputs[0]
            if hasattr(first_pool, 'tables') and hasattr(first_pool.tables, '_tables'):
                for name, info in first_pool.tables._tables.items():
                    if name not in tables:  # Don't overwrite our result table
                        # Is this table the content-bearing map of a stream?
                        stream = output_streams.get(name)
                        if stream is not None and getattr(stream, "map_table", None):
                            dest_path = stream.map_table
                        else:
                            dest_path = self.table_path(name)
                        tables[name] = TableInfo(
                            name=name,
                            path=dest_path,
                            columns=info.info.columns,
                            description=info.info.description,
                        )

            # Propagate rendering_parameters from first pool that has them
            for pool in self.pool_outputs:
                rp = getattr(pool, 'rendering_parameters', None)
                if rp:
                    output_streams["rendering_parameters"] = rp
                    break

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
