# Data Management

[← Back to Tool Reference](../ToolReference.md)

---

### Filter

Filters structures or sequences based on metric criteria. Uses pandas query expressions to select items meeting specified thresholds.

**Environment**: `ProteinEnv`

**Parameters**:
- `data`: Union[ToolOutput, StandardizedOutput] (required) - Table input to filter
- `pool`: Union[ToolOutput, StandardizedOutput] = None - Structure/sequence pool for copying filtered items
- `expression`: str (required) - Pandas query-style filter expression (e.g., "distance < 3.5 and confidence > 0.8")
- `max_items`: Optional[int] = None - Maximum items to keep after filtering
- `sort_by`: Optional[str] = None - Column name to sort by before applying max_items
- `sort_ascending`: bool = True - Sort order (True = ascending, False = descending)

**Outputs**:
- Filtered pool with same structure as input
- All the tables of the upstream tool given as input will be copied and fitered based on the expression, and maintain the same table name. For example, after filtering an output of MergeTables, you can access filtered.tables.merged.
- `tables.missing`:

  | id | structure | msa |
  |----|-----------|-----|

**Example**:
```python
from PipelineScripts.filter import Filter

filtered = Filter(
    data=distances.tables.analysis,
    pool=boltz,
    expression="distance < 3.5 and confidence_score > 0.85",
    max_items=10,
    sort_by="distance"
)
```

---

### Rank

Ranks entries based on a metric (column or computed expression), renames IDs to sequential format, and optionally copies structures/compounds in ranked order. Useful for generating ranked lists with standardized naming.

**Environment**: `ProteinEnv`

**Parameters**:
- `data`: Union[ToolOutput, StandardizedOutput, TableInfo, tuple] (required) - Table input to rank. Supports tuple notation: `tool.tables.table_name.column_name` (metric auto-extracted)
- `pool`: Union[ToolOutput, StandardizedOutput] = None - Structure/sequence pool for copying ranked items
- `metric`: str = None - Column name or expression for ranking (e.g., "pLDDT" or "0.8*pLDDT+0.2*affinity"). Optional if using tuple notation for data.
- `ascending`: bool = False - Sort order (False = descending/best first, True = ascending)
- `prefix`: str = "rank" - Prefix for renamed IDs (e.g., "rank" produces rank_1, rank_2, ...)
- `top`: Optional[int] = None - Limit to top N entries after ranking

**Outputs**:
- Ranked pool with same structure as input (when pool provided)
- `tables.ranked`: Full ranked table with all columns

  | id | source_id | metric | {variable_columns} | {original_columns} |
  |----|-----------|--------|-------------------|-------------------|

- `tables.metrics`: Summary table with only ranking info

  | id | source_id | {metric_column} |
  |----|-----------|-----------------|

**Output Columns**:
- `id`: Renamed IDs (e.g., rank_1, rank_2, ...)
- `source_id`: Original IDs
- `metric`: Computed metric column (if expression used)
- Individual variable columns if metric is an expression with multiple variables (e.g., pLDDT, affinity)
- All other original columns preserved

**Metric Types**:
- **Column reference**: Simple column name (e.g., "pLDDT", "affinity")
- **Expression**: Computed metric using pandas eval syntax (e.g., "0.8*pLDDT + 0.2*affinity", "pLDDT - 2*rmsd")

**Example**:
```python
from PipelineScripts.rank import Rank

# Rank using tuple notation (metric auto-extracted from column reference)
ranked = Rank(
    data=boltz.tables.structures.pLDDT,  # metric="pLDDT" inferred
    pool=boltz,
    prefix="model",
    top=10
)

# Rank by single column (explicit metric)
ranked = Rank(
    data=analysis.tables.merged,
    pool=boltz,
    metric="pLDDT",
    ascending=False,  # Higher is better
    prefix="model",
    top=10
)

# Rank by computed expression with pool mode
ranked = Rank(
    data=merged.tables.merged,
    pool=boltz,
    metric="0.8*pLDDT + 0.2*binding_affinity",
    prefix="design",
    top=20
)
```

---

### SelectBest

Selects the single best structure or sequence based on optimization criteria. Supports single or multi-objective optimization with configurable weights.

**Environment**: `ProteinEnv`

**Parameters**:
- `pool`: Union[ToolOutput, StandardizedOutput, List[Union[ToolOutput, StandardizedOutput]]] (required) - Single or list of tool outputs to select from
- `tables`: Union[List[Union[ToolOutput, StandardizedOutput, TableInfo, str]], List[str]] (required) - Tables to evaluate for selection
- `metric`: str (required) - Primary metric to optimize
- `mode`: str = "max" - Optimization direction ("max" or "min")
- `weights`: Optional[Dict[str, float]] = None - Dictionary of {metric_name: weight} for multi-metric selection
- `tie_breaker`: str = "first" - How to break ties ("first", "random", or metric name)
- `composite_function`: str = "weighted_sum" - How to combine metrics (weighted_sum, product, min, max)
- `name`: str = "best" - Name for output structure file

**Outputs**:
- Single best structure/sequence with same format as input pool

**Example**:
```python
from PipelineScripts.select_best import SelectBest

best = SelectBest(
    pool=boltz,
    tables=[distances.tables.analysis],
    metric="distance",
    mode="min"
)

# Multi-objective selection
best_multi = SelectBest(
    pool=boltz,
    tables=[analysis.tables.merged],
    metric="composite_score",
    weights={"binding_affinity": 0.6, "pLDDT": 0.4},
    mode="max"
)
```

---

### RemoveDuplicates

Removes duplicate structures or sequences from a pool. Supports deduplication by sequence, structure similarity, or ID matching.

**Environment**: `ProteinEnv`

**Parameters**:
- `pool`: Union[ToolOutput, StandardizedOutput] (required) - Items to deduplicate
- `history`: Optional[Union[ToolOutput, StandardizedOutput, List]] = None - Previous tables for cross-cycle deduplication
- `compare`: str = "sequence" - Comparison method (sequence, structure, id)
- `similarity_threshold`: float = 1.0 - Similarity threshold for structure comparison (1.0 = exact match)

**Outputs**:
- Deduplicated pool with same structure as input
- `tables.removed`:

  | id | reason |
  |----|--------|

**Example**:
```python
from PipelineScripts.remove_duplicates import RemoveDuplicates

unique = RemoveDuplicates(
    pool=lmpnn,
    compare="sequence"
)
```

---

### MergeTables

Combines multiple tables by joining on a common key column. Enables integration of metrics from different analysis tools.

**Environment**: `ProteinEnv`

**Parameters**:
- `tables`: List[Union[ToolOutput, StandardizedOutput, TableInfo, str]] (required) - List of tables to merge
- `key`: str = "id" - Join column name
- `prefixes`: Optional[List[str]] = None - Prefixes for columns from each table
- `suffixes`: Optional[List[str]] = None - Suffixes for columns from each table
- `how`: str = "inner" - Join type (inner, outer, left, right)
- `calculate`: Optional[Dict[str, str]] = None - Derived column expressions {new_col: expression}

**Outputs**:
- `tables.merged`: Combined table with columns from all inputs

**Example**:
```python
from PipelineScripts.merge_tables import MergeTables

merged = MergeTables(
    tables=[distances.tables.analysis, plip.tables.interactions],
    prefixes=["dist_", "plip_"],
    key="id",
    calculate={"score": "dist_distance + plip_energy"}
)
```

---

### ConcatenateTables

Stacks multiple tables vertically (row-wise). Useful for combining results from multiple cycles or parallel runs.

**Environment**: `ProteinEnv`

**Parameters**:
- `tables`: List[Union[ToolOutput, StandardizedOutput, TableInfo, str]] (required) - List of tables to concatenate
- `fill`: str = "N/A" - Value for missing columns
- `ignore_index`: bool = True - Reset index in concatenated output

**Outputs**:
- `tables.concatenated`: Row-wise concatenation of all input tables

**Example**:
```python
from PipelineScripts.concatenate_tables import ConcatenateTables

concat = ConcatenateTables(
    tables=[cycle1_results, cycle2_results, cycle3_results],
    fill="N/A"
)
```

---

### SliceTable

Extracts a subset of rows and/or columns from a table. Enables data sampling and column selection.

**Environment**: `ProteinEnv`

**Parameters**:
- `table`: Union[ToolOutput, StandardizedOutput, TableInfo, str] (required) - Input table to slice
- `start`: int = 0 - Starting row index
- `end`: Optional[int] = None - Ending row index (None = to end)
- `step`: int = 1 - Step size for slicing
- `columns`: Optional[List[str]] = None - Specific columns to keep (None = all columns)

**Outputs**:
- `tables.sliced`: Sliced table

**Example**:
```python
from PipelineScripts.slice_table import SliceTable

sliced = SliceTable(
    table=results.tables.analysis,
    start=0,
    end=100,
    columns=["id", "distance", "confidence"]
)
```

---

### ExtractMetrics

Extracts and aggregates specific metrics from tables. Supports grouping and various aggregation functions for data summarization.

**Environment**: `ProteinEnv`

**Parameters**:
- `tables`: List[Union[ToolOutput, StandardizedOutput, TableInfo, str]] (required) - Input tables
- `metrics`: List[str] (required) - Metric column names to extract
- `group_by`: Optional[str] = None - Column to group by for aggregation
- `aggregation`: str = "mean" - Aggregation function (mean, median, min, max, sum, std)
- `pivot`: bool = False - Pivot metrics to columns

**Outputs**:
- `tables.extracted`: Extracted metrics table

**Example**:
```python
from PipelineScripts.extract_metrics import ExtractMetrics

metrics = ExtractMetrics(
    tables=[boltz.tables.confidence],
    metrics=["complex_plddt", "ptm"],
    group_by="input_file",
    aggregation="mean"
)
```

---

### AverageByTable

Computes averages of metrics grouped by a specified column. Useful for summarizing results across multiple structures or cycles.

**Environment**: `ProteinEnv`

**Parameters**:
- `tables`: List[Union[ToolOutput, StandardizedOutput, TableInfo, str]] (required) - Input tables
- `group_by`: str (required) - Column to group by
- `metrics`: List[str] (required) - Metric columns to average
- `weights`: Optional[Dict[str, float]] = None - Weights for each metric

**Outputs**:
- `tables.averaged`: Averaged metrics by group

**Example**:
```python
from PipelineScripts.average_by_table import AverageByTable

averaged = AverageByTable(
    tables=[cycle1.tables.analysis, cycle2.tables.analysis],
    group_by="structure_id",
    metrics=["distance", "confidence"]
)
```

---
