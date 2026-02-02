# Data Management

[← Back to Tool Reference](../ToolReference.md)

---

### Panda

Unified tool for pandas-style table transformations. Supports filtering, sorting, ranking, merging, concatenation, grouping, and calculated columns through a declarative operation-based API.

**Replaces**: Filter, Rank, SelectBest, MergeTables, ConcatenateTables, SliceTable

**Does NOT replace**:
- **RemoveDuplicates**: Use for cross-table deduplication (filter current pool against history)
- **ExtractMetrics**: Use for creating separate CSV files per metric for statistical software

**Environment**: `biopipelines`

**Parameters**:
- `table`: Union[TableInfo, StandardizedOutput, str] - Single table input (mutually exclusive with `tables`)
- `tables`: List[Union[TableInfo, StandardizedOutput, str]] - Multiple tables for merge/concat operations
- `operations`: List[Operation] (required) - Sequence of operations to apply
- `pool`: Optional[StandardizedOutput] = None - Tool output for pool mode (copies structures matching filtered IDs)
- `rename`: Optional[str] = None - If provided, output IDs will be renamed to `{rename}_1`, `{rename}_2`, etc.

**Available Operations**:

| Operation | Description | Example |
|-----------|-------------|---------|
| `Panda.filter(expr)` | Filter rows using pandas query | `Panda.filter("pLDDT > 80 and distance < 5.0")` |
| `Panda.sort(by, ascending)` | Sort by column(s) | `Panda.sort("affinity", ascending=False)` |
| `Panda.head(n)` | Keep first N rows | `Panda.head(10)` |
| `Panda.tail(n)` | Keep last N rows | `Panda.tail(5)` |
| `Panda.sample(n, frac)` | Random sample | `Panda.sample(n=100, random_state=42)` |
| `Panda.rank(by, prefix, ascending)` | Add rank column | `Panda.rank(by="score", ascending=False)` |
| `Panda.drop_duplicates(subset, keep)` | Remove duplicates | `Panda.drop_duplicates(subset="sequence")` |
| `Panda.merge(on, how, prefixes)` | Join tables horizontally | `Panda.merge(on="id", prefixes=["apo_", "holo_"])` |
| `Panda.concat(fill, add_source)` | Stack tables vertically | `Panda.concat(fill="", add_source=True)` |
| `Panda.calculate(exprs)` | Add computed columns | `Panda.calculate({"delta": "holo - apo"})` |
| `Panda.groupby(by, agg)` | Group and aggregate | `Panda.groupby("category", {"score": "mean"})` |
| `Panda.select_columns(cols)` | Keep specified columns | `Panda.select_columns(["id", "score"])` |
| `Panda.drop_columns(cols)` | Remove columns | `Panda.drop_columns(["temp"])` |
| `Panda.rename(mapping)` | Rename columns | `Panda.rename({"old": "new"})` |
| `Panda.fillna(value, column)` | Fill missing values | `Panda.fillna(0)` |
| `Panda.pivot(index, columns, values)` | Wide format | `Panda.pivot("id", "metric", "value")` |
| `Panda.melt(id_vars, value_vars)` | Long format | `Panda.melt(id_vars="id")` |
| `Panda.average_by_source()` | Average per source table | `Panda.average_by_source()` |

**Outputs**:
- `tables.result`: Transformed table with all operations applied
- `tables.missing`: IDs filtered out (when using pool mode)
- Pool mode: Structures/compounds/sequences matching filtered IDs

**Examples**:

```python
from PipelineScripts.panda import Panda

# Filter (replaces Filter tool)
filtered = Panda(
    table=boltz.tables.confidence,
    operations=[
        Panda.filter("confidence_score > 0.8")
    ]
)

# Sort + head (replaces SelectBest)
best = Panda(
    table=boltz.tables.confidence,
    operations=[
        Panda.sort("confidence_score", ascending=False),
        Panda.head(5)
    ]
)

# Rank with renamed IDs (replaces Rank tool)
ranked = Panda(
    table=boltz.tables.confidence,
    operations=[
        Panda.sort("confidence_score", ascending=False)
    ],
    rename="best",  # Output: best_1, best_2, ...
    pool=boltz
)

# Merge tables (replaces MergeTables)
merged = Panda(
    tables=[apo.tables.affinity, holo.tables.affinity],
    operations=[
        Panda.merge(on="id", prefixes=["apo_", "holo_"]),
        Panda.calculate({"delta": "holo_affinity - apo_affinity"})
    ]
)

# Concatenate tables (replaces ConcatenateTables)
combined = Panda(
    tables=[cycle0.tables.results, cycle1.tables.results],
    operations=[
        Panda.concat(fill="", add_source=True)
    ]
)

# Complex pipeline with pool mode
complex_result = Panda(
    tables=[boltz1.tables.confidence, boltz2.tables.confidence],
    operations=[
        Panda.concat(fill="", add_source=True),
        Panda.filter("confidence_score > 0.6"),
        Panda.sort("confidence_score", ascending=False),
        Panda.head(10)
    ],
    rename="top",
    pool=boltz1
)
```

---

### RemoveDuplicates

**Cross-table deduplication** tool for iterative design cycles. Filters sequences from a current pool that already exist in historical data, preventing recomputation of identical sequences across cycles.

**When to use RemoveDuplicates vs Panda.drop_duplicates()**:
- Use **RemoveDuplicates** for cross-table deduplication (filter current pool against history from previous cycles)
- Use **Panda.drop_duplicates()** for within-table deduplication (remove duplicates from a single table)

**Environment**: `biopipelines`

**Parameters**:
- `pool`: Union[StandardizedOutput, TableInfo, str] (required) - Current cycle sequences to check for duplicates
- `history`: Optional[Union[StandardizedOutput, TableInfo, str]] = None - Historical sequences from previous cycles (None for first cycle)
- `compare`: str = "sequence" - Column name to compare for duplicates

**Outputs**:
- `tables.sequences`: Filtered unique sequences not present in history
- `tables.missing`: Sequences filtered out due to duplication

  | id | structure | msa |
  |----|-----------|-----|

**Example**:
```python
from PipelineScripts.remove_duplicates import RemoveDuplicates

# First cycle - remove self-duplicates only
unique_sequences = RemoveDuplicates(
    pool=composer,
    history=None,
    compare="sequence"
)

# Subsequent cycles - filter against accumulated history
unique_sequences = RemoveDuplicates(
    pool=composer_current,
    history=all_sequences_seen,
    compare="sequence"
)
```

---

### ExtractMetrics

**Specialized metrics extraction** tool for statistical analysis software. Extracts specific metric columns from multiple tables and creates **one CSV file per metric**, with each column representing a different table/condition. Output format is optimized for tools like GraphPad Prism.

**When to use ExtractMetrics vs Panda**:
- Use **ExtractMetrics** when you need separate CSV files for each metric across multiple conditions (for Prism, etc.)
- Use **Panda** for standard table transformations within the pipeline

**Environment**: `biopipelines`

**Parameters**:
- `tables`: List[Union[TableInfo, str]] (required) - Input tables (typically one per cycle/condition)
- `metrics`: List[str] (required) - Metric column names to extract
- `table_names`: Optional[List[str]] = None - Custom column names for output (defaults to Table_0, Table_1, ...)

**Outputs**:
- `tables.{metric_name}`: One CSV file per metric

  | Table_0 | Table_1 | Table_2 | ... |
  |---------|---------|---------|-----|
  | value   | value   | value   | ... |

**Example**:
```python
from PipelineScripts.extract_metrics import ExtractMetrics

# Extract multiple metrics across cycles for Prism analysis
metrics_extract = ExtractMetrics(
    tables=[cycle0.tables.merged,
            cycle1.tables.merged,
            cycle2.tables.merged],
    metrics=["affinity_delta", "affinity_delta_R", "affinity_delta_S"],
    table_names=["Cycle0", "Cycle1", "Cycle2"]
)

# Output files:
# - affinity_delta.csv (columns: Cycle0, Cycle1, Cycle2)
# - affinity_delta_R.csv (columns: Cycle0, Cycle1, Cycle2)
# - affinity_delta_S.csv (columns: Cycle0, Cycle1, Cycle2)
```

---

### SelectionEditor

Modifies PyMOL-formatted selection strings (e.g., "3-45+58-60") with structure-aware operations. Validates all operations against actual PDB residue numbering and automatically merges overlapping/adjacent ranges.

**Installation**: Requires an environment containing pandas (e.g. biopipelines).

**Parameters**:
- `selection`: tuple (required) - Table column reference (e.g., `tool.tables.structures.designed`)
- `structures`: Optional[Union[ToolOutput, List[str]]] = None - Input structures (auto-detected from selection source if not provided)
- `expand`: int = 0 - Number of residues to add on each side of intervals
- `shrink`: int = 0 - Number of residues to remove from each side of intervals
- `shift`: int = 0 - Number of residues to shift all intervals (+/-)
- `invert`: bool = False - Whether to invert the selection (select complement)

**Outputs**:
- `tables.selections`:

  | id | pdb | {column_name} | original_{column_name} |
  |----|-----|---------------|------------------------|

**Example**:
```python
from PipelineScripts.selection_editor import SelectionEditor
from PipelineScripts.distance_selector import DistanceSelector
from PipelineScripts.ligand_mpnn import LigandMPNN

# Expand binding site selection by 2 residues
distances = DistanceSelector(
    structures=rfdaa,
    ligand="LIG",
    distance=5
)

expanded = SelectionEditor(
    selection=distances.tables.selections.within,
    expand=2
)

# Use expanded selection with LigandMPNN
lmpnn = LigandMPNN(
    structures=rfdaa,
    ligand="LIG",
    redesigned=expanded.tables.selections.within
)

# Invert selection to get everything except binding site
fixed_region = SelectionEditor(
    selection=distances.tables.selections.within,
    invert=True
)

# Shrink selection
tighter = SelectionEditor(
    selection=distances.tables.selections.within,
    shrink=1
)
```

**Key Features**:
- **Structure-aware**: Validates against actual PDB residue numbers (handles gaps, non-sequential numbering)
- **Range merging**: Automatically merges overlapping/adjacent intervals (e.g., "1-4+6-10" with expand=1 becomes "1-11")
- **Auto-detection**: Structures parameter is optional; can auto-detect from selection source
