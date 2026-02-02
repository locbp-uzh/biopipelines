# Data Management

[← Back to Tool Reference](../ToolReference.md)

---

### Panda

Unified tool for pandas-style table transformations. Supports filtering, sorting, ranking, merging, concatenation, grouping, and calculated columns through a declarative operation-based API.

**Replaces**: Filter, Rank, SelectBest, MergeTables, ConcatenateTables, SliceTable

**Environment**: `ProteinEnv`

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
