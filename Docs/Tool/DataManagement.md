# Data Management

[← Back to Tool Reference](../ToolReference.md)

---

## Panda

Unified pandas-style table transformations. Replaces Filter, Rank, SelectBest, MergeTables, ConcatenateTables, SliceTable.

**Environment**: `biopipelines`

**Parameters**:
- `table`: TableInfo | StandardizedOutput | str - Single table input
- `tables`: List[...] - Multiple tables (for merge/concat)
- `operations`: List[Operation] - Sequence of operations
- `pool`: StandardizedOutput - Copy files matching filtered IDs
- `rename`: str - Rename output IDs to `{rename}_1`, `{rename}_2`, ...

**Operations**:

| Operation | Example |
|-----------|---------|
| `filter(expr)` | `Panda.filter("pLDDT > 80")` |
| `sort(by, ascending)` | `Panda.sort("score", ascending=False)` |
| `head(n)` | `Panda.head(10)` |
| `tail(n)` | `Panda.tail(5)` |
| `sample(n, frac)` | `Panda.sample(n=100)` |
| `rank(by, prefix)` | `Panda.rank(by="score")` |
| `drop_duplicates(subset)` | `Panda.drop_duplicates(subset="sequence")` |
| `merge(on, prefixes)` | `Panda.merge(on="id", prefixes=["a_", "b_"])` |
| `concat(fill, add_source)` | `Panda.concat(fill="")` |
| `calculate(exprs)` | `Panda.calculate({"delta": "a - b"})` |
| `groupby(by, agg)` | `Panda.groupby("cat", {"score": "mean"})` |
| `select_columns(cols)` | `Panda.select_columns(["id", "score"])` |
| `drop_columns(cols)` | `Panda.drop_columns(["temp"])` |
| `rename(mapping)` | `Panda.rename({"old": "new"})` |
| `fillna(value)` | `Panda.fillna(0)` |
| `pivot(index, columns, values)` | `Panda.pivot("id", "metric", "value")` |
| `melt(id_vars)` | `Panda.melt(id_vars="id")` |
| `average_by_source()` | `Panda.average_by_source()` |

**Outputs**:
- `tables.result` - Transformed table
- `tables.missing` - Filtered out IDs (pool mode)
- Pool mode: structures/compounds/sequences matching IDs

**Examples**:

```python
from PipelineScripts.panda import Panda

# Filter
filtered = Panda(
    table=boltz.tables.confidence,
    operations=[Panda.filter("confidence_score > 0.8")]
)

# Sort + head (replaces SelectBest)
best = Panda(
    table=boltz.tables.confidence,
    operations=[
        Panda.sort("confidence_score", ascending=False),
        Panda.head(5)
    ]
)

# Rank with renamed IDs
ranked = Panda(
    table=boltz.tables.confidence,
    operations=[Panda.sort("score", ascending=False)],
    rename="best",  # Output: best_1, best_2, ...
    pool=boltz
)

# Merge tables
merged = Panda(
    tables=[apo.tables.affinity, holo.tables.affinity],
    operations=[
        Panda.merge(on="id", prefixes=["apo_", "holo_"]),
        Panda.calculate({"delta": "holo_affinity - apo_affinity"})
    ]
)

# Concatenate tables
combined = Panda(
    tables=[cycle0.tables.results, cycle1.tables.results],
    operations=[Panda.concat(fill="", add_source=True)]
)

# Multi-pool selection (select best from multiple sources)
best = Panda(
    tables=[cycle1.tables.result, cycle2.tables.result],
    operations=[
        Panda.concat(add_source=True),
        Panda.sort("metric", ascending=True),
        Panda.head(1)
    ],
    pool=[cycle1, cycle2],  # Pools match tables
    rename="best"
)
```

---

## ExtractMetrics

Creates separate CSV files per metric for statistical software (GraphPad Prism).

**Environment**: `biopipelines`

**Parameters**:
- `tables`: List[TableInfo | str] - Input tables (one per condition)
- `metrics`: List[str] - Column names to extract
- `table_names`: List[str] = None - Custom column names

**Outputs**:
- `tables.{metric}` - One CSV per metric with columns for each table

**Example**:

```python
from PipelineScripts.extract_metrics import ExtractMetrics

metrics = ExtractMetrics(
    tables=[cycle0.tables.merged, cycle1.tables.merged, cycle2.tables.merged],
    metrics=["affinity_delta", "pLDDT"],
    table_names=["Cycle0", "Cycle1", "Cycle2"]
)
# Output: affinity_delta.csv, pLDDT.csv
```

---

## SelectionEditor

Modifies PyMOL selection strings (e.g., "3-45+58-60") with structure-aware operations.

**Environment**: `biopipelines`

**Parameters**:
- `selection`: tuple - Table column reference (e.g., `tool.tables.structures.designed`)
- `structures`: StandardizedOutput = None - Auto-detected from selection source
- `expand`: int = 0 - Residues to add on each side
- `shrink`: int = 0 - Residues to remove from each side
- `shift`: int = 0 - Shift all intervals (+/-)
- `invert`: bool = False - Select complement

**Outputs**:
- `tables.selections`: | id | pdb | {column} | original_{column} |

**Example**:

```python
from PipelineScripts.selection_editor import SelectionEditor
from PipelineScripts.distance_selector import DistanceSelector

distances = DistanceSelector(structures=rfdaa, ligand="LIG", distance=5)

# Expand by 2 residues
expanded = SelectionEditor(
    selection=distances.tables.selections.within,
    expand=2
)

# Invert selection
fixed = SelectionEditor(
    selection=distances.tables.selections.within,
    invert=True
)
```
