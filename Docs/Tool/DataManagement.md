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
| `calculate(exprs)` | `Panda.calculate({"delta": "a - b", "k2": "cos(angle) ** 2"})` |
| `groupby(by, agg)` | `Panda.groupby("cat", {"score": "mean"})` |
| `select_columns(cols)` | `Panda.select_columns(["id", "score"])` |
| `drop_columns(cols)` | `Panda.drop_columns(["temp"])` |
| `rename(mapping)` | `Panda.rename({"old": "new"})` |
| `fillna(value)` | `Panda.fillna(0)` |
| `pivot(index, columns, values)` | `Panda.pivot("id", "metric", "value")` |
| `melt(id_vars)` | `Panda.melt(id_vars="id")` |
| `average_by_source()` | `Panda.average_by_source()` |

**Streams** (pool mode): `structures`, `sequences`, `compounds` (matching filtered IDs)

**Tables**:
- `result` - Transformed table
- `missing` - Filtered out IDs (pool mode)

**Examples**:

```python
from biopipelines.panda import Panda

# Filter
filtered = Panda(
    tables=boltz.tables.confidence,
    operations=[Panda.filter("confidence_score > 0.8")]
)

# Sort + head (replaces SelectBest)
best = Panda(
    tables=boltz.tables.confidence,
    operations=[
        Panda.sort("confidence_score", ascending=False),
        Panda.head(5)
    ]
)

# Rank with renamed IDs
ranked = Panda(
    tables=boltz.tables.confidence,
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

# Calculate with math functions (cos, sin, sqrt, log, exp, radians, degrees, pi, ...)
# Expressions can reference columns defined earlier in the same calculate call
fret = Panda(
    tables=[distances.tables.result, angles.tables.angles],
    operations=[
        Panda.merge(on="id"),
        Panda.calculate({
            "kappa2": "cos(orientation) ** 2",
            "R0_eff": "49.0 * (kappa2 / 0.6667) ** (1.0 / 6.0)",
            "efficiency": "1 / (1 + (distance / R0_eff) ** 6)"
        })
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

## ReMap

Renames IDs across all streams and tables from a source tool output. At execution time, files are symlinked and CSV tables are rewritten with new IDs.

**Environment**: `biopipelines`

**Parameters**:
- `source`: StandardizedOutput - Tool output whose IDs will be renamed
- `onto`: str | list | dict | list[tuple] | DataStream | StandardizedOutput - Target ID specification
- `map`: StandardizedOutput = None - Intermediate tool for provenance bridging

**`onto` specification**:
| Type | Behavior |
|------|----------|
| `str` | Auto-number: `"design"` → `design_1`, `design_2`, ... |
| `list[str]` | Explicit new IDs (matched to streams with same length) |
| `dict` | Selective: `{"old_id": "new_id"}` |
| `list[tuple]` | Same as dict: `[("old_id", "new_id")]` |
| `DataStream` | Align onto stream's IDs |
| `StandardizedOutput` | Align onto tool's IDs (all streams must have same IDs) |

**Streams**: All source streams whose IDs can be fully remapped (others are discarded with a message)

**Tables**: Remapped copies of all source tables

**Examples**:

```python
from biopipelines.remap import ReMap

# Auto-numbered
remapped = ReMap(source=tool_a, onto="design")

# Explicit list
remapped = ReMap(source=tool_a, onto=["kinase_apo", "kinase_holo"])

# Dict mapping
remapped = ReMap(source=tool_a, onto={"prot1": "complex_A", "prot2": "complex_B"})

# List of tuples
remapped = ReMap(source=tool_a, onto=[("prot1", "complex_A"), ("prot2", "complex_B")])

# Align onto another tool's IDs
remapped = ReMap(source=tool_a, onto=tool_b)

# Align onto a specific stream
remapped = ReMap(source=tool_a, onto=tool_b.streams.structures)

# Use intermediate tool as provenance bridge
remapped = ReMap(source=tool_a, onto=tool_c, map=tool_b)
```

---

## ExtractMetrics

Creates separate CSV files per metric for statistical software (GraphPad Prism).

**Environment**: `biopipelines`

**Parameters**:
- `tables`: List[TableInfo | str] - Input tables (one per condition)
- `metrics`: List[str] - Column names to extract
- `table_names`: List[str] = None - Custom column names

**Tables**:
- `{metric}` - One CSV per metric with columns for each table

**Example**:

```python
from biopipelines.extract_metrics import ExtractMetrics

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

**Tables**:
- `selections`: | id | pdb | {column} | original_{column} |

**Example**:

```python
from biopipelines.selection_editor import SelectionEditor
from biopipelines.distance_selector import DistanceSelector

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
