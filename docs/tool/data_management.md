# Data Management

[← Back to Tool Reference](../tool_reference.md)

These tools reshape and route the **tables and streams** flowing through a pipeline — filtering and merging metric tables, gathering parallel runs, renaming IDs, and reformatting for external software. To load an *existing* CSV/Excel file into a pipeline, see [Table](inputs_io.md#table) under Inputs & I/O.

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

## Panda

Unified pandas-style table transformations. Replaces Filter, Rank, SelectBest, MergeTables, ConcatenateTables, SliceTable.

**Environment**: `biopipelines`

**Parameters**:
- `tables`: TableInfo | StandardizedOutput | str | List[...] = None - One table or a list of tables (a list enables per-frame ops and merge/concat)
- `operations`: List[Operation] = None - Sequence of operations
- `pool`: StandardizedOutput | List[StandardizedOutput] = None - Copy files matching filtered IDs (list = one pool per input table)
- `rename`: str = None - Rename output IDs to `{rename}_1`, `{rename}_2`, ...
- `ignore_missing`: bool = True - Tolerate missing columns/tables instead of raising
- `prune_redundant_provenance`: bool = True - Drop redundant `<axis>.id` provenance columns from the result

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
| `merge(on, how, prefixes)` | `Panda.merge(prefixes=["a_", "b_"])` |
| `concat(fill)` | `Panda.concat(fill="")` (auto-tags rows with `Panda.SOURCE` when >1 inputs) |
| `calculate(exprs)` | `Panda.calculate({"delta": "a - b", "k2": "cos(angle) ** 2"})` |
| `zscore(columns, by, sign)` | `Panda.zscore(["plddt","aggr"], sign={"aggr":-1})` — standardize to `<col>_z` for scale-fair combining (optional per-group `by=`, sign flip for lower-is-better) |
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

# Merge tables (on=None uses biopipelines ID matching by default)
merged = Panda(
    tables=[apo.tables.affinity, holo.tables.affinity],
    operations=[
        Panda.merge(prefixes=["apo_", "holo_"]),
        Panda.calculate({"delta": "holo_affinity - apo_affinity"})
    ]
)

# Tip: when both tables share a literal `id` column with no auto-rename
# in between, prefer Panda.merge(on="id") — it skips the suffix / provenance
# matcher and uses pandas.merge directly. ~20× faster on large tables and
# the result is identical when the IDs already match exactly.
fast_merged = Panda(
    tables=[boltz.tables.confidence, boltz.tables.affinity],
    operations=[Panda.merge(on="id")]
)

# Calculate with math functions (cos, sin, sqrt, log, exp, radians, degrees, pi, ...)
# Expressions can reference columns defined earlier in the same calculate call
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

# Concatenate tables. With >1 inputs, every row is implicitly tagged with
# its origin index in the internal `Panda.SOURCE` column for the rest of
# the chain (stripped from the final CSV). Reference it as a column name
# in groupby / filter / etc.
combined = Panda(
    tables=[cycle0.tables.results, cycle1.tables.results],
    operations=[Panda.concat(fill="")]
)

# Multi-pool selection (select best from multiple sources). The implicit
# Panda.SOURCE column drives multi-pool file routing too: each surviving
# row's structure is copied from its origin pool.
best = Panda(
    tables=[cycle1.tables.result, cycle2.tables.result],
    operations=[
        Panda.concat(),
        Panda.sort("metric", ascending=True),
        Panda.head(1)
    ],
    pool=[cycle1, cycle2],  # Pools match tables
    rename="best"
)

# Best per source — natural form: per-frame ops broadcast over a list
# of input tables (sort/head/tail/sample/filter/etc. apply independently
# to each table), then concat stacks the results.
best_per_source = Panda(
    tables=[a.tables.result, b.tables.result, c.tables.result],
    operations=[
        Panda.sort("score", ascending=False),  # applied per table
        Panda.head(1),                         # best row per table
        Panda.concat(),                        # stack into one frame
    ],
)

# Renaming the survivors. `rename=` runs once at the very end of the chain
# on whatever rows are left, so just pass it on the Panda(...) call — no
# extra operation step needed. Concat preserves input order, so the
# default below gives best_1 = a's best, best_2 = b's best, best_3 = c's.
ranked_per_source = Panda(
    tables=[a.tables.result, b.tables.result, c.tables.result],
    operations=[
        Panda.sort("score", ascending=False),
        Panda.head(1),
        Panda.concat(),
    ],
    rename="best",  # → best_1, best_2, best_3 (one per source)
)

# Want best_1 to be the GLOBALLY top-scoring row instead of the first
# source's best? Add a cross-source sort after concat so rename ranks
# across sources rather than mirroring input order.
ranked_globally = Panda(
    tables=[a.tables.result, b.tables.result, c.tables.result],
    operations=[
        Panda.sort("score", ascending=False),  # per table
        Panda.head(1),                         # per table
        Panda.concat(),                        # stack
        Panda.sort("score", ascending=False),  # cross-source rank
    ],
    rename="best",
)

# Same outcome, post-concat form: useful when you also want a cross-source
# aggregate (mean/std/etc.) — groupby/pivot/melt/average_by_source must
# come AFTER concat because they aggregate across sources.
mean_per_source = Panda(
    tables=[a.tables.result, b.tables.result, c.tables.result],
    operations=[
        Panda.concat(),
        Panda.groupby(Panda.SOURCE, {"score": "mean"}),
    ],
)
```

---

## Pool

Gathers N `StandardizedOutput`s from parallel runs of the **same upstream tool** into one combined `StandardizedOutput`. Designed to pair with `with Parallel():` (see *Parallel batches* under [Resources](../user_manual.md#resources) in the user manual) for the canonical fan-out / fan-in pattern.

**Environment**: `biopipelines`

**Parameters**:
- `runs`: List of two or more `StandardizedOutput` objects, all from runs of the same upstream tool. Must expose identical stream-name sets and identical table-name sets, with matching `format` per shared stream. Pool raises with a descriptive message on mismatch.
- `recount_prefix` (optional): if set, replace the default per-run id suffix (`<orig_id>_<pool_idx>`) with a flat 1-based renumber across all rows of the pool, producing ids `<recount_prefix>_1`, `<recount_prefix>_2`, …, `<recount_prefix>_N`. Original ids are preserved in an `original.id` column and `pool.path` still records the source-run index. Counted at config time when every input has fully-resolved ids; otherwise the framework emits a lazy `<recount_prefix>_[<N>]` pattern that the runtime resolves.

**Type-agnostic semantics**: Pool iterates `runs[0].streams.items()` and treats every `DataStream` uniformly — there is no special-casing for structures vs sequences vs compounds. The same applies to tables. As long as all inputs share the same stream / table names and formats, Pool works.

**ID renumbering**: by default every output id is `<orig_id>_<pool_idx>` where `pool_idx` is the 1-based source-run position in `Pool(runs=runs)`. The typical case — same upstream tool, same parameters, identical original ids across runs — therefore no longer collides. Use `recount_prefix=` for a flat renumber instead (see above).

**Provenance**: every emitted map_table carries an extra column `pool.path` whose value is the source-run index (`1`, `2`, …, `N`) for each row. This matches the spirit of the existing `<axis>.id` provenance convention but tracks the *parallel-run* axis specifically.

**Streams / Tables**: every shared stream and table on `runs[0]` appears on the pooled output with concatenated rows.

**Examples**:

```python
from biopipelines.entities import PDB
from biopipelines.pipeline import Parallel, Resources
from biopipelines.rfdiffusion import RFdiffusion
from biopipelines.pool import Pool

# Run RFdiffusion 10x in parallel, gather all 100 designs into one stream
seed = PDB("4AKE")
runs = []
with Parallel():
    for _ in range(10):
        Resources(gpu="A100", time="6:00:00")
        runs.append(RFdiffusion(pdb=seed, num_designs=10))

Resources(gpu="A100", time="12:00:00")
combined = Pool(runs=runs)
# combined.streams.structures has 100 ids of the form
# design_1_1, design_2_1, ..., design_10_1, design_1_2, ..., design_10_10
# combined map_table also carries a pool.path column.
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

## Selection

Combines and modifies PyMOL-formatted selection strings using composable operations applied left-to-right.

**Environment**: `biopipelines`

**Parameters**:
- `*ops`: Sequence of `SelectionOp` objects (from `Selection.add`, `Selection.subtract`, `Selection.expand`, `Selection.shrink`, `Selection.shift`, `Selection.invert`)
- `structures`: StandardizedOutput = None - Required for structure-aware ops (expand, shrink, shift, invert)

**Operations**:
- `Selection.add(*refs)` - Union of one or more column references
- `Selection.subtract(*refs)` - Remove residues from running selection
- `Selection.expand(n)` - Add n residues on each side
- `Selection.shrink(n)` - Remove n residues from each side
- `Selection.shift(n)` - Shift all intervals by n
- `Selection.invert()` - Select complement

**Tables**:
- `selections`: `id | selection | n_residues` — `selection` is the chain-aware PyMOL selection string (e.g. `"A12+A45-47"`); `n_residues` is the number of residues it contains, so a Selection result can be used directly as a size metric (e.g. to measure how many pocket residues survive a set-difference).

**Example**:

```python
from biopipelines.selection import Selection
from biopipelines.distance_selector import DistanceSelector

distances = DistanceSelector(structures=rfdaa, ligand=rfdaa, distance=5)

# Expand by 2 residues
expanded = Selection(
    Selection.add(distances.tables.selections.within),
    Selection.expand(2),
    structures=rfdaa,
)

# Union two columns then invert
fixed = Selection(
    Selection.add(fuse.tables.sequences.L1, fuse.tables.sequences.L2),
    Selection.invert(),
    structures=rfdaa,
)
```
