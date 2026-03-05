# Lazy IDs: Compact ID Representation for DataStreams

**Status**: Design document
**Date**: 2026-03-05
**Authors**: G. Quargnali

---

## 1. Motivation

Currently, DataStreams store **fully expanded** ID lists at config time. For example, RFdiffusion with 50 designs on input `5HG6` produces:

```python
ids = ["5HG6_0", "5HG6_1", "5HG6_2", ..., "5HG6_49"]  # 50 strings
files = ["/path/5HG6_0.pdb", "/path/5HG6_1.pdb", ..., "/path/5HG6_49.pdb"]  # 50 strings
```

This has two problems:

1. **Verbosity**: For large pipelines (e.g., 10 PDBs x 50 RFdiffusion designs x 3 ProteinMPNN sequences = 1,500 IDs), the fully expanded lists are bulky and make config output, map_tables, and JSON files unnecessarily large at config time.

2. **Unpredictable outputs**: Some tools produce outputs whose count or naming depends on runtime data. The current Mutagenesis tool requires a fixed `position` integer precisely because it must predict exact output IDs at config time. If Mutagenesis accepted a per-row selection (e.g., "mutate all linker residues"), the number of output variants per input would depend on linker length — which varies per Fuse variant and is only known at runtime.

The `files_contain_wildcards` flag on DataStream is a partial workaround for the file side of this problem, but IDs themselves remain fully expanded.

---

## 2. Proposed Design: Lazy ID Patterns

### 2.1 Core Idea

Replace fully-expanded ID lists with **compact pattern expressions** that describe the ID space declaratively. These patterns are expanded lazily — fully at config time when possible, or partially/at-runtime when the expansion depends on runtime data.

### 2.2 Pattern Syntax

| Pattern | Meaning | Example |
|---------|---------|---------|
| `<1..N>` | Numeric range from 1 to N | `5HG6_<1..50>` = 50 IDs |
| `<0..N>` | Numeric range from 0 to N | `5HG6_<0..49>` = 50 IDs |
| `<A B C>` | Explicit set of values | `5HG6_42<S A L K>` = 4 IDs |
| `[_suffix]` | Optional suffix (0 or 1) | `5HG6[_1]` = `5HG6` or `5HG6_1` |
| `[_<pattern>]` | Arbitrary-count suffix | `5HG6[_<N><S A L K>]` = unknown count at config time |
| nested | Compositions | `5HG6_<1..50>_<1..3>` = 150 IDs |

**Key properties:**
- `<1..50>` is always expandable at config time (deterministic count).
- `<S A L K>` is always expandable (explicit set).
- `[_<pattern>]` means "zero or more suffixes matching this pattern" — the count depends on runtime data (e.g., how many positions match a selection). This is the mechanism that enables tools like Mutagenesis to work with variable-length selections.

### 2.3 The `[...]` Bracket: Runtime-Dependent Multiplicity

The square bracket `[...]` is the key innovation. It marks a suffix segment whose **repetition count is unknown at config time** but will be resolved at runtime.

Example — Mutagenesis on Fuse linker residues:

```python
fusions = Fuse(sequences=[donor, cam, acceptor],
               linker="GSG", linker_lengths=["0-3", "0-3"])
# fusions.ids pattern: "EBFP_<0..3>_CaM_<0..3>_EYFP"
# (16 variants with different linker lengths)

mut = Mutagenesis(original=fusions,
                  selection=fusions.tables.sequences.L1,  # per-row positions
                  mutate_to="SALK")
# mut.ids pattern: "EBFP_<0..3>_CaM_<0..3>_EYFP[_<N><S A L K>]"
# The [_<N><S A L K>] part means: for each linker position N,
# substitute with S, A, L, or K. Count depends on linker length per variant.
```

At config time, we know the input has 16 IDs and each will produce `len(linker) * 4` mutants — but `len(linker)` varies per variant (0 to 3 residues), so total output count is not a simple formula.

At **runtime**, the pipe script reads actual sequences, resolves selections, and produces the fully expanded CSV with all IDs.

### 2.4 File Path Patterns

File paths follow the same principle, using the ID pattern as the variable part:

```python
# Current (50 explicit paths):
files = ["/path/5HG6_0.pdb", "/path/5HG6_1.pdb", ..., "/path/5HG6_49.pdb"]

# Proposed (one pattern):
file_pattern = "/path/<id>.pdb"
```

The `<id>` token in file paths means "substitute the actual ID here." This replaces `files_contain_wildcards` entirely.

For tools where output files have a different naming convention (not `<id>.ext`), the pattern can include tool-specific templates:

```python
file_pattern = "/path/boltz_results/<id>/predictions/model_0.cif"
```

### 2.5 DataStream Changes

```python
@dataclass
class DataStream:
    name: str = ""
    id_pattern: str = ""           # NEW: compact pattern like "5HG6_<1..50>"
    file_pattern: str = ""         # NEW: "/path/<id>.pdb"
    ids: List[str] = field(...)    # Kept: lazily expanded from id_pattern
    files: List[str] = field(...)  # Kept: lazily expanded from file_pattern
    map_table: str = ""
    format: str = "pdb"
    is_lazy: bool = False          # NEW: True when [..] brackets present (runtime-dependent)
```

**Expansion rules:**
- `id_pattern` with only `<..>` tokens (no `[..]`): **fully expandable** at config time. `ids` list is populated immediately. `len()` returns exact count. Iteration and indexing work normally.
- `id_pattern` with `[..]` tokens: **lazy**. `is_lazy = True`. `ids` list may be partially expanded (the deterministic prefix part). `len()` returns the count of the deterministic prefix. Full expansion happens at runtime when the map_table CSV is read.

### 2.6 Downstream Tool Wiring

When a downstream tool receives a lazy DataStream:
- It can still access `map_table` (the CSV path) — this is the primary data contract.
- It can read `id_pattern` to understand the ID structure.
- For combinatorics (`predict_output_ids_with_provenance`), lazy streams are treated as single-axis inputs where the axis IDs are the **deterministic prefix** (the part before `[..]`). The `[..]` suffix is propagated.

Example:
```python
# Mutagenesis output: "EBFP_<0..3>_CaM_<0..3>_EYFP[_<N><S A L K>]"
# AlphaFold on mutagenesis output:
af = AlphaFold(sequences=mut)
# af.ids pattern: "EBFP_<0..3>_CaM_<0..3>_EYFP[_<N><S A L K>]"
# (same pattern — AlphaFold is 1:1 input→output)
```

For tools that multiply IDs (e.g., ProteinMPNN with `num_sequences=3`):
```python
pmpnn = ProteinMPNN(structures=rfd, num_sequences=3)
# Input pattern: "5HG6_<1..50>"
# Output pattern: "5HG6_<1..50>_<1..3>"  (appends a new range)
```

---

## 3. Implementation Impact

### 3.1 What Changes

| Component | Change |
|-----------|--------|
| `DataStream` | Add `id_pattern`, `file_pattern`, `is_lazy` fields. Lazy expansion logic in `__iter__`, `__getitem__`, `__len__`. |
| `DataStream.ids` | Still a list. Populated eagerly when pattern is deterministic, partially when lazy. |
| `DataStream.files` | Still a list. Populated from `file_pattern` + expanded IDs. Empty when lazy. |
| `files_contain_wildcards` | **Removed**. Replaced by `file_pattern` with `<id>` token. |
| `create_map_table()` | Only called at config time for deterministic patterns. For lazy patterns, the map_table CSV is written at runtime by the pipe script. Config time writes a "template" map_table (or nothing). |
| `generate_multiplied_ids()` | Returns a **pattern** (e.g., `"parent_<1..3>"`) instead of expanded list when called from `get_output_files()`. |
| `predict_output_ids_with_provenance()` | Works with patterns. For cartesian products, composes patterns (e.g., `"prot_<1..3>+lig_<1..5>"`). |
| Tool `get_output_files()` | Returns patterns in DataStream. Each tool defines its suffix pattern. |
| `biopipelines_io.py` (runtime) | `iterate_files()` and `iterate_values()` read the map_table CSV (which is fully expanded at runtime). No change needed — they already work from CSV. |
| `id_map_utils.py` | Pattern-aware matching. `<..>` segments in patterns can match any value during parent/child resolution. |
| `base_config.py` pipeline summary | Display patterns instead of listing all IDs. More readable. |
| Mutagenesis | Now accepts `selection` (string or table column ref) instead of just `position`. Output pattern includes `[..]` suffix for variable-count outputs. |
| HelpScripts `pipe_*.py` | No change for most. They already write the runtime CSV. They just need to write expanded IDs to map_table. |

### 3.2 What Does NOT Change

- **Runtime behavior**: pipe scripts already produce fully-expanded CSVs. They continue to do so.
- **Map table CSV format**: still `id, file, value, ...` columns. Always fully expanded at runtime.
- **Tool bash scripts**: no change. They call pipe scripts which handle data.
- **biopipelines_io runtime functions**: `iterate_files()`, `iterate_values()`, `resolve_file()` all read CSVs — they work on expanded data.
- **The `+` separator for multi-axis IDs**: unchanged.
- **The `_` separator for parent→child suffixes**: unchanged.
- **id_map_utils matching logic**: the 5-priority matching (exact, provenance, child, parent, sibling) still works — patterns just add a way to describe the space.

### 3.3 Migration Path

1. Add `id_pattern` and `file_pattern` to DataStream (backward-compatible: default to `""`).
2. Add pattern parser module (`id_patterns.py`) with expand/compose/is_lazy functions.
3. Update `generate_multiplied_ids()` to return patterns when possible.
4. Update each tool's `get_output_files()` to use patterns (can be done incrementally, tool by tool).
5. Remove `files_contain_wildcards` once all tools use `file_pattern`.
6. Update Mutagenesis to accept `selection` parameter with lazy output.

---

## 4. Detailed Examples

### 4.1 Simple Pipeline: RFdiffusion → ProteinMPNN → AlphaFold

```
RFdiffusion(pdb="5HG6", num_designs=50)
  ids pattern: "5HG6_<0..49>"
  file pattern: "/path/RFdiffusion_1/5HG6_<0..49>.pdb"   (= "/path/RFdiffusion_1/<id>.pdb")
  len: 50 (deterministic)

ProteinMPNN(structures=rfd, num_sequences=3)
  ids pattern: "5HG6_<0..49>_<1..3>"
  file pattern: none (CSV-based)
  map_table: "/path/ProteinMPNN_2/sequences.csv"
  len: 150 (deterministic)

AlphaFold(sequences=pmpnn)
  ids pattern: "5HG6_<0..49>_<1..3>"
  file pattern: "/path/AlphaFold_3/<id>.pdb"
  len: 150 (deterministic)
```

All deterministic — patterns expand fully at config time, exactly as today but more compact.

### 4.2 FRET Pipeline with Mutagenesis (the motivating case)

```
Fuse(sequences=[donor, cam, acceptor], linker="GSG", linker_lengths=["0-3", "0-3"])
  ids pattern: "EBFP_<0 1 2 3>_CaM_<0 1 2 3>_EYFP"
  len: 16 (deterministic)
  tables.sequences columns: id, sequence, lengths, S1, L1, S2, L2, S3

Mutagenesis(original=fusions, selection=fusions.tables.sequences.L1, mutate_to="SALK")
  ids pattern: "EBFP_<0 1 2 3>_CaM_<0 1 2 3>_EYFP[_<N><S A L K>]"
  len: 16 (deterministic prefix; actual total unknown at config time)
  is_lazy: True

  At runtime, for variant with L1="74-76" (3 positions):
    EBFP_3_CaM_2_EYFP_74S, EBFP_3_CaM_2_EYFP_74A, EBFP_3_CaM_2_EYFP_74L, ...
    EBFP_3_CaM_2_EYFP_75S, EBFP_3_CaM_2_EYFP_75A, ...
    EBFP_3_CaM_2_EYFP_76S, EBFP_3_CaM_2_EYFP_76A, ...
    = 12 mutants for this variant

  For variant with L1="" (0-length linker, no positions):
    = 0 mutants (this variant may be in missing table or have original only)

Boltz2(proteins=mut)
  ids pattern: "EBFP_<0 1 2 3>_CaM_<0 1 2 3>_EYFP[_<N><S A L K>]"
  is_lazy: True (inherits laziness)
  file pattern: "/path/Boltz2_4/<id>.cif"
```

### 4.3 Multi-input Combinatorics (Gnina docking)

```
ProteinMPNN → sequences: "5HG6_<0..49>_<1..3>"  (150 IDs)
CompoundLibrary → compounds: "aspirin", "ibuprofen", "caffeine"  (3 IDs)

Gnina(structures=af, ligands=compounds)
  ids pattern: "5HG6_<0..49>_<1..3>+<aspirin ibuprofen caffeine>"
  len: 450 (deterministic, cartesian product)
  file pattern: "/path/Gnina_5/<id>.sdf"
```

The `+` multi-axis separator composes naturally with `<..>` patterns.

---

## 5. Pattern Parsing Specification

### 5.1 Grammar

```
pattern     := segment+
segment     := literal | range | set | optional
literal     := <any characters except < > [ ]>
range       := '<' start '..' end '>'
set         := '<' value (' ' value)* '>'
optional    := '[' pattern ']'
```

### 5.2 Core Functions

```python
def parse_pattern(pattern: str) -> PatternNode:
    """Parse pattern string into AST."""

def expand_pattern(pattern: str) -> List[str]:
    """Fully expand pattern to list of strings.
    Raises LazyPatternError if pattern contains [...] brackets."""

def try_expand_pattern(pattern: str) -> Tuple[List[str], bool]:
    """Expand what we can. Returns (expanded_ids, is_complete).
    For lazy patterns, expands the deterministic prefix."""

def is_lazy(pattern: str) -> bool:
    """True if pattern contains [...] brackets."""

def compose_patterns(parent: str, suffix: str) -> str:
    """Compose parent_suffix (e.g., 'A_<1..3>' + '<1..5>' → 'A_<1..3>_<1..5>')."""

def pattern_count(pattern: str) -> Optional[int]:
    """Return total count of expanded IDs, or None if lazy."""

def pattern_to_regex(pattern: str) -> str:
    """Convert pattern to regex for matching actual IDs against the pattern."""
```

---

## 6. Pros

1. **Compact representation**: `"5HG6_<0..49>_<1..3>"` instead of 150 explicit strings. Config output, JSON files, and pipeline summaries become more readable.

2. **Enables runtime-dependent tools**: Mutagenesis with per-row selections, future tools that produce variable outputs per input. The `[..]` mechanism cleanly separates "known at config time" from "known at runtime."

3. **Unifies wildcard handling**: `file_pattern = "/path/<id>.pdb"` replaces the ad-hoc `files_contain_wildcards` flag. One mechanism for all cases.

4. **Backward compatible**: Tools can continue to produce explicit ID lists. Patterns are opt-in. Migration can happen tool-by-tool.

5. **Better pipeline display**: Instead of printing 1,500 IDs, the pipeline summary shows `"5HG6_<0..49>_<1..3>"` — immediately communicating the structure.

6. **Natural composition**: Patterns compose well. ProteinMPNN on RFdiffusion output naturally becomes `parent_pattern + "_<1..N>"`. No special logic needed.

7. **Consistent with biopipelines philosophy**: Config time = declare intent; runtime = resolve data. Lazy patterns are a declaration of intent that gets resolved at runtime.

## 7. Cons and Risks

1. **Complexity in DataStream**: `__len__`, `__iter__`, `__getitem__`, `filter_by_ids` all need pattern awareness. Lazy streams can't be indexed by integer (which integer? the pattern hasn't been expanded). Need clear error messages.

2. **Downstream tool confusion**: A tool receiving a lazy DataStream needs to know it can't enumerate IDs. Most tools don't need to — they just pass `map_table` to their pipe script. But tools that do per-ID config-time logic (like Boltz2 writing per-sequence YAML files) need adaptation.

3. **Map table at config time**: Currently, `create_map_table()` writes a CSV with all IDs at config time. For lazy patterns, this CSV either doesn't exist or is partial. Downstream tools that read map_tables at config time (rare, but possible) would break.

4. **Pattern nesting complexity**: Deep nesting like `A_<1..3>[_<N><S A L K>[_<1..5>]]` could get confusing. Need clear documentation and limits on nesting depth.

5. **id_map_utils matching**: Pattern matching adds complexity. When a downstream tool asks "which IDs from stream A match stream B?", and one or both use patterns, the matching logic needs to reason about patterns — not just string comparison.

6. **Debugging**: When something goes wrong, patterns hide the actual IDs. The user sees `"5HG6_<0..49>_<1..3>"` but the actual failing ID is `5HG6_23_2`. Need good error messages that show both the pattern and the specific failing ID.

7. **Two representations**: Having both `id_pattern` and `ids` list is a source of truth duplication. Need clear rules about which is authoritative (answer: `id_pattern` is the source of truth; `ids` is a cache).

8. **Config-time file writing**: Tools like `generate_script()` currently reference specific file paths (e.g., checking `if [ ! -f "/path/5HG6_0.pdb" ]`). With patterns, they'd need to iterate at runtime instead. Most tools already do this through pipe scripts, but some have inline bash checks.

---

## 8. Edge Cases and Concerns

### 8.1 What about Panda pool?

Panda with `pool=` copies structure files based on filtered IDs. With lazy patterns, the pool stream might be lazy. At runtime, this is fine — the pipe script reads the CSV and copies files. At config time, Panda just needs to know the pool's `file_pattern` to generate the copy commands (or defer to its pipe script).

### 8.2 What about Bundle/Each combinatorics?

Patterns compose with the `+` separator for multi-axis cartesian products. For Bundle, all IDs appear together — the pattern represents the full set. No fundamental conflict.

### 8.3 What about iteration in pipeline scripts?

```python
for item in rfd:  # Iterates over 50 single-item DataStreams
    tool = SomeTool(structures=item)
```

With deterministic patterns, iteration expands the pattern and yields items as before. With lazy patterns, iteration yields the **deterministic prefix** items (each carrying the lazy suffix). This means the `SomeTool` instance per item would still have a lazy DataStream — which is correct, because each prefix item may produce a variable number of outputs.

### 8.4 What about len()?

- Deterministic pattern: `len()` returns exact count (e.g., 150).
- Lazy pattern: `len()` returns the deterministic prefix count (e.g., 16 for the Fuse output before Mutagenesis). Could also return `None` or raise, but returning the prefix count seems most useful — it tells you "at least this many distinct input lineages."

### 8.5 Zero-output variants

In the Mutagenesis example, a variant with a 0-length linker has no positions to mutate. The runtime pipe script would produce 0 rows for that input. The `missing` table would capture it. This is already how tools handle filtered-out inputs.

---

## 9. Ideas for Improvement

### 9.1 Pattern Metadata

Store metadata alongside patterns to help downstream tools:

```python
id_pattern_info = {
    "pattern": "5HG6_<0..49>_<1..3>",
    "prefix_count": 50,      # RFdiffusion outputs
    "suffix_count": 3,       # ProteinMPNN sequences per structure
    "total_count": 150,      # Total (None if lazy)
    "lazy_segments": []       # List of [...] segment descriptions
}
```

### 9.2 Named Segments

Allow naming segments for provenance:

```python
"5HG6_<design:0..49>_<seq:1..3>"
```

This way, when composing, we know what each segment represents. Provenance columns can be auto-generated: `design.id`, `seq.id`.

### 9.3 Validate at Runtime

Add a runtime validation step in pipe scripts:

```python
# In pipe script, after generating all IDs:
validate_ids_against_pattern(actual_ids, expected_pattern)
```

This catches bugs where the runtime output doesn't match the config-time pattern.

### 9.4 Consider NOT storing `ids` list at all for pattern-based streams

If the pattern is the source of truth, maybe `ids` should only exist as a computed property:

```python
@property
def ids(self) -> List[str]:
    if self.is_lazy:
        raise LazyPatternError("Cannot enumerate lazy pattern IDs at config time")
    return expand_pattern(self.id_pattern)
```

This avoids duplication but is a bigger change. Could be phase 2.

---

## 10. Implementation Priority

### Phase 1: Pattern Infrastructure
- `id_patterns.py` module with parsing, expansion, composition
- `DataStream` gains `id_pattern` and `file_pattern` fields
- Patterns are optional: tools can still pass explicit `ids` lists

### Phase 2: Tool Migration (deterministic patterns only)
- Update `generate_multiplied_ids()` to produce patterns
- Update RFdiffusion, ProteinMPNN, AlphaFold, Boltz2, etc.
- Remove `files_contain_wildcards`, replace with `file_pattern`

### Phase 3: Lazy Patterns
- Implement `[..]` bracket parsing and `is_lazy` flag
- Update Mutagenesis to accept `selection` parameter
- Update `predict_output_ids_with_provenance()` for lazy composition
- Test with FRET + Mutagenesis pipeline

### Phase 4: Cleanup
- Remove explicit `ids` list from DataStream for pattern-based streams (make it computed)
- Update pipeline display to show patterns
- Update documentation
