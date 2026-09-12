# BioPipelines Developer Manual

## Index

- [Architecture](#architecture)
  - [Two-Phase Execution](#two-phase-execution)
  - [Core Classes](#core-classes)
  - [Data Flow](#data-flow)
  - [Internal Conventions](#internal-conventions)
    - [Value-Based `csv` Streams](#value-based-csv-streams)
    - [`resi-csv` Streams](#resi-csv-streams)
    - [The Ligand Contract: compounds = chemistry, structures = coordinates](#the-ligand-contract-compounds--chemistry-structures--coordinates)
    - [Contract Enforcement](#contract-enforcement)
    - [Internal Tools and Input Shorthands](#internal-tools-and-input-shorthands)
- [IDs: Configuration Time vs Execution Time](#ids-configuration-time-vs-execution-time)
  - [ID Patterns](#id-patterns)
  - [Lazy IDs](#lazy-ids)
  - [The Rule](#the-rule)
  - [How to Iterate Structures in Generated Bash](#how-to-iterate-structures-in-generated-bash)
  - [How to Resolve a Single File](#how-to-resolve-a-single-file)
  - [Passing IDs to pipe scripts](#passing-ids-to-pipe-scripts)
- [Tool Development](#tool-development)
  - [Creating a Tool](#creating-a-tool)
  - [Required Methods](#required-methods)
  - [Renaming a Parameter: `PARAMETER_ALIASES`](#renaming-a-parameter-parameter_aliases)
  - [Forwarding Untyped Kwargs: `FORWARD_UNKNOWN_KWARGS`](#forwarding-untyped-kwargs-forward_unknown_kwargs)
  - [Path Descriptors](#path-descriptors)
  - [Script Generation](#script-generation)
  - [Output Prediction](#output-prediction)
    - [File Templates with \<id\>](#file-templates-with-id)
  - [Map Table Contract](#map-table-contract)
  - [Install Scripts and the `$INSTALL_SUCCESS` Contract](#install-scripts-and-the-install_success-contract)
  - [Environments and `ENV_NAME`](#environments-and-env_name)
- [Shell Safety](#shell-safety)
  - [Why](#why)
  - [Where Validation Lives](#where-validation-lives)
  - [Wiring It Up in a New Tool](#wiring-it-up-in-a-new-tool)
  - [What Not to Validate This Way](#what-not-to-validate-this-way)
- [pipe scripts Development](#pipe-scripts-development)
  - [biopipelines_io Module](#biopipelines_io-module)
  - [pdb_parser Module](#pdb_parser-module)
  - [Table References](#table-references)
  - [Error Handling](#error-handling)
- [ID Generation and Provenance](#id-generation-and-provenance)
  - [Output ID Rules](#output-id-rules)
  - [Provenance Columns](#provenance-columns)
  - [Shared Utilities](#shared-utilities)
  - [Pipeline vs SLURM Agreement](#pipeline-vs-slurm-agreement)
- [Testing](#testing)
  - [Scope](#testing-scope)
  - [The Mock Tool](#the-mock-tool)
  - [Running the Suite](#running-the-suite)
- [Code Principles](#code-principles)
- [Working with Git](#working-with-git)
- [Working with Claude Code](#working-with-claude-code)

---

## Architecture

### Two-Phase Execution

| Phase | Location | What Happens |
|-------|----------|--------------|
| **Configuration time** | `biopipelines/` | Python generates bash scripts, predicts outputs |
| **Execution time** | `pipe_scripts/` | Bash scripts execute, `pipe_*.py` scripts run |

Tools never execute computations directly. They:
1. Validate parameters
2. Predict output file paths
3. Generate bash scripts

### Core Classes

**BaseConfig** (`base_config.py`) — Abstract base for all tools:
```python
class MyTool(BaseConfig):
    TOOL_NAME = "MyTool"

    def validate_params(self): ...
    def configure_inputs(self, pipeline_folders): ...
    def generate_script(self, script_path): ...
    def get_output_files(self): ...
```

**DataStream** (`datastream.py`) — Unified container for files and IDs:
```python
DataStream(
    name="structures",
    ids=["prot_1", "prot_2"],
    files=["/path/prot_1.pdb", "/path/prot_2.pdb"],
    map_table="/path/to/map.csv",
    format="pdb"
)
```

DataStream attributes:
- `ids` — List of identifiers (may contain compact patterns)
- `files` — Either a list of file paths or a single string path. Three forms:
  - empty list `[]` — value-based stream; content lives in `map_table`
  - list `["<id>.pdb"]` or `["a.pdb", "b.pdb", ...]` — per-ID files (template or explicit)
  - string `"path/to/shared.fasta"` — *shared-file* form: a single artifact (e.g. multi-record FASTA) covering all ids. `is_shared_file` returns True.
- `map_table` — CSV with additional metadata
- `format` — Data format (pdb, cif, fasta, csv, smiles, etc.)

Streams may advertise more than one possible file format with a pipe-delimited
`format`, e.g. `"pdb|cif"` when upstream intentionally leaves structures in
mixed PDB/mmCIF form. Keep the serialized field as that string for compatibility,
but do not hand-roll checks like `stream.format in ("pdb", "cif", "pdb|cif")`.
Use the helpers on `DataStream`:

```python
stream.formats                         # ("pdb", "cif")
stream.has_format("pdb")               # membership
stream.has_only_formats("pdb", "cif")  # True for "pdb", "cif", or "pdb|cif"
```

`has_only_formats()` is a whitelist/subset check: it means there is no advertised
format outside the allowed set. Thus `"pdb"` passes `has_only_formats("pdb", "cif")`,
while `"pdb|sdf"` fails. Tools that require homogeneous PDB input should use
`has_only_formats("pdb")`, not `has_format("pdb")`.

**Shared-file streams and slicers.** When a tool (e.g. Panda) filters a shared-file stream down to a subset of ids, the underlying artifact must also be sliced — copying the whole file would leak stale records past the filter. Format-aware slicers live in `biopipelines/stream_slicers.py` keyed by stream `format`. Built-in slicers cover `fasta`/`fa` and `csv`. To add another format, decorate a function with `@register("sdf")` (or the relevant format key) — Panda picks it up automatically. No silent fallback to "copy whole file": an unregistered format raises `ValueError`.

**TableInfo** (`data_containers.py`, re-exported from `base_config.py`) — Metadata for CSV outputs:
```python
TableInfo(
    name="results",
    path="/path/to/results.csv",
    columns=["id", "score", "pLDDT"],
    description="Analysis results"
)
```

Reading any public attribute off a `TableInfo` gives you a column reference — `results.score` is `TableReference(path, "score")` — whether or not the name appears in `columns`, because a tool often knows a table's path long before its schema. Names starting with `_` are the one exception: they raise `AttributeError`, so `hasattr(table, "_entries")` and the copy/pickle protocol's dunder lookups answer honestly instead of every duck-type probe against a table coming back true. Do not duck-type a `TableReference` with `hasattr(x, "path") and hasattr(x, "column")` — a `TableInfo` still satisfies that, since `path` and `column` are legal column names; use `isinstance(x, TableReference)`.

**StandardizedOutput** (`outputs.py`, re-exported from `base_config.py`) — Tool output wrapper:
```python
output.structures  # DataStream
output.sequences   # DataStream
output.compounds   # DataStream
output.msas        # DataStream
output.tables      # TableContainer
output.output_folder  # str
```

### Data Flow

```
Entity/Tool → DataStream → Downstream Tool
     ↓
   Tables (TableInfo)
     ↓
  Column References → (TableInfo, "column")
```

Tools communicate via StandardizedOutput. Each tool predicts outputs that downstream tools consume:

```python
rfd = RFdiffusion(contigs="50-100", num_designs=5)
# rfd.structures is a DataStream with predicted paths

mpnn = ProteinMPNN(structures=rfd, num_sequences=2)
# mpnn receives rfd.structures DataStream
```

### Internal Conventions

These are framework-wide rules that aren't obvious from any single tool but that every tool must respect. They're invariants, not suggestions — a tool that violates them will silently mislead downstream tools.

#### Value-Based `csv` Streams

A DataStream is *value-based* when its content lives entirely in its `map_table` CSV rather than in per-id files. Such a stream has `files=[]` and **`format="csv"`**. The rule is sharp:

- A value-based stream's `format` is **`"csv"`**. What the stream *carries* — SMILES, a residue code, a protein sequence — is expressed as **columns of the map_table** (the stream's metadata), and is **never** encoded in the `format` field. There is no `format="smiles"` / `format="code"` / `format="sequence"`; those are all `format="csv"` with the relevant column present.
- The signal that a stream is value-based is **`files == []`**, not the format string. `format` is purely descriptive (`pdb`, `sdf`, `csv`, `a3m`, `resi-csv`, …) and a file-based stream may legitimately use a csv-shaped format (e.g. a per-id MSA). Consumers that need to know whether to collect per-id files (e.g. `panda.py`) test `stream.files` directly.
- The map_table is still written **only at runtime by the pipe script** (see [Map Table Contract](#map-table-contract)); `get_output_files()` just declares `files=[]`, `map_table=<path>`, `format="csv"`.

Read values from a value-based stream at runtime with `get_value(ds, id, column=...)` or `iterate_values(ds, columns=[...])` from `biopipelines_io` — never by parsing the CSV by hand.

#### `resi-csv` Streams

A `resi-csv` stream is the per-residue counterpart: it carries **one CSV file per id** (file-based), and each CSV has multiple rows, one per residue, with a `resi` column plus one or more value columns. CABSflex's RMSF output is the canonical producer (`biopipelines/cabsflex.py`, `format="resi-csv"`, columns `id, chain, resi, rmsf`), and `Selection` consumes it stream-based (reading the `resi` column against a `"column op value"` threshold). Use `resi-csv` whenever a tool emits a per-residue numeric profile that another tool will threshold or select on.

#### The Ligand Contract: compounds = chemistry, structures = coordinates

Ligands split cleanly across two streams, and every tool must honour the split:

- The **`compounds`** stream is always a value-based `csv` stream (`files=[]`, `format="csv"`). It carries the ligand *chemistry and identity* — `id, format, code, lookup, source, ccd, cid, cas, smiles, name, formula, file_path` — in its map_table. **compounds never carries coordinate files.**
- The **`structures`** stream carries the ligand *coordinates* — the `.sdf` / `.pdb` / `.cif` / `.mol2` files (`format="sdf"`, etc.). When a tool needs a 3-D ligand, the user produces it explicitly with `OpenBabel(compounds=lig, convert_3d="sdf")`, whose `structures` output is the sdf and whose `compounds` output is the chemistry passthrough.

Two consequences:

1. **Tools never take a ligand SDF path or a bare 3-letter `ligand_code` string.** They take a `Ligand` (or any tool's compounds/structures output) and read the residue `code` (and `smiles`) from the compounds stream's map_table at runtime. A ligand that only names an existing HETATM code is constructed as `Ligand(codes="ZIT")` — a one-row `compounds` csv (`format="csv"`, `code="ZIT"`, empty `smiles`) and no structures stream. `codes` is one parameter serving three modes and what else the caller passes selects between them (code-only, labelled chemistry, or carving a residue out of `structures=`), so a consumer must not assume a `compounds` stream carries `smiles` just because its `code` is set — see [Which mode `codes` selects](tool/inputs_io.md#ligand) for the mode table.
2. **Producers that rename a ligand emit an updated compounds stream.** A tool like Boltz2 assigns its own residue codes (`LIG`, `LIG01`, …) when it writes the complex. It must emit a fresh `compounds` stream carrying the **same ids and SMILES** as its input but with the `code` column overwritten with the codes it actually assigned. Because the ids are unchanged, `{compounds}.id` provenance columns (see [Provenance Columns](#provenance-columns)) stay joinable — a rename changes a *cell*, not an *id*. Passing that output downstream then yields the correct code automatically, with nobody having to guess.

#### Contract Enforcement

The conventions above used to exist only as prose, which is how the manual's own tool template came to declare `compounds` as `format="sdf"` for years while five tools copied it. `biopipelines/contract_enforcement.py` is now the framework's single home for the checks that hold them, and a violation announces itself on stderr:

```
[contract:compounds_format] compounds stream declared format='sdf', expected 'csv' (in MyTool)
    compounds carries ligand chemistry and identity in its map_table; coordinates belong on a structures stream. See the Ligand Contract in docs/developer_manual.md.
```

**Every check is a pure function.** A `check_*` function takes the facts it judges and returns a `Violation` — a frozen dataclass of `check`, `message` and an optional `hint`, rendered as the two lines above — or `None`, which is a pass. It never prints, never raises, and never looks at a severity. Acting on the result is `report(violation)`'s job, and it is the only thing that does: `None` returns immediately, `off` returns, `warn` prints `violation.render()` to stderr, `raise` raises `ContractViolation` carrying the same text. `check_all(*violations)` reports several in order and stops at the first `raise`-severity one. Keeping the judgement and the reaction apart is what makes a check testable without capturing output and reusable from a call site that wants a different reaction.

**The severity ladder is `off` / `warn` / `raise`**, resolved per check by `severity_for(check)` in this order:

1. `BIOPIPELINES_ENFORCE_<CHECK>` in the environment — e.g. `BIOPIPELINES_ENFORCE_COMPOUNDS_FORMAT=raise` to promote one check for one run.
2. `SEVERITY[<check>]`, the table at the top of the module.
3. `DEFAULT_SEVERITY`, which is `warn`, for a check with no table entry.

**There is no global switch.** `BIOPIPELINES_ENFORCE` was removed: it promised two things and delivered neither. `=off` did not silence everything, because the rows reporting a working construction rather than a fault were skipped in both directions; and `=raise`, the setting this manual named for CI, refused every pipeline using `FORWARD_UNKNOWN_KWARGS` — a supported feature — until an exemption set was added to patch that, which is what broke the other direction. Name the checks you want fatal, one variable each. An unrecognised value is ignored rather than rejected, and says so once, so a typo falls through to the next rung instead of turning enforcement off.

**Checks ship at `warn`.** A new check goes into `SEVERITY` as `WARN`, so a violation is reported without failing a pipeline that would otherwise run — the framework tells you about a contract breach, it does not hold your job hostage to one. Promote a row to `RAISE` only once a warn phase has shown the tree clean for it; until then, `BIOPIPELINES_ENFORCE_<CHECK>=raise` is how you make one fatal for your own run or in CI.

**Reports are deduplicated.** `report()` keeps the rendered text of every warning it has printed in a module-level set and prints each distinct line once. Without it a violation inside a per-id loop prints hundreds of identical lines and buries everything else. `reset_reported()` clears the cache — tests need it, and so does a long-lived session (a notebook) where you want to see the same warning again after a fix.

**Adding a check.** Four steps, no more:

1. Write a `check_<name>(...) -> Optional[Violation]` that returns `Violation(check="<name>", message=..., hint=...)` or `None`. Put the actionable instruction in the `hint` — what to change, and what goes wrong if it is not changed.
2. Add `"<name>": WARN` to `SEVERITY`.
3. Call it from **the one place** the contract applies, through `report()` or `check_all()`. `check_stream()` is the model: `DataStream.__post_init__` is the single place a stream is built, so every stream-level contract is checked there — above its early returns, so pattern and shared-file streams are covered too — and nowhere else.
4. Add a case to `tests/test_contract_enforcement.py`.

**Naming the offender.** A `DataStream` does not know which tool built it, and must stay constructible standalone because pipe scripts import it directly at runtime — so the tool name cannot travel through its signature. Instead `BaseConfig.__init_subclass__` wraps each subclass's `get_output_files()` in `contract_enforcement.tool_context(TOOL_NAME)`, a `contextvars` context that the stream checks read `where` from. Wrapping once in `__init_subclass__` rather than at ~90 call sites is what keeps that free.

**Three class attributes a tool may set:**

- `USER_STREAM_NAMES` (default `False`) — set it to `True` when the tool's *user* names its streams, so no registry can know them in advance. `Mock(streams={...})` and `Scripting` are the two, and `tool_context` carries the flag so `check_stream_name` is skipped for that tool only. The format checks still apply: a user-named stream is still a stream, so `compounds`-as-`sdf` and a value-based stream with a non-csv format are still reported. Only the "is this name in `KNOWN_STREAM_NAMES`?" question is dropped, and only because it is unanswerable.
- `FORWARD_UNKNOWN_KWARGS` (default `""`, off) — set it to `"argparse"` or `"hydra"` to declare that this wrapper renders a constructor keyword it does not type onto the upstream command line, which switches `check_kwargs` from `check_no_unknown_kwargs` to the `check_probable_typo` + `check_forwarded_kwargs` pair. It only works on a wrapper whose own bash contains that command, and `assert_can_forward` refuses it at class creation otherwise. See [Forwarding Untyped Kwargs](#forwarding-untyped-kwargs-forward_unknown_kwargs).
- `ENV_NAME` — the tool's own default environment name, declared beside `TOOL_NAME` and `TOOL_VERSION`. See [Environments and `ENV_NAME`](#environments-and-env_name); it is not a contract-enforcement attribute, but it is the other class attribute a new tool has to remember to declare.

A new stream name is declared by adding it to `KNOWN_STREAM_NAMES` in the module. That is not a permission list — nothing is forbidden, and an unregistered name is a `warn`, never an error. It exists so a typo (`scors` for `scores`) is caught at construction rather than discovered later as a consumer that reads the stream back by name and silently gets nothing.

Leftover constructor kwargs are checked through the same machinery: `BaseConfig.__init__` calls `check_kwargs(...)` once, after the subclass has bound its own named parameters, so anything left is a framework key from `RESERVED_KWARGS`, a typo, or an option the tool forwards upstream.

**The roster.** `SEVERITY` in `contract_enforcement.py` is the authoritative list; these are its rows today, all at `warn`:

| Check | Fires when | Silence with |
|---|---|---|
| `stream_name` | a stream's name is not in `KNOWN_STREAM_NAMES` | `BIOPIPELINES_ENFORCE_STREAM_NAME=off` |
| `compounds_format` | a `compounds` stream declares a format other than `csv` | `…_COMPOUNDS_FORMAT=off` |
| `value_based_format` | a `files=[]` stream declares a non-csv format | `…_VALUE_BASED_FORMAT=off` |
| `unknown_kwargs` | a constructor kwarg matches no parameter, on a tool that forwards nothing (`check_no_unknown_kwargs`), or is a near-miss on a real one (`check_probable_typo`) | `…_UNKNOWN_KWARGS=off` |
| `forwarded_kwargs` | an untyped kwarg is about to be rendered onto a forwarding tool's upstream command line | `…_FORWARDED_KWARGS=off` |
| `pattern_selection` | a pattern set selects 0 of a non-empty row set, so the selection silently did nothing | `…_PATTERN_SELECTION=off` |
| `deprecated_alias` | a `DEPRECATED_ALIASES` spelling is used — the value still binds | `…_DEPRECATED_ALIAS=off` |
| `code_only_ligand` | a `Ligand` is built from `codes=` alone, so it carries a residue code and no chemistry | `…_CODE_ONLY_LIGAND=off` |
| `id_match_consistency` | a remap's id-match score is low enough to be worth seeing, though nothing gates on it | `…_ID_MATCH_CONSISTENCY=off` |

`forwarded_kwargs`, `deprecated_alias`, `code_only_ligand` and `id_match_consistency` report a construction that *works*, so they are the rows a blanket "make everything fatal" should never have touched — which is why there is no longer a switch that could. Promote any of them by name when you want it fatal.

Two of them came in with the alias mechanism below. Both report something that *works* rather than something broken: a deprecated spelling still binds its parameter, and a code-only `Ligand` is a legitimate object — the notice exists because a forgotten `smiles=` now produces one silently instead of raising, and the failure would otherwise surface in whatever downstream tool needed the molecule. `check_code_only_ligand` is deliberately not called for the `ligand="LIG"` shorthand, where naming a residue code is the entire intent.

#### Internal Tools and Input Shorthands

A tool constructed with the reserved kwarg `_internal=True` is a **real** tool: it auto-registers in the active pipeline, generates a script, and executes in normal execution order. It is only *hidden from the public layout*. This is how an ergonomic shorthand like `ligand="LIG"` materializes the entity it stands for without cluttering the user's numbered steps.

The framework keeps three counters on the pipeline:

- `execution_order` — every tool, the true run order (drives nothing user-visible directly).
- `public_step_order` — public tools only; drives public output-folder names and the displayed `Step NNN`.
- `internal_order` — internal tools only; drives `.internal/` folder names.

Each tool gets a `script_basename` (set in `set_pipeline_context`) that is the **single source of truth** for its `RunTime/<basename>.sh`, `Logs/<basename>.log`, and `ToolOutputs/<basename>.json`. Public tools use `NNN_<Tool>` (the public step); internal tools use `.internal/NNN_<Tool>` (the internal order, nested in a `.internal/` subdir of each). Never re-derive these names from a list index — read `tool.script_basename`. Output folders follow the same split: public tools go under `<Folder stack>/NNN_<Tool>/`, internal tools under `.internal/NNN_<Tool>/` (internal placement ignores the `Folder()` stack).

Because a public tool's *folder* number and its *script/manifest* number both come from `public_step`, they always agree (`002_Tool/` ↔ `RunTime/002_Tool.sh`). The numbering you read off filenames is therefore the public step order, **not** the raw `execution_order` — those differ whenever internal tools are present (an internal tool bumps `execution_order` but not `public_step`). True run order is expressed only by the sequence of invocations in `pipeline.sh` (which iterates `self.tools`), where an internal tool appears right before the consumer that created it. Don't infer run order from public filenames.

**The shared resolver.** Tools that accept the StandardizedOutput / DataStream / shorthand triad normalize the input with `resolve_basic_input` (in `input_standardization.py`) instead of hand-rolling an `isinstance` ladder:

```python
from .input_standardization import resolve_basic_input
from .ligand import Ligand

# StandardizedOutput -> streams.compounds; DataStream -> itself;
# str "LIG" -> Ligand(codes="LIG", _internal=True) -> its compounds stream
self.ligand_stream = resolve_basic_input(ligand, Ligand, "compounds", "codes", allow_none=False)
```

The signature is `resolve_basic_input(obj, cls, stream, argument, *, allow_none=True)`. A bare string is promoted to `cls(**{argument: obj}, _internal=True)` and its `<stream>` is returned. The promotion works both inside and outside a pipeline: inside, the entity auto-registers and `cls(...)` already returns a `StandardizedOutput`; standalone it returns the raw tool instance, which the resolver wraps with `StandardizedOutput(entity.get_output_files())` to reach the same stream. Inside a pipeline the auto-registration also means the entity is constructed *during the consuming tool's `__init__`* and therefore runs **before** the consumer — ordering is correct for free. Only a non-string, non-DataStream, non-StandardizedOutput value raises. Apply the shorthand **only to the parameter whose stream the tool actually consumes** — e.g. a tool's compounds-reading `ligand` gets it, but a coordinate-reading `reference_ligand` (which needs a real 3-D structure) does not.

---

## IDs: Configuration Time vs Execution Time

### ID Patterns

DataStream IDs can be stored in compact pattern form instead of listing every ID explicitly:

| Pattern | Meaning | Expansion |
|---------|---------|-----------|
| `prot_<0..2>` | Numeric range | `prot_0`, `prot_1`, `prot_2` |
| `<A B C>` | Enumeration | `A`, `B`, `C` |
| `<0..1>_<A B>` | Multi-slot — cartesian product of the slots | `0_A`, `0_B`, `1_A`, `1_B` |
| `prot_<N><S A L K>` | Multi-slot, first slot a one-element set | `prot_NS`, `prot_NA`, `prot_NL`, `prot_NK` |

A slot's content is parsed by `_parse_slot` in `id_patterns.py`: it is a numeric range only when it contains `..`, and otherwise a **space-separated enumeration**. So `<N>` is not a placeholder or a wildcard — it is the one-element set containing the literal letter `N`, and `expand_pattern('prot_<N>')` returns `['prot_N']`. Multiple slots in one string multiply out as a cartesian product in left-to-right, row-major order.

The `<..>` angle-bracket patterns are **deterministic** — they can be fully expanded at configuration time, however many slots they have. Laziness comes only from `[...]` square brackets (next section); `is_lazy('prot_<N><S A L K>')` is `False`.

### Lazy IDs

Bracket `[...]` segments mark parts of an ID that depend on runtime data and **cannot be expanded at configuration time**:

```
prot[_<N><S A L K>]
```

Here `prot` is the deterministic prefix, and `[_<N><S A L K>]` is a lazy suffix whose actual values come from an upstream tool's output (which doesn't exist yet at config time). Lazy IDs resolve at runtime by **selecting** the rows of the DataStream's `map_table` CSV that the pattern covers — the map_table is the source of truth, the pattern is the selector.

A pattern is just an id with unresolved slots. There is no separate "partially expanded" concept: `ids` already holds patterns that may mix deterministic `<..>` slots and lazy `[...]` brackets. At configuration time, `ids_expanded` on a lazy DataStream expands the deterministic slots **and keeps the brackets**, so each element is still a valid lazy pattern:
```python
# ds.ids == ["prot_<1 2>[_<N>]"]
ds.ids_expanded  # → ["prot_1[_<N>]", "prot_2[_<N>]"] (still lazy, not concrete)
```

At execution time (with `_runtime_mode=True`, as set by `load_datastream()`), `ids_expanded` reads the full set of concrete IDs straight from the map_table:
```python
ds = load_datastream("structures.json")
ds.ids_expanded  # → ["prot_1_1S", "prot_2_2A", ...] (concrete rows)
```

To select map_table rows against a pattern at runtime (the shared, one-map helper used by combinatorics, stitch, and any bulk consumer), use `id_patterns.select_ids(patterns, row_ids)` or `id_patterns.resolve_pattern_ids(patterns, map_table)`. Matching is exact everywhere except the unresolved slots: a deterministic slot matches its own values, and inside a `[...]` bracket the *literal* text stays literal while only the `<..>` slots become wildcards. So `design_1[_<N>]` covers `design_1_7` but not `design_10_7` — the bracket's leading `_` is what separates them. (`glob_from_lazy` is the separate, genuinely glob-shaped helper, used for matching real filenames on disk rather than map_table rows.) Because the rows are the source of truth, this honors any upstream filter automatically — a filtered stream simply has fewer rows.

#### Where the delimiter sits decides what the bracket means

A `[...]` bracket is a **zero-or-more repetition** of its content, and each slot value is bounded by the `_` suffix delimiter — a value never crosses one. Those two rules are the whole contract; the two spellings below are not special cases but consequences of where the caller puts the delimiter, so no producer has to know about either.

| Pattern | Compiles to | Covers | Does not cover |
|---------|-------------|--------|----------------|
| `parent[_<?>]` | `parent(?:_[^_]+)*` | `parent`, `parent_4E`, `parent_4E_6U`, `parent_4E_6U_9Z` | `parentX` |
| `parent_[<?>]` | `parent_[^_]*` | `parent_4E`, `parent_1` | `parent`, `parent_4E_6U` |

Note the `<?>`. A slot is a literal value set, so `parent[_<N>]` compiles to `parent(?:_N)*` and covers `parent`, `parent_N`, `parent_N_N` — the letter N, nothing else. Write `<?>` for one segment and `<#>` for a number.

With the delimiter **inside** the bracket, each repetition carries its own, so the suffix is genuinely optional and genuinely repeatable — this is the spelling for a chain of appended suffixes.

With the delimiter **outside**, the value cannot cross a delimiter, so the group collapses to at most one suffix. Zero repetitions would leave the bare `parent_`, which is never a real id, so in practice this spelling selects exactly the ids that are `parent_` followed by one value — which is what a flat renumber (`name_1`, `name_2`, …) wants. `[^_]+` means "one or more characters that are not `_`"; `[^_]*` is the same with zero allowed.

Bounding the value at the delimiter is also what keeps selection exact: `design_1[_<?>]` cannot reach into `design_10_5`.

### IdSet — the type that owns these rules

An id is a string and an id set was a `List[str]`, so every rule above lived as a free function in `id_patterns` and had to be remembered at each call site. `IdSet` (`biopipelines/idset.py`) gives them one home: a function that takes an `IdSet` cannot be handed a bare list that skipped the rule.

It is a pure value type — ids and nothing else, no map_table and no file handles — so it is comparable, hashable and testable on its own, and equality is over the ids alone.

**Reading it.** `len(ids)` is how many ids are *declared*, so a pattern counts once; `ids.count()` is how many they expand to, counting a lazy bracket's prefix only. Those differ deliberately: at configuration time an expansion is often not knowable, and a count you cannot compute is worse than the one you can. `ids.enumerated()` gives every id the set stands for — all of them when the set is expandable, and ids that still carry their unexpandable parts when it is not. `ids.select(rows)` picks the map_table rows the patterns cover.

**The four operations**, which are the only ways ids are combined:

| Operation | Call | Example |
|---|---|---|
| bundling | `ids.bundled()` | `['l1','l2']` → `['l1+l2']` |
| cartesian product | `a.product(b)` | `['p1','p2'] × ['l1','l2']` → `['p1+l1','p1+l2','p2+l1','p2+l2']` |
| suffix multiplication | `ids.multiplied_by_suffix('<1..3>')` | `['5HG6_<0..4>']` → `['5HG6_<0..4>_<1..3>']` |
| renaming | `ids.renamed(mapping)` | cardinality-preserving; a rename that would collapse two ids raises |

`+` composes axes and `_` separates a parent from its child, which is what lets `prot1+lig1` be told apart from a suffix pattern. A suffix carrying a `+` is refused for that reason.

**Composing axes.** `compose_axes([(ids, mode), ...])` is the axis-aware entry point, and it is what tools should use: a `"bundle"` axis contributes one `+`-joined prefix however many ids it holds, a `"each"` axis contributes one id per row. **Bundle prefixes come first, ahead of the iterated axes and regardless of declared order.** That convention has a consequence worth internalizing: an id cannot be decomposed back by splitting on `+` and pairing positions with declared axes, because a bundle occupies as many positions as it has members and not the position it was declared in. Record what each axis contributed while composing instead. Getting this wrong silently mis-assigns provenance columns rather than failing.

`product` on its own is the plain cartesian product and respects declared order; it does not know about modes, because `bundled()` has already collapsed a bundled axis to an ordinary one-element set. The convention lives in `compose_axes`, not in `product`.

#### Two bracket vocabularies, and why they must stay apart

`<...>` and `[...]` answer different questions, and a third notation in `id_map_utils` answers a third. Confusing them is silent, so they are spelled so they cannot be.

| Notation | Module | Means | Example |
|---|---|---|---|
| `<...>` | `id_patterns` | a finite set of **literal values**, compacting an id that is already **predictable** | `<0..2>` is three ids; `<A V>` is two; `<N>` is the single id `N` |
| `[...]` | `id_patterns` | an **unpredictable** value, not knowable until the upstream tool has run | `prot[_<N>]` covers `prot`, `prot_4E`, `prot_4E_6U` |
| `<#>` `<?>` | `id_map_utils` | a **match class** for walking up a parent chain: digits, or one segment | `{"*": "*_<?>"}` strips `protein_1_19A` to `protein` |

**`<N>` is the letter N, not a number.** These brackets exist to compact ids that are already known, so a slot is a set of values and `expand_pattern('prot_<N>')` is `['prot_N']`. Write `<0..9>` or `<1 2 3>` when you mean digits.

**A `[...]` bracket cannot be *enumerated*, but it can be *constrained*.** Its values arrive at runtime, so nothing there expands at configuration time — yet a slot can still say what shape the value will take, and saying so is what stops `design_1[_<#>]` from covering `design_1_ZZZ`:

| Slot inside a bracket | Matches | Example |
|---|---|---|
| `<#>` | one or more digits | `design_1[_<#>]` covers `design_1_99`, not `design_1_1A` |
| `<?>` | any one segment | `design_1[_<?>]` covers both |
| `<A V>`, `<0..9>` | exactly those values | `prot[_<#><A V>]` covers `prot_5A`, not `prot_5X` |
| a single bare word | **that literal text** | `[<hits>]` covers `hits`; `[_<N>]` covers `_N` |

`<#>` and `<?>` are the only two classes. Everything else inside a slot is a value set, including a single bare word — so `[<hits>]` matches the literal `hits` and nothing else, and `[_<N>]` matches `_N`, not a number. Both spellings existed before the classes did and used to fall through to a wildcard; `rcsb.py` was migrated to `[<?>]` in RCSB 1.4 for exactly this reason. If a pattern of yours reads like a placeholder, it is now a literal.

**Adjacent classes are rewritten, not refused.** `[<#><?>]` compiles to `(?:\d[^_]+)?` rather than the naive `(?:\d+[^_]+)*`. Both quantifiers would otherwise compete for the same characters — `<?>` is `[^_]+`, which includes digits — and rejecting an id would cost one attempt per ordered partition of its run, 395 ms at 22 characters and doubling with each one. The rewrite matches exactly the same strings: a piece with a strictly wider neighbour gives up its quantifier because that neighbour absorbs what it left, one quantifier covers a run of equal-width pieces, and the outer `*` becomes `?` because two iterations are then always subsumed by one. `[<#><#>]` still means two or more digits, which no single class can say.

**`id_map_utils` uses punctuation, because a letter cannot survive here.** Its classes used to be `<N>` for digits and `<S>` for a segment, which read as the exact opposite of `id_patterns` and are ambiguous besides: `<N><S E>` reads as serine or glutamate at position N, while `<15 16><N>` reads as asparagine at positions 15 and 16. The classes are now `<#>` and `<?>`. **The letter spellings no longer do anything**: `{"*": "*_<S>"}` leaves `protein_1_19A` unstripped rather than walking up its parent chain, and it does so silently — there is no deprecation warning on this path. Rewrite any stored id_map to `<#>` or `<?>`.

The delimiter and the "one segment" rule have a single definition, `id_patterns.SUFFIX_DELIMITER` and `id_patterns.SEGMENT_REGEX`, which `id_map_utils` imports — so changing what separates a parent from a child moves both modules at once.

### The Rule

**Never call `ids_expanded` or `files_expanded` at configuration time to generate per-ID bash commands.** This breaks when the DataStream has lazy patterns, because the full ID list doesn't exist yet.

Instead:
- **Serialize the DataStream** to JSON via `save_json()` at config time
- **Expand IDs at runtime** inside the generated bash script using `Resolve.stream_ids()` or inside pipe scripts using `load_datastream()` + `ids_expanded` / `iterate_files()`

`DataStream.records()` is the exception by intent: it is a notebook/result
inspection helper for streams whose outputs have already been materialized. It
returns `StreamRecord` objects with `record.id`, `record.file` for file-backed
streams, and requested map-table columns via `record.values` or dot access for
value-backed streams:

```python
for record in output.streams.structures.records():
    print(record.id, record.file)

for record in ligand.streams.compounds.records(columns=["smiles"]):
    print(record.id, record.smiles)
```

Do not use `records()` from tool wrappers to generate scripts or predict
outputs. Tool wrappers still serialize streams and defer expansion to generated
bash or pipe scripts as described below.

### How to Iterate Structures in Generated Bash

Use `Resolve.stream_ids()` to get a bash expression that prints all expanded IDs at runtime, and `resolve_stream_item` (sourced by `activate_environment()`) to resolve each ID to its file path:

```python
from .biopipelines_io import Resolve

def generate_script(self, script_path):
    self.structures_stream.save_json(self.structures_json)

    script = "#!/bin/bash\n"
    script += self.activate_environment()  # sources resolve_stream_item.sh
    script += f"""
for struct_id in {Resolve.stream_ids(self.structures_json)}; do
    PDB_FILE=$(resolve_stream_item "{self.structures_json}" "$struct_id")
    echo "Processing $struct_id: $PDB_FILE"
    python run.py --pdb "$PDB_FILE"
done
"""
    return script
```

This generates bash like:
```bash
for struct_id in $(python "/path/to/resolve_stream_ids.py" "/path/to/structures.json"); do
    PDB_FILE=$(resolve_stream_item "/path/to/structures.json" "$struct_id")
    ...
done
```

The loop works correctly whether the DataStream has literal IDs, deterministic patterns, or lazy patterns.

`stream_ids(ds_json, index=None, valid_set=True)` takes a third argument worth knowing about. With `valid_set=True` (the default) the printed ids are restricted to those whose file is actually present on disk, which is what you want for a file-backed stream: a filtered `Pool`/`Panda` declares every original id but materializes only the survivors, so looping over the declared list would hand the tool paths that do not exist. The restriction is skipped automatically for a stream that declares no files at all (a value-based `compounds` stream, say), where file presence says nothing — so the default is safe there too and still prints every declared id. Pass `valid_set=False` only when you deliberately want the declared ids of a file-backed stream, including ones whose files are absent.

#### Using `eval` for pre-formatted argument variables

When a pipe script outputs a pre-formatted command-line fragment with embedded quotes (e.g. `--fixed_residues "A10 A11 A12"`), storing it in a bash variable and expanding it with `$VAR` does **not** interpret the embedded quotes — bash treats them as literal characters and word-splits on spaces. Use `eval` so the shell re-parses the line and respects the embedded quotes:

```python
script += f"""
    OPTIONS=$(python helper.py "{self.config_json}" "$struct_id")
    eval python run.py --pdb '"$PDB_FILE"' $OPTIONS
"""
```

Note that `$PDB_FILE` is wrapped in `'"..."'` (single-quoted double quotes) so `eval` preserves the quoting around the variable expansion. Only use `eval` when you need embedded quotes interpreted; for simple arguments, plain `python run.py --flag "$VAR"` is preferred.

#### Per-ID table values inside the loop (the Colab-safe pattern)

A common need is a *different* value per id read from an upstream table — e.g. one contig string per input PDB, supplied as a `(TableInfo, "col")` column reference. There is a `Resolve.table_column(ref, "$ID")` helper that emits an inline lookup, but it wraps the call in `<env_manager> run -n biopipelines`, which **does not exist on Colab** (there tools run in base Python, with no `biopipelines` conda env) and spawns a fresh Python process *per id*. Use it only for a genuine one-off single lookup.

For a per-id value across a loop, follow the same shape every other tool uses: a **pipe script resolves all ids in one pass** (`load_table` + `lookup_table_value` from `biopipelines_io`), writes a small `{id: value}` JSON, and the bash loop reads each id's value from that JSON via a tiny reader script invoked as **plain `python`**. Plain `python` runs under the already-activated tool env, where the pipe script's `sys.path.insert(0, repo_root)` makes `biopipelines` importable — so it works identically on cluster and Colab, with no second env and no per-id process spawn. The RFdiffusion family is the reference (`pipe_rfdiffusion_contigs.py` builds the JSON, `resolve_rfdiffusion_contigs.py` reads one id). The general rule: **never wrap a helper in `<mgr> run -n biopipelines` inside generated bash** — it breaks Colab; resolve via a pipe script under the activated env instead.

### How to Resolve a Single File

When a tool processes one fixed input (not iterating over all IDs), use `Resolve.stream_item()`:

```python
# At config time — produces a bash expression, not the actual path
first_id = self.structures_stream.ids[0]
script += f'INPUT_PDB={Resolve.stream_item(self.structures_json, first_id)}\n'
```

This generates:
```bash
INPUT_PDB=$(resolve_stream_item "/path/structures.json" "prot_1")
```

**Warning:** `ids[0]` is only safe when the stream has deterministic IDs (literal or patterns like `<1..5>`, which expand on indexing). If the stream has lazy IDs (e.g. `<N>`), `ids[0]` returns the unexpanded pattern, which may represent multiple structures. For lazy streams, resolve all IDs at runtime and take the first:

```python
script += f'FIRST_ID={Resolve.stream_ids(self.structures_json, index=0)}\n'
script += f'INPUT_PDB={Resolve.stream_item(self.structures_json, "$FIRST_ID")}\n'
```

### Passing IDs to pipe scripts

Don't build per-ID data structures at config time. Instead, pass the DataStream JSON path to pipe scripts and let them expand IDs at runtime:

```python
# BAD — breaks with lazy IDs
id_map = {}
for sid, path in zip(ds.ids_expanded, ds.files_expanded):
    id_map[os.path.basename(path)] = sid

# GOOD — pipe_script builds the map at runtime
script += f'python {self.helper_py} --ds-json "{self.structures_json}"\n'
```

In the pipe script:
```python
from biopipelines.biopipelines_io import load_datastream, iterate_files

ds = load_datastream(sys.argv[1])
for struct_id, pdb_path in iterate_files(ds):
    process(struct_id, pdb_path)
```

### Filtered streams vs. the raw map_table

A stream's `map_table` is the full **warehouse** of rows an upstream tool produced. A stream's `ids` are a **view** over it. Filtering a stream (`output["CP1"]`, slicing, `filter_by_ids()`) narrows the `ids` but leaves the `map_table` untouched — the filtered-out rows are still physically in the CSV.

So the hazard is: **a consumer that reads the raw `map_table` directly (`cp map_table`, `pd.read_csv(map_table)`, `os.listdir` of an upstream folder) gets the whole warehouse and silently ignores the filter** — it reprocesses every row, not just the selected ids.

The rule has one branch:

- **If the consumer resolves through the DataStream at runtime** — `load_datastream()` + `iterate_files()` / `iterate_values()` / `Resolve.stream_ids()` — it iterates `ids_expanded`, so it is **filtered by construction**. This is the default transport; prefer it whenever the runtime script can consume a DataStream.
- **If and only if the consumer bypasses datastream resolution** — its binary needs a concrete CSV / FASTA / native file and reads the table directly — it **must materialize the filtered table first**. Never feed the raw `map_table` to such a consumer.

Materialize the filtered CSV with the shared `BaseConfig` block, emitted *before* the consumer step:

```python
# In generate_script(), before activating the tool env:
script += self.generate_filtered_map_table_block(
    self.sequences_json,            # stream JSON saved at config time
    self.filtered_sequences_csv,    # runtime-materialized, filtered CSV
    required_columns=["id", "sequence"],
)
# ... then activate the tool env and pass self.filtered_sequences_csv to the binary
```

The block runs `materialize_filtered_map_table.py` (which calls `write_filtered_map_table()`) under the **biopipelines** env — it imports `biopipelines`, so it cannot run in the tool env. It projects the `map_table` to `ids_expanded` in stream order; an unfiltered stream gets a byte-for-byte copy (fast path). `write_filtered_map_table` is representation-agnostic (ids + columns only) — any FASTA/format conversion the binary needs is derived from this CSV in the tool's own pipe script, not in `biopipelines_io`.

Combinatorics consumers (the unified Boltz/AlphaFold config builder) are the same rule seen from the other side: they don't use datastream resolution either, so they carry the filter as `ids` in the combinatorics config and the runtime loader intersects the source rows with them. Either way, *some* explicit filter step is applied — the only consumer that needs none is the datastream-resolving one.

### Summary Table

| Operation | Config time | Execution time |
|-----------|------------|----------------|
| Store DataStream | `ds.save_json(path)` | — |
| Load DataStream | — | `load_datastream(path)` |
| Get all IDs | `ds.ids` (compact patterns) | `ds.ids_expanded` (full list) |
| Iterate in bash | `Resolve.stream_ids(json)` | — |
| Resolve one file in bash | `Resolve.stream_item(json, id)` | — |
| Iterate in Python | — | `iterate_files(ds)` |
| Resolve one file in Python | — | `resolve_file(ds, id)` |
| Build per-ID data | **Don't** | Do it in pipe scripts |
| Feed a raw CSV to a binary | `generate_filtered_map_table_block(...)` | materializes filtered CSV |

---

## Tool Development

### Creating a Tool

1. Create `biopipelines/my_tool.py`
2. Inherit from `BaseConfig`
3. Implement required methods
4. Create helper script `pipe_scripts/pipe_my_tool.py` if needed

```python
"""MyTool - brief description."""

import os
from typing import Dict, List, Any, Union

try:
    from .base_config import BaseConfig, StandardizedOutput, TableInfo
    from .file_paths import Path
    from .datastream import DataStream
    from .biopipelines_io import Resolve
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo
    from file_paths import Path
    from datastream import DataStream
    from biopipelines_io import Resolve
```

> **Why the dual import?** Tools are normally imported as part of the `biopipelines` package (relative imports via `.base_config`, etc.). The `except ImportError` fallback adds the module's own directory to `sys.path` so the same file can also be imported standalone — useful for debugging or running a tool file directly. Always include this pattern in new tool files.

```python
class MyTool(BaseConfig):
    TOOL_NAME = "MyTool"
    TOOL_VERSION = "1.0"
    # The env this tool builds and runs in. Nothing mechanical catches its
    # absence; the tool just fails to install. See "Environments and ENV_NAME".
    ENV_NAME = "mytool"

    # Path descriptors — route each artefact into its canonical sub-folder.
    # Never use ``self.output_folder`` directly; go through:
    #   self.configuration_path(*parts) — config-time input JSONs/YAMLs
    #   self.execution_path(*parts)     — raw runtime dumps
    #   self.stream_folder(name)        — per-stream files + map_table
    #   self.stream_map_path(name)      — <stream>/<name>_map.csv
    #   self.table_path(name)           — tables/<name>.csv (standalone TableInfo)
    #   self.extras_path(*parts)        — catch-all for ancillary files
    results_csv = Path(lambda self: self.table_path("results"))
    structures_json = Path(lambda self: self.configuration_path(".input_structures.json"))
    helper_py = Path(lambda self: os.path.join(self.folders["pipe_scripts"], "pipe_my_tool.py"))

    def __init__(self,
                 structures: Union[DataStream, StandardizedOutput],
                 param1: str,
                 param2: int = 10,
                 **kwargs):
        # Resolve input to DataStream
        if isinstance(structures, StandardizedOutput):
            self.structures_stream = structures.streams.structures
        elif isinstance(structures, DataStream):
            self.structures_stream = structures
        else:
            raise ValueError(f"structures must be DataStream or StandardizedOutput")

        self.param1 = param1
        self.param2 = param2

        super().__init__(**kwargs)

    def validate_params(self):
        if not self.structures_stream or len(self.structures_stream) == 0:
            raise ValueError("structures cannot be empty")
        if self.param2 <= 0:
            raise ValueError("param2 must be positive")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        self.folders = pipeline_folders

    def generate_script(self, script_path: str) -> str:
        # Serialize DataStream for runtime access
        self.structures_stream.save_json(self.structures_json)

        script_content = "#!/bin/bash\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()
        script_content += f"""
echo "Running MyTool"
python {self.helper_py} \\
    --ds-json "{self.structures_json}" \\
    --param1 "{self.param1}" \\
    --param2 {self.param2} \\
    --output "{self.results_csv}"
"""
        script_content += self.generate_completion_check_footer()
        return script_content

    def get_output_files(self) -> Dict[str, Any]:
        tables = {
            "results": TableInfo(
                name="results",
                path=self.results_csv,
                columns=["id", "param1", "param2", "result"],
                description="MyTool results"
            )
        }

        return {
            "structures": DataStream.empty("structures", "pdb"),
            "sequences": DataStream.empty("sequences", "fasta"),
            "compounds": DataStream.empty("compounds", "csv"),  # always "csv" — see the Ligand Contract
            "tables": tables,
            "output_folder": self.output_folder
        }
```

### Required Methods

| Method | Purpose |
|--------|---------|
| `validate_params()` | Validate parameters, raise ValueError on failure |
| `configure_inputs(pipeline_folders)` | Set `self.folders`, configure input sources |
| `generate_script(script_path)` | Return bash script as string |
| `get_output_files()` | Return dict with DataStreams and tables |

### Renaming a Parameter: `PARAMETER_ALIASES`

Every tool constructor ends in `**kwargs`, so an unknown key is silently accepted and only the leftover-kwargs contract notices it. That makes renaming a parameter dangerous in a specific way: a pipeline written against the old spelling keeps running, the value never binds, and the tool does whatever it does with the parameter at its default. Four near-misses between our own sibling tools all failed exactly like that, which is why `BaseConfig` grew two class attributes:

```python
class RFdiffusion3(BaseConfig):
    TOOL_NAME = "RFdiffusion3"
    # RFdiffusion, RFdiffusion2 and RFdiffusionAllAtom all spell this `contigs`.
    PARAMETER_ALIASES = {"contigs": "contig"}

class Ligand(BaseConfig):
    TOOL_NAME = "Ligand"
    PARAMETER_ALIASES = {"code": "codes"}
    DEPRECATED_ALIASES = ("code",)
```

- **`PARAMETER_ALIASES`** is `{alternative spelling: current parameter name}`. Declare a renamed parameter's old spelling here so a pipeline written against it still binds the new parameter, and declare a synonym a user reasonably expects the same way — `RFdiffusion3`'s `contigs` is the plural its three sibling wrappers all take, so the plural has to bind rather than be forwarded.
- **`DEPRECATED_ALIASES`** is the subset of those keys that are a *retired* spelling rather than a first-class synonym. Only these report a `deprecated_alias` contract line, so a synonym stays silent instead of warning on every legitimate use.

**Aliases resolve before the constructor binds.** `__init_subclass__` wraps the tool's `__init__` in `_alias_resolving_init`, which rewrites the kwargs dict on the way in. That ordering is the whole point: an aliased key is never a leftover kwarg, so it is neither reported as a probable typo nor — on a tool with [`FORWARD_UNKNOWN_KWARGS`](#forwarding-untyped-kwargs-forward_unknown_kwargs) — rendered onto the upstream command line. `RFdiffusion3(contigs=…)` used to do both: reported as a misspelling, then passed to hydra as a bogus override, leaving `self.contig` empty while the job ran to completion producing unconditioned de-novo backbones.

**Two guards fire at class creation**, not at call time, because an alias that does not work is exactly the silence the mechanism exists to remove:

1. An alias whose *target* no constructor parameter names — the alias would land back in `**kwargs` and change nothing. `ValueError`, naming the unbound targets.
2. A `DEPRECATED_ALIASES` entry absent from `PARAMETER_ALIASES` — a deprecation for a spelling that binds nothing. `ValueError`.

**A reserved key is copied, not moved.** For a key in `contract_enforcement.RESERVED_KWARGS` — in practice `name`, the job name — the value is written under the tool's new parameter name *and left in place*, so the tool and the framework both see it. That is what let `Table` keep `name=` as a synonym for `table_name=` while the job name finally reaches `BaseConfig` again. For any other key, passing both spellings raises rather than picking one:

```
ValueError: Ligand: code= is another spelling of codes=; pass one or the other, not both.
```

And for a reserved key, passing both is not a clash: `Table(path, table_name="metrics", name="jobname")` names the table `metrics` and the job `jobname`.

**When to add an alias, and when not to.** A permanent synonym is worth it when the value is a *handle other code reaches for* — `Table`'s `name` becomes `tool.tables.<name>`, so silently redirecting it breaks a consumer further down the pipeline. A clean break is better when the old spelling only affects this tool's own output names — `Fuse`'s label became `prefix` with **no** alias, so `name` on `Fuse` means the job name exactly as it does everywhere else. Whole-tool aliases (`Structure = PDB`, `Compound = Ligand`) are a different thing entirely: a bare module-level assignment, no `PARAMETER_ALIASES` involved, and `TOOL_NAME` unchanged so the registry, the docs index and the config's `environments:` / `folders:` keys stay single-entry per tool.

### Forwarding Untyped Kwargs: `FORWARD_UNKNOWN_KWARGS`

`TOOL_NAME`, `TOOL_VERSION`, `ENV_NAME` and `USER_STREAM_NAMES` say what a tool *is*. `FORWARD_UNKNOWN_KWARGS` says something narrower and rarer: that this wrapper hands a constructor keyword it does not recognise straight to the program it runs.

```python
class RFdiffusion(BaseConfig):
    TOOL_NAME = "RFdiffusion"
    FORWARD_UNKNOWN_KWARGS = "hydra"      # "" (default, off) | "argparse" | "hydra"
```

**What it is for.** Some tools are a wrapper around an inference script, a binary, or a container that takes its parameters in free form, and the framework has not typed every one of those arguments — there are hundreds, they change between upstream releases, and most are never used. Without the marker, `RFdiffusion(denoiser_noise_scale=0.5)` lands in `**kwargs`, is parked in `self.params`, and is never read: the pipeline runs to completion with the option at its default. With it, the leftover keys are rendered onto the upstream command line, so an advanced user can reach an option the wrapper never curated instead of waiting for someone to type it.

**The dialect** is how the upstream program spells an option, and there are two:

| Value | Renders `omit_AAs="CX"` as | Use when the upstream parser is |
|---|---|---|
| `"argparse"` | `--omit_AAs CX` | Python `argparse`/`click`, or any `--flag value` binary |
| `"hydra"` | `omit_AAs=CX` | hydra / OmegaConf, taking dotted `key=value` overrides |

Read it off the upstream command, not off the language: `python inference.py --num_seq 8` is argparse, `python run_inference.py inference.num_designs=8` is hydra. Under `argparse` a `True` value becomes a bare flag (`--score_only`) and a list repeats the flag; under `hydra` `True` has to be spelled out (`inference.cyclic=True`) because a bare token is not a valid override, and a list becomes one bracketed override (`ppi.hotspot_res=[A30,A33]`) because a repeated hydra override keeps only its last value. `False` and `None` render as nothing at all under both — "off" is the upstream default the wrapper never touched, and emitting a flag would change it.

**The architectural fact that decides whether the marker works at all.** Most wrappers in this tree do not build an upstream command line. They serialize their parameters into a config JSON and emit one line of bash:

```bash
python pipe_mytool.py --config /path/to/_configuration/mytool_config.json
```

The *pipe script* then assembles argv on the compute node. For such a wrapper there is no wrapper-written command for a forwarded token to join, so the marker is accepted and the tokens are dropped — the exact silence it exists to remove. **Only a wrapper whose own generated bash contains the upstream command can take the marker as a one-liner.** Thirteen of the ninety tools do today: the RFdiffusion family, the MPNN family, AlphaFold, Boltz2, NeuralPLexer, LASErMPNN, Aggrescan3D, CABSflex and DynamicBind. Before setting it, open your `generate_script()` and look for the upstream program's name. If you find it, interpolate the tokens into that command; if you find only `pipe_<tool>.py --config`, you have two honest options — give the pipe script a typed parameter instead, or move the upstream command into the wrapper first.

Rendering is one call, in whichever shape the command needs:

```python
script += self.extra_args_echo()                       # announce them in the step log
...
f"python {self.inference_py} {options} {self.extra_args_bash()}"   # one bash fragment
opts += self.extra_args_tokens()                       # raw argv tokens, to join yourself
self.extra_args_bash_tokens()                          # each token double-quoted
```

`extra_args_echo()` alone is not forwarding — it prints the tokens and nothing else — and `contract_enforcement.assert_can_forward` refuses a marker on a class whose wrapper source never calls one of the three rendering accessors. That check runs at class creation, so a marker that cannot work is an import error naming the tool and the fix, rather than a flag that is accepted and then silently ignored. It is a necessary condition, not a sufficient one; `tests/test_forwarded_kwargs.py` discovers every tool carrying the marker and requires the token to appear on a real command line in the generated `NNN_<Tool>.sh`, with `echo` lines excluded so the announcement cannot answer for the command. A tool that gains the marker is tested by having gained it — there is no list to remember to edit.

**The trade you are making.** Forwarding turns a configuration-time report into a runtime failure. A non-forwarding tool answers `MyTool(num_desgins=8)` at the user's desk with `[contract:unknown_kwargs] … did you mean 'num_designs'?`; a forwarding tool renders `--num_desgins 8` as written, and the complaint arrives from upstream, on a compute node, after the queue. Only the typo check survives the switch: `check_probable_typo` still names a key within the similarity cutoff of a real parameter name, because a near-miss is almost certainly a misspelling either way — but a key that resembles nothing is passed through unvalidated, since the wrapper has no model of the upstream parser to check it against. So the marker is a good trade on a large, stable option surface where the curated subset will always be a fraction of what upstream accepts, and a bad one on a tool with two or three flags, where typing the third flag is less work than the failure mode.

**Shell safety.** Every forwarded key and value goes through `_validate_extra_arg` at construction. Keys are checked against `_SAFE_EXTRA_KEY_RE` in `base_config.py`, so a key that is not a plain Python identifier — a dotted hydra override such as `denoiser.noise_scale_ca` — reaches the constructor only as `Tool(**{"denoiser.noise_scale_ca": 0.5})`; check that regex before assuming a given upstream flag has a keyword spelling at all. Values refuse the four characters that break double-quoted interpolation — `"`, `` ` ``, `$`, backslash — and also `;`, `|`, `&`, `<`, `>`, `(`, `)`, CR and LF. That second list is wider than the freeform denylist elsewhere in the framework for a specific reason: LigandMPNN emits its command through `eval`, which re-parses the line and strips one layer of quoting, so the double quotes the renderer puts around a forwarded value do not survive to contain a `;`. See [Shell Safety](#shell-safety); this is the same fourth layer, with the `eval` case accounted for.

### Path Descriptors and the Canonical Layout

Every tool's ``output_folder`` has a predictable sub-layout that the framework creates automatically — but only *after* ``get_output_files()`` returns, via ``_materialize_output_layout()``. **Never invent a directory name or ``mkdir`` a path you derived yourself**: route every output path through one of the helpers below, and in the generated bash rely on the layout the framework has already created rather than adding your own ``mkdir -p``.

Two places do have to create their directory themselves, because they run before that layout exists:

- **``configure_inputs()``**, which the pipeline calls *before* ``_materialize_output_layout()``. A tool that writes a config-time input file there must ``os.makedirs(os.path.dirname(path), exist_ok=True)`` first — `mmseqs2.py`, `sequence.py` and `table.py` all do. The path itself still comes from a helper.
- **A pipe script writing per-stream files**, which may run before its stream folder exists; ``os.makedirs(self.stream_folder(name), exist_ok=True)`` before the first write is correct. (``create_map_table()`` already does this for the map itself.)

Note that ``BaseConfig.stream_folder``'s own docstring states the folder "is not created here — the tool (or its pipe script) is responsible for calling ``os.makedirs``", which reads as a broader licence than the rule above. Follow the rule: create nothing you can get from a helper, and makedirs only in the two cases listed.

Instead of hand-built paths, use these helpers:

| Artefact | Helper | Resolves to |
|---|---|---|
| Config-time inputs (JSONs, YAMLs passed to the CLI) | `self.configuration_path(*parts)` | `<output_folder>/_configuration/<parts>` |
| Raw model dumps (Boltz's `boltz_results_*`, etc.) | `self.execution_path(*parts)` | `<output_folder>/_execution/<parts>` |
| Per-stream files (`<id>.pdb`, …) | `self.stream_folder(name)` | `<output_folder>/<name>/` |
| Stream's map_table CSV | `self.stream_map_path(name)` | `<output_folder>/<name>/<name>_map.csv` |
| Standalone TableInfo CSV | `self.table_path(name)` | `<output_folder>/tables/<name>.csv` |
| Ancillary files (session files, info dumps) | `self.extras_path(*parts)` | `<output_folder>/_extras/<parts>` |

Use `Path` descriptors for lazy path evaluation:

```python
from .file_paths import Path

class MyTool(BaseConfig):
    # Evaluated on first access, after output_folder is set.
    structures_json = Path(lambda self: self.configuration_path(".input_structures.json"))
    results_csv = Path(lambda self: self.table_path("results"))
    helper_py = Path(lambda self: os.path.join(self.folders["pipe_scripts"], "pipe_my_tool.py"))
```

**Content-bearing streams.** When a stream's map_table IS its content
table (e.g. Sequence's `sequences.csv` holds `id, sequence` and doubles
as the map), point both the `DataStream.map_table` and the `TableInfo.path`
at `<stream>/<stream>.csv` — one file, one source of truth, no duplicate
CSV. `stream_map_path(name)` gives `<stream>_map.csv`, which is the
right default for lineage-only streams (most producers) but NOT for
content-bearing ones.

### Script Generation

Use provided helpers and runtime resolution:

```python
def generate_script(self, script_path: str) -> str:
    # Save DataStream for runtime access
    self.structures_stream.save_json(self.structures_json)

    script_content = "#!/bin/bash\n"
    script_content += self.generate_completion_check_header()  # Skip if completed
    script_content += self.activate_environment()              # Activate conda + source resolve_stream_item.sh

    # Option A: Pass DataStream to pipe_script (preferred for complex processing)
    script_content += f"""
python {self.helper_py} --ds-json "{self.structures_json}" --output "{self.results_csv}"
"""

    # Option B: Bash for-loop (for simple per-structure commands)
    script_content += f"""
for struct_id in {Resolve.stream_ids(self.structures_json)}; do
    PDB_FILE=$(resolve_stream_item "{self.structures_json}" "$struct_id")
    python external_tool.py --input "$PDB_FILE" --output "{self.output_folder}/$struct_id.out"
done
"""
    script_content += self.generate_completion_check_footer()  # Check outputs
    return script_content
```

### Containerization (develop with it in mind)

Wherever possible, write a tool so it can switch between **environment execution** (the binary runs in an activated conda env) and **container execution** (the binary runs inside an Apptainer/Singularity `.sif`) without any code change — the choice is purely config-driven. The seam is `BaseConfig.container_prefix()`: it returns `"<executor> exec --nv -B <binds> <image> "` when a `containers: <ToolName>: <image>` entry is configured, and `""` otherwise. A tool opts in by prefixing the line that runs the heavy binary:

```python
# Direct-in-bash heavy command (preferred shape):
script += f'{self.container_prefix()}python "{self.inference_script}" --in "{x}" ...'
```

When a container is configured the command runs in the `.sif`; when not, the prefix is empty and it runs in the activated host env. The `.sif` only needs the tool's own binary/model — not biopipelines — so keep host-side orchestration (stream IO, CSV assembly, pre/post helpers that `import biopipelines`) **outside** the prefix, in the host env.

**Activation still happens in container mode.** `container_prefix()` wraps only the binary command — it does not replace `activate_environment()`. A container-mode script both activates the host env and prefixes the binary: `activate host env → {container_prefix}binary`. This is deliberate, because the surrounding helper/orchestration Python (which imports biopipelines, pandas, BioPython) runs in the activated host env while only the binary runs in the image. For a tool that runs fully in a container, set its env-map entry to a lightweight host env that just satisfies the helpers — e.g. `RFdiffusion2: "biopipelines"` and `PLIP: "biopipelines"` rather than the heavy model env, since inference happens in the `.sif`.

**Two architectures:**

- **Heavy command emitted directly in bash** (RFdiffusion, ESMFold, NeuralPLexer, DiffDock): the wrapper's `generate_script` writes the model/binary invocation as its own bash line. Prefix that line — done. This is the shape to aim for in new tools.
- **Binary invoked inside a host helper** (`python "{self.helper_py}"` where the helper imports biopipelines and `subprocess.run`s the binary): you cannot prefix the `python "{helper}"` line — that would force the biopipelines-importing orchestrator into the image. Instead pass the prefix down to the helper (a `--container-prefix "{self.container_prefix()}"` CLI arg or a config-JSON key) and have the helper prepend it to the binary command with `container_argv_prefix()` from `biopipelines_io`:

```python
from biopipelines.biopipelines_io import container_argv_prefix
cmd = container_argv_prefix(args.container_prefix) + ["mytool", "-i", pdb, ...]
subprocess.run(cmd, ...)
```

If the helper resolves the binary off the host PATH (`shutil.which`) or sets host-env-specific variables (e.g. `BABEL_LIBDIR`), skip that resolution when a prefix is present — inside the image the binary is on the container PATH and the host paths/vars are wrong.

**When containerization is impractical:** a tool that runs its model **in-process** (`import torch; model.run(...)`) in the same script that imports biopipelines, or that dispatches multiple host-env sub-Pythons (each resolved from a conda env path), has no single binary boundary to wrap. Either split it into a standalone, biopipelines-free model-runner the wrapper calls with `{self.container_prefix()}` directly (the ESMFold two-script pattern), or — if that refactor isn't worth it — emit a warning when `self.uses_container()` is true stating container execution isn't supported for this tool, and run in env mode. BP-native pure-Python tools (Selection, Panda, Pool, …) run in the biopipelines env and never containerize.

**Per-iteration overhead (known limitation).** A helper that loops over a stream and prefixes the binary inside the loop pays one `apptainer exec` startup per iteration (SIF mount + namespace + `--nv` GPU probe, ~0.3–2 s each). For heavy models run a handful of times this is noise; for a light binary (mkdssp, fpocket, reduce, xtb) over hundreds of inputs it can dominate. The current helpers accept this cost — none has a container image wired up today, so it is hypothetical. If a light looping tool ever gets a real image and the overhead bites, the escape hatch is a **persistent container instance** rather than restructuring the host loop: `apptainer instance start <image> <name>` once before the loop, run each iteration as `apptainer exec instance://<name> <binary> ...` (namespaces already set up — cheap), and `apptainer instance stop <name>` after. This keeps the biopipelines-importing loop on the host while amortizing container startup across all iterations.

### Output Prediction

Return standardized output structure. Use compact `ids` patterns (not expanded) and `<id>` file templates:

#### File Templates with \<id\>

`<id>` is a placeholder in the `files` list that represents all output files at once. Instead of listing one file per ID, you provide a single-element list like `["<id>.pdb"]`. At expansion time, `<id>` is replaced with each expanded ID. For example:

- `files=["<id>.pdb"]` + `ids=["prot_<0..2>"]` → `prot_0.pdb`, `prot_1.pdb`, `prot_2.pdb`

`files` carries meaning in its **type**: a list is per-ID, a bare string is the *shared-file* form (one artifact covering every ID). So `files="<id>.pdb"` as a string is not a template — every ID resolves to that same literal path.

This prevents a length-mismatch validation error (1 file vs N ids) and is **required** when inputs may carry lazy IDs — building per-ID file paths with f-strings embeds bracket patterns into paths, causing `LazyPatternError`. The implementation lives in `datastream.py:_has_file_template()` (detects the pattern) and `id_patterns.py:expand_file_pattern()` (performs the substitution).

```python
def get_output_files(self) -> Dict[str, Any]:
    # Keep IDs compact — don't expand. Route per-ID files into the
    # structures/ stream folder; its map_table lives alongside them.
    structure_ids = self.structures_stream.ids
    structure_files = [self.stream_path("structures", "<id>.pdb")]
    structures_map = self.stream_map_path("structures")

    structures = DataStream(
        name="structures",
        ids=structure_ids,
        files=structure_files,
        map_table=structures_map,
        format="pdb"
    )

    tables = {
        "results": TableInfo(
            name="results",
            path=self.results_csv,
            columns=["id", "score"],
            description="Analysis results"
        )
    }

    return {
        "structures": structures,
        "sequences": DataStream.empty("sequences", "fasta"),
        "compounds": DataStream.empty("compounds", "csv"),  # always "csv" — see the Ligand Contract
        "tables": tables,
        "output_folder": self.output_folder
    }
```

`get_output_files()` returns a plain dict, **not** a `StandardizedOutput`. The wrapping into `StandardizedOutput` happens automatically in the `ToolOutput.output` property; both classes live in `biopipelines/outputs.py` and are re-exported from `base_config.py`, so either import path works. This is by design:

- Pipeline infrastructure (`pipeline.py`, `base_config.py:get_id_provenance()`) iterates the raw dict with `.items()` and `isinstance()` checks before any user accesses it — dicts are natural for this.
- `StandardizedOutput.__init__` takes a dict and destructures it into `.streams`, `.tables`, etc. — returning `StandardizedOutput` from tools would just add an object whose constructor immediately unpacks it back.
- Keeps tools simple: tool authors build plain dicts, the framework handles the user-facing API.

**Rule:** tool authors return dicts; consumers get `StandardizedOutput` with dot-notation (e.g., `tool.output.structures`).

### Map Table Contract

Every non-empty output DataStream must have a valid `map_table` CSV at runtime. Downstream tools rely on `load_datastream()` → `ids_expanded` / `iterate_files()` to discover what the upstream tool actually produced, and this reads from the map_table.

**Rule: `create_map_table()` is a runtime-only helper.** Never call it from `get_output_files()` (or anywhere else that runs at config time). `get_output_files()` is purely declarative: it returns a `DataStream` carrying the predicted ids, the `<id>`-templated file pattern, and the `map_table` path — nothing on disk. The pipe script the tool launches is the sole writer of that CSV, and it writes only the rows whose files actually exist. This keeps the map honest under partial failure and keeps `Pipeline.save()` O(#tools), independent of design count.

Declarative `get_output_files()` example:

```python
def get_output_files(self):
    ids = self.structures_stream.ids
    files = [self.stream_path("structures", "<id>.pdb")]
    structures = DataStream(
        name="structures",
        ids=ids,
        files=files,
        map_table=self.stream_map_path("structures"),
        format="pdb",
    )
    return {"structures": structures, "tables": {}, "output_folder": self.output_folder}
```

Runtime writer (in the pipe script) drops failed items and emits provenance columns when relevant:

```python
rows = []
for sid, in_path in iterate_files(ds):
    out_path = os.path.join(out_dir, f"{sid}.pdb")
    if run_tool(in_path, out_path):
        rows.append({"id": sid, "file": out_path})
pd.DataFrame(rows, columns=["id", "file"]).to_csv(args.map_csv, index=False)
```

File templates with `<id>` are the preferred form for the declared stream — compact, lazy-id friendly, and at runtime `iterate_files()` expands `<id>` against the ids the pipe script actually wrote into the map.

**Rule:** if a tool returns a DataStream with a `map_table` path, that CSV **must exist and be accurate** by the time the tool's script finishes. Downstream tools will read it.

**When config time is allowed.** The config-time ban exists because a normal tool's outputs do not exist yet — its ids are predictions, so writing a map at config time would record rows that have not been produced. The discriminator is therefore not the phase but the question *"does the data this map describes exist at the moment the map is written?"* For a producing tool the answer is no until its run finishes (→ runtime only). Where the answer is yes, or where a config-time map is only a placeholder the runtime pass overwrites, a config-time map is legitimate. Three cases in tree:

- **`Load`** imports data that **already exists on disk** when the pipeline is authored, so its ids are genuinely known at config time. Note that `Load` does not call `create_map_table` at all: it reads the upstream producer's map_table and reuses it, writing a `.rebased_<map filename>` copy beside it only when the run folder has moved and the recorded paths need rewriting.
- **`BoltzGen`** calls `create_map_table` at config time inside `get_output_files()` (twice — for its final-ranked and refolded structure streams) because the rank ids and their glob file patterns are fully determined by `budget`, not by the run.
- **`Mock`** does so behind its `map_table_strategy` flag (`"config"` or `"both"`), where a config-time map is deliberately a test fixture; for lazy outputs it writes only the deterministic prefix and the pipe script overwrites it at runtime.

Outside these, write the map from the pipe script.

`Load` resolves against the upstream **map_table**, not by globbing: the map lists one row per id the producer actually wrote, so it is the authoritative answer to what exists. Two things follow. First, a `files` entry containing `<id>` is a *template*, not a path — statting it always fails, so it must be expanded via the map before any existence check. Second, `Load` narrows the stream's ids to the map's rows, so a stream declaring 300 ids whose producer emitted 290 loads as 290 and downstream tools never inherit the 10 phantom ids. Globbing is the fallback only when a stream has no readable map_table, and that case is reported as *unresolved* rather than as missing files. With `validate_files=False` none of this runs — ids propagate as declared, and downstream tools must consume `missing.csv` themselves.

### Install Scripts and the `$INSTALL_SUCCESS` Contract

Every tool exposes a classmethod `_install_script(cls, folders, env_manager, force_reinstall, **kwargs)` that returns the bash for `Tool.install()`. The framework wraps this script in a `_Installer` step that exports a single env var, `$INSTALL_SUCCESS`, pointing at a sentinel file under the install step's output folder. The completion check then surfaces the install as **COMPLETED** or **FAILED** based on whether that file exists.

The contract is one-line: **the install script must `touch "$INSTALL_SUCCESS"` only after a real verification succeeds — typically an import test in the env it just created.** The framework deletes the sentinel before the script runs, so a stale file from a prior failed run cannot fake success.

For tools that create their own conda env, end the install with a verification block:

```python
return f"""echo "=== Installing MyTool ==="
{skip_block}
{env_install_block}

# Verify installation
if {env_manager} run -n MyToolEnv python -c "import mytool" >/dev/null 2>&1; then
    touch "$INSTALL_SUCCESS"
    echo "=== MyTool installation complete ==="
else
    echo "ERROR: MyTool verification failed (cannot import mytool)"
    exit 1
fi
"""
```

For tools that do nothing at install time (utility tools that piggyback on the `biopipelines` env, or tools that re-use another tool's env like `ProteinEnv` / `MutationEnv`), unconditionally touch the marker:

```python
return """echo "=== MyTool ==="
echo "Uses biopipelines environment (no additional installation needed)."
touch "$INSTALL_SUCCESS"
echo "=== MyTool ready ==="
"""
```

For the **skip path** (an existing-install detector at the top of the script that calls `exit 0`), also touch the marker before exiting — the env-existence check is the verification in that branch:

```bash
if {env_manager} env list 2>/dev/null | grep -q "MyToolEnv"; then
    echo "MyTool already installed, skipping. Use force_reinstall=True to reinstall."
    touch "$INSTALL_SUCCESS"
    exit 0
fi
```

**Pick a real verification.** Prefer an import that exercises the heaviest binary dependency (`import torch`, `import boltz`, `import pymol`) over a shallow check like "directory exists". For binary tools (GNINA), check that the executable is present and `-x`. For tools whose install only downloads weights, check the weight file exists.

### Environments and `ENV_NAME`

A tool that builds its own conda/mamba environment declares that environment's default name as a class attribute, beside `TOOL_NAME` and `TOOL_VERSION`:

```python
class MyTool(BaseConfig):
    TOOL_NAME = "MyTool"
    TOOL_VERSION = "1.0"
    ENV_NAME = "mytool"
```

**Nothing mechanical will remind you.** The pre-commit gate (`versions/check_tool_edits.py`) and `tests/test_registry_consistency.py` between them check the changelog, the tool index, the reference, the category doc, the README row and the `environments:` entries — none of them checks that `ENV_NAME` exists. A tool that omits it imports fine, passes the whole suite, and then fails to install. Declare it.

**Where the installed env name actually comes from.** `_install_script` must never spell its environment out; it calls `cls._install_env(env_manager)` and installs into whatever that returns. The resolution order in `BaseConfig._install_env` is:

1. An explicit `env_name=` argument from the caller, if one was passed.
2. The tool's `environments:` entry in the active `config.<variant>.yaml`, keyed by `TOOL_NAME`. A list-valued entry contributes its first element.
3. Under `env_manager: pip` only — `cls.ENV_NAME`, because a laptop variant declares no `environments:` block at all.
4. Otherwise **it raises**. There is no hardcoded fallback, and in particular no fallback to `ENV_NAME` outside pip mode.

The config is consulted **before** the pip fallback, not after. It used to be the other way round, which made pip the one manager where the wrapper's own literal beat an explicit `environments:` entry — so a laptop config that did name an environment was ignored.

The config entry is the one `_load_environments` activates at run time, so reading the same entry at install time is what keeps the two pointed at the same environment. An install script that hard-codes its own name builds one env while the pipeline activates another, and that mismatch does not surface until a job reaches a compute node and cannot find its interpreter.

**The naming standard: a spec file is named after the environment it defines.** An `environments:` entry names an *environment*, and an environment is defined by `environments/<env>.<variant>.yaml`, falling back to `environments/<env>.yaml` when the env resolves identically across variants. Nothing is keyed by tool name, which is why several tools can share one environment — 23 share `biopipelines`, six share `ProteinEnv`, three share `SE3nv` — and why redirecting a tool is a one-line config change:

```yaml
environments:
  DSSP: MyDSSPEnv      # then environments/MyDSSPEnv.yaml must exist, or be built elsewhere
```

Across `config.cluster.yaml`'s 77 entries there are 42 distinct environment names and 45 spec files, and 41 of the 42 have a spec. **The exception is deliberate and worth knowing:** `RF3: modelforge` names an environment that the upstream installer builds, so no biopipelines spec exists for it on cluster, colab or container. So a missing spec is *not* an error — the emitted script tries the variant spec, then the shared one, and otherwise says it is assuming the environment already exists instead of running `env create -f` on a file that is not there. That keeps a hand-built or vendor-built env usable.

That test lives in the generated bash rather than in Python, and it has to: the install script is generated on your machine and executed on the cluster, so the repo path it names is remote and no config-time `os.path.isfile` can see it.

**`env_manager: pip` cannot build an environment at all.** It has no environment manager, so `_env_install_block` refuses with a message naming the variant, rather than emitting `pip env create -f ...` and `pip run -n ...` — neither of which is a pip command, so the old script could only fail later and less legibly.

**Why a raise is the right answer.** On `daint`, 29 tools that declare an `ENV_NAME` have no `environments:` entry at all, because they genuinely are not set up there — aarch64 has no x86-64 binaries and no PyG wheels, and `config.daint.yaml` is explicitly exempt from the three-variant parity test for that reason. Falling back to `ENV_NAME` for those tools would happily build `mytool` on a login node and then activate nothing, or activate an env the runtime side never agreed to. The raise instead names the tool, the variant and the config file:

```
MyTool: no environment is configured for it in config.daint.yaml (variant 'daint'), so
there is no name to install into. Add an environments: entry --
environments.MyTool: "<env name>" -- or install MyTool under a variant that configures it.
```

An install-time refusal that names all three beats the same failure arriving later, on a compute node, as an activation that finds nothing.

An `environments:` entry that is present but **empty** (`MyTool: ""`) means something different and gets its own message: the tool's image supplies the interpreter, so there is no environment to build and installing under that variant is unsupported. `BioEmu` on `daint` is the one such entry today.

So `ENV_NAME` is the tool's declared default, consulted only in pip mode — and the reason to declare it anyway is that it is the single place a reader can find which env the tool means, instead of reading it out of an install script's bash.

**A multi-environment tool's secondary envs are not config-driven.** `DynamicBind` needs two: an inference env and a relax env. The first comes from `environments.DynamicBind` like everyone else's; the second is a class attribute, `RELAX_ENV_NAME = "dynamicbind_relax"`, with no `environments:` entry in any config variant. Install and runtime both address it through that attribute (`cls.RELAX_ENV_NAME` in `_install_script`, `ConfigManager().get_env_python_command(self.RELAX_ENV_NAME)` in `generate_script`), so the two agree and this is not a bug — but it does mean a site cannot redirect a secondary environment the way it can redirect a primary one. `_load_environments` already accepts a list-valued `environments:` entry (`["dynamicbind", "dynamicbind_relax"]`, addressed by `activate_environment(index=N)`) and `_install_env` takes its first element, so the plumbing for config-driven secondaries exists; no config uses it today. If you write a second multi-env tool, know that you are copying a hardcoded pattern.

**Pip-mode installs are non-functional.** Under `env_manager: pip` the shared helpers emit `pip env create -f <yaml>` and `pip run -n <env> ...`, which are not pip commands — pip has no `env` or `run` subcommand, and neither has ever installed anything. This is also why `_install_env` can safely keep returning `ENV_NAME` in pip mode: that path never built an environment it could get wrong. Colab, the variant that uses `env_manager: pip`, runs tools in the pre-existing base Python and installs them by other means. Do not spend an afternoon debugging why a pip-mode install produced nothing.

---

## Shell Safety

### Why

BioPipelines emits bash scripts that are executed on the user's account (laptop, HPC). Several user-supplied strings — pipeline identifiers, tool parameters like `ligand` or `contigs`, config file paths — are interpolated directly into those scripts. The risk model is *footgun, not external attacker*: a space or quote in a project name silently produces a broken script; a `` ` `` or `$(...)` produces confusing execution; and these errors surface far from the Python code the user wrote.

The rule is: **validate user-supplied strings at Python-construction time, before any bash is written.** A Python `ValueError` naming the offending parameter is the behaviour we want — never a cryptic bash syntax error at runtime.

### Where Validation Lives

Four layers, all enforced before any bash is written:

1. **`Pipeline(project=..., job=...)`** — `biopipelines/pipeline.py` (`_validate_identifier`). Restricted to `[A-Za-z0-9._-]+`, plus explicit rejection of `.`, `..`, and leading `-`. These become folder names and `sbatch --job-name` values.

2. **Pipeline description** — escaped at emission (`_escape_for_double_quotes` in `pipeline.py`) so a stray `"`, `` ` ``, `$`, or `\` can't break the surrounding `echo` literal. Validation happens at emission rather than at `__init__`, so the stored `self.description` stays the user's original text for logs and metadata.

3. **Config file (`config.<variant>.yaml`)** — `biopipelines/config_manager.py` (`_validate_shell_safety`, called from `_load_config`). Covers `folders.*.*`, `containers.*`, `machine.{username, scheduler.modules}`. Denylists `"`, `` ` ``, `$`, `\`. The three enum-shaped machine fields (`env_manager`, `scheduler`, `container_executor`) use an **allowlist**, not a denylist, because they are interpolated *unquoted* inside `eval "$(<mgr> shell hook ...)"` where a `;` would suffice to smuggle a command.

4. **Per-tool free-form strings** — each tool's `validate_params()` calls `_validate_freeform_string` (defined in `biopipelines/base_config.py`) on every user-supplied string param that reaches bash. Same `" \` $ \\` denylist.

### Wiring It Up in a New Tool

When adding a new tool, for each `__init__` parameter that is a user-supplied string and ends up in the generated bash (via `echo`, CLI flag, `get_config_display()`, filename built from it), add one line to `validate_params()`:

```python
from .base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string

class MyTool(BaseConfig):
    def __init__(self, ligand: str, positions: Optional[str] = None, ...):
        self.ligand = ligand
        self.positions = positions
        super().__init__(**kwargs)

    def validate_params(self):
        # ...existing checks...
        _validate_freeform_string("ligand", self.ligand)
        _validate_freeform_string("positions", self.positions)
```

For list-typed params, iterate:

```python
for i, spec in enumerate(self.extra_args):
    _validate_freeform_string(f"extra_args[{i}]", spec)
```

For `Union[str, ...]` params, guard the call:

```python
if isinstance(self.position, str):
    _validate_freeform_string("position", self.position)
```

### What Not to Validate This Way

Some parameters are **expression-shaped by design** and legitimately contain characters the denylist would reject:

- Panda `filter` / `calculate` expressions (`pandas.eval`/`query` syntax).
- `Load.filter_input` pandas queries.

These have their own context-specific validators (e.g. `Panda._validate_expression`) or produce a clear error at use time. Do **not** route them through `_validate_freeform_string` — it would reject legitimate input.

Enum-shaped string parameters (e.g. `mode: str` constrained to a fixed set) also don't need the freeform helper; a direct `if value not in {...}: raise` in `validate_params()` is both stricter and clearer.

Tests in `tests/test_shell_safety.py` cover all four validation layers with positive and negative cases, including wiring spot-checks for representative tools. Add a spot-check there when a new tool introduces a non-obvious user-string surface.

---

## pipe scripts Development

pipe scripts (`pipe_scripts/pipe_*.py`) execute at execution time. They process data, generate outputs, and communicate results back to the pipeline.

**Key rule**: pipe scripts must not generate bash code. They process data and write output files (CSV, JSON, FASTA, etc.).

### biopipelines_io Module

The `biopipelines.biopipelines_io` module provides utilities for reading DataStreams and tables at execution time:

```python
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from biopipelines.biopipelines_io import (
    # DataStream utilities
    load_datastream,      # Load DataStream from JSON or dict
    iterate_files,        # Iterate (id, file_path) pairs
    iterate_values,       # Iterate (id, value_dict) pairs from map_table
    resolve_file,         # Get single file for an ID
    get_value,            # Get single value from map_table
    get_all_values,       # Get all values for an ID

    # Table reference utilities
    load_table,           # Load table from path or TABLE_REFERENCE
    lookup_table_value,   # Look up value for an ID
    iterate_table_values, # Iterate (id, value) pairs
)
```

> **Why `sys.path.insert`?** pipe scripts run on SLURM nodes where `biopipelines` is not an installed package. The `sys.path.insert(0, ...)` line adds the repository root so that `from biopipelines.biopipelines_io import ...` resolves correctly. This boilerplate is required in every pipe script that imports from `biopipelines`.

#### DataStream Iteration

For file-based streams (structures, sequences):

```python
ds = load_datastream("/path/to/structures.json")

# Iterate over all files (handles wildcards and lazy patterns)
for struct_id, struct_file in iterate_files(ds):
    process_structure(struct_id, struct_file)

# Get single file
file_path = resolve_file(ds, "protein_1")
```

For value-based streams (SMILES, sequences in map_table):

```python
ds = load_datastream("/path/to/compounds.json")

# Iterate with specific columns
for comp_id, values in iterate_values(ds, columns=['smiles', 'name']):
    smiles = values['smiles']

# Get single value
smiles = get_value(ds, "ligand_001", column="smiles")
```

`load_datastream()` sets `_runtime_mode=True`, which means `ids_expanded` reads the full set of IDs from the map_table CSV. This is how lazy patterns get resolved at runtime.

### pdb_parser Module

The `pdb_parser.py` module provides PDB parsing and selection utilities for pipe scripts:

```python
from biopipelines.pdb_parser import (
    # Data
    Atom,               # NamedTuple: x, y, z, atom_name, res_name, res_num, chain, element
    STANDARD_RESIDUES,  # Set of 20 standard amino acid 3-letter codes

    # Parsing
    parse_pdb_file,     # PDB path → List[Atom]
    get_protein_sequence, # List[Atom] → Dict[chain, sequence]

    # Selection — single entry point
    resolve_selection,  # Selection string + atoms → List[Atom]

    # PyMOL range helpers (pure string operations, no atoms needed)
    parse_pymol_ranges, # "3-45+58-60" → [(3, 45), (58, 60)]
    format_pymol_ranges, # [3, 4, 5, 10, 11] → "3-5+10-11"

    # Distance
    calculate_distance,   # Atom, Atom → float
    calculate_distances,  # List[Atom], List[Atom], metric → float

    # Fixed-column field accessors — read ATOM/HETATM fields through these
    field_atom_name,    # line → "CA"
    field_res_name,     # line → "ALA" or "A1EI4"
    field_chain,        # line → "A", or "" when the record leaves it blank
    field_res_seq,      # line → "201" as a string, insertion code excluded
    field_coords,       # line → (x, y, z) floats
    field_element,      # line → "C", falling back to the atom name's first letter
)
```

#### Field accessors — never inline PDB column offsets

Fourteen pipe scripts plus `converters.py` read ATOM/HETATM fields, and they all go through the `field_*` accessors rather than slicing columns themselves. That is not only DRY: the PDB layout is no longer fixed. The residue-name field is 3 wide, but the PDB now issues **5-character CCD codes** (`A1EI4`), and a record carrying one shifts every field after the code right by its excess width — 2 columns, measured against RCSB's own output.

`_residue_shift(line)` derives that offset from the code's width and every accessor applies it, so one code path handles both layouts. Do not reintroduce a tokenizing (`line.split()`) shortcut: it looks equivalent and fails on three real inputs — a blank chain id (routine in tool-generated ligand PDBs) shifts the token list by one so `field_res_seq` returns a coordinate; an insertion code glues onto the residue number; and two adjacent full-width `%8.3f` coordinates (`-100.123-100.123`) are a single token. Each case made `int()` raise, and the atom was **silently dropped** — an entire ligand disappearing from a parse while every distance and box computed from it looked plausible.

Both parser loops now report what they skip (`_report_dropped`). If you add a third, report too: the silence is what kept that bug invisible.

#### resolve_selection

`resolve_selection(selection, atoms)` is the single function for all atom/residue selection. It handles:

| Syntax | Example | Meaning |
|--------|---------|---------|
| Residue number | `145`, `-1` | Select all atoms of residue 145; last residue |
| Range | `10-20` | Residues 10 through 20 |
| Multiple | `10+15+20` | Residues 10, 15, and 20 |
| Residue.atom | `10.CA`, `-1.C` | Alpha-carbon of residue 10; C of last residue |
| Ligand.atom | `LIG.Cl` | Chlorine atom of ligand LIG |
| Sequence context | `D in IGDWG` | Aspartate within the IGDWG motif |

#### parse_pymol_ranges / format_pymol_ranges

Pure string ↔ tuple conversion for PyMOL-style range strings. No atoms or structures needed:

```python
# String → tuples
parse_pymol_ranges("3-45+58-60")  # [(3, 45), (58, 60)]

# Numbers → string
format_pymol_ranges([3, 4, 5, 10, 11])  # "3-5+10-11"
```

#### Selection format convention

Inter-tool selections (table columns like `within`, `beyond`, `designed`) must use **chain-aware format** (`"A1-50+B10"`) so downstream tools can correctly identify residues. `sele_utils.py` provides the shared conversion functions:

- `sele_to_list(s)` — parse any format (chain-aware, chainless, legacy) → sorted `(chain, resnum)` tuples
- `chain_aware_sele(residues)` / `list_to_sele(a)` — tuples → compact chain-aware string

Tools that accept user-provided position strings (e.g. `fixed="10-20"`) should include a `chain` parameter (default `"A"`) to fill in missing chain info at runtime. Chainless format (`"1-50+10"`) is acceptable only for tool-internal use (e.g. ProteinMPNN's native jsonl format).

### Table References

Tools can pass per-structure data (e.g., fixed positions) to pipe scripts via table references. The format is:

```
TABLE_REFERENCE:/path/to/table.csv:column_name
```

Use `biopipelines_io` to resolve these:

```python
from biopipelines.biopipelines_io import load_table, lookup_table_value, iterate_table_values

# Parse reference and load table
table, column = load_table("TABLE_REFERENCE:/path/to/positions.csv:within")

# Look up value for a specific structure. ALWAYS pass the consumed stream's
# map_table — see "Always pass map_table_paths" below.
ds = load_datastream(structures_json)
maps = [ds.map_table] if ds.map_table else None
positions = lookup_table_value(table, "protein_1", column, map_table_paths=maps)

# Or iterate over all structures
for struct_id, positions in iterate_table_values(table, structure_ids, column,
                                                 map_table_paths=maps):
    print(f"{struct_id}: {positions}")
```

The lookup handles ID matching automatically:
1. Try `pdb` column with `.pdb` extension
2. Try `id` column exact match
3. Try ID mapping (strip suffixes like `_1`, `_2`)
4. Try `<stream>.id` provenance — **only when `map_table_paths` is supplied**

#### Always pass `map_table_paths`

Tiers 1–3 are all string operations on the id. They break the moment an upstream tool renames ids: `Panda_1` has no string relationship to `design_3`, so the lookup raises even though the rename is fully recorded. Tier 4 is what recovers it, and it is unreachable unless the caller supplies the map_tables.

The map to pass is the **map_table of the stream the caller is iterating** — the CSV whose `<stream>.id` columns map the caller's ids back to the ids the referenced table is keyed on. Every pipe script that resolves a table reference already has this: it called `load_datastream()` to get the ids in the first place, so `ds.map_table` is in scope. Scripts handed a map_table CSV directly (`--sequences`, `--sequences-csv`) pass that path.

Do **not** pass the referenced table's own path. A table's `<stream>.id` column records where *its* rows came from, not where its ids were later renamed to — it points the wrong way and will not match.

#### Example: Processing Structures with Per-Structure Data

```python
#!/usr/bin/env python3
"""Example pipe script using biopipelines_io utilities."""

import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from biopipelines.biopipelines_io import (
    load_datastream, iterate_files,
    load_table, lookup_table_value
)

def main():
    structures_json = sys.argv[1]
    positions_ref = sys.argv[2]  # TABLE_REFERENCE:path:column
    output_csv = sys.argv[3]

    # Load structures — ids_expanded works here (runtime mode)
    ds = load_datastream(structures_json)

    # Load positions table
    table, column = load_table(positions_ref)

    # The stream's map_table routes ids an upstream tool renamed.
    maps = [ds.map_table] if ds.map_table else None

    results = []
    for struct_id, struct_file in iterate_files(ds):
        # Get per-structure positions
        positions = lookup_table_value(table, struct_id, column, map_table_paths=maps)

        # Process structure with its positions
        result = process(struct_file, positions)
        results.append({"id": struct_id, "result": result})

    # Write output
    pd.DataFrame(results).to_csv(output_csv, index=False)

if __name__ == "__main__":
    main()
```

### Error Handling

pipe scripts should handle per-item failures gracefully: skip failures, collect partial results, and report at the end.

**Pattern:** wrap the per-item loop body in `try/except`, accumulate failures, write whatever succeeded, and exit with an error only if everything failed.

```python
def main():
    ds = load_datastream(sys.argv[1])
    output_csv = sys.argv[2]

    results = []
    failed = []

    for struct_id, struct_file in iterate_files(ds):
        try:
            result = process(struct_file)
            results.append({"id": struct_id, "score": result})
        except Exception as e:
            print(f"WARNING: {struct_id} failed: {e}", file=sys.stderr)
            failed.append(struct_id)

    # Always write partial results (even if some items failed)
    if results:
        pd.DataFrame(results).to_csv(output_csv, index=False)

    # Report summary
    if failed:
        print(f"Failed {len(failed)}/{len(failed)+len(results)}: {failed}", file=sys.stderr)

    # Exit with error only if ALL items failed
    if not results:
        sys.exit(1)
```

Key points:
- **Always write partial results** — downstream tools can work with a subset.
- **Print failure summary to stderr** — it appears in SLURM logs for debugging.
- **`sys.exit(1)` only if nothing succeeded** — partial success is still useful in a pipeline.

---

## ID Generation and Provenance

All output ID generation and provenance tracking is centralized in `combinatorics.py`. Tool configs and pipe scripts must never implement their own ID logic.

### Output ID Rules

When a tool combines multiple input axes (e.g., proteins × ligands), the output ID is always the full cartesian product of all iterated axes joined with `+`:

| Inputs | Output IDs |
|--------|-----------|
| 1 protein × 3 ligands | `prot1+lig1`, `prot1+lig2`, `prot1+lig3` |
| 2 proteins × 3 ligands | `prot1+lig1`, `prot1+lig2`, ..., `prot2+lig3` |
| 2 proteins × 1 ligand | `prot1+lig1`, `prot2+lig1` |

There are no shortcuts (e.g., dropping single-element axes from the ID). This keeps the ID format predictable and eliminates case-specific logic.

When a tool multiplies inputs by a suffix (e.g., ProteinMPNN generates N sequences per structure), the output ID is `{parent_id}_{suffix}`:

| Tool | Suffix pattern | Example |
|------|---------------|---------|
| ProteinMPNN | `_{seq_num}` | `prot1_1`, `prot1_2` |
| LigandMPNN | `_{seq_num}` | `prot1_1`, `prot1_2` |
| Mutagenesis | `_{position}{aa}` | `prot1_50A`, `prot1_50V` |

### Provenance Columns

Every map_table CSV includes provenance columns named `{alias}.id` that track which input from each axis produced each output row. The alias matches the **tool parameter name** (e.g., `proteins` not `sequences` for Boltz2), making provenance columns unambiguous when a tool accepts multiple inputs of the same stream type.

Example `structures_map.csv` from Boltz2 with 1 protein × 3 ligands:

```csv
id,file,value,proteins.id,ligands.id
prot1+lig1,/path/prot1+lig1.pdb,,prot1,lig1
prot1+lig2,/path/prot1+lig2.pdb,,prot1,lig2
prot1+lig3,/path/prot1+lig3.pdb,,prot1,lig3
```

For multiplier tools, the provenance column tracks the parent:

```csv
id,structures.id,sequence,score,...
struct1_1,struct1,MKTVRQ...,0.95,...
struct1_2,struct1,AETGFT...,0.91,...
struct2_1,struct2,MKTVRQ...,0.88,...
```

Downstream tools can join on any provenance column:
```python
# Find all outputs for a specific protein
df[df['proteins.id'] == 'prot1']

# Join two tables sharing a ligand provenance
pd.merge(confidence, affinity, on=['id', 'ligands.id'])
```

### Shared Utilities

All ID generation uses functions from `combinatorics.py`:

| Function | Use case |
|----------|----------|
| `predict_output_ids_with_provenance()` | Multi-axis tools (Boltz2) at configuration time |
| `predict_output_ids()` | Same, without provenance |
| `predict_single_output_id()` | Single-row ID at execution time (pipe scripts) |
| `generate_multiplied_ids()` | Multiplier tools (ProteinMPNN, Mutagenesis) |
| `generate_multiplied_ids_pattern()` | Same, but keeps compact ID patterns |

#### Using `predict_output_ids_with_provenance` (multi-axis tools)

Each kwarg is a `(value, stream_name)` tuple — the key becomes the provenance column name, `stream_name` tells which stream to extract IDs from. Bare values (no tuple) use the key as both alias and stream name.

```python
from .combinatorics import predict_output_ids_with_provenance

predicted_ids, provenance = predict_output_ids_with_provenance(
    proteins=(self.proteins, "sequences"),
    ligands=(self.ligands, "compounds")
)
# provenance = {"proteins": ["prot1", "prot1", ...], "ligands": ["lig1", "lig2", ...]}
```

`get_output_files()` only declares `DataStream(ids=predicted_ids, files=structure_files, map_table=map_path)`; the pipe script reads `predicted_ids` / `provenance` from the `CombinatoricsConfig` JSON at runtime and writes the map_table CSV with `{stream}.id` provenance columns.

#### Using `generate_multiplied_ids_pattern` (multiplier tools)

Prefer `generate_multiplied_ids_pattern` over `generate_multiplied_ids` to keep IDs compact:

```python
from .combinatorics import generate_multiplied_ids_pattern

suffix_pattern = f"<1..{self.num_sequences}>"
sequence_ids = generate_multiplied_ids_pattern(
    self.structures_stream.ids, suffix_pattern,
    input_stream_name="structures"
)
```

### Pipeline vs SLURM Agreement

The `CombinatoricsConfig` JSON file stores pre-computed `predicted_ids` and `provenance` at configuration time. Pipe scripts read these stored values instead of re-computing:

```json
{
  "axes": { ... },
  "predicted_ids": ["prot1+lig1", "prot1+lig2", "prot1+lig3"],
  "provenance": {
    "sequences": ["prot1", "prot1", "prot1"],
    "compounds": ["lig1", "lig2", "lig3"]
  }
}
```

For pipe scripts that iterate and need single-row IDs, use `predict_single_output_id()` which mirrors the config-time logic exactly:

```python
from combinatorics import predict_single_output_id

config_id = predict_single_output_id(
    sequences=("each", protein_ids, prot_idx, [], False),
    compounds=("each", ligand_ids, lig_idx, static_ids, static_first)
)
```

---

## Testing

The repository ships with an automated pytest suite under `tests/` and a GitHub Actions workflow (`.github/workflows/tests.yml`) that runs it on every push/PR across Python 3.10 / 3.11 / 3.12.

### Testing Scope

The suite is deliberately scoped to **pure-logic correctness and configuration-time wiring**, plus Mock-tool-driven runtime checks. It does **not** exercise any external ML tool (RFdiffusion, AlphaFold, Boltz2, ProteinMPNN, CABSflex, …), any GPU code, or any real model weights — those require dedicated environments, driver-coupled stacks, and hours of compute that are outside what CI can reasonably reproduce, and are outside the framework's own contract.

What the suite covers:

| Module                        | Coverage                                                                                                                                                                               |
| ----------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `test_id_patterns.py`         | Compact `<a..b>` / `<A B>` expansion, deterministic vs lazy classification, `expand_at` (including the single-pattern integer-indexing regression), `append_suffix`, file-template `<id>`. |
| `test_datastream.py`          | `__len__`, integer and slice `__getitem__`, iteration, `StandardizedOutput` ID-based selection, out-of-range behavior.                                                                 |
| `test_combinatorics.py`       | `Bundle` / `Each` axis resolution, `predict_output_ids_with_provenance`, `{alias}.id` provenance columns, `AxisConfig` / `CombinatoricsConfig` round-trip.                             |
| `test_pipeline_generation.py` | End-to-end `Pipeline.save()` emitting a runnable `pipeline.sh`, expected-outputs JSON, ToolOutputs metadata; full bash execution of a 3-Mock chain on POSIX hosts.                     |
| `test_mock.py`                | Every documented Mock pattern: explicit / source IDs, `Bundle` / `Each`, deterministic + lazy `children` with `produce`, `map_table_strategy`, `missing`, table `fill`, multi-stream. |
| `test_provenance.py`          | Multi-hop parent→child provenance through chained Mocks and Mock → Panda → Mock cycles; `<stream>.id` / `<stream>.parent` provenance columns on the materialized `<stream>_map.csv` files. |
| `test_folders.py`             | Filesystem-layout prediction and folder resolution.                                                                                                                                    |
| `test_remap.py`               | ID remapping rules.                                                                                                                                                                    |
| `test_panda.py`               | Table transformations, filter/sort/head/tail/sample rename paths.                                                                                                                      |
| `test_id_filter_audit.py`     | Filtered-stream consumers: `write_filtered_map_table` (subset, order, byte-identical fast path, missing/required columns), combinatorics carrying filtered ids, the runtime materializer pipe script, and stitch/Boltz loaders honoring an id filter. |

All tests use the pip-mode fixture `tests/fixtures/config.local.yaml` (no conda, no containers, no SLURM), so CI and local runs need only `pip install -e ".[test]"`.

### The Mock Tool

The Mock tool (`biopipelines/mock.py` + `pipe_scripts/pipe_mock.py`) is the primary driver for runtime tests. It is a **stub-output generator** that implements the full BaseConfig contract — streams, tables, map tables, `children` / `produce`, `missing`, `source` — but instead of running any model it just creates empty files and well-formed CSV map tables at the paths the framework predicts.

This makes it the right tool for testing framework plumbing end-to-end:

- It exercises the same code paths a real tool does: config serialization to `mock_config.json`, script emission with `generate_completion_check_header`, activation, and the runtime reading of its own config.
- Tests can compose it arbitrarily (Mock → Mock, Mock → Panda → Mock, multi-axis `Each(a) × Each(b)` into a fan-out Mock) without any external dependency.
- Because the runtime is deterministic and cheap (~3 s for the whole suite), provenance and ID-expansion invariants can be asserted against real materialized `<stream>_map.csv` files rather than against mocked-in-memory state.

When adding a framework-level feature (a new ID pattern, a new provenance column, a new cycle pattern), prefer writing the test against Mock rather than against a real tool. Only reach for a real tool when the behavior under test is specific to that tool's payload script.

### Running the Suite

```bash
pip install -e ".[test]"
pytest tests/ -v
```

`tests/conftest.py` writes a per-run report to `tests/test_results.csv` / `.xlsx` with one row per test — `input`, `expected`, `actual`, `matched`, `duration_s`, and any error details — so the suite doubles as a legible artifact of what was tested, not just a pass/fail count.

The end-to-end test that actually invokes `pipeline.sh` via `bash` (`test_generated_pipeline_sh_executes_end_to_end`) is skipped on Windows because MSYS/Git-bash reinterprets the backslashes in embedded Windows paths as escape characters. It runs on every push/PR on the Linux CI runner.

---

## Code Principles

### No Fallbacks

Code must crash explicitly. Never guess values or use defaults for missing data.

```python
# BAD
def get_config(path):
    if not os.path.exists(path):
        return {"default": "value"}  # Fallback!
    return read_config(path)

# GOOD
def get_config(path):
    if not os.path.exists(path):
        raise ValueError(f"Config not found: {path}")
    return read_config(path)
```

### Single Definition

Define paths once, use everywhere:

```python
# BAD
def generate_script(self):
    csv = os.path.join(self.output_folder, "results.csv")
    # ... later ...
    another_csv = os.path.join(self.output_folder, "results.csv")  # Duplicate!

# GOOD
results_csv = Path(lambda self: os.path.join(self.output_folder, "results.csv"))

def generate_script(self):
    # Use self.results_csv consistently
```

### Tool Agnostic

Tools work without knowing upstream tool identities:

```python
# BAD
if isinstance(input_tool, Boltz2):
    structures = input_tool.boltz2_specific_output

# GOOD
if hasattr(input_tool, 'structures'):
    structures = input_tool.structures
else:
    raise ValueError("Input must have structures attribute")
```

### Prediction-Based

Never check file existence during configuration:

```python
# BAD
def configure_inputs(self):
    if os.path.exists(predicted_file):
        self.input_file = predicted_file

# GOOD
def configure_inputs(self):
    # Predict path without checking existence
    self.input_file = os.path.join(self.folders["data"], "file.txt")
    # File will exist at execution time
```

### pipe scripts Don't Write Bash

pipe scripts run Python at execution time. They must produce data outputs (CSV, JSON, FASTA) — never bash scripts. If a tool needs per-structure bash commands, use a `for` loop in the generated script with `Resolve.stream_ids()`.

### Validate User Strings at Construction Time

If a user-supplied string reaches the generated bash, validate it in Python before any script is written — never rely on bash to catch the problem. See [Shell Safety](#shell-safety).

---

## Working with Git

### Daily Workflow

```bash
git status
git pull origin main
git checkout -b feature/my-tool
# Make changes
git add biopipelines/my_tool.py
git commit -m "Add MyTool"
git push origin feature/my-tool
```

### What to Commit

- Tool implementations (`biopipelines/`)
- Helper scripts (`pipe_scripts/`)
- Documentation (`docs/`)

### What to Ignore

- Generated scripts (`RunTime/`)
- Job outputs (`*_[0-9]*/`)
- Cache (`__pycache__/`)

---

## Working with Claude Code

### Start Sessions with Context

Open biopipelines folder and start with `prompt.txt` to give Claude repository context.

### Use Plan Mode

```
/plan Create a new tool for binding analysis
```

### Reference Existing Patterns

```
"Before creating MyTool, read Distance to understand the pattern"
"Check how Boltz2 handles multiple input types and use the same pattern"
```

### Be Explicit

```
# BAD
"Update the tool to handle MSAs"

# GOOD
"Update Boltz2 to:
1. Accept msas parameter of type Union[str, ToolOutput]
2. Extract MSA files from tables.msas
3. Pass MSA folder path to boltz predict with --msa-cache flag
4. Do NOT create fallback MSA generation"
```

### Verify Changes

```bash
# Quick syntax check
python -c "from biopipelines.my_tool import MyTool; print('OK')"

# Test with local output (writes to ./BioPipelines/)
with Pipeline("Test", "Debug", "Testing", local_output=True):
    result = MyTool(...)
    print(result)
```

#### Testing a pipe script in Isolation

Create a mock DataStream JSON and run the script directly:

```bash
# Create mock input
echo '{"name":"structures","ids":["prot_1","prot_2"],"files":["/tmp/prot_1.pdb","/tmp/prot_2.pdb"],"map_table":"","format":"pdb"}' > /tmp/test_ds.json

# Run the pipe_script
python pipe_scripts/pipe_my_tool.py /tmp/test_ds.json /tmp/output.csv

# Inspect output
cat /tmp/output.csv
```

#### Testing Config-Time Logic

Instantiate the tool with `local_output=True` and inspect the generated outputs:

```python
from biopipelines.pipeline import Pipeline
from biopipelines.my_tool import MyTool

with Pipeline("Test", "Debug", "Testing", local_output=True):
    result = MyTool(structures=..., param1="test")
    # Check predicted output IDs
    print(result.output.structures.ids)
    # Read the generated bash script
    with open(result.script_path) as f:
        print(f.read())
```

> **Note:** The `tests/` directory holds the automated pytest suite (run `pytest tests/`). Local pipeline runs (via `local_output=True`) write into `outputs/`, which is gitignored.

### Critical Files (rarely need changes)

- `base_config.py` (`BaseConfig` itself; `ToolOutput`/`StandardizedOutput` live in `outputs.py` and the table and stream containers in `data_containers.py`, both re-exported here)
- `pipeline.py`
- `datastream.py`
- `combinatorics.py`

If Claude suggests changing these, question why. Usually the tool should adapt to the base class.
