# BioPipelines Developer Manual

## Index

- [Architecture](#architecture)
  - [Two-Phase Execution](#two-phase-execution)
  - [Core Classes](#core-classes)
  - [Data Flow](#data-flow)
- [IDs: Configuration Time vs Execution Time](#ids-configuration-time-vs-execution-time)
  - [ID Patterns](#id-patterns)
  - [Lazy IDs](#lazy-ids)
  - [The Rule](#the-rule)
  - [How to Iterate Structures in Generated Bash](#how-to-iterate-structures-in-generated-bash)
  - [How to Resolve a Single File](#how-to-resolve-a-single-file)
  - [Passing IDs to HelpScripts](#passing-ids-to-helpscripts)
- [Tool Development](#tool-development)
  - [Creating a Tool](#creating-a-tool)
  - [Required Methods](#required-methods)
  - [Path Descriptors](#path-descriptors)
  - [Script Generation](#script-generation)
  - [Output Prediction](#output-prediction)
    - [File Templates with \<id\>](#file-templates-with-id)
- [HelpScript Development](#helpscript-development)
  - [biopipelines_io Module](#biopipelines_io-module)
  - [pdb_parser Module](#pdb_parser-module)
  - [Table References](#table-references)
  - [Error Handling](#error-handling)
- [ID Generation and Provenance](#id-generation-and-provenance)
  - [Output ID Rules](#output-id-rules)
  - [Provenance Columns](#provenance-columns)
  - [Shared Utilities](#shared-utilities)
  - [Pipeline vs SLURM Agreement](#pipeline-vs-slurm-agreement)
- [Code Principles](#code-principles)
- [Working with Git](#working-with-git)
- [Working with Claude Code](#working-with-claude-code)

---

## Architecture

### Two-Phase Execution

| Phase | Location | What Happens |
|-------|----------|--------------|
| **Configuration time** | `biopipelines/` | Python generates bash scripts, predicts outputs |
| **Execution time** | `HelpScripts/` | Bash scripts execute, `pipe_*.py` scripts run |

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
- `files` — List of file paths (may be empty for value-based formats)
- `map_table` — CSV with additional metadata
- `format` — Data format (pdb, cif, fasta, csv, smiles, etc.)

**TableInfo** (`base_config.py`) — Metadata for CSV outputs:
```python
TableInfo(
    name="results",
    path="/path/to/results.csv",
    columns=["id", "score", "pLDDT"],
    description="Analysis results",
    count=100
)
```

**StandardizedOutput** (`base_config.py`) — Tool output wrapper:
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

---

## IDs: Configuration Time vs Execution Time

### ID Patterns

DataStream IDs can be stored in compact pattern form instead of listing every ID explicitly:

| Pattern | Meaning | Expansion |
|---------|---------|-----------|
| `prot_<0..2>` | Numeric range | `prot_0`, `prot_1`, `prot_2` |
| `<A B C>` | Enumeration | `A`, `B`, `C` |
| `prot_<N><S A L K>` | Deterministic prefix + lazy suffix | See below |

The `<..>` angle-bracket patterns are **deterministic** — they can be fully expanded at configuration time.

### Lazy IDs

Bracket `[...]` segments mark parts of an ID that depend on runtime data and **cannot be expanded at configuration time**:

```
prot[_<N><S A L K>]
```

Here `prot` is the deterministic prefix, and `[_<N><S A L K>]` is a lazy suffix whose actual values come from an upstream tool's output (which doesn't exist yet at config time). Lazy IDs expand at runtime by matching patterns against IDs found in the DataStream's `map_table` CSV.

At configuration time, `ids_expanded` on a lazy DataStream returns only the deterministic prefix:
```python
ds.ids_expanded  # → ["prot"] (incomplete)
```

At execution time (with `_runtime_mode=True`, as set by `load_datastream()`), `ids_expanded` reads the full set of IDs from the map_table:
```python
ds = load_datastream("structures.json")
ds.ids_expanded  # → ["prot_1S", "prot_2A", "prot_3L", ...] (complete)
```

### The Rule

**Never call `ids_expanded` or `files_expanded` at configuration time to generate per-ID bash commands.** This breaks when the DataStream has lazy patterns, because the full ID list doesn't exist yet.

Instead:
- **Serialize the DataStream** to JSON via `save_json()` at config time
- **Expand IDs at runtime** inside the generated bash script using `Resolve.stream_ids()` or inside HelpScripts using `load_datastream()` + `ids_expanded` / `iterate_files()`

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

Note: `ids[0]` (not `ids_expanded[0]`) is safe here because we only need the compact pattern, not the expanded list.

### Passing IDs to HelpScripts

Don't build per-ID data structures at config time. Instead, pass the DataStream JSON path to HelpScripts and let them expand IDs at runtime:

```python
# BAD — breaks with lazy IDs
id_map = {}
for sid, path in zip(ds.ids_expanded, ds.files_expanded):
    id_map[os.path.basename(path)] = sid

# GOOD — HelpScript builds the map at runtime
script += f'python {self.helper_py} --ds-json "{self.structures_json}"\n'
```

In the HelpScript:
```python
from biopipelines.biopipelines_io import load_datastream, iterate_files

ds = load_datastream(sys.argv[1])
for struct_id, pdb_path in iterate_files(ds):
    process(struct_id, pdb_path)
```

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
| Build per-ID data | **Don't** | Do it in HelpScripts |

---

## Tool Development

### Creating a Tool

1. Create `biopipelines/my_tool.py`
2. Inherit from `BaseConfig`
3. Implement required methods
4. Create helper script `HelpScripts/pipe_my_tool.py` if needed

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

> **Why the dual import?** Tools are normally imported as part of the `biopipelines` package (relative imports via `.base_config`, etc.). The `except ImportError` fallback adds the module's own directory to `sys.path` so the same file can also be imported standalone — useful for debugging or running a tool file directly. Always include this pattern in new tool files.


class MyTool(BaseConfig):
    TOOL_NAME = "MyTool"

    # Path descriptors (lazy evaluation)
    results_csv = Path(lambda self: os.path.join(self.output_folder, "results.csv"))
    structures_json = Path(lambda self: os.path.join(self.output_folder, ".input_structures.json"))
    helper_py = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_my_tool.py"))

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
                description="MyTool results",
                count=len(self.structures_stream)
            )
        }

        return {
            "structures": DataStream.empty("structures", "pdb"),
            "sequences": DataStream.empty("sequences", "fasta"),
            "compounds": DataStream.empty("compounds", "sdf"),
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

### Path Descriptors

Use `Path` descriptors for lazy path evaluation:

```python
from .file_paths import Path

class MyTool(BaseConfig):
    # Evaluated when accessed, after output_folder is set
    results_csv = Path(lambda self: os.path.join(self.output_folder, "results.csv"))
    helper_py = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_my_tool.py"))
```

### Script Generation

Use provided helpers and runtime resolution:

```python
def generate_script(self, script_path: str) -> str:
    # Save DataStream for runtime access
    self.structures_stream.save_json(self.structures_json)

    script_content = "#!/bin/bash\n"
    script_content += self.generate_completion_check_header()  # Skip if completed
    script_content += self.activate_environment()              # Activate conda + source resolve_stream_item.sh

    # Option A: Pass DataStream to HelpScript (preferred for complex processing)
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

### Output Prediction

Return standardized output structure. Use compact `ids` patterns (not expanded) and `<id>` file templates:

#### File Templates with \<id\>

`<id>` is a placeholder in the `files` list that represents all output files at once. Instead of listing one file per ID, you provide a single-element list like `["<id>.pdb"]`. At expansion time, `<id>` is replaced with each expanded ID. For example:

- `files=["<id>.pdb"]` + `ids=["prot_<0..2>"]` → `prot_0.pdb`, `prot_1.pdb`, `prot_2.pdb`

This prevents a length-mismatch validation error (1 file vs N ids) and is **required** when inputs may carry lazy IDs — building per-ID file paths with f-strings embeds bracket patterns into paths, causing `LazyPatternError`. The implementation lives in `datastream.py:_has_file_template()` (detects the pattern) and `id_patterns.py:expand_file_pattern()` (performs the substitution).

```python
def get_output_files(self) -> Dict[str, Any]:
    # Keep IDs compact — don't expand
    structure_ids = self.structures_stream.ids
    structure_files = [os.path.join(self.output_folder, "<id>.pdb")]

    structures = DataStream(
        name="structures",
        ids=structure_ids,
        files=structure_files,
        map_table=self.structures_csv,
        format="pdb"
    )

    tables = {
        "results": TableInfo(
            name="results",
            path=self.results_csv,
            columns=["id", "score"],
            description="Analysis results",
            count=len(self.structures_stream)
        )
    }

    return {
        "structures": structures,
        "sequences": DataStream.empty("sequences", "fasta"),
        "compounds": DataStream.empty("compounds", "sdf"),
        "tables": tables,
        "output_folder": self.output_folder
    }
```

`get_output_files()` returns a plain dict, **not** a `StandardizedOutput`. The wrapping into `StandardizedOutput` happens automatically in the `ToolOutput.output` property (`base_config.py:1768`). This is by design:

- Pipeline infrastructure (`pipeline.py`, `base_config.py:get_id_provenance()`) iterates the raw dict with `.items()` and `isinstance()` checks before any user accesses it — dicts are natural for this.
- `StandardizedOutput.__init__` takes a dict and destructures it into `.streams`, `.tables`, etc. — returning `StandardizedOutput` from tools would just add an object whose constructor immediately unpacks it back.
- Keeps tools simple: tool authors build plain dicts, the framework handles the user-facing API.

**Rule:** tool authors return dicts; consumers get `StandardizedOutput` with dot-notation (e.g., `tool.output.structures`).

---

## HelpScript Development

HelpScripts (`HelpScripts/pipe_*.py`) execute at execution time. They process data, generate outputs, and communicate results back to the pipeline.

**Key rule**: HelpScripts must not generate bash code. They process data and write output files (CSV, JSON, FASTA, etc.).

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

> **Why `sys.path.insert`?** HelpScripts run on SLURM nodes where `biopipelines` is not an installed package. The `sys.path.insert(0, ...)` line adds the repository root so that `from biopipelines.biopipelines_io import ...` resolves correctly. This boilerplate is required in every HelpScript that imports from `biopipelines`.

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

The `pdb_parser.py` module provides PDB parsing and selection utilities for HelpScripts:

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
)
```

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

### Table References

Tools can pass per-structure data (e.g., fixed positions) to HelpScripts via table references. The format is:

```
TABLE_REFERENCE:/path/to/table.csv:column_name
```

Use `biopipelines_io` to resolve these:

```python
from biopipelines.biopipelines_io import load_table, lookup_table_value, iterate_table_values

# Parse reference and load table
table, column = load_table("TABLE_REFERENCE:/path/to/positions.csv:within")

# Look up value for a specific structure
positions = lookup_table_value(table, "protein_1", column)

# Or iterate over all structures
for struct_id, positions in iterate_table_values(table, structure_ids, column):
    print(f"{struct_id}: {positions}")
```

The lookup handles ID matching automatically:
1. Try `pdb` column with `.pdb` extension
2. Try `id` column exact match
3. Try ID mapping (strip suffixes like `_1`, `_2`)

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

    results = []
    for struct_id, struct_file in iterate_files(ds):
        # Get per-structure positions
        positions = lookup_table_value(table, struct_id, column)

        # Process structure with its positions
        result = process(struct_file, positions)
        results.append({"id": struct_id, "result": result})

    # Write output
    pd.DataFrame(results).to_csv(output_csv, index=False)

if __name__ == "__main__":
    main()
```

### Error Handling

HelpScripts should handle per-item failures gracefully: skip failures, collect partial results, and report at the end.

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

When a tool combines multiple input axes (e.g., proteins × ligands), the output ID is always the full cartesian product of all iterated axes joined with `_`:

| Inputs | Output IDs |
|--------|-----------|
| 1 protein × 3 ligands | `prot1_lig1`, `prot1_lig2`, `prot1_lig3` |
| 2 proteins × 3 ligands | `prot1_lig1`, `prot1_lig2`, ..., `prot2_lig3` |
| 2 proteins × 1 ligand | `prot1_lig1`, `prot2_lig1` |
| All bundled | `bundled_complex` |
| Single iterated axis | IDs from that axis directly |

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
prot1_lig1,/path/prot1_lig1.pdb,,prot1,lig1
prot1_lig2,/path/prot1_lig2.pdb,,prot1,lig2
prot1_lig3,/path/prot1_lig3.pdb,,prot1,lig3
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
    bundled_name="bundled_complex",
    proteins=(self.proteins, "sequences"),
    ligands=(self.ligands, "compounds")
)
# provenance = {"proteins": ["prot1", "prot1", ...], "ligands": ["lig1", "lig2", ...]}

create_map_table(map_path, predicted_ids, files=structure_files, provenance=provenance)
```

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

#### Passing provenance to `create_map_table`

```python
from .datastream import create_map_table

create_map_table(
    output_path, ids=predicted_ids, files=structure_files,
    provenance=provenance  # Dict keys become {key}.id columns
)
```

### Pipeline vs SLURM Agreement

The `CombinatoricsConfig` JSON file stores pre-computed `predicted_ids` and `provenance` at configuration time. Pipe scripts read these stored values instead of re-computing:

```json
{
  "axes": { ... },
  "predicted_ids": ["prot1_lig1", "prot1_lig2", "prot1_lig3"],
  "provenance": {
    "sequences": ["prot1", "prot1", "prot1"],
    "compounds": ["lig1", "lig2", "lig3"]
  }
}
```

For pipe scripts that iterate and need single-row IDs, use `predict_single_output_id()` which mirrors the pipeline-time logic exactly:

```python
from combinatorics import predict_single_output_id

config_id = predict_single_output_id(
    bundled_name="bundled_complex",
    sequences=("each", protein_ids, prot_idx),
    compounds=("each", ligand_ids, lig_idx)
)
```

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

### HelpScripts Don't Write Bash

HelpScripts run Python at execution time. They must produce data outputs (CSV, JSON, FASTA) — never bash scripts. If a tool needs per-structure bash commands, use a `for` loop in the generated script with `Resolve.stream_ids()`.

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
- Helper scripts (`HelpScripts/`)
- Documentation (`Docs/`)

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

#### Testing a HelpScript in Isolation

Create a mock DataStream JSON and run the script directly:

```bash
# Create mock input
echo '{"name":"structures","ids":["prot_1","prot_2"],"files":["/tmp/prot_1.pdb","/tmp/prot_2.pdb"],"map_table":"","format":"pdb"}' > /tmp/test_ds.json

# Run the HelpScript
python HelpScripts/pipe_my_tool.py /tmp/test_ds.json /tmp/output.csv

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

> **Note:** No automated test suite exists yet; the `tests/` directory is for manual output inspection.

### Critical Files (rarely need changes)

- `base_config.py`
- `pipeline.py`
- `datastream.py`
- `standardized_output.py`

If Claude suggests changing these, question why. Usually the tool should adapt to the base class.
