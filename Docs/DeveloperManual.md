# BioPipelines Developer Manual

## Index

- [Architecture](#architecture)
  - [Two-Phase Execution](#two-phase-execution)
  - [Core Classes](#core-classes)
  - [Data Flow](#data-flow)
- [Tool Development](#tool-development)
  - [Creating a Tool](#creating-a-tool)
  - [Required Methods](#required-methods)
  - [Path Descriptors](#path-descriptors)
  - [Script Generation](#script-generation)
  - [Output Prediction](#output-prediction)
- [HelpScript Development](#helpscript-development)
  - [pipe_biopipelines_io Module](#pipe_biopipelines_io-module)
  - [pdb_parser Module](#pdb_parser-module)
  - [Table References](#table-references)
- [Code Principles](#code-principles)
- [Working with Git](#working-with-git)
- [Working with Claude Code](#working-with-claude-code)

---

## Architecture

### Two-Phase Execution

| Phase | Location | What Happens |
|-------|----------|--------------|
| **Pipeline time** | biopipelines/ | Python generates bash scripts, predicts outputs |
| **SLURM time** | HelpScripts/ | Bash scripts execute, `pipe_*.py` scripts run |

Tools never execute computations directly. They:
1. Validate parameters
2. Predict output file paths
3. Generate bash scripts

### Core Classes

**BaseConfig** (`base_config.py`) - Abstract base for all tools:
```python
class MyTool(BaseConfig):
    TOOL_NAME = "MyTool"

    def validate_params(self): ...
    def configure_inputs(self, pipeline_folders): ...
    def generate_script(self, script_path): ...
    def get_output_files(self): ...
```

**DataStream** (`datastream.py`) - Unified container for files and IDs:
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
- `ids` - List of identifiers
- `files` - List of file paths (may be empty for value-based formats)
- `map_table` - CSV with additional metadata
- `format` - Data format (pdb, cif, fasta, csv, smiles, etc.)

**TableInfo** (`base_config.py`) - Metadata for CSV outputs:
```python
TableInfo(
    name="results",
    path="/path/to/results.csv",
    columns=["id", "score", "pLDDT"],
    description="Analysis results",
    count=100
)
```

**StandardizedOutput** (`base_config.py`) - Tool output wrapper:
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
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo
    from file_paths import Path
    from datastream import DataStream


class MyTool(BaseConfig):
    TOOL_NAME = "MyTool"

    # Path descriptors (lazy evaluation)
    results_csv = Path(lambda self: os.path.join(self.output_folder, "results.csv"))
    helper_py = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_my_tool.py"))

    def __init__(self,
                 structures: Union[DataStream, StandardizedOutput],
                 param1: str,
                 param2: int = 10,
                 **kwargs):
        # Resolve input to DataStream
        if isinstance(structures, StandardizedOutput):
            self.structures_stream = structures.structures
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
        script_content = "#!/bin/bash\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()
        script_content += f"""
echo "Running MyTool"
python {self.helper_py} \\
    --input "{','.join(self.structures_stream.files)}" \\
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

Use provided helpers:

```python
def generate_script(self, script_path: str) -> str:
    script_content = "#!/bin/bash\n"
    script_content += self.generate_completion_check_header()  # Skip if completed
    script_content += self.activate_environment()              # Activate conda
    script_content += """
# Your commands here
python script.py --args
"""
    script_content += self.generate_completion_check_footer()  # Check outputs
    return script_content
```

### Output Prediction

Return standardized output structure:

```python
def get_output_files(self) -> Dict[str, Any]:
    # Predict file paths
    structure_files = [
        os.path.join(self.output_folder, f"{sid}.pdb")
        for sid in self.predicted_ids
    ]

    # Create DataStreams
    structures = DataStream(
        name="structures",
        ids=self.predicted_ids,
        files=structure_files,
        map_table=self.structures_csv,
        format="pdb"
    )

    # Create tables
    tables = {
        "results": TableInfo(
            name="results",
            path=self.results_csv,
            columns=["id", "score"],
            description="Analysis results",
            count=len(self.predicted_ids)
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

---

## HelpScript Development

HelpScripts (`HelpScripts/pipe_*.py`) execute at SLURM runtime. They process data, generate outputs, and communicate results back to the pipeline.

### pipe_biopipelines_io Module

The `pipe_biopipelines_io.py` module provides utilities for reading DataStreams and tables at SLURM runtime:

```python
from pipe_biopipelines_io import (
    # DataStream utilities
    load_datastream,      # Load DataStream from JSON or dict
    iterate_files,        # Iterate (id, file_path) pairs
    iterate_values,       # Iterate (id, value_dict) pairs from map_table
    resolve_file,         # Get single file for an ID
    get_value,            # Get single value from map_table
    get_all_values,       # Get all values for an ID

    # Table reference utilities
    load_table,           # Load table from path or DATASHEET_REFERENCE
    lookup_table_value,   # Look up value for an ID
    iterate_table_values, # Iterate (id, value) pairs
)
```

#### DataStream Iteration

For file-based streams (structures, sequences):

```python
ds = load_datastream("/path/to/structures.json")

# Iterate over all files (handles wildcards automatically)
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

### pdb_parser Module

The `pdb_parser.py` module provides PDB parsing and selection utilities for HelpScripts:

```python
from pdb_parser import (
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
DATASHEET_REFERENCE:/path/to/table.csv:column_name
```

Use `pipe_biopipelines_io` to resolve these:

```python
from pipe_biopipelines_io import load_table, lookup_table_value, iterate_table_values

# Parse reference and load table
table, column = load_table("DATASHEET_REFERENCE:/path/to/positions.csv:within")

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
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from pipe_biopipelines_io import (
    load_datastream, iterate_files,
    load_table, lookup_table_value
)

def main():
    structures_json = sys.argv[1]
    positions_ref = sys.argv[2]  # DATASHEET_REFERENCE:path:column
    output_csv = sys.argv[3]

    # Load structures
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
    # File will exist at SLURM runtime
```

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

### Critical Files (rarely need changes)

- `base_config.py`
- `pipeline.py`
- `datastream.py`
- `standardized_output.py`

If Claude suggests changing these, question why. Usually the tool should adapt to the base class.
