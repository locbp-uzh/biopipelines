# Utilities

[← Back to Tool Reference](../ToolReference.md)

---

## Entity Types

Basic input types exported from `PipelineScripts.entities`:

```python
from PipelineScripts.entities import PDB, Sequence, Ligand, CompoundLibrary
```

---

### PDB

Fetches protein structures with priority: `local_folder` → `PDBs/` → RCSB download.

**Environment**: `biopipelines`

**Parameters**:
- `pdbs`: str | List[str] - PDB IDs or folder path
- `ids`: str | List[str] = None - Custom IDs (defaults to pdbs)
- `format`: str = "pdb" - Output format ("pdb" or "cif")
- `local_folder`: str = None - Check first before PDBs/
- `biological_assembly`: bool = False - Download biological assembly
- `remove_waters`: bool = True - Remove water molecules

**Operations** (optional positional args):
- `PDB.Rename(old, new)` - Rename residues (e.g., for RFdiffusion3 compatibility)

**Outputs**:
- `structures` - DataStream of PDB files
- `sequences` - DataStream of sequences
- `compounds` - DataStream of ligands extracted from structures
- `tables.structures`: | id | pdb_id | file_path | format | source |
- `tables.sequences`: | id | sequence |
- `tables.compounds`: | id | code | smiles | ccd |

**Examples**:

```python
from PipelineScripts.entities import PDB

# Simple fetch
protein = PDB("4ufc")

# Multiple with custom IDs
proteins = PDB(["4ufc", "1aki"], ids=["POI1", "POI2"])

# From folder
proteins = PDB("/path/to/structures")

# With ligand renaming
pdb = PDB("structure.pdb", PDB.Rename("LIG", ":L:"))

# Biological assembly
pdb = PDB("4ufc", biological_assembly=True, remove_waters=False)
```

---

### Sequence

Creates sequences from strings with auto-detection (protein/DNA/RNA).

**Environment**: `biopipelines`

**Parameters**:
- `seq`: str | List[str] - Sequence string(s)
- `type`: str = "auto" - Sequence type ("auto", "protein", "dna", "rna")
- `ids`: str | List[str] = None - Custom IDs (defaults to "seq_N")

**Outputs**:
- `sequences` - DataStream of sequences
- `tables.sequences`: | id | sequence | type | length |

**Examples**:

```python
from PipelineScripts.entities import Sequence

# Single sequence
seq = Sequence("MKTVRQERLKSIVRILERSKEPVSGAQ", ids="my_protein")

# Multiple sequences
seqs = Sequence(["MKTVRQ...", "AETGFT..."], ids=["p1", "p2"])

# DNA sequence
dna = Sequence("ACGTACGT", type="dna", ids="my_gene")
```

---

### Ligand

Fetches small molecules from RCSB (CCD) or PubChem (name, CID, CAS) or generates from SMILES.

**Environment**: `biopipelines`

**Parameters**:
- `lookup`: str | List[str] = None - Lookup values (CCD, CID, CAS, or name)
- `ids`: str | List[str] = None - Custom IDs (defaults to lookup)
- `codes`: str | List[str] = None - 3-letter PDB residue codes (defaults to lookup[:3])
- `source`: str = None - Force "rcsb" or "pubchem" (auto-detects if None)
- `local_folder`: str = None - Check first before Ligands/
- `output_format`: str = "pdb" - Output format ("pdb" or "cif")
- `smiles`: str | List[str] = None - Direct SMILES input (bypasses lookup)

**Auto-detection** (when source=None):
- 1-3 uppercase alphanumeric → RCSB (CCD)
- Purely numeric → PubChem (CID)
- XX-XX-X format → PubChem (CAS)
- Otherwise → PubChem (name)

**Outputs**:
- `structures` - DataStream of ligand files
- `compounds` - DataStream with SMILES data
- `tables.compounds`: | id | code | lookup | source | smiles | formula |

**Examples**:

```python
from PipelineScripts.entities import Ligand

# RCSB by CCD code
atp = Ligand("ATP")

# PubChem by name
aspirin = Ligand("aspirin", codes="ASP")

# PubChem by CID
caffeine = Ligand("2157", ids="caffeine", codes="CAF")

# PubChem by CAS
ibuprofen = Ligand("15687-27-1", ids="ibuprofen", codes="IBU")

# Direct SMILES
ethanol = Ligand(smiles="CCO", ids="ethanol", codes="ETH")

# Multiple SMILES
ligands = Ligand(smiles=["CCO", "CC(=O)O"], ids=["ethanol", "acetic"], codes=["ETH", "ACE"])
```

---

### CompoundLibrary

Creates compound collections from dictionaries with optional combinatorial expansion.

**Environment**: `biopipelines`

**Parameters**:
- `library`: str | Dict - Dictionary with SMILES or path to CSV
- `primary_key`: str = None - Root key for expansion (enables `<key>` placeholders)
- `covalent`: bool = False - Generate CCD/PKL files for covalent binding
- `validate_smiles`: bool = True - Validate SMILES during expansion
- `conformer_method`: str = "UFF" - Conformer method (UFF, OpenFF, DFT)

**Outputs**:
- `compounds` - DataStream with library compounds
- `tables.compounds`: | id | format | smiles | ccd | {branching_keys} |

**Examples**:

```python
from PipelineScripts.entities import CompoundLibrary

# Simple dictionary (no expansion)
library = CompoundLibrary({
    "aspirin": "CC(=O)OC1=CC=CC=C1C(=O)O",
    "caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
})

# With expansion using <key> placeholders
library = CompoundLibrary(
    library={
        "scaffold": "<linker><fluorophore>",
        "linker": ["CCOCC", "CCOCCOCC"],
        "fluorophore": ["c1ccc(N)cc1"]
    },
    primary_key="scaffold"
)
# Generates 2 compounds: linker×fluorophore combinations
```

---

## Loading Outputs

### LoadOutput

Loads a single tool output JSON file.

**Parameters**:
- `output_json`: str - Path to JSON file
- `filter`: str = None - Pandas query filter
- `validate_files`: bool = True - Check file existence

**Example**:

```python
from PipelineScripts.load import LoadOutput

prev = LoadOutput(
    output_json="/path/to/ToolOutputs/003_Boltz2.json",
    filter="confidence_score > 0.8"
)
```

### LoadOutputs

Loads multiple outputs from a ToolOutputs folder.

**Parameters**:
- `path`: str - Path to job folder or ToolOutputs folder
- `tool`: str = None - Filter by tool name
- `suffix`: str = None - Filter by filename suffix
- `ascending`: bool = True - Sort order

**Returns**: Dictionary mapping identifiers to LoadOutput objects

**Example**:

```python
from PipelineScripts.load import LoadOutputs

# All Boltz2 outputs
all_boltz = LoadOutputs("/path/to/job/", tool="Boltz2")

# Access: all_boltz["003_Boltz2"], all_boltz["007_Boltz2"], etc.

# Filter by suffix
cycle10 = LoadOutputs("/path/to/job/", suffix="Cycle10")
```

---

## MSA Generation

### MMseqs2

Generates multiple sequence alignments for structure prediction.

**Environment**: `ProteinEnv`

**Parameters**:
- `sequences`: str | List[str] | StandardizedOutput - Input sequences
- `output_format`: str = "csv" - Output format (csv, a3m)
- `timeout`: int = 3600 - Server timeout in seconds

**Outputs**:
- `msas` - DataStream of MSA files
- `tables.msas`: | id | sequence_id | sequence | msa_file |

**Example**:

```python
from PipelineScripts.mmseqs2 import MMseqs2

msas = MMseqs2(sequences=lmpnn, timeout=7200)
```

---

## Visualization

### PyMOL

Creates PyMOL sessions using a declarative operation-based API. 

**Resources**: Requires a GPU node.

**Environment**: `ProteinEnv`

**Default Settings** (applied automatically):
```
show cartoon
set seq_view, 1
set cartoon_gap_cutoff, 0
set sphere_scale, 0.2
set ray_trace_mode, 1
set ray_shadows, 0
set spec_reflect, 0
set ray_trace_frames, 1
set ray_trace_color, gray20
```

To override defaults, use `PyMOL.Set()` before other operations.

**Operations**:
- `Names(prefix, basename, suffix)` - Set ID → name mapping
- `Load(structures)` - Load structures
- `Color(structures, selection, color)` - Color selection
- `ColorAF(structures, upper=100)` - Color by pLDDT (AlphaFold scheme)
- `Align(method, target)` - Align objects (methods: "align", "super", "cealign")
- `Show(structures, representation, selection)` - Show representation
- `Hide(structures, representation, selection)` - Hide representation
- `Set(setting, value, selection)` - Set PyMOL setting
- `Save(filename)` - Save session
- `Render(structures, orient_selection, width, height, filename, dpi)` - Render single PNG
- `RenderEach(structures, ...)` - Render each structure individually as PNG

**Outputs**:
- `renders` - DataStream of PNG files (when using Render/RenderEach)
- `session_file` - Path to .pse session file

**RenderEach Parameters**:
- `structures`: Structures to render
- `orient_selection`: Selection to orient towards (default: "resn LIG")
- `color_protein`: "plddt" for confidence coloring, or color name (default: "plddt")
- `color_ligand`: "byatom" for element colors, or color name (default: "byatom")
- `ligand_selection`: Selection for ligand (default: "resn LIG")
- `plddt_upper`: Upper bound for pLDDT (default: 1 for Boltz2, use 100 for AlphaFold)
- `width`, `height`: Image dimensions (default: 1920x1080)
- `dpi`: Image DPI (default: 300)
- `background`: Background color (default: "white")

**Examples**:

```python
from PipelineScripts.pymol import PyMOL

# Basic session with alignment
PyMOL(session="my_session",
    PyMOL.Load(boltz),
    PyMOL.ColorAF(boltz),
    PyMOL.Align("align")
)

# Override default settings
PyMOL(session="custom",
    PyMOL.Set("ray_shadows", 1),
    PyMOL.Set("cartoon_gap_cutoff", 10),
    PyMOL.Load(structures),
    PyMOL.ColorAF(structures)
)

# Render each structure individually
PyMOL(session="renders",
    PyMOL.RenderEach(
        structures=boltz,
        orient_selection="resn LIG",
        color_protein="plddt",
        plddt_upper=1,  # Boltz2 uses 0-1 confidence
        width=1920,
        height=1080
    )
)

# AlphaFold structures (pLDDT 0-100)
PyMOL(session="af_renders",
    PyMOL.RenderEach(
        structures=alphafold,
        color_protein="plddt",
        plddt_upper=100
    )
)
```

### Plot

Creates publication-ready plots from CSV data.

**Environment**: None (uses matplotlib at SLURM time)

**Operations**:
- `Scatter(data, x, y, color)` - Scatter plot
- `Histogram(data, x)` - Distribution histogram
- `Bar(data, x, y)` - Bar chart
- `Column(data, y, labels, style)` - Prism-style column plot
- `HeatMap(data, columns)` - Correlation heatmap

**Column styles**: "column", "simple_bar", "scatter", "box", "floating_bar"

**Example**:

```python
from PipelineScripts.plot import Plot

Plot(
    Plot.Scatter(
        data=results.tables.analysis,
        x="pLDDT",
        y="affinity",
        color="source"
    )
)

Plot(
    Plot.Column(
        data=[cycle0.tables.conf, cycle1.tables.conf],
        y="plddt",
        labels=["Cycle 0", "Cycle 1"],
        style="box"
    )
)
```
