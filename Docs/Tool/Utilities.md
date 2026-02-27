# Utilities

[← Back to Tool Reference](../ToolReference.md)

---

## Entity Types

Basic input types exported from `biopipelines.entities`:

```python
from biopipelines.entities import PDB, RCSB, Sequence, Ligand, CompoundLibrary
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
- `chain`: str = "longest" - Which chain to extract sequence from. "longest" (default) picks the longest chain. Specify a chain letter (e.g. "A", "B") to select that chain.

**Operations** (optional positional args):
- `PDB.Rename(old, new)` - Rename residues (e.g., for RFdiffusion3 compatibility)

**Streams**: `structures`, `sequences`, `compounds`

**Tables**:
- `structures`: | id | pdb_id | file_path | format | source |
- `sequences`: | id | sequence |
- `compounds`: | id | code | smiles | ccd |

**Examples**:

```python
from biopipelines.entities import PDB

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

### RCSB

Searches the RCSB PDB Search API v2 and downloads matching structures.

**Environment**: `biopipelines`

**Parameters**:
- `*queries`: One or more query objects (combined with logical_operator)
- `max_results`: int = 10 - Maximum PDB entries to return
- `return_type`: str = "entry" - Result granularity ("entry", "assembly", "polymer_entity", "polymer_instance")
- `sort`: str = "score" - Sort field ("score", "resolution", "release_date", or RCSB attribute)
- `format`: str = "pdb" - Output format ("pdb" or "cif")
- `ids`: str | List[str] = None - Custom IDs (defaults to PDB IDs)
- `remove_waters`: bool = True - Remove water molecules
- `chain`: str = "longest" - Which chain to extract sequence from
- `logical_operator`: str = "and" - How to combine queries ("and" or "or")

**Query Types**:
- `RCSB.Text(value)` - Full-text search (supports +, |, -, quoted phrases)
- `RCSB.Sequence(sequence, identity_cutoff, evalue_cutoff, sequence_type)` - BLAST-like similarity
- `RCSB.SeqMotif(pattern, pattern_type, sequence_type)` - Motif search (prosite/simple/regex)
- `RCSB.Structure(entry_id, assembly_id)` - 3D structure similarity
- `RCSB.StrucMotif(residues, pdb_id, assembly_id)` - Structure motif search
- `RCSB.Chemical(value, match_type, descriptor_type)` - Chemical similarity (by SMILES/InChI)
- `StructureAttribute.*` - Typed attribute search (see below)

**StructureAttribute** — typed descriptors organised by RCSB category. Each attribute exposes comparison methods that return a query object:

| Method | Description |
|--------|-------------|
| `.equals(value)` | Exact match. Accepts an enum member or plain value. |
| `.is_any(*values)` | Value is one of the options (OR). |
| `.contains(text)` | Contains words (case-insensitive). |
| `.contains_phrase(text)` | Contains exact phrase. |
| `.between(low, high)` | Numeric or date range (inclusive by default). |
| `.greater_than(value, or_equal=False)` | Greater than (or equal). |
| `.less_than(value, or_equal=False)` | Less than (or equal). |
| `.exists()` | Attribute is not null. |
| `.not_equals(value)` | Negated exact match. |
| `.not_any(*values)` | Negated "in" operator. |

**StructureAttribute categories and attributes**:

`StructureAttribute.IDsAndKeywords`
- `EntryId`, `StructureKeywords`, `AdditionalKeywords`

`StructureAttribute.StructureDetails`
- `StructureTitle`, `DepositDate`, `ReleaseDate`, `RevisionDate`

`StructureAttribute.EntryFeatures`
- `MolecularWeight`, `NumberOfProteinEntities`, `NumberOfRnaEntities`, `NumberOfDnaEntities`
- `TotalPolymerResidues`, `EntryPolymerComposition`, `EntryPolymerTypes`

`StructureAttribute.PolymerMolecularFeatures`
- `PolymerEntityDescription`, `PolymerEntityType`, `PolymerEntitySequenceLength`
- `ScientificName`, `TaxonomyId`, `GeneName`, `EnzymeClassification`

`StructureAttribute.NonpolymerFeatures`
- `ComponentIdentifier`, `FormulaWeight`, `LigandQscore`

`StructureAttribute.AssemblyFeatures`
- `SymmetryType`, `SymmetrySymbol`, `OligomericState`

`StructureAttribute.Methods`
- `ExperimentalMethod`, `Resolution`

**Value enums** (for use with `.equals()` / `.is_any()`):
- `PolymerEntityType` — `PROTEIN`, `DNA`, `RNA`, `NA_HYBRID`, `OTHER`
- `SymmetryType` — `ASYMMETRIC`, `HOMO_2_MER`, `HOMO_3_MER`, `HOMO_4_MER`, `HOMO_5_MER`, `HOMO_6_MER`, `HOMO_7_MER`, `PSEUDO_SYMMETRIC`, `HETEROMERIC`
- `ExperimentalMethod` — `X_RAY_DIFFRACTION`, `ELECTRON_MICROSCOPY`, `SOLUTION_NMR`, `SOLID_STATE_NMR`, `NEUTRON_DIFFRACTION`, `ELECTRON_CRYSTALLOGRAPHY`, `FIBER_DIFFRACTION`, `EPR`, `FLUORESCENCE_TRANSFER`, `INFRARED_SPECTROSCOPY`

**Streams**: `structures`, `sequences`, `compounds`

**Tables**:
- `structures`: | id | pdb_id | file_path | format | source |
- `sequences`: | id | sequence |
- `compounds`: | id | code | smiles | ccd |
- `search_results`: | id | pdb_id | result_id | score |
- `entry_info`: | id | pdb_id | title | resolution | method | molecular_weight_kda | organism | entity_description | protein_entity_count | residue_count | citation_title | citation_journal | citation_year | citation_authors | release_date | deposit_date |

**Examples**:

```python
from biopipelines.entities import RCSB
from biopipelines.rcsb import StructureAttribute, PolymerEntityType, SymmetryType, ExperimentalMethod

# Typed attribute search (recommended)
results = RCSB(
    StructureAttribute.NonpolymerFeatures.ComponentIdentifier.is_any("FMN"),
    StructureAttribute.PolymerMolecularFeatures.PolymerEntityType.equals(PolymerEntityType.PROTEIN),
    StructureAttribute.AssemblyFeatures.SymmetryType.equals(SymmetryType.ASYMMETRIC),
    StructureAttribute.PolymerMolecularFeatures.PolymerEntitySequenceLength.between(80, 300),
    max_results=100,
)

# Text search
results = RCSB(
    RCSB.Text("insulin receptor"),
    max_results=10
)

# Resolution + organism filter
results = RCSB(
    StructureAttribute.PolymerMolecularFeatures.ScientificName.equals("Homo sapiens"),
    StructureAttribute.Methods.Resolution.less_than(2.0),
    max_results=20
)

# Sequence similarity search
results = RCSB(
    RCSB.Sequence("MKTVRQERLKSIVRILERSKEPVSGAQ", identity_cutoff=0.9),
    max_results=50
)

# Structure similarity
results = RCSB(
    RCSB.Structure("4HHB", assembly_id=1),
    max_results=10
)

# Sequence motif (PROSITE format)
results = RCSB(
    RCSB.SeqMotif("C-x(2,4)-C-x(3)-[LIVMFYWC]-x(8)-H-x(3,5)-H", pattern_type="prosite"),
    max_results=100
)

# Chemical similarity (by SMILES)
results = RCSB(
    RCSB.Chemical("c1ccc(cc1)C(=O)O", match_type="fingerprint-similarity"),
    max_results=10
)

# Mixed: text + typed attributes
results = RCSB(
    RCSB.Text("kinase"),
    StructureAttribute.Methods.Resolution.less_than(2.5),
    StructureAttribute.PolymerMolecularFeatures.ScientificName.equals("Homo sapiens"),
    max_results=50,
    sort="resolution"
)

# Use downstream
af = AlphaFold(proteins=results)
```

---

### Sequence

Creates sequences from strings with auto-detection (protein/DNA/RNA).

**Environment**: `biopipelines`

**Parameters**:
- `seq`: str | List[str] - Sequence string(s)
- `type`: str = "auto" - Sequence type ("auto", "protein", "dna", "rna")
- `ids`: str | List[str] = None - Custom IDs (defaults to "seq_N")

**Streams**: `sequences`

**Tables**:
- `sequences`: | id | sequence | type | length |

**Examples**:

```python
from biopipelines.entities import Sequence

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

**Streams**: `structures`, `compounds`

**Tables**:
- `compounds`: | id | code | lookup | source | smiles | formula |

**Examples**:

```python
from biopipelines.entities import Ligand

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

**Streams**: `compounds`

**Tables**:
- `compounds`: | id | format | smiles | ccd | {branching_keys} |

**Examples**:

```python
from biopipelines.entities import CompoundLibrary

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

### Load

Loads a single tool output JSON file.

**Parameters**:
- `output_json`: str - Path to JSON file
- `filter`: str = None - Pandas query filter
- `validate_files`: bool = True - Check file existence

**Example**:

```python
from biopipelines.load import Load

prev = Load(
    output_json="/path/to/ToolOutputs/003_Boltz2.json",
    filter="confidence_score > 0.8"
)
```

### LoadMultiple

Loads multiple outputs from a ToolOutputs folder.

**Parameters**:
- `path`: str - Path to job folder or ToolOutputs folder
- `tool`: str = None - Filter by tool name
- `suffix`: str = None - Filter by filename suffix
- `ascending`: bool = True - Sort order

**Returns**: Dictionary mapping identifiers to Load objects

**Example**:

```python
from biopipelines.load import LoadMultiple

# All Boltz2 outputs
all_boltz = LoadMultiple("/path/to/job/", tool="Boltz2")

# Access: all_boltz["003_Boltz2"], all_boltz["007_Boltz2"], etc.

# Filter by suffix
cycle10 = LoadMultiple("/path/to/job/", suffix="Cycle10")
```

---

## MSA Generation

### MMseqs2

Generates multiple sequence alignments for structure prediction.

**Environment**: `biopipelines`

**Parameters**:
- `sequences`: str | List[str] | StandardizedOutput - Input sequences
- `output_format`: str = "csv" - Output format (csv, a3m)
- `timeout`: int = 3600 - Server timeout in seconds

**Streams**: `msas`

**Tables**:
- `msas`: | id | sequences.id | sequence | msa_file |

**Example**:

```python
from biopipelines.mmseqs2 import MMseqs2

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
- `ColorAlign(reference, targets, ...)` - Color by sequence alignment against reference
- `Align(method, target)` - Align objects (methods: "align", "super", "cealign")
- `Show(structures, representation, selection)` - Show representation
- `Hide(structures, representation, selection)` - Hide representation
- `Set(setting, value, selection)` - Set PyMOL setting
- `Center(selection="all")` - Center view on selection without changing orientation
- `Zoom(selection="all", buffer=5.0)` - Zoom to fit selection with buffer (Angstroms)
- `Orient(selection="all")` - Orient view to show selection from best angle
- `Save(filename)` - Save session
- `Render(structures, orient_selection, width, height, filename, dpi)` - Render single PNG
- `RenderEach(structures, ...)` - Render each structure individually as PNG

**Streams**: `images` (when using Render/RenderEach)

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

**ColorAlign Parameters**:
- `reference`: Reference structure for sequence alignment
- `targets`: Target structure(s) to color based on alignment
- `identical`: Color for identical residues (default: "white")
- `similar`: Color for similar amino acid groups (default: "wheat")
- `different`: Color for non-similar mismatches (default: "wheat")
- `notcovered`: Color for gaps/unaligned regions (default: "gray50")
- `show_mutations`: Show sticks for mutated residues (default: True)

**Amino Acid Groups** (for similar coloring):
- Aliphatic: A, V, L, I, M
- Aromatic: F, W, Y
- Polar: S, T, N, Q
- Positive: K, R, H
- Negative: D, E
- Special: C, G, P

**Examples**:

```python
from biopipelines.pymol import PyMOL

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

# Color by sequence alignment against reference
PyMOL(session="alignment_view",
    PyMOL.ColorAlign(
        reference=xrc_structure,
        targets=designed_sequences,
        identical="white",
        similar="palegreen",
        different="salmon",
        notcovered="gray50",
        show_mutations=True
    )
)
```

### Plot

Creates plots from CSV data.

**Environment**: None (uses matplotlib at execution time)

**Operations**:
- `Scatter(data, x, y, color)` - Scatter plot
- `Histogram(data, x)` - Distribution histogram
- `Bar(data, x, y)` - Bar chart
- `Column(data, y, labels, style)` - Prism-style column plot
- `HeatMap(data, columns)` - Correlation heatmap

**Column styles**: "column", "simple_bar", "scatter", "box", "floating_bar"

**Example**:

```python
from biopipelines.plot import Plot

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
