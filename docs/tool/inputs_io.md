# Inputs & I/O

[← Back to Tool Reference](../tool_reference.md)

These tools bring data **into** a pipeline (structures, sequences, ligands, compound libraries, existing tables) and take results **out** (plots, PyMOL sessions, reloaded previous runs). The basic input entities are exported from `biopipelines.entities`:

```python
from biopipelines.entities import PDB, RCSB, Sequence, Ligand, CompoundLibrary
```

For MSA generation/conversion see [MSAs](msas.md); for molecule format conversion and descriptors see [Cheminformatics](cheminformatics.md).

---

### CompoundLibrary

Creates compound collections from dictionaries with optional combinatorial expansion.

**Environment**: `biopipelines`

**Parameters**:
- `library`: str | Dict - Dictionary with SMILES, path to CSV, or path to `.cdxml` file (ChemDraw R-group enumeration)
- `primary_key`: str = None - Root key for expansion. Auto-detected by scanning keys in insertion order for `<placeholder>` patterns.
- `covalent`: bool = False - Generate CCD/PKL files for covalent binding
- `validate_smiles`: bool = True - Validate SMILES during expansion
- `generate_images`: bool = False - Generate a PNG image per compound using RDKit (no extra dependencies)

**Streams**: `compounds`, `images` (if `generate_images=True`)

**Tables**:
- `compounds`: | id | format | smiles | ccd | {branching_keys} |

**Examples**:

```python
from biopipelines.entities import CompoundLibrary

# Simple dictionary (no expansion): keys are compound IDs, values are SMILES
library = CompoundLibrary({
    "aspirin": "CC(=O)OC1=CC=CC=C1C(=O)O",
    "caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
})

# With expansion using <key> placeholders — primary key auto-detected by order of appearance
library = CompoundLibrary({
    "scaffold": "<linker><fluorophore>",
    "linker": ["CCOCC", "CCOCCOCC"],
    "fluorophore": ["c1ccc(N)cc1"]
})
# Generates 2 compounds; branching columns 'linker' and 'fluorophore' track substituents

# From CDXML file (ChemDraw R-group enumeration)
# Fragment names defined in ChemDraw are used in branching columns
library = CompoundLibrary("my_library.cdxml")

# From CSV file (expansion supported if SMILES column contains <placeholders>)
library = CompoundLibrary("my_library.csv")

# With 2D molecule images
library = CompoundLibrary({...}, generate_images=True)
```

**CDXML R-Group Enumeration**:

Draw the following in a single ChemDraw `.cdxml` file:

- **Core scaffold**: the main molecule with **R1**, **R2**, etc. at substitution points (using ChemDraw's R-group substitution table)
- **R-group fragments**: defined in the R-group table, one per substitution point

To assign display names to substituents (shown in branching columns instead of SMILES), add chemical property labels to the fragment structures in ChemDraw (Structure > Chemical Properties, or use text labels linked to fragments via the ChemDraw property tool).

The output CSV includes branching columns (R1, R2, ...) showing the substituent name or SMILES at each position.

Requires RDKit (`conda install -c conda-forge rdkit`).

---

### Ligand

Fetches small molecules from RCSB (CCD) or PubChem (name, CID, CAS) or generates from SMILES.

**Environment**: `biopipelines`

**Parameters**:
- `lookup`: str | List[str] | Dict[str, str] = None - Lookup values (CCD, CID, CAS, name), or path to a `.txt` file (one SMILES per line) or `.cdxml` file (ChemDraw molecules). File types are auto-detected by extension.
- `ids`: str | List[str] = None - Custom IDs (defaults to lookup / file-derived names / "smilesN")
- `codes`: str | List[str] = None - Residue codes carried on the compounds stream, 1-5 alphanumeric (extended CCD); defaults to the lookup value or "LIG"
- `code`: str | List[str] = None - Name an existing HETATM residue by its code (1-5 alphanumeric) **without** downloading or generating any structure (compounds-only). Hand the result to tools that only need the residue code (LigandMPNN, PoseBusters, PLIP, …).
- `source`: str = None - Force "rcsb" or "pubchem" (auto-detects if None)
- `local_folder`: str = None - Check first before ligands/
- `smiles`: str | List[str] | Dict[str, str] = None - Direct SMILES input (bypasses lookup)
- `structures`: DataStream | StandardizedOutput = None - Extract the ligand's bound coordinates from these structures (keeps the crystal pose) while taking chemistry from `code`/`smiles`. See [The Ligand Contract](../developer_manual.md#the-ligand-contract-compounds--chemistry-structures--coordinates).
- `generate_images`: bool = False - Generate PNG images per ligand using RDKit

**Auto-detection** (when source=None):
- 1-5 uppercase alphanumeric → RCSB (CCD)
- Purely numeric → PubChem (CID)
- XX-XX-X format → PubChem (CAS)
- Otherwise → PubChem (name)

**Streams**: `structures`, `compounds`, `images` (if `generate_images=True`)

The `structures` stream is **SDF** for fetched/generated ligands (download, PubChem, SMILES) and the input structure's own format (pdb/cif) for the `structures=` HETATM-extract path; it is absent in `code`-only mode. `Ligand` does not write PDB/CIF coordinates itself — for those, convert afterwards with `OpenBabel(compounds=lig, convert_3d="pdb"|"cif")`. The residue code (1-5 chars, extended-CCD capable) lives on the `compounds` stream and is independent of the coordinate file's format.

**Tables**:
- `compounds`: | id | code | lookup | source | smiles | formula |

**Examples**:

```python
from biopipelines.entities import Ligand

# RCSB by CCD code
atp = Ligand("ATP")

# PubChem by name
aspirin = Ligand("aspirin", codes="ASP")

# Direct SMILES
ethanol = Ligand(smiles="CCO", ids="ethanol", codes="ETH")

# From a .txt file (one SMILES per line): IDs become myligands1, myligands2, ...
ligands = Ligand("myligands.txt")

# From a .cdxml file: each molecule becomes a ligand; ChemDraw names used as IDs
ligands = Ligand("my_ligands.cdxml")
```

A `code`-only Ligand (`Ligand(code="ZIT")`) names an existing HETATM residue without downloading or generating any structure — hand it to tools that only need the ligand's residue code (LigandMPNN, PoseBusters, PLIP, …). To turn a SMILES-carrying Ligand into a 3-D file, run [OpenBabel](cheminformatics.md#openbabel).

---

### Load / LoadMultiple {#load}

Reload previously produced pipeline outputs back into a run.

**Load** — loads a single tool output JSON file.

**Parameters**:
- `path`: str - Path to the tool output JSON file
- `filter`: str = None - Pandas query filter
- `validate_files`: bool = True - Check file existence

```python
from biopipelines.load import Load

prev = Load(
    path="/path/to/ToolOutputs/003_Boltz2.json",
    filter="confidence_score > 0.8"
)
```

**LoadMultiple** — loads multiple outputs from a ToolOutputs folder. Returns a dict mapping identifiers to `Load` objects.

**Parameters**:
- `path`: str - Path to job folder or ToolOutputs folder
- `tool`: str = None - Filter by tool name
- `suffix`: str = None - Filter by filename suffix
- `in_suffix`: str | List[str] = None - Keep only outputs whose suffix contains this (any of these) substring(s)
- `not_in_suffix`: str | List[str] = None - Drop outputs whose suffix contains this (any of these) substring(s)
- `ascending`: bool = True - Sort order

```python
from biopipelines.load import LoadMultiple

# All Boltz2 outputs
all_boltz = LoadMultiple("/path/to/job/", tool="Boltz2")
# Access: all_boltz["003_Boltz2"], all_boltz["007_Boltz2"], etc.

# Filter by suffix
cycle10 = LoadMultiple("/path/to/job/", suffix="Cycle10")
```

---

### PDB

Fetches protein structures with priority: `local_folder` → `pdbs/` → RCSB download.

**Environment**: `biopipelines`

**Parameters**:
- `pdbs`: str | List[str] - PDB IDs or folder path
- `ids`: str | List[str] = None - Custom IDs (defaults to pdbs)
- `convert`: str = None - Target format to convert to ("pdb", "cif", or None). When None, no conversion is performed: structures are kept in whatever format they are found locally or downloaded from RCSB. The structures stream will have format "pdb|cif".
- `local_folder`: str = None - Check first before pdbs/
- `biological_assembly`: bool = False - Download biological assembly
- `remove_waters`: bool = True - Remove water molecules
- `chain`: str | List[str] = "auto" - Which chain(s) to extract the sequence from.
  - `"auto"` (default): emit one row `<id>,<sequence>` per input — the longest chain of the structure (or the longest entity returned by RCSB FASTA when the input is a PDB code). The structure file on disk keeps every chain.
  - `"all"`: emit one row `<id>_<chain_letter>,<sequence>` per chain present in the structure (sequences pulled from the RCSB FASTA endpoint when the input is a PDB code, otherwise from the structure file). No aggregate longest row — each id identifies a single chain unambiguously.
  - `["A", "C", ...]`: explicit list of chain letters. Emits one row `<id>_<chain_letter>,<sequence>` per listed chain (cardinality known at config time). Useful when the structure has many chains but only a few are biologically relevant.
  - `"A"` / `"B"` / ... : single chain. The structure file on disk is filtered to that chain only and a single row `<id>,<sequence_for_that_chain>` is emitted.
- `split_chains`: bool = False - When True (and `chain` is `"all"` or a list), additionally split the structure file into one `.pdb` (or `.cif`) per chain letter at `<output>/<custom_id>_<letter>.<ext>`, with the structures stream IDs becoming `<custom_id>_<letter>`. Mutually exclusive with the single-chain forms (`"auto"` / explicit chain letter).
- `fetch_compounds`: bool = True - Fetch SMILES for ligands from RCSB. Set to False when compounds are not needed to avoid per-ligand API calls.

**Operations** (optional positional args):
- `PDB.rename(old, new)` - Rename residues (e.g., for RFdiffusion3 compatibility)
- `PDB.remove(selection, remove_hetatm=True)` - Remove residues by selection. `selection` is a PyMOL-style string (`"1-83"`, `"1-83+90-95"`, `"A1-83"`), a table column reference (`tool.tables.structures.col`), or an equivalent `(TableInfo, "col")` tuple. A column reference is resolved per-structure by ID match (single-row table broadcasts to all; unmatched structures are left unchanged) — so each structure can be truncated by its own value from an upstream table. By default (`remove_hetatm=True`) a HETATM (ligand/ion) numbered inside the selected range is removed along with the protein residues. A bound ligand survives truncation when it sits *outside* the removed range — that's the usual "dock against the full protein, then drop a segment" case. **If a ligand shares a residue number with the removed span (e.g. a dye is residue 1 on chain B while you drop protein residues 1-83), chain-qualify the selection — `PDB.remove("A1-83")` — so only the protein chain is touched and the ligand is left intact.** Alternatively set `remove_hetatm=False` to drop only protein ATOM records and keep all HETATM.
- `PDB.break_bond(atom1, atom2)` - Sever a bond between two atoms (PyMOL `unbond` style), removing any CONECT/LINK records joining them; coordinates untouched. Atoms use `<residue>.<atom>` syntax (`"A145.SG"`, `"LIG.C12"`). Use it to detach a covalently-tethered ligand after it has served its purpose — e.g. break the BG–Cys145 thioether on PLACER's covalent conformers so downstream design (RFdiffusion) sees the dye as a separate non-covalent HETATM.
- `PDB.rotate_bond(atom1, atom2, angle)` - Rotate the fragment on `atom2`'s side of the `atom1`–`atom2` bond by `angle` degrees about the bond axis; `atom1`'s side stays fixed. The moving fragment is found by graph connectivity *within `atom2`'s residue* (atoms reachable from `atom2` without crossing the bond), so the protein backbone is never disturbed and only the chosen rotatable bond's downstream atoms move. If the bond lies in a ring (`atom1` reachable from `atom2` another way) the op refuses and leaves coordinates untouched. Atoms use `<residue>.<atom>` syntax. Use it to re-orient a flexible ligand moiety while keeping a covalent anchor fixed — e.g. swing a dye's fluorophore toward a target region by rotating ~180° about a linker bond.

**Streams**: `structures`, `sequences`, `compounds`

**Tables**:
- `structures`: | id | pdb_id | file_path | format | source |
- `sequences`: | id | sequence |
- `compounds`: | id | code | smiles | ccd |

**Stream cardinality vs `chain` / `split_chains`**:

| Setting | `structures` rows / input | `sequences` rows / input |
|---|---|---|
| `chain="auto"` (default) | 1 | 1 (longest) |
| `chain="A"` (or any letter) | 1 (filtered to that chain) | 1 (that chain) |
| `chain="all"`, `split_chains=False` | 1 (full structure, all chains) | one per chain in the structure (lazy ids) |
| `chain="all"`, `split_chains=True` | one per chain in the structure (lazy ids) | one per chain in the structure (lazy ids) |
| `chain=["A","C"]`, `split_chains=False` | 1 (full structure, all chains) | 2 (`<id>_A`, `<id>_C`, literal ids) |
| `chain=["A","C"]`, `split_chains=True` | 2 (`<id>_A`, `<id>_C`, literal ids) | 2 (`<id>_A`, `<id>_C`, literal ids) |

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
pdb = PDB("structure.pdb", PDB.rename("LIG", ":L:"))

# Biological assembly
pdb = PDB("4ufc", biological_assembly=True, remove_waters=False)

# Per-chain sequences, single structure file (e.g. for ProteinMPNN over a complex)
hba = PDB("1A3N", chain="all")           # 4 sequence rows: 1A3N_A, 1A3N_B, 1A3N_C, 1A3N_D

# Per-chain sequences AND per-chain structure files
ake_chains = PDB("4AKE", chain="all", split_chains=True)
                                          # structures: 4AKE_A.pdb, 4AKE_B.pdb
                                          # sequences:  4AKE_A, 4AKE_B

# Subset of chains by explicit list (e.g. drop a non-biological copy)
hba_ab = PDB("1A3N", chain=["A", "B"])    # sequences: 1A3N_A, 1A3N_B (literal ids)
```

---

### Plot

Creates plots from CSV data.

**Environment**: None (uses matplotlib at execution time)

**Operations**:
- `Scatter(data, x, y, color)` - Scatter plot
- `Histogram(data, x)` - Distribution histogram
- `Bar(data, x, y)` - Bar chart
- `Column(data, y, labels, style)` - Prism-style column plot
- `Line(data, x, y, labels, color)` - Line plot for evolution/trends
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

Plot(
    Plot.Line(
        data=[best1.tables.result, best2.tables.result, best3.tables.result],
        x="iteration",
        y="affinity_probability_binary",
        labels=["best1", "best2", "best3"],
        title="Affinity Evolution"
    )
)
```

---

### PyMOL

Creates PyMOL sessions using a declarative operation-based API.

> **⚠ Ray-traced rendering is VERY SLOW.** With the default settings (`ray_trace_mode, 1`) every PNG is ray-traced at ~5–10 s each. `RenderEach` renders one image *per structure*, so on a few hundred–thousand structures it takes **hours** and can silently blow a short batch walltime. Only render a small, curated set (e.g. the top-N winners after scoring), never a full design pool. To inspect many structures cheaply, save a session (`Save`) and open it interactively instead.

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
- `RenderEach(structures, ...)` - Render each structure individually as PNG. **SLOW** (~5–10 s ray-traced per structure); pass only a small curated set, not a full design pool — see the warning above.

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

---

### RCSB

Searches the RCSB PDB Search API v2 and downloads matching structures.

**Environment**: `biopipelines`

**Parameters**:
- `*queries`: One or more query objects (combined with logical_operator)
- `max_results`: int = 10 - Maximum PDB entries to return
- `return_type`: str = "entry" - Result granularity ("entry", "assembly", "polymer_entity", "polymer_instance")
- `sort`: str = "score" - Sort field ("score", "resolution", "release_date", or RCSB attribute)
- `convert`: str = None - Target format to convert to ("pdb", "cif", or None). When None, no conversion is performed: each structure is kept in whatever format is downloaded from RCSB (PDB if available, CIF as fallback). The structures stream will have format "pdb|cif".
- `ids`: str | List[str] = None - Custom IDs (defaults to PDB IDs)
- `remove_waters`: bool = True - Remove water molecules
- `chain`: str = "auto" - Which chain to extract the sequence from. Same semantics as `PDB.chain` (`"auto"` longest only, `"all"` per-chain rows, or an explicit letter).
- `fetch_compounds`: bool = True - Fetch SMILES for ligands from RCSB. Set to False when compounds are not needed to avoid per-ligand API calls.
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
- `search_results`: | id | pdb_id | result_id | score | title | resolution | method | molecular_weight_kda | organism | entity_description | protein_entity_count | residue_count | citation_title | citation_journal | citation_year | citation_authors | release_date | deposit_date |
- `missing`: | pdb_id | error_message | source | attempted_path |

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

### Table

Direct table construction from an existing CSV or Excel file. Loads the file at config time, exposes it to downstream tools through the standard `tables.<name>` interface, and (for `.xlsx`/`.xls` inputs) writes a sibling CSV so the run stays CSV-native.

**Environment**: `biopipelines` (no extra installation).

**Parameters**:
- `path`: str (required) — Path to a `.csv`, `.xlsx`, or `.xls` file.
- `name`: str = `"data"` — Logical table name; downstream tools reference it as `<step>.tables.<name>`.
- `description`: str = `""` — Free-form description captured in provenance.

**Tables**:
- `<name>`: the loaded table, columns preserved verbatim.

**Example**:
```python
from biopipelines.table import Table

# CSV input → directly available as tbl.tables.metrics
tbl = Table("metrics.csv", name="metrics", description="Per-design metrics")

# Use a column as a per-structure reference for a downstream tool
ProteinMPNN(structures=proteins, redesigned=(tbl.tables.metrics, "designed_positions"))
```

---

### UniProt

Fetches sequences and annotations from the UniProt REST API. Takes one or more UniProt accessions and produces a `sequences` stream (from the canonical FASTA endpoint) plus an annotations table with organism, length, GO terms, Pfam domains, EC number, and review status — ready for downstream filtering or design.

**Environment**: `biopipelines`

**Parameters**:
- `accessions`: str | List[str] (required) — One or more UniProt accessions (e.g. `"P69905"`).

**Streams**: `sequences`

**Tables**:
- `sequences`: | id | sequence | type | length |
- `annotations`: | id | name | organism | taxonomy_id | length | go_terms | pfam | ec_number | reviewed |

**Example**:
```python
from biopipelines.uniprot import UniProt

# Fetch human hemoglobin alpha and feed straight into a predictor
hba = UniProt("P69905")
af = AlphaFold(proteins=hba)

# Multiple accessions
kinases = UniProt(["P00533", "P04626", "P06239"])
```
