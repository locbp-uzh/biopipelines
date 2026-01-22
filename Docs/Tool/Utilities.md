# Utilities

[← Back to Tool Reference](../ToolReference.md)

---

### LoadOutput / LoadOutputs

Load previously saved tool outputs for reuse in new pipelines. Enables incremental development and filtering of existing results at pipeline runtime.

**Environment**: `ProteinEnv`

#### LoadOutput (single file)

Loads a single tool output JSON file.

**Parameters**:
- `output_json`: str (required) - Path to tool output JSON file (in ToolOutputs folder)
- `filter`: Optional[str] = None - Pandas query-style filter expression
- `validate_files`: bool = True - Check file existence when loading

**Outputs**:
- Same structure as original tool that created the output

**Example**:
```python
from PipelineScripts.load_output import LoadOutput

previous_boltz = LoadOutput(
    output_json="/path/to/job/ToolOutputs/003_Boltz2.json",
    filter="confidence_score > 0.8"
)
```

#### LoadOutputs (multiple files)

Loads multiple tool outputs from a ToolOutputs folder with filtering and sorting.

**Parameters**:
- `path`: str (required) - Path to the job folder or ToolOutputs folder
  - If job folder: automatically appends "ToolOutputs"
  - If ToolOutputs folder: uses directly
- `tool`: Optional[str] = None - Filter by tool name (e.g., "MergeTables", "Filter", "Boltz2")
- `suffix`: Optional[str] = None - Filter by filename suffix (e.g., "Cycle10", "Affinity")
  - Matches pattern: NNN_ToolName_Suffix.json
- `ascending`: bool = True - Sort order (True = ascending by name, False = descending)
- `**load_output_kwargs`: Additional parameters passed to LoadOutput (e.g., `filter`, `validate_files`)

**Returns**:
- Dictionary mapping identifiers to LoadOutput objects
  - Keys: `"{execution_order}_{tool_name}_{suffix}"` or `"{execution_order}_{tool_name}"`
  - Values: LoadOutput instances

**Example**:
```python
from PipelineScripts.load_output import LoadOutputs

# Load all MergeTables outputs from a job
data = LoadOutputs(
    path="/shares/user/BioPipelines/Project/Job_001",
    tool="MergeTables"
)
# Access outputs: data["003_MergeTables"], data["007_MergeTables"], etc.

# Load all outputs from Cycle10 with validation disabled
cycle10_data = LoadOutputs(
    path="/shares/user/BioPipelines/Project/Job_001/ToolOutputs",
    suffix="Cycle10",
    validate_files=False
)

# Load all Filter outputs in descending order
filters = LoadOutputs(
    path="/shares/user/BioPipelines/Project/Job_001",
    tool="Filter",
    ascending=False
)

# Combine filters: specific tool with suffix
merged_cycle5 = LoadOutputs(
    path="/shares/user/BioPipelines/Project/Job_001",
    tool="MergeTables",
    suffix="Cycle5"
)
```

**Practical Tips**:
- Instead of copy-pasting ligand SMILES across pipelines, you can create a compound library, and then load the smiles passing an id filter:
```python
from PipelineScripts.load_output import LoadOutput
from PipelineScripts.compound_library import CompoundLibrary
from PipelineScripts.boltz2 import Boltz2

# Pipeline 1 to store the library
# imports, pipeline instantiation, ...
CompoundLibrary({
  "Compound1": "CCNCNNCC(=O)C",
  "Compound2": "CCNCNNCC(=O)C",
  ...
})
# submit

# Pipeline 2 running calculations with one of the compounds
# imports, pipeline instantiation, ...
compound1 = LoadOutput(
    output_json="/path/to/job/ToolOutputs/001_CompoundLibrary.json",
    filter='id == "Compound1"'  # quotes are important for proper pandas query here: x is a column name; "x" is a string.
)
boltz = Boltz2(
    proteins=HaloTag,
    ligands=compound1
)
# submit
```

### MMseqs2

Generates multiple sequence alignments (MSAs) for protein sequences. Used for improving structure prediction quality by providing evolutionary information.

**Environment**: `ProteinEnv`

**Parameters**:
- `sequences`: Union[str, List[str], ToolOutput, StandardizedOutput] (required) - Input sequences
- `output_format`: str = "csv" - Output format (csv, a3m)
- `timeout`: int = 3600 - Timeout in seconds for server response

**Outputs**:
- `tables.msas`:

  | id | sequence_id | sequence | msa_file |
  |----|-------------|----------|----------|

**Example**:
```python
from PipelineScripts.mmseqs2 import MMseqs2

msas = MMseqs2(
    sequences=lmpnn,
    timeout=7200
)
```

---

### MMseqs2Server

Runs an MMseqs2 server for local MSA generation. Automatically started by MMseqs2 client when needed; manual setup typically not required.

**Environment**: None (doesn't require mamba)

**Parameters**:
- `port`: int = 8000 - Server port
- `host`: str = "0.0.0.0" - Server host
- `workers`: int = 4 - Number of worker processes

**Note**: This server tool must be run separately to provide MMseqs2 as a service. However, the MMseqs2 client automatically runs it when needed, so manual server setup is typically not required.

---

### CompoundLibrary

Creates and manages chemical compound libraries. Supports combinatorial SMILES expansion and optional generation of covalent binding files.

**Environment**: `ProteinEnv`

**Parameters**:
- `library`: Union[str, Dict[str, Union[str, List[str]]]] (required) - Dictionary with expansion keys or path to CSV file
- `primary_key`: Optional[str] = None - Root key for expansion when library is dictionary
- `covalent`: bool = False - Generate CCD/PKL files for covalent ligand binding
- `validate_smiles`: bool = True - Validate SMILES strings during expansion
- `conformer_method`: str = "UFF" - Method for conformer generation (UFF, OpenFF, DFT)

**Outputs**:
- `compounds`: CSV file with compound library
- `tables.compounds`:

  | id | format | smiles | ccd | {branching_keys} |
  |----|--------|--------|-----|------------------|

**Examples**:
```python
from PipelineScripts.compound_library import CompoundLibrary

# With primary_key for combinatorial expansion
library = CompoundLibrary(
    library={
        "scaffold": "<linker><fluorophore>",
        "linker": ["CCOCC", "CCOCCOCC"],
        "fluorophore": ["c1ccc(N)cc1"]
    },
    primary_key="scaffold"
)

# Without primary_key - direct SMILES list
library = CompoundLibrary(
    library={
        "compounds": ["CCO", "CCCO", "CCCCO", "c1ccccc1"]
    }
)
```

---

### PDB

Fetches protein structures with automatic fallback: checks local folders first, then downloads from RCSB PDB if not found locally. Downloads are automatically cached in PDBs/ folder for reuse.

**Environment**: `ProteinEnv`

**Parameters**:
- `pdbs`: Union[str, List[str]] (required) - PDB IDs to fetch (e.g., "4ufc" or ["4ufc","1abc"])
- `ids`: Optional[Union[str, List[str]]] - Custom IDs for renaming (defaults to pdbs)
- `format`: str = "pdb" - File format ("pdb" or "cif")
- `local_folder`: Optional[str] = None - Custom local folder to check first (before PDBs/)
- `biological_assembly`: bool = False - Download biological assembly from RCSB
- `remove_waters`: bool = True - Remove water molecules

**Fetch Priority**:
For each PDB ID, searches in order:
1. `local_folder` (if parameter provided)
2. `./PDBs/` folder in repository
3. Download from RCSB PDB (cached to PDBs/ for reuse)

**Outputs**:
- `structures`: List of structure files
- `tables.structures`:

  | id | pdb_id | file_path | format | file_size | source |
  |----|--------|-----------|--------|-----------|--------|

- `tables.sequences`:

  | id | sequence |
  |----|----------|

- `tables.failed`:

  | pdb_id | error_message | source | attempted_path |
  |--------|---------------|--------|----------------|

**Examples**:
```python
from PipelineScripts.pdb import PDB

# Automatic fallback: checks PDBs/, then downloads
pdb = PDB(
    pdbs=["4ufc", "1abc"],
    ids=["POI1", "POI2"]
)

# Check custom folder first, then PDBs/, then download
pdb = PDB(
    pdbs=["my_structure"],
    local_folder="/path/to/my/pdbs"
)

# Download biological assembly
pdb = PDB(
    pdbs="4ufc",
    ids="MyProtein",
    biological_assembly=True,
    remove_waters=False
)
```

---

### Ligand

Fetches small molecule ligands from RCSB PDB or PubChem. Downloads SDF files and converts to PDB format with proper atom numbering (sequential serial numbers, chain A, residue number 1). Supports lookup by CCD code (RCSB), compound name, CID, or CAS number (PubChem).

**Environment**: `biopipelines`

**Parameters**:
- `ids`: Union[str, List[str]] (required) - Output identifier(s) for filenames (e.g., "my_ligand" -> my_ligand.pdb)
- `codes`: Optional[Union[str, List[str]]] - 3-letter PDB residue code(s) to use in PDB file (e.g., "LIG"). Defaults to lookup value (truncated to 3 chars, uppercased).
- `lookup`: Optional[Union[str, List[str]]] - Lookup value(s) for fetching. Can be:
  - RCSB CCD codes: "ATP", "GDP", "HEM"
  - PubChem CID: "2244"
  - PubChem CAS: "50-78-2"
  - PubChem name: "aspirin", "caffeine"
  Defaults to codes if not provided (backward compatibility).
- `source`: Optional[str] = None - Force source ("rcsb" or "pubchem"). If None, auto-detects based on lookup format.
- `local_folder`: Optional[str] = None - Custom folder to check before Ligands/ cache

**Auto-Detection** (when source=None):
- 1-3 uppercase alphanumeric -> RCSB (CCD)
- Purely numeric -> PubChem (CID)
- XX-XX-X format -> PubChem (CAS)
- Otherwise -> PubChem (name)

**Fetch Priority**:
1. `local_folder` (if provided)
2. `Ligands/` folder (cached downloads)
3. Download from RCSB or PubChem (cached to Ligands/ for reuse)

**PDB Normalization**:
Downloaded ligands are normalized with:
- Sequential atom serial numbers starting from 1
- Chain ID = 'A'
- Residue number = 1
- Residue name = specified `codes` parameter

**Outputs**:
- `structures`: List of ligand PDB files
- `tables.compounds`:

  | id | code | lookup | source | ccd | cid | cas | smiles | name | formula | file_path |
  |----|------|--------|--------|-----|-----|-----|--------|------|---------|-----------|

- `tables.failed`:

  | lookup | error_message | source | attempted_path |
  |--------|---------------|--------|----------------|

**Examples**:
```python
from PipelineScripts.ligand import Ligand
from PipelineScripts.rfdiffusion3 import RFdiffusion3

# RCSB ligand by CCD code (auto-detect)
lig = Ligand(ids="atp_ligand", lookup="ATP")

# PubChem ligand by compound name
lig = Ligand(ids="aspirin_lig", lookup="aspirin", codes="ASP")

# PubChem by CID
lig = Ligand(ids="caffeine", lookup="2157", codes="CAF")

# PubChem by CAS number
lig = Ligand(ids="ibuprofen", lookup="15687-27-1", codes="IBU")

# Multiple ligands with explicit codes
lig = Ligand(
    ids=["lig1", "lig2"],
    lookup=["ATP", "aspirin"],
    codes=["ATP", "LIG"]
)

# Force PubChem source for a CCD-like code
lig = Ligand(ids="test", lookup="ATP", source="pubchem")

# Use in RFdiffusion3 workflow
ligand = Ligand(ids="my_ligand", lookup="2244", codes="LIG")
rfd3 = RFdiffusion3(input=ligand, contig="50-80,A1-100")
```

---

### PyMOL

Creates PyMOL sessions using a declarative operation-based API. Supports ID-based matching between structures and table columns for per-structure selections and naming.

**Environment**: `ProteinEnv` (or any with pymol-open-source installed)

**Operations**:
- `Names(prefix, basename, suffix)` - Set up ID → PyMOL name mapping
- `Load(structures)` - Load structures with current naming
- `Color(structures, selection, color)` - Color selection on structures
- `ColorAF(structures, upper=100)` - Color by AlphaFold pLDDT (B-factor), upper scales thresholds
- `Align(method, target)` - Align all loaded objects
- `Show(structures, representation, selection)` - Show representation (structures optional)
- `Hide(structures, representation, selection)` - Hide representation (structures optional)
- `Set(setting, value, selection)` - Set PyMOL setting
- `Save(filename)` - Save session (auto-saves if omitted)

These operations are accessible as static methods with  `PyMOL.<operation>` (see examples).

**Parameters**:
- `session`: str = "session" - Output session filename (without .pse)
- `*operations`: Sequence of PyMOL operations

**ID-based matching**:
- `Names()` creates an ID → PyMOL name mapping from a table column
- `Load()` uses this mapping to name objects (or uses ID if no mapping)
- `Color()` looks up selections per-structure from table columns by ID

**Outputs**:
- PyMOL session file (.pse)

**Example**:
```python
from PipelineScripts.pymol import PyMOL
from PipelineScripts.fuse import Fuse
from PipelineScripts.alphafold import AlphaFold

# Fuse domains with linkers
fused = Fuse(proteins=[A, B, C], linker="GSGAG", linker_lengths=["2-4", "2-4"])
folded = AlphaFold(sequences=fused)

# Create PyMOL session with domain coloring
PyMOL("domain_visualization",
    # Set naming: ID "sensor_2_3" (from the basename table) → PyMOL object "s_2_3_parts"
    PyMOL.Names(prefix="s", basename=fused.tables.sequences.lengths, suffix="parts"),
    PyMOL.Load(folded),  # folded IDs match fused table IDs
    # Color domains and linkers using per-structure selections from table
    PyMOL.Color(folded, selection=fused.tables.sequences.D1, color="white"),
    PyMOL.Color(folded, selection=fused.tables.sequences.L1, color="orange"),
    PyMOL.Color(folded, selection=fused.tables.sequences.D2, color="marine"),
    PyMOL.Color(folded, selection=fused.tables.sequences.L2, color="wheat"),
    PyMOL.Color(folded, selection=fused.tables.sequences.D3, color="white"),
    # Load again with pLDDT coloring under different names
    PyMOL.Names(prefix="sensor", basename=fused.tables.sequences.lengths, suffix="pLDDT"),
    PyMOL.Load(folded),
    PyMOL.ColorAF(folded),
    # Align all to first loaded object using align algorithm (sequence-based)
    PyMOL.Align("align")
)
```

---

### Plot

Creates publication-ready plots from CSV tables using a declarative operation-based API. Supports multiple plot types including scatter, histogram, bar, column (Prism-style), and heatmap.

**Environment**: None (uses matplotlib/seaborn at SLURM runtime)

**Operations**:
- `Scatter(data, x, y, ...)` - X vs Y scatter plot with optional color grouping
- `Histogram(data, x, ...)` - Distribution histogram with statistics
- `Bar(data, x, y, ...)` - Categorical bar chart
- `Column(data, y, ...)` - Prism-style column plot comparing multiple data sources
- `HeatMap(data, ...)` - Correlation matrix or pivot table heatmap

These operations are accessible as static methods with `Plot.<operation>` (see examples).

**Common Parameters** (all operations):
- `data`: Data source (ToolOutput, TableInfo, or list for Column)
- `title`: str = None - Plot title
- `xlabel`, `ylabel`: str = None - Axis labels
- `figsize`: Tuple[float, float] = (8, 6) - Figure size in inches

**Plot.Scatter Parameters**:
- `x`: str (required) - Column name for X axis
- `y`: str (required) - Column name for Y axis
- `color`: str = None - Column name for color grouping
- `x_tick_rotation`: float = 0 - X-axis label rotation in degrees
- `y_tick_rotation`: float = 0 - Y-axis label rotation in degrees
- `grid`: bool = True - Show grid lines
- `color_legend_loc`: str = "upper right" - Legend location
- `color_legend_outside`: bool = False - Place legend outside plot to the right

**Plot.Column Parameters**:
- `y`: str (required) - Column name to compare across sources
- `labels`: List[str] = None - X-axis labels for each data source
- `style`: str = "column" - Plot style:
  - `"column"` - Bars with error bars and overlaid scatter points (default)
  - `"simple_bar"` - Just bars with error bars
  - `"scatter"` - Just scattered points with mean line
  - `"box"` - Box and whiskers plot
  - `"floating_bar"` - Horizontal bars showing mean ± error
- `show_error`: str = "sd" - Error bar type ("sd", "sem", "ci", or None)
- `color_groups`: List[str] = None - Group names for coloring (same group = same color)
- `colors`: List[str] = None - Explicit color list (hex or named)
- `color_legend_title`: str = None - Title for color legend
- `color_legend_loc`: str = "upper right" - Legend location
- `color_legend_outside`: bool = False - Place legend outside plot to the right
- `x_tick_rotation`: float = 0 - X-axis label rotation in degrees
- `y_tick_rotation`: float = 0 - Y-axis label rotation in degrees
- `grid`: bool = True - Show grid lines

**Plot.HeatMap Parameters**:
- Pivot mode: `x`, `y`, `value` columns for pivot table
- Correlation mode: `columns` list for correlation matrix
- `cmap`: str = "viridis" - Colormap name
- `annotate`: bool = True - Show values in cells

**Outputs**:
- PNG plot files in output folder
- `tables.metadata`: Summary of generated plots

**Examples**:
```python
from PipelineScripts.plot import Plot
from PipelineScripts.load_output import LoadOutputs

# Load multiple AlphaFold outputs
folded = LoadOutputs(path="/path/to/job/", tool="AlphaFold")

# Extract metadata from output names
exps = [(k.split('_')[2], k.split('_')[3]) for k in folded.keys()]  # (protein, analyte)
proteins = [protein for protein, analyte in exps]
analytes = [analyte for protein, analyte in exps]
confidences = [v.tables.confidence for k, v in folded.items()]

# Column plot with color grouping and legend outside
Plot(
    Plot.Column(
        data=confidences,
        y="plddt",
        labels=analytes,              # X-axis shows analyte names
        color_groups=proteins,        # Color by protein
        color_legend_title="Protein",
        color_legend_outside=True,    # Legend outside plot area
        xlabel="Analyte",
        ylabel="pLDDT",
        x_tick_rotation=45,
        grid=False
    )
)

# Box plot style
Plot(
    Plot.Column(
        data=confidences,
        y="plddt",
        labels=analytes,
        style="box",                  # Box and whiskers
        color_groups=proteins,
        color_legend_title="Protein",
        color_legend_outside=True
    )
)

# Simple scatter style (points only)
Plot(
    Plot.Column(
        data=confidences,
        y="plddt",
        labels=analytes,
        style="scatter",
        color_groups=proteins
    )
)

# Scatter plot with color grouping from table column
Plot(
    Plot.Scatter(
        data=analysis.tables.results,
        x="pLDDT",
        y="affinity",
        color="source",
        title="pLDDT vs Affinity",
        grid=True
    )
)

# Correlation heatmap
Plot(
    Plot.HeatMap(
        data=results,
        columns=["pLDDT", "affinity", "contacts"],
        title="Metric Correlations"
    )
)
```

---
