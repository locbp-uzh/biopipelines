# BioPipelines User Manual

## Index

- [What is BioPipelines?](#what-is-biopipelines)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Core Concepts](#core-concepts)
  - [Pipeline and SLURM Runtime](#pipeline-and-slurm-runtime)
  - [Entity Types](#entity-types)
  - [DataStream vs Tables](#datastream-vs-tables)
  - [Combinatorics: Bundle and Each](#combinatorics-bundle-and-each)
  - [Table Column References](#table-column-references)
- [Resources](#resources)
- [On-the-fly Execution](#on-the-fly-execution)
- [Job Submission](#job-submission)
- [Data Management with Panda](#data-management-with-panda)
- [Filesystem Structure](#filesystem-structure)
- [Troubleshooting](#troubleshooting)

---

## What is BioPipelines?

BioPipelines is a Python framework for writing bioinformatics workflows that  encompasses bash and script orchestration and slurm submission. It does not execute computations directly - instead, it predicts the filesystem structure and creates scripts that will be executed on SLURM clusters. 

BioPipelines was designed to maximize pipeline clarity and conciseness, as shown in the following example:

```python
#imports omitted
with Pipeline(project="Examples",
              job="RFD-ProteinMPNN-AlphaFold2",
              description="Redesign of N terminus domain of lysozyme"):
    Resources(gpu="A100", 
              time="4:00:00",
              memory="16GB")
    lysozyme = PDB("168L")
    rfd = RFdiffusion(pdb=lysozyme,
                        contigs='50-70/A81-140', #redesign N terminus
                        num_designs=3)
    pmpnn = ProteinMPNN(structures=rfd, 
                      num_sequences=2)
    af = AlphaFold(proteins=pmpnn)
```

---

## Installation

Clone the repository and set up the environment:

```bash
# 1. Clone the repository
git clone https://github.com/locbp-uzh/biopipelines
cd biopipelines

# 2. Create the biopipelines conda environment
mamba env create -f Environments/biopipelines.yaml
mamba activate biopipelines

# 3. Install the package (editable mode, so updates via git pull take effect immediately)
pip install -e .

# 4. Optional, for usage in Jupyter notebooks
ipython kernel install --user --name biopipelines
```

Edit config.yaml to fit your cluster configuration.

---

## Quick Start

Navigate to the ExamplePipelines folder and, after activating the biopipelines environment, run:

```bash
biopipelines-submit <example pipeline name>.py    # from a .py script
biopipelines-submit <example pipeline name>.ipynb  # from a Jupyter notebook
```

Note: you must have installed the relevant environments and configured them in config.yaml

---

## Core Concepts

### Pipeline and SLURM Runtime

BioPipelines operates alternatively in one or two phases:

One phase (Jupyter/Colab notebooks):
| Phase | What Happens | Executer | Where |
|-------|--------------|-------|-------|
| **Notebook cell** | Generation and running of bash scripts, prediction of output paths and files | Python + bash | Cluster |

Two phases (biopipelines-submit):
| Phase | What Happens | Executer | Where |
|-------|--------------|-------|-------|
| **Pipeline** | Generation of bash scripts, prediction of output paths and files | Python | Cluster |
| **SLURM** | Bash scripts execute, files are created | Slurm | Cluster |

### DataStream vs Tables

Tools output two types of data containers:

**DataStream** - Unified container supporting ID tracking, and association of IDs to files (e.g. .pdb, .cif, .sdf) or values (e.g. protein/dna sequences):

```python
# Access structures from a tool via the streams container
for struct_id, pdb_path in boltz.streams.structures:
    print(f"{struct_id}: {pdb_path}")

# Count items
print(f"Generated {len(boltz.streams.structures)} structures")
print(f"Expected ids: {boltz.streams.structures.ids}")
```

DataStream types accessed via `tool.streams.<name>`:
- `streams.structures` - PDB/CIF files
- `streams.sequences` - ID-tracked table with id, sequence columns
- `streams.compounds` - CSV with SMILES column, or SDF/PDB files
- `streams.msas` - A3M or CSV files
- `streams.images` - PNG files

**Tables (TableInfo)** - Rich metadata about CSV files. They do not track IDs.

```python
# Access table metadata via .info
info = tool.tables.confidence.info
print(info.path)        # /path/to/confidence.csv
print(info.columns)     # ["id", "pTM", "complex_plddt", ...]
print(info.description) # "Confidence scores"
print(info.count)       # Number of rows

# Access column references for downstream tools
tool.tables.confidence.plddt  # Returns (TableInfo, "plddt") tuple
```

In the ToolReference, one can find for each tool what is the expected output in terms of streams and tables, and use this information to write pipelines.

### Common inputs

Basic input types can be imported from `biopipelines/entities.py`. Importantly, for models having entities such as PDB paths, proteins sequences or ligand smiles as parameters, we always pass an entity object rather than a string to ensure representation coherence across the repository.

| Entity | Purpose |
|--------|---------|
| `PDB` | Fetch protein structures |
| `Sequence` | Defines proteins and polynucleotides from strings |
| `Ligand` | Fetch small molecules |
| `CompoundLibrary` | Create compound collections |
| `Table` | Load existing CSV files |

**PDB** - Fetches from local folders or RCSB with priority: `local_folder` → `<biopipelines>/PDBs/` → RCSB download. 
It also generates protein sequences for each of the proteins. 
If an RCSB code is provided, ligands will also be downloaded and will be available with their smiles/ccd.

```python
# Simple fetch
protein = PDB("4ufc")

# Multiple with custom IDs
proteins = PDB(["4ufc", "1aki"], ids=["POI1", "POI2"])

# From folder
proteins = PDB("/path/to/structures") # if the folder contains cif files, add format="cif"
```

**Sequence** - Creates sequences with auto-detection (protein/DNA/RNA):

```python
# Single sequence
seq = Sequence("MKTVRQERLKSIVRILERSKEPVSGAQ", ids="my_protein")

# Multiple sequences
seqs = Sequence(["MKTVRQ...", "AETGFT..."], ids=["p1", "p2"])

# Multiple DNA sequences from a file
seqs = Sequence("/path/to/sequences.csv", type="dna") # must have columns id, sequence
```

**Ligand** - Fetches from RCSB (CCD codes) or PubChem (names, CID, CAS):

```python
# RCSB by CCD code
atp = Ligand("ATP")

# PubChem by name
aspirin = Ligand("aspirin", codes="ASP")

# Direct SMILES
ethanol = Ligand(smiles="CCO", ids="ethanol", codes="ETH")
```

**CompoundLibrary** - Creates compound collections:

```python
# Simple dictionary (no expansion)
library = CompoundLibrary({
    "aspirin": "CC(=O)OC1=CC=CC=C1C(=O)O",
    "caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
})

# With expansion using <key> placeholders
library = CompoundLibrary(
    library={
        "scaffold": "<aryl><amide>",
        "aryl": ["C1(=CC(F)=CC=C1)", "C1(=CC(O)=CC=C1)"],
        "amide": ["C(=O)N","C(=O)NC","C(=O)NCC(F)(F)F"]
    },
    primary_key="scaffold"
)
```

**Table** - Loads existing CSV files:

```python
# Load a CSV file (columns auto-detected)
metrics = Table("/path/to/metrics.csv")

# Access via tables.data (default name)
Panda(
    tables=[metrics.tables.data],
    operations=[Panda.sort("score", ascending=False)]
)

# Use column reference for per-structure data
ProteinMPNN(
    structures=proteins,
    redesigned=(metrics.tables.data, "designed_positions")
)

# Custom table name
previous = Table("/path/to/results.csv", name="previous_run")
# Access via: previous.tables.previous_run
```

### Combinatorics: Bundle and Each

Control how multiple inputs combine in tools like Boltz2:

| Mode | Behavior | Example |
|------|----------|---------|
| `Each` (default) | Cartesian product | 2 proteins × 3 ligands = 6 predictions |
| `Bundle` | Group as one entity | 2 proteins bundled + 3 ligands = 3 predictions |

```python
from biopipelines.combinatorics import Bundle, Each

# Default: Each protein with each ligand (6 predictions)
boltz = Boltz2(
    proteins=Each(protein_a, protein_b),
    ligands=ligand_library  # 3 ligands
)

# Bundle ligands: Each protein with all ligands together (2 predictions)
boltz = Boltz2(
    proteins=Each(protein_a, protein_b),
    ligands=Bundle(ligand_library)
)

# Nested: Each ligand bundled with a cofactor (3 predictions per protein)
# Affinity calculated for library ligand (first in bundle)
boltz = Boltz2(
    proteins=protein_a,
    ligands=Bundle(Each(ligand_library), cofactor)
)
```


### Table Column References

Reference columns from upstream tables using tuple syntax:

```python
# RFdiffusion outputs a table with 'designed' column
rfd = RFdiffusion(contigs="50-100", num_designs=5)

# Pass column reference to downstream tool
lmpnn = LigandMPNN(
    structures=rfd,
    ligand="LIG",
    redesigned=rfd.tables.structures.designed  # Tuple: (TableInfo, "designed")
)
```

Hint: if you don't remember the table or column name, you can look it up in the ToolReference.

At SLURM time, the column value is resolved per-structure by ID matching.

---

## Resources

Set compute resources before tools:

```python
with Pipeline("Project", "Job", "Description"):
    Resources(gpu="A100", memory="32GB", time="24:00:00")    # Specific GPU
    Resources(gpu="32GB|80GB|96GB", memory="32GB", time="24:00:00")  # Memory-based
    Resources(memory="128GB", time="24:00:00", cpus=32)      # CPU-only
```

GPU options: `"T4"`, `"L4"`, `"V100"`, `"A100"`, `"H100"`, `"H200"`, `"24GB"`, `"32GB"`, `"80GB"`, `"96GB"`, `"any"`, `"high-memory"`

**Batch dependencies**: Multiple `Resources()` calls create sequential batches:

```python
with Pipeline("Project", "Job", "Description"):
    Resources(gpu="V100", time="4:00:00")    # Batch 1
    tool1 = RFdiffusion(...)

    Resources(time="2:00:00")                # Batch 2 (waits for Batch 1)
    tool2 = Panda(...)
```

---

## On-the-fly Execution

For interactive prototyping in Jupyter notebooks or Google Colab, on-the-fly mode is **enabled automatically** — no extra arguments needed. Each tool's bash script is executed immediately when the tool is added, so you see results step by step. The pipeline stays active across notebook cells, so each cell can add new tools:

```python
# Cell 1: Create pipeline and run first tools
Pipeline("Examples", "interactive_test",
         description="Quick test run")

lysozyme = PDB("168L")
rfd = RFdiffusion(pdb=lysozyme,
                  contigs='50-70/A81-140',
                  num_designs=3)
# rfd has already finished running at this point
```

```python
# Cell 2: Continue adding tools to the same pipeline
pmpnn = ProteinMPNN(structures=rfd, num_sequences=2)
# pmpnn has already finished running at this point
```

You can also force on-the-fly mode explicitly with `on_the_fly=True` (e.g., when running locally with plain Python outside a notebook).

Key differences from normal mode:
- The `with` statement is **optional** — the pipeline context stays active across cells
- `Resources()` is **optional** and ignored for execution purposes
- Tools run sequentially as they are added — each tool finishes before the next one starts
- stdout/stderr is streamed in real-time (visible in notebooks and terminals)
- SLURM submission is skipped
- Output is written to `./BioPipelines/` (current directory) instead of shared storage
- The completion check mechanism is preserved, so re-running a notebook skips already-completed steps
- Re-running the cell that creates the `Pipeline()` starts a new pipeline

---

## Job Submission

**Submit to SLURM**:
```bash
biopipelines-submit /path/to/pipeline.py
biopipelines-submit /path/to/pipeline.ipynb   # directly from a notebook
```

**Run directly**:
```bash
biopipelines-run /path/to/pipeline.py
biopipelines-run /path/to/pipeline.ipynb       # directly from a notebook
```

Both `.py` scripts and `.ipynb` notebooks are supported. When a notebook is provided, code cells are automatically extracted (skipping shell commands and IPython magics) and executed as a script.

These commands work from any directory as long as biopipelines environment is activated. Alternatively, you can run the scripts directly from the biopipelines root:
```bash
cd biopipelines
./submit /path/to/pipeline.py
./run /path/to/pipeline.py
```

**Resubmit** existing job:

```bash
./resubmit /path/to/job/RunTime/slurm.sh
```

**External dependencies** - Wait for other SLURM jobs:

```python
with Pipeline("Project", "Job", "Description"):
    Dependencies("12345678")  # Wait for job ID. Also accepts lists
    Resources(gpu="V100", time="4:00:00")
    ...
```

---

## Data Management with Panda

Panda provides pandas-style table transformations:

```python
from biopipelines.panda import Panda

# Filter rows
filtered = Panda(
    tables=boltz.tables.confidence,
    operations=[Panda.filter("pLDDT > 80")]
)

# Sort and take top N
best = Panda(
    tables=boltz.tables.confidence,
    operations=[
        Panda.sort("confidence_score", ascending=False),
        Panda.head(5)
    ]
)

# Merge tables horizontally
merged = Panda(
    tables=[apo.tables.affinity, holo.tables.affinity],
    operations=[
        Panda.merge(on="id", prefixes=["apo_", "holo_"]),
        Panda.calculate({"delta": "holo_affinity - apo_affinity"})
    ]
)

# Concatenate tables vertically
combined = Panda(
    tables=[cycle0.tables.results, cycle1.tables.results],
    operations=[Panda.concat(fill="", add_source=True)]
)

# Pool mode: copy structures matching filtered IDs
filtered_with_files = Panda(
    tables=boltz.tables.confidence,
    operations=[Panda.filter("pLDDT > 80")],
    pool=boltz  # Copy matching structures
)

# Rename output IDs after sorting
ranked = Panda(
    tables=boltz.tables.confidence,
    operations=[Panda.sort("score", ascending=False)],
    rename="best",  # Output: best_1, best_2, ...
    pool=boltz
)
```

Available operations: `filter`, `sort`, `head`, `tail`, `sample`, `rank`, `drop_duplicates`, `merge`, `concat`, `calculate`, `groupby`, `select_columns`, `drop_columns`, `rename`, `fillna`, `pivot`, `melt`, `average_by_source`

---

## Filesystem Structure

```
/shares/<user>/BioPipelines/<project>/<job>_<NNN>/
├── RunTime/                    # Execution scripts
│   ├── pipeline.sh
│   ├── 001_<tool>.sh
│   └── 002_<tool>.sh
├── Logs/                       # Execution logs
│   ├── 001_<tool>.log
│   └── 002_<tool>.log
├── ToolOutputs/                # Tool output predictions (JSON)
│   ├── 001_<tool>.json
│   └── 002_<tool>.json
├── 001_<Tool>/                 # Tool outputs
│   ├── structures/
│   └── results.csv
└── 002_<Tool>/
    └── ...
```

Configure paths and environments in `config.yaml` at repository root.

---

## Troubleshooting

**Path errors**: Run from BioPipelines root directory.

**Environment issues**: Check conda environments in `config.yaml`.

**Resource limits**: Adjust GPU/memory in `Resources()`.

**Missing files**: Check `Logs/<NNN>_<tool>.log`.

**Local output**: Write results to the current directory instead of the config-defined path:

```python
with Pipeline("Test", "Debug", "Testing", local_output=True):
    ...
```

Note: `local_output` defaults to `True` automatically when `on_the_fly` is enabled (i.e., in Jupyter notebooks).

**Load previous outputs**:

```python
from biopipelines.load import LoadOutput, LoadOutputs

# Single output
prev = LoadOutput("/path/to/ToolOutputs/001_Boltz2.json")

# Multiple outputs by tool name
all_boltz = LoadOutputs("/path/to/job/", tool="Boltz2")
```
