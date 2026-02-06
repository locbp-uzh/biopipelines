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
- [Job Submission](#job-submission)
- [Data Management with Panda](#data-management-with-panda)
- [Filesystem Structure](#filesystem-structure)
- [Troubleshooting](#troubleshooting)

---

## What is BioPipelines?

BioPipelines is a Python framework that generates bash scripts for bioinformatics workflows. It does not execute computations directly - instead, it predicts the filesystem structure and creates scripts that will be executed on SLURM clusters.

---

## Installation
Login to your cluster via terminal or the website ([S3IT Apps](https://apps.s3it.uzh.ch) > Clusters > Shell access). In your home directory, clone the biopipelines repository:
```bash
git clone https://gitlab.uzh.ch/locbp/public/biopipelines
```

---

## Quick Start

```python
from PipelineScripts.pipeline import *
from PipelineScripts.entities import *
from PipelineScripts.boltz2 import Boltz2

with Pipeline("MyProject", "JobName", "Description"):
    Resources(gpu="A100", time="4:00:00", memory="16GB")

    # Fetch a protein structure
    protein = PDB("4ufc", ids="LYZ")

    # Predict structure with a ligand
    prediction = Boltz2(
        proteins=protein,
        ligands=Ligand("ATP")
    )

    print(prediction)
```

---

## Core Concepts

### Pipeline and SLURM Runtime

BioPipelines operates in two phases:

| Phase | What Happens | Executer | Where |
|-------|--------------|-------|-------|
| **Pipeline time** | Generation of bash scripts, prediction of output paths and files | Python | Cluster |
| **SLURM time** | Bash scripts execute, files are created | Slurm | Cluster |

Tools predict their outputs before execution. These predictions enable chaining:

```python
rfd = RFdiffusion(contigs="50-100", num_designs=5)
# rfd contains predicted paths (files don't exist yet)
mpnn = ProteinMPNN(structures=rfd, num_sequences=2)
# mpnn uses rfd's predicted outputs
```

### Entity Types

Basic input types can be imported conveniently from `PipelineScripts/entities.py`:

| Entity | Purpose |
|--------|---------|
| `PDB` | Fetch protein structures |
| `Sequence` | Defines proteins and polynucleotides from strings |
| `Ligand` | Fetch small molecules |
| `CompoundLibrary` | Create compound collections |
| `Table` | Load existing CSV files |

**PDB** - Fetches from local folders or RCSB with priority: `local_folder` → `<biopipelines>/PDBs/` → RCSB download. 
It also generate protein sequences for each of the proteins. 
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
seqs = Sequence("/path/to/sequences.csv", type="dna")
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

### DataStream vs Tables

Tools output two types of data containers:

**DataStream** - Unified container supporting ID tracking, and association of IDs to files (e.g. PDBs) or values (e.g. protein sequence):

```python
# Access structures from a tool via the streams container
for struct_id, pdb_path in boltz.streams.structures:
    print(f"{struct_id}: {pdb_path}")

# Count items
print(f"Generated {len(boltz.streams.structures)} structures")
print(f"Expected ids: {boltz.streams.structures.ids}")
```

DataStream types accessed via `tool.streams.<type>`:
- `streams.structures` - PDB/CIF files
- `streams.sequences` - Table with id, sequence columns
- `streams.compounds` - CSV with SMILES column, or SDF/PDB files
- `streams.msas` - A3M or CSV files

**Tables (TableInfo)** - Rich metadata about CSV files. They do not track IDs.

```python
# Access table path
path = tool.tables.confidence  # Returns path string

# Access table metadata
info = tool.tables.info("confidence")
print(info.path)        # /path/to/confidence.csv
print(info.columns)     # ["id", ""pTM", "complex_plddt", ...]
print(info.description) # "Confidence scores"
print(info.count)       # Number of rows
```

### Combinatorics: Bundle and Each

Control how multiple inputs combine in tools like Boltz2:

| Mode | Behavior | Example |
|------|----------|---------|
| `Each` (default) | Cartesian product | 2 proteins × 3 ligands = 6 predictions |
| `Bundle` | Group as one entity | 2 proteins bundled + 3 ligands = 3 predictions |

```python
from PipelineScripts.combinatorics import Bundle, Each

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

## Job Submission

```bash
cd biopipelines
./submit /path/to/pipeline.py
```

For memory errors during submission (This usually happens during first submission as the environment biopipelines is created):

```bash
srun --mem=8G --time=1:00:00 --pty bash
./submit /path/to/pipeline.py
```

Resubmit existing job:

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
from PipelineScripts.panda import Panda

# Filter rows
filtered = Panda(
    table=boltz.tables.confidence,
    operations=[Panda.filter("pLDDT > 80")]
)

# Sort and take top N
best = Panda(
    table=boltz.tables.confidence,
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
    table=boltz.tables.confidence,
    operations=[Panda.filter("pLDDT > 80")],
    pool=boltz  # Copy matching structures
)

# Rename output IDs after sorting
ranked = Panda(
    table=boltz.tables.confidence,
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

**Debug mode**: Test locally without SLURM:

```python
with Pipeline("Test", "Debug", "Testing", debug=True):
    ...
```

**Load previous outputs**:

```python
from PipelineScripts.load import LoadOutput, LoadOutputs

# Single output
prev = LoadOutput("/path/to/ToolOutputs/001_Boltz2.json")

# Multiple outputs by tool name
all_boltz = LoadOutputs("/path/to/job/", tool="Boltz2")
```
