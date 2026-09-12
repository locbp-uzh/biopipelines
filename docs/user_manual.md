# BioPipelines User Manual

## Index

- [What is BioPipelines?](#what-is-biopipelines)
- [Ways to run BioPipelines](#ways-to-run-biopipelines)
- [Usage with an AI coding assistant](#usage-with-an-ai-coding-assistant)
- [Installation](#installation)
- [Google Colab](#google-colab)
- [CSCS Alps / Daint](#installation-cscs-alps--daint)
- [Local](#installation-local--linux--macos--windows)
- [Quick Start](#quick-start)
- [Core Concepts](#core-concepts)
  - [Configuration and Execution Time](#configuration-and-execution-time)
  - [Entity Types](#entity-types)
  - [DataStream vs Tables](#datastream-vs-tables)
  - [Combinatorics: Bundle and Each](#combinatorics-bundle-and-each)
  - [Table Column References](#table-column-references)
- [Resources](#resources)
- [Naming Runs with Suffix](#naming-runs-with-suffix)
- [Grouping Outputs with Folder](#grouping-outputs-with-folder)
- [On-the-fly Execution](#on-the-fly-execution)
- [Job Submission](#job-submission)
- [Saving Without Submitting, and Background Services](#saving-without-submitting-and-background-services)
- [Data Management with Panda](#data-management-with-panda)
- [When IDs Disappear: the `missing` Table](#when-ids-disappear-the-missing-table)
- [Reading Tables Back in Python](#reading-tables-back-in-python)
- [Filesystem Structure](#filesystem-structure)
- [Troubleshooting](#troubleshooting)

---

## What is BioPipelines?

BioPipelines is a Python framework for writing and executing protein engineering workflows that run with the same syntax on HPC clusters (SLURM, LSF, PBS), Jupyter, and Google Colab.

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

## Ways to run BioPipelines

There are three ways to use BioPipelines, in increasing order of hands-on effort:

1. **With an AI coding assistant (recommended for non-programmers).** You clone the repo and on your computer and run an AI coding assistant (Claude Code, Codex, …) inside it. The assistant reads the framework's prompts under `llm/`, interviews you about your biological problem, and writes and runs the pipeline for you. You never write Python yourself — you describe the protocol you want in plain language. This is the easiest entry point and the one most users should start with. See [Usage with an AI coding assistant](#usage-with-an-ai-coding-assistant).
2. **On an HPC cluster (SLURM, LSF, or PBS/Torque).** You write pipeline scripts yourself and submit them with `biopipelines-submit`. Best for large production runs on institutional compute. See [Installation (Slurm HPC)](#installation-slurm-hpc).
3. **On Google Colab.** You write pipeline cells in a notebook and run them inline with a Colab GPU, no SLURM needed. Best for quick interactive prototyping. See [Installation (Google Colab)](#installation-google-colab).
4. **Locally (Linux / macOS / Windows).** You run pipelines inline on your own machine if they only comprise lightweight base-environment tools. On Windows, use WSL. See [Installation (Local)](#installation-local--linux--macos--windows).

These are not mutually exclusive: a common pattern is to let an AI assistant author a pipeline (mode 1) and then submit the resulting script to a cluster (mode 2) or run it on Colab (mode 3).

---

## Usage with an AI coding assistant

If you are not comfortable writing Python, the simplest way to use BioPipelines is to let an AI coding assistant drive it for you. You describe the computational protocol you want in plain language; the assistant translates that into a BioPipelines pipeline, runs it (on a cluster or Colab), inspects the outputs, and iterates.

The repository ships with prompt files under `llm/` that configure the assistant with the framework's contract — typed streams, the `Pipeline`/`Resources` API, the available tools — so it grounds its suggestions in what BioPipelines actually provides instead of guessing.

### 1. Clone the repository

```bash
git clone https://github.com/locbp-uzh/biopipelines
cd biopipelines
```

### 2. Open an AI coding assistant inside the repo

Start your assistant from the repository root so it can read the `llm/` prompts and the rest of the codebase. For example, with [Claude Code](https://claude.com/claude-code):

```bash
claude
```

or with [Codex](https://openai.com/codex/):

```bash
codex
```

Any repository-aware assistant works (Claude Code, Codex, Cursor, Copilot Chat, …). The only requirement is that it can read files in the working directory.

### 3. Point the assistant at the right prompt and state your goal

The `llm/` folder holds two session prompts, depending on what you want to do, backed by three per-backend references (`cluster.md`, `colab.md`, `daint.md`):

- **`llm/pipelines.md`** — to *use* the framework: design and run a pipeline for a specific biological problem. This is what most users want.
- **`llm/development.md`** — to *change* the framework itself: add a tool wrapper, fix a bug, refactor internals.

Open your session with a message like:

> Read and follow `llm/pipelines.md`. I want to redesign the N-terminal domain of lysozyme with RFdiffusion, then inverse-fold and validate with AlphaFold.

The assistant reads the prompt, loads the framework's documentation, interviews you about any open choices (which tools, how many designs, which execution mode), and then writes and runs the pipeline.

### 4. One-time setup for where the pipeline will run

The assistant authors the pipeline locally, but the heavy compute runs elsewhere — on a cluster or on Colab. Set that up once:

- **Cluster:** follow `llm/cluster.md` to add an ssh alias and the `llm/log.sh` command logger, then copy `llm/resources.md.template` to `llm/resources.md` and let the assistant fill in your cluster's partitions, GPU types, and walltime policy. This gives the assistant honest defaults instead of generic guesses.
- **Colab:** no setup needed here — the assistant hands you notebook cells to run, and you paste back the outputs.

After that, the assistant can push the pipeline, submit it, tail the logs, and report back the results — all from within the same session.

---

## Installation (Slurm HPC)

This section walks through a SLURM cluster as the worked example. LSF and PBS/Torque clusters install identically — only the `machine.scheduler` block differs (set `name` to `lsf` or `pbs`, or let `bp-config auto` detect it). See [Other batch schedulers (LSF, PBS/Torque)](#other-batch-schedulers-lsf-pbstorque) under Job Submission.

### 1. Setting up BioPipelines

```bash
# Clone the repository
git clone https://github.com/locbp-uzh/biopipelines
cd biopipelines

# Create the biopipelines conda environment (use mamba, conda, or micromamba)
mamba env create -f environments/biopipelines.yaml
mamba activate biopipelines

# Install the package (editable mode, so updates via git pull take effect immediately)
pip install -e .

# Optional: register a Jupyter kernel if you plan to work with Jupyter
pip install ipython
ipython kernel install --user --name biopipelines
```

### 2. Configuring your machine

Each machine you run on (cluster, laptop, colab) is described by a `config.<variant>.yaml` at the repo root. The `bp-config` CLI manages those files.

**Your edits go into an overlay, not the committed file.** `config.<variant>.yaml` ships committed defaults and stays pristine. Anything you change with `bp-config` (auto / edit / set) is written to a gitignored `.config.<variant>.yaml` — the *overlay* — alongside it. At load time the overlay is deep-merged on top of the base, so the active config is `base ⊕ your-edits`. Two consequences: a `git pull` that adds or changes repo defaults never clobbers your local settings, and your overlay carries only the keys you actually changed (the rest tracks the repo). The overlay is disposable — delete `.config.<variant>.yaml` to reset to the committed defaults. `bp-config show` prints the merged result.

```bash
# Probe the host (env manager, scheduler, modules, container runtime, username, git email) and write the discovered values under `machine:` into the overlay of the variant of your choice. Pop a picker, choose `cluster`, hit Enter.
bp-config auto
```

Once `machine.username` matches your Unix user, that variant becomes the auto-detected default — so subsequent `bp-config` commands hit `cluster.yaml` without needing `--variant cluster`.

```bash
# Open the active config in an interactive TUI to fill in the rest (folders + per-tool conda environments)
bp-config edit
```

The editor shows the merged config (base + overlay), so you navigate the whole tree; on save it writes only your changes to the overlay. In the editor:

- **Folders highlighted in red** still hold placeholder values — point
  them at real paths on the cluster (e.g. shared scratch, output dirs,
  per-tool repo locations).
- **Environments** — for any tool you've already installed under a
  custom env name, set `environments.<Tool>` to that name so
  BioPipelines activates the right env at execution time. Tools you
  haven't installed yet can stay at the default; `Tool.install()` will
  create them later.
- **`folders.infrastructure.scripts`** — the folder the `Scripting` tool
  searches for a bare script filename (default `<biopipelines>/my_scripts`).
  Keep your `Scripting` scripts there and call `Scripting("my_step.py", …)`
  without a path.

Other useful `bp-config` subcommands:

```bash
bp-config path             # path of the overlay (created if absent); --base for the committed file
bp-config list             # list every config.<variant>.yaml in the repo
bp-config show             # print the resolved config (base merged with overlay)
bp-config set <key> <val>  # set one dotted key non-interactively (writes the overlay)
bp-config folder <key>     # resolve one folder path
bp-config env <Tool>       # the conda env configured for a tool
```

### 3. Installing tools

Each external tool (RFdiffusion, ProteinMPNN, AlphaFold, …) lives in its own conda env. Install them on demand:

```python
from biopipelines.pipeline import *
from biopipelines import RFdiffusion, ProteinMPNN, AlphaFold

with Pipeline("Setup", "install", description="Install tools"):
    RFdiffusion.install()
    ProteinMPNN.install()
    AlphaFold.install()
```

`.install()` clones the upstream repo (where applicable) and creates the per-tool conda env. Re-running is a no-op once the env exists; pass `force_reinstall=True` to rebuild.

*Note: installation was optimized on the S3IT UZH HPC (Ubuntu 24.04, SLURM 25.05) and might require adjusting e.g. due to CUDA version mismatch.*

### 4. Submitting a pipeline

After activating the biopipelines environment:

```bash
# Cluster: generate scripts and submit them to the scheduler (SLURM/LSF/PBS)
bp-submit example_pipelines/<pipeline>.py

# Laptop / interactive: run inline (no scheduler); each tool's bash script
# executes immediately as the tool is added (need resources on node)
bp-run example_pipelines/<pipeline>.py
```

Note: you must have installed the relevant tool environments
(`Tool.install()`) and configured them in your `config.<variant>.yaml`.

---

---

## Installation (Google Colab)

BioPipelines runs on Google Colab with GPU support out of the box. No SLURM needed — tools are installed via `micromamba` into isolated environments, matching the cluster behavior.

### 1. Setting up BioPipelines

Run this cell at the top of your Colab notebook:

```python
# Cell 1: Install BioPipelines and micromamba
!git clone https://github.com/locbp-uzh/biopipelines
%cd biopipelines
!pip install -e ".[colab]"
!wget -q https://github.com/mamba-org/micromamba-releases/releases/latest/download/micromamba-linux-64 -O /usr/local/bin/micromamba && chmod +x /usr/local/bin/micromamba
```

On Colab the `biopipelines` env is **not** created: tools mapped to it run in base Python, where `pip install -e ".[colab]"` installed the deps (the `colab` extra adds conda-only packages like `openbabel-wheel`). That env is never activated on Colab, so creating it just wastes setup time. The `micromamba` binary is still needed — each per-tool `.install()` creates its own env with it.

BioPipelines automatically detects the Colab environment and loads `config.colab.yaml` instead of the site variant it would otherwise pick. This sets `env_manager: "micromamba"`, which means:

- Each tool gets its own isolated conda environment (same as on the cluster)
- Tool installation uses `micromamba` to create environments from YAML specs
- Scheduler-related settings are disabled

### 2. Installing tools

Install the tools you need using `.install()`. This only needs to run once per Colab session (completed steps are skipped on re-run):

```python
from biopipelines.pipeline import *
from biopipelines import RFdiffusion, ProteinMPNN, AlphaFold

with Pipeline("Setup", "install", description="Install tools"):
    RFdiffusion.install()
    ProteinMPNN.install()
    AlphaFold.install()
```

### 3. Running Pipelines

After installation, pipelines work exactly as described in the rest of this manual. On-the-fly execution is enabled automatically in notebooks:

```python
Pipeline("Examples", "RFD-ProteinMPNN-AF2",
         description="Redesign of N terminus domain of lysozyme")

lysozyme = PDB("168L")
rfd = RFdiffusion(pdb=lysozyme,
                  contigs='50-70/A81-140',
                  num_designs=3)
pmpnn = ProteinMPNN(structures=rfd, num_sequences=2)
af = AlphaFold(proteins=pmpnn)
```

### Key Differences from Cluster

| | Cluster | Google Colab |
|---|---|---|
| **Environment manager** | mamba/conda/micromamba | micromamba (hardcoded) |
| **Execution** | Batch scheduler (SLURM/LSF/PBS) or on-the-fly | On-the-fly only |
| **GPU** | Configured via `Resources()` | Colab's assigned GPU |

Colab sessions are ephemeral. Installed tools and generated outputs are lost when the runtime disconnects. Mount Google Drive or download results before the session ends.

---

## Installation (CSCS Alps / Daint)

Daint needs its own variant because three of its properties break the assumptions the SLURM install above makes:

- **aarch64 (GH200)** — x86-64 conda builds and `.sif` images do not apply.
- **No lmod modules, no apptainer** — GPU tools run through the CSCS Container Engine.
- **Per-project billing** — a job without an account is rejected.

### 1. Setting up BioPipelines

Use a venv, not conda, and put it on `$SCRATCH`. `$STORE` is backed up and shared with your project, but its quota allows only 150,000 files, and Python environments exhaust that long before the 1 TB of space — 13 tools' venvs come to ~101,000 files for 15 GB. Scratch has 1M inodes; the cost is its 30-day access-time purge, so a tool left unused for a month needs reinstalling.

```bash
git clone https://github.com/locbp-uzh/biopipelines
cd biopipelines-locbp

# There is no `python` on the login nodes — only python3 (3.6) and python3.11.
/usr/bin/python3.11 -m venv $SCRATCH/venvs/biopipelines
source $SCRATCH/venvs/biopipelines/bin/activate

# aarch64 wheels; the conda yaml's rdkit/openbabel/py3dmol are pip packages here.
pip install -r environments/biopipelines.pip.daint.txt
pip install -e .
```

Activating the venv is what supplies a `python` on PATH, which the generated tool scripts call.

### 2. Configuring your machine

```bash
bp-config set machine.billing_account <your-project> --variant daint
```

This becomes `#SBATCH --account=`, and CSCS rejects jobs without it. It is the *project* charged for compute, not your username. As with any variant, your edit lands in the gitignored overlay `.config.daint.yaml`, not the committed file.

### 3. Running

```bash
BIOPIPELINES_CONFIG_VARIANT=daint biopipelines-submit my_pipeline.py
```

### GPU tools

GPU tools do not use a conda env. Each is mapped in `config.daint.yaml`'s `edf:` block to an Environment Definition File, and its whole script runs via `srun --environment=<edf>`. The image comes from NVIDIA NGC, which publishes GH200-native PyTorch, so torch and CUDA arrive prebuilt rather than being ported. A tool's venv is layered on the image with `--system-site-packages` and must be built inside the same container.

Not every tool is available. Tools needing an x86-64-only container image do not run on Daint — RFdiffusion is the clearest case, since the RosettaCommons images have no arm64 build. PyMOL, mkdssp and ProteinMPNN *do* run there, supplied by their EDF image. `config.daint.yaml`'s `environments:` block is a routing table, not a record of what has been verified: it names an env for every tool that has one, whether or not the tool has been run. `llm/daint.md` carries the list actually verified by running, and the blocked list with a reason for each.

### Nodes are billed whole

Every GPU partition is `OverSubscribe=EXCLUSIVE`, and every node is 4× GH200 / 288 CPUs / 870 GB. A job requesting one GPU is allocated and billed for the entire node, so splitting a sweep into many single-GPU jobs wastes 4× on GPUs and 288× on CPUs. There are no A100s — `Resources(gpu="A100")` has no meaning here.

---

## Installation (Local — Linux / macOS / Windows)

You can run pipelines with only lightweight tools (input preparation, table transforms, and small analyses) locally. On Windows, do this inside [WSL](https://learn.microsoft.com/windows/wsl/) (`wsl --install -d Ubuntu`).

Install into a virtualenv with Python `>=3.10,<3.13` (only the base deps — pandas, numpy, biopython, rdkit; no conda):

```bash
cd /path/to/biopipelines
python3 -m venv .venv && source .venv/bin/activate
pip install -e .
```

The `local` config variant (`env_manager: pip`, `scheduler: none`) is the default when no batch scheduler is present; check with `bp-config machine`. Then run a script with `bp-run`:

```bash
bp-run my_pipeline.py # or: BIOPIPELINES_OTF=1 python my_pipeline.py
```

---

## Core Concepts

### Configuration and Execution Time

BioPipelines operates alternatively in one or two phases:

One phase (Jupyter/Colab notebooks):
| Phase | What Happens | Executer |
|-------|--------------|-------|
| **Notebook cell** | Generation and running of bash scripts, prediction of output paths and files | Python + bash |

Two phases (biopipelines-submit):
| Phase | What Happens | Executer |
|-------|--------------|-------|
| **Configuration** | Generation of bash scripts, prediction of output paths and files | Python |
| **Execution** | Bash scripts execute, files are created | Bash + Python |

### DataStream vs Tables

Tools output two types of data containers:

**DataStream** - Unified container supporting ID tracking, and association of IDs to files (e.g. .pdb, .cif, .sdf) or values (e.g. protein/dna sequences):

DataStream types accessed via `tool.streams.<name>`:
- `streams.structures` - PDB/CIF/SDF coordinate files (a ligand's 3-D coordinates live here)
- `streams.sequences` - ID-tracked table with id, sequence columns
- `streams.compounds` - value-based CSV (the ligand's chemistry/identity: `smiles`, `code`)
- `streams.msas` - A3M or CSV files
- `streams.images` - PNG files

In notebooks and on-the-fly runs, a tool has finished by the time the next
Python line runs. Use `records()` to inspect the materialized items in a stream:

```python
for pdb in af.streams.structures.records():
    print(pdb.id, pdb.file)

for compound in lig.streams.compounds.records(columns=["smiles", "code"]):
    print(compound.id, compound.smiles, compound.code)
```

`records()` is for inspecting results after they exist. To pass stream items to
another BioPipelines tool during pipeline construction, pass the stream itself
or iterate the stream as single-item DataStreams.

**Tables (TableInfo)** - Rich metadata about CSV files. They do not track IDs.

```python
# Access table metadata via .info
info = tool.tables.confidence.info
print(info.path)        # /path/to/confidence.csv
print(info.columns)     # ["id", "pTM", "complex_plddt", ...]
print(info.description) # "Confidence scores"

# Access column references for downstream tools
tool.tables.confidence.plddt  # Returns (TableInfo, "plddt") tuple
```

In the ToolReference, one can find for each tool what is the expected output in terms of streams and tables, and use this information to write pipelines.

### Common inputs

Basic input types can be imported from `biopipelines/entities.py`. Importantly, for models having entities such as PDB paths, proteins sequences or ligand smiles or codes as parameters, we always pass an entity object rather than a string to ensure representation coherence across the repository. The same is true for ligand codes (some tools alter them in output structures).

| Entity | Purpose |
|--------|---------|
| `PDB` | Fetch protein structures |
| `Sequence` | Defines proteins and polynucleotides from strings |
| `Ligand` | Fetch small molecules |
| `CompoundLibrary` | Create compound collections |
| `Table` | Load existing CSV files |

#### Structure and Compound

`Structure` is another name for `PDB` and `Compound` for `Ligand` — the same tools, named for the streams they emit (`structures`, `compounds`), and `PDB` accepts mmCIF as readily as PDB. Both spellings are first-class and neither is deprecated; `Structure is PDB` and `Compound is Ligand` are literally true, because each is a bare alias rather than a subclass.

Outputs keep the original name either way. `TOOL_NAME` stays `"PDB"` and `"Ligand"`, and that single value drives the step's output folder name, the config's `environments:` and `folders:` keys, the log lines, and any `Load()` pointing at an existing folder. So a `Structure(...)` step still writes `001_PDB/` and still logs `PDB`, and a `Compound(...)` step still writes `001_Ligand/`. That is the cost of keeping one registry entry per tool: if you rename a step in your script, the paths on disk do not follow.

**PDB** - Fetches from local folders or RCSB with priority: `local_folder` → `<biopipelines>/pdbs/` → RCSB download.
It also generates protein sequences for each of the proteins.
If an RCSB code is provided, ligands will also be downloaded and will be available with their smiles/ccd.

```python
# Simple fetch
protein = PDB("4ufc")

# Multiple with custom IDs
proteins = PDB(["4ufc", "1aki"], ids=["POI1", "POI2"])

# From folder
proteins = PDB("/path/to/structures")  # convert defaults to None (pdb|cif, no conversion); pass convert="pdb" to convert all to PDB
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

**Ligand** - Fetches from RCSB (CCD codes), PubChem (names, CID, CAS), SMILES, or CDXML:

```python
# RCSB by CCD code
atp = Ligand("ATP")

# PubChem by name
aspirin = Ligand("aspirin", codes="ASP")

# Direct SMILES
ethanol = Ligand(smiles="CCO", ids="ethanol", codes="ETH")

# From CDXML file: each molecule is a separate ligand; ChemDraw names used as IDs.
# There is no cdxml= parameter — pass the path as the lookup value and the
# .cdxml suffix is recognized automatically (a .txt of SMILES works the same way).
ligands = Ligand("my_ligands.cdxml")

# Code-only: name an existing HETATM residue code (no download, no SMILES).
# Produces a compounds stream you hand to tools that read a ligand's code
# (LigandMPNN, PoseBusters, PLIP, RFdiffusionAllAtom, RFdiffusion3, …).
lig = Ligand(codes="ZIT")

# Carve a bound ligand out of complexes, keeping the crystal coordinates.
posed = Ligand(codes="STI", structures=complexes)
```

**One `codes` parameter, and the rest of the call selects the mode.** `codes` is the residue code, and what you pass *alongside* it decides what the `Ligand` is:

| Call | Mode |
|---|---|
| `Ligand(codes="UNL")` | code-only — names an existing HETATM residue: no chemistry, no `structures` stream |
| `Ligand(smiles=…, ids=…, codes="BGD")` | chemistry you supplied, plus the residue label it carries downstream (the dominant real usage) |
| `Ligand(codes="STI", structures=complexes)` | carves that residue out of the given structures, keeping their coordinates |
| `Ligand("aspirin")` | `lookup` is the first positional, so this fetches chemistry — **not** code-only |

This used to be two parameters: `codes=` for the label on a ligand with chemistry, `code=` for code-only construction. One letter apart, mutually exclusive, and confusing them silently produced a chemistry-free stub. `code=` still binds `codes=` so old scripts keep running, but it will soon be deprecated and each use prints one line:

```
[contract:deprecated_alias] Ligand: code= is a synonym for codes= and will soon be deprecated.
```

Because a forgotten `smiles=` now yields a stub instead of an error, the code-only path announces itself too:

```
[contract:code_only_ligand] Ligand(codes='UNL') carries no chemistry, only a residue code
    Pass lookup or smiles to retrieve chemical information.
```

Both lines are severity-gated like every other contract check — `BIOPIPELINES_ENFORCE_CODE_ONLY_LIGAND=off` or `BIOPIPELINES_ENFORCE_DEPRECATED_ALIAS=off` silences one for a run.

**Two calls that used to raise and now work.** `Ligand(codes=…, structures=…)` was rejected as "not compatible" and is now the carve path; `Ligand(code=…, smiles=…)` was rejected as "mutually exclusive" and is now the labelled-chemistry case. A script written against the old API does not fail on these — it does something, so check it against the table above.

A code-only Ligand has no SMILES, so it cannot be converted to 3-D. To get a 3-D ligand (e.g. an SDF for docking-adjacent tools), start from a Ligand that carries a SMILES and run OpenBabel:

```python
aspirin = Ligand("aspirin")                       # has SMILES
sdf = OpenBabel(compounds=aspirin, convert_3d="sdf")
# sdf.streams.structures -> the SDF; sdf.streams.compounds -> chemistry passthrough
```

**Ligand string shorthand.** Tools that read a ligand by its residue code (LigandMPNN, PLIP, PoseBusters, RFdiffusionAllAtom, RFdiffusion3) accept a bare string in place of a `Ligand`: `ligand="LIG"` is shorthand for `ligand=Ligand(codes="LIG")`. It creates a code-only ligand — no chemistry, no SMILES — and is exactly equivalent to constructing the `Ligand` yourself. Filtering by residue code is the whole point of the shorthand, so the `code_only_ligand` notice is deliberately suppressed for it; write `Ligand(codes="LIG")` out by hand and you get the notice, because there the omission of `smiles=` might be an accident. For a ligand with chemistry (to fetch SMILES, or to convert to 3-D), pass an explicit `Ligand("ATP")` / `Ligand(smiles=...)` instead; the string shorthand never fetches or downloads. The auto-created ligand is registered as an internal step (see [Filesystem Structure](#filesystem-structure)).

```python
# These two are equivalent:
LigandMPNN(structures=rfd, ligand="STI")
LigandMPNN(structures=rfd, ligand=Ligand(codes="STI"))
```

**CompoundLibrary** - Creates compound collections:

```python
# Simple dictionary (no expansion): keys are compound IDs, values are SMILES
library = CompoundLibrary({
    "aspirin": "CC(=O)OC1=CC=CC=C1C(=O)O",
    "caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
})

# With expansion using <key> placeholders — primary key auto-detected by order of appearance
library = CompoundLibrary({
    "scaffold": "<aryl><amide>",
    "aryl": ["C1(=CC(F)=CC=C1)", "C1(=CC(O)=CC=C1)"],
    "amide": ["C(=O)N","C(=O)NC","C(=O)NCC(F)(F)F"]
})
# Generates 2×3=6 compounds; branching columns 'aryl' and 'amide' track substituents

# From CDXML file (ChemDraw R-group enumeration); names defined in ChemDraw are used
library = CompoundLibrary("my_library.cdxml")

# From CSV file (expansion supported if SMILES column contains <placeholders>)
library = CompoundLibrary("my_library.csv")

# With 2D molecule images (PNG per compound, uses RDKit — no extra dependencies)
library = CompoundLibrary({...}, generate_images=True)
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
    redesigned=metrics.tables.data.designed_positions
)

# Custom table name
previous = Table("/path/to/results.csv", table_name="previous_run")
# Access via: previous.tables.previous_run
```

The parameter is `table_name`, not `name`: `name` is framework-reserved as the job name, and `Table` used to capture it, so the job name of every `Table` step silently vanished. `name=` still works as a deprecated synonym — it names the handle a downstream step reaches the table by (`previous.tables.previous_run`), and silently changing that would have broken every consumer — but it prints a `[contract:deprecated_alias]` line and now also sets the job name. Pass both and `table_name` names the table while `name` is only the job name.

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

**Output ID naming**: Output IDs are always the full cartesian product of all iterated axis IDs joined with `+`. For example, 1 protein (`prot1`) × 3 ligands (`lig1`, `lig2`, `lig3`) produces IDs `prot1+lig1`, `prot1+lig2`, `prot1+lig3`. There are no shortcuts — even with a single protein, the protein ID is always included. The `+` separator is deliberately distinct from `_`, which is reserved for parent→child suffixes (`protein_1`, `protein_2`), so a multi-axis ID is never mistaken for a suffixed one.

**Provenance columns**: A tool's postprocessing step can append one `{stream_name}.id` column per input axis (e.g. `sequences.id`, `compounds.id`), tracking which input produced each output row. These are added at runtime, so they are not in the table's declared columns — and which tables get them is per-tool: Boltz2 adds them to `structures`, `confidence` and `affinity`, while ESMFold2 adds them to `structures` only. Check the tool's page. Where present they make filtering and joining direct:

```python
# Filter Boltz2 results for a specific protein
df = pd.read_csv(boltz.tables.confidence.info.path)
prot1_results = df[df['sequences.id'] == 'prot1']

# Filter for a specific ligand
lig2_results = df[df['compounds.id'] == 'lig2']
```


### Table Column References

Reference columns from upstream tables using tuple syntax:

```python
# RFdiffusion outputs a table with 'designed' column
rfd = RFdiffusion(contigs="50-100", num_designs=5)

# Pass column reference to downstream tool
lmpnn = LigandMPNN(
    structures=rfd,
    ligand=Ligand(codes="LIG"),  # code read from the compounds stream at runtime
    redesigned=rfd.tables.structures.designed  # Tuple: (TableInfo, "designed")
)
```

Hint: if you don't remember the table or column name, you can look it up in the ToolReference.

At execution time, the column value is resolved per-structure by ID matching.

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

**Parallel batches**: Wrap `Resources()` calls in `with Parallel():` to run them as siblings instead of sequentially. The next batch opened after the block fans in on all of them:

```python
with Pipeline("Project", "Job"):
    Resources(gpu="A100", time="2:00:00")
    seed = PDB("4AKE")

    runs = []
    with Parallel():
        for i in range(10):
            Resources(gpu="A100", time="6:00:00")   # sibling batch
            runs.append(RFdiffusion(pdb=seed, num_designs=10))

    Resources(gpu="A100", time="12:00:00")          # fan-in: waits for all 10
    sequences = ProteinMPNN(structures=Pool(runs=runs))
```

Each iteration must call `Resources()` to open its own sibling batch. `Dependencies()` and nested `Parallel()` are disallowed inside the block.

**Packing a node**: on a machine that allocates whole nodes (`machine.node_exclusive`, e.g. CSCS Daint, where every job gets all 4 GPUs and 288 CPUs whether or not it asks), one sibling job per iteration wastes most of each node. `Parallel(pack=N)` instead submits a single job whose siblings are concurrent job steps, N per node:

```python
with Pipeline("Project", "Job"):
    with Parallel(pack=4):
        Resources(gpu="gh", time="6:00:00")   # the allocation, shared by all tasks
        for lig in ligands:
            with Run():                       # one task
                poses = Boltz2(ligand=lig, ...)
                PoseBusters(structures=poses)
```

`Resources()` is called once and describes the whole allocation; `Run()` delimits one task. `pack` is tasks per node, so the node count is derived — 10 tasks at `pack=4` requests 3 nodes. Tools inside a `Run()` run in order and may mix containerized and plain tools freely; tasks run concurrently. Each task gets a proportional share of the node's cores and memory unless `Run(cpus=..., memory=...)` says otherwise. Explicit memory is useful for asymmetric layouts such as a large in-memory database server beside smaller GPU workers.

`Run()` outside a packed block raises, as does a second `Resources()` inside one.

Packing requires the machine to declare its node geometry, and refuses to engage otherwise: `machine.node_exclusive: true`, a positive `machine.cores_per_node`, and — when any task requests GPUs — `machine.gpus_per_node`. The scheduler must be SLURM, since the tasks are emitted as `srun` job steps. Each requirement is checked when the block opens, so a misconfigured machine fails immediately rather than after the allocation is granted. `cores_per_node` is what divides the node between the tasks; without it every step would silently receive a single core.

**Splitting a stream across tasks**: iterating a DataStream yields one stream per id, which is too fine when the work is meant to be shared between a few workers. `chunks()` groups instead:

```python
with Parallel(pack=4):
    Resources(gpu="gh", time="6:00:00")
    for chunk in designs.chunks(4):        # 4 chunks; sizes derived
        with Run():
            Boltz2(structures=chunk, ...)
```

`chunks(4)` gives 4 chunks with any remainder spread over the leading ones (10 items → 3/3/2/2); `chunks(size=50)` fixes the items per chunk and derives the count. Chunks carry the stream's name, format, and files, and a shared-file stream keeps every chunk pointing at the same artifact.

Lazy ids split on their deterministic outer axis and keep their `[...]` suffix, so a stream of `prot_<0..9>[_<N><A V>]` chunks into groups of `prot_N[_<N><A V>]` whose runtime fan-out is still deferred — you do not need to wait for a stream to be fully expanded to split it. An id that is lazy at the top level has no such axis and raises.

---

## Naming Runs with Suffix

`Suffix("label")` sets a label that is appended to the folder and script names of every tool created after it, until you change or clear it. It is the only way to tell two runs of the *same* tool apart: without it, `002_Boltz2/` and `005_Boltz2/` differ only by step number, and after you insert a step upstream those numbers shift.

```python
with Pipeline("Project", "Job"):
    Resources(gpu="A100")

    Suffix("apo")
    apo = Boltz2(proteins=target)                    # 001_Boltz2_apo/

    Suffix("holo")
    holo = Boltz2(proteins=target, ligands=lig)      # 002_Boltz2_holo/

    Suffix()                                          # clears it
    Panda(tables=[apo, holo])                        # 003_Panda/
```

The label lands on the output folder (`NNN_<Tool>_<suffix>/`) and on the matching `RunTime/NNN_<Tool>_<suffix>.sh`, `Logs/NNN_<Tool>_<suffix>.log`, and `ToolOutputs/NNN_<Tool>_<suffix>.json`, so a run is identifiable from any of them. `LoadMultiple(..., suffix="holo")` (see [Troubleshooting](#troubleshooting)) then reloads exactly that run.

`Suffix()` must be called inside a `Pipeline` block, and it applies from that point forward — it is a cursor, not a per-tool argument. Set it again before each group of tools you want labelled, and call `Suffix()` with no argument to go back to unlabelled names. Prefer a short, meaningful label (`apo`, `cycle3`, `batch1`) over a number, since the step number is already in the name.

---

## Grouping Outputs with Folder

`Folder("name")` is a context manager that nests the output folders of the tools created inside it under a named subdirectory. It is purely organizational — it does not change execution order, resources, batching, or dependencies. The global step counter keeps running, so the numbers still reflect true execution order:

```python
with Pipeline("Project", "Job"):
    Resources()
    Tool1()                 # 001_Tool1/
    with Folder("group"):
        Tool2()             # group/002_Tool2/
    Tool3()                 # 003_Tool3/
```

Blocks nest (`with Folder("a"): with Folder("b"):` produces `a/b/...`). Folder names must be filesystem-safe (`[A-Za-z0-9._-]`, no path separators); `.internal` is reserved for the framework.

**Downloading a group (Colab).** Bind the block and call `.download()` after it to zip the folder and trigger a browser download. This works in on-the-fly / Colab runs (where tools execute as they are added, so the files exist by the time the block closes); in submit mode it raises, since the outputs don't exist yet at that point.

```python
with Folder("Results") as results:
    af = AlphaFold(proteins=pmpnn)
    PoseBusters(structures=af, ligand="LIG")
results.download()   # zips Results/ -> Results.zip and downloads it in Colab
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
- Scheduler submission is skipped
- Output is written to `./BioPipelines/` (current directory) instead of shared storage
- The completion check mechanism is preserved, so re-running a notebook skips already-completed steps
- Re-running the cell that creates the `Pipeline()` starts a new pipeline

---

## Job Submission

**Submit to SLURM, LSF or PBS/Torque**:
```bash
biopipelines-submit /path/to/pipeline.py
biopipelines-submit /path/to/pipeline.ipynb   # directly from a notebook
```

Both `.py` scripts and `.ipynb` notebooks are supported. When a notebook is provided, code cells are automatically extracted (skipping shell commands and IPython magics) and executed as a script.

This command works from any directory as long as the biopipelines environment is activated. Alternatively, you can run the script directly from the biopipelines root:
```bash
cd biopipelines
./submit /path/to/pipeline.py
```

For interactive / inline execution (no SLURM), construct the
``Pipeline(...)`` with ``on_the_fly=True`` and run the Python file
directly. Each tool's bash script is generated and executed inline as
the corresponding wrapper is called.

### Notes on LSF and PBS/Torque schedulers

`Resources(gpu=, memory=, time=, cpus=)` is portable across schedulers; the active backend translates it best-effort into native directives. Because clusters vary, two caveats:

- **LSF** — `memory` is converted to MB and emitted as `#BSUB -M <MB>` plus a matching `#BSUB -R "rusage[mem=<MB>]"`; `time` is converted to LSF's `[HH:]MM` form (`"24:00:00"` → `"24:00"`). If your site uses different memory units (`LSF_UNIT_FOR_LIMITS`) or queue conventions, override with `lsf_options` (e.g. `Resources(q="normal")` → `#BSUB -q normal`).
- **PBS/Torque** — `mem`, `walltime`, `ncpus`, and `ngpus` are emitted as separate `#PBS -l` requests, which Torque and OpenPBS accept. Sites (often PBS Pro) that require a single chunk like `select=1:ncpus=4:ngpus=1:mem=16gb` should pass it via `pbs_options` (`Resources(select="1:ncpus=4:ngpus=1:mem=16gb")`).

A constraint-style GPU spec with no native equivalent (e.g. `high-memory`, `80GB|96GB`) falls back to a plain GPU-count request and prints a warning. Any other `Resources(**options)` kwargs are stored under the active scheduler's options key (`slurm_options` / `lsf_options` / `pbs_options`) and emitted verbatim.

Generated scripts are stemmed by scheduler: `slurm_batch*.sh`, `lsf_batch*.sh`, `pbs_batch*.sh`.

**Resubmit** existing job:

```bash
./resubmit /path/to/job/RunTime/slurm_batch1.sh   # or lsf_batch1.sh / pbs_batch1.sh
```

`resubmit` **strips the script's dependency directives** and says which ones it removed. It has to: `submit` writes the real job ids into the batch scripts when it resolves their placeholders, so the script on disk waits on the *original* run's jobs — ids that have since completed, failed, or aged out of the scheduler. Left in place they leave the resubmitted job pending on a dependency that can never be satisfied, or get it rejected outright. Pass `--keep-dependencies` for the one case where they still mean something: the parent batch is still queued or running and you want the resubmission chained behind it.

The script itself is never modified — it is piped to the scheduler — so `RunTime/` stays an accurate record of what the original run submitted.

**External dependencies** - Wait for other scheduler jobs:

```python
with Pipeline("Project", "Job", "Description"):
    Dependencies("12345678")  # Wait for job ID. Also accepts lists
    Resources(gpu="V100", time="4:00:00")
    ...
```

---

## Saving Without Submitting, and Background Services

### `Save()`

By default a pipeline submits itself when its `with` block exits. Calling `Save()` inside the block writes all the scripts and manifests but **suppresses that auto-submission**, so you can inspect (or hand-edit) what was generated before anything runs:

```python
with Pipeline("Project", "Job", description="Dry run"):
    Resources(gpu="A100")
    rfd = RFdiffusion(...)
    seqs = LigandMPNN(structures=rfd)
    Save()          # writes RunTime/, Logs/, pipeline.sh — submits nothing
```

You can then submit the generated `pipeline.sh` by hand later, or just read the scripts to check the commands. `Save()` must be called inside a `Pipeline` block; where you put it in the block does not matter, since it is the exit that would otherwise submit.

### `Service()`

Some tools are **servers**: they start up, stay running, and are consumed by a later step while still alive. `MMseqs2Server` is the canonical case. Submitting such a server as an ordinary step does not work, because the chain would wait for it to *finish* — and a server only exits on its own idle timeout.

`Service()` is a context manager for exactly this. The batch inside it keeps its normal dependency on whatever precedes it, but the chain does not wait for it to complete: the first batch *after* the block waits only for the server to be **running**. That also removes the scheduling race you get from submitting a server separately, where the server sits queued at low priority while the client's allocation burns wall-clock.

```python
with Pipeline("Project", "Job"):
    Resources(gpu="A100")
    seqs = LigandMPNN(structures=rfd)          # runs first

    with Service():
        Resources(memory="900GB", cpus=32, time="24:00:00")
        MMseqs2Server(mode="cpu")              # starts, then stays up

    Resources(memory="16GB")
    msas = MMseqs2(sequences=seqs)             # starts once the server is up

    Resources(gpu="A100")
    Boltz2(proteins=seqs, msas=msas)           # waits for the MSAs normally
```

Three rules: a `Service()` block must contain **exactly one** batch (one `Resources()` plus the daemon tool); it cannot be nested or placed inside a `Parallel()` block, and it cannot contain `Dependencies()`. Because nothing waits for the daemon's exit code, the pipeline's success never hinges on it — the server self-terminates on its idle timeout once its consumer has drained the queue.

---

## Data Management with Panda

Panda provides pandas-style table transformations:

```python
from biopipelines import Panda

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
        Panda.merge(prefixes=["apo_", "holo_"]),
        Panda.calculate({"delta": "holo_affinity - apo_affinity"})
    ]
)

# Concatenate tables vertically
combined = Panda(
    tables=[cycle0.tables.results, cycle1.tables.results],
    operations=[Panda.concat(fill="")]
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

Available operations: `filter`, `sort`, `head`, `tail`, `sample`, `rank`, `drop_duplicates`, `merge`, `concat`, `calculate`, `zscore`, `groupby`, `select_columns`, `drop_columns`, `rename`, `fillna`, `pivot`, `melt`, `average_by_source`

---

## When IDs Disappear: the `missing` Table

A pipeline step can finish successfully with **fewer ids than it started with**. A filter drops rows, a structure fails to converge, a compound has no SMILES — and the ids concerned simply stop appearing downstream. This is normal and intended, but it is also the usual answer to "my run finished, so why do I have 260 designs instead of 300?"

The mechanism is a standard table called `missing`, written to `<tool folder>/tables/missing.csv` by over forty tools. Its schema is fixed:

| Column | Meaning |
|---|---|
| `id` | The id that was removed |
| `removed_by` | The step that removed it, as `<NNN>_<Tool>` (e.g. `005_Panda`) |
| `kind` | `filter` if the tool dropped it on purpose, `failure` if the tool tried and failed |
| `cause` | A short human-readable reason |

Two things follow from this being a real table rather than a log message.

**It propagates.** Each tool merges its inputs' `missing` tables into its own, de-duplicating by id, so the `missing.csv` of the *last* step in a chain is a complete account of everything lost anywhere upstream, with the step that lost it named. You do not have to walk back through the pipeline — read the last one:

```python
import pandas as pd

if "missing" in final_step.tables:                       # not declared when nothing upstream lost an id
    lost = pd.read_csv(final_step.tables.missing.info.path)
    print(lost.groupby(["removed_by", "kind"]).size())   # who dropped how many, and why
    print(lost[lost["kind"] == "failure"])               # the ones that are actually errors
```

**It is what keeps the completion check honest.** Without it, a tool that declared 300 outputs and wrote 260 files would be reported FAILED. The check reads `missing.csv` and excuses exactly the ids it accounts for: rows propagated from an upstream step, and rows this step marked `kind="filter"`. A row this step itself wrote with `kind="failure"` is **not** excused — that is a genuine failure and the step is still reported FAILED. So a green run with a populated `missing.csv` means "these ids were deliberately dropped", not "these ids broke silently".

When you want to know where your ids went, read `missing.csv` first and the logs second.

---

## Reading Tables Back in Python

Every tool exposes its output tables through `.tables`, and every table carries its own metadata. The path is the thing you actually want — hand it to `pandas` and analyze the run like any other CSV:

```python
import pandas as pd

df = pd.read_csv(boltz.tables.confidence.info.path)
best = df.sort_values("confidence_score", ascending=False).head(10)
```

`.tables` supports the operations you would expect of a mapping, which is how you explore a tool's output without looking it up in the reference:

```python
list(boltz.tables.keys())              # ['structures', 'confidence', 'sequences', 'msas', 'affinity', 'compounds']
"affinity" in boltz.tables             # True
boltz.tables["confidence"]             # the path, as a string
boltz.tables.confidence.info.columns   # ['id', 'input_file', 'confidence_score', 'ptm', 'iptm', ...]
boltz.tables.confidence.info.name      # 'confidence'
```

Note the two access styles: `boltz.tables["confidence"]` returns the **path string**, while `boltz.tables.confidence` returns the **table object** whose `.info` holds `path`, `columns`, `name`, and `description`. Streams work the same way through `.streams` (`list(boltz.streams.keys())` → `['structures', 'sequences', 'compounds', 'msas']`).

For a run that has already finished, `Load()` reaches the same tables on disk — but how you get at them depends on where you call it. **Inside** a `Pipeline` block it returns the usual object, so `.tables` works exactly as above:

```python
with Pipeline("Project", "Analysis"):
    Resources()
    prev = Load("/path/to/job/003_Boltz2")
    df = pd.read_csv(prev.tables.confidence.info.path)
```

**Outside** a pipeline — the common case in a notebook or a one-off analysis script — `Load()` returns the tool object itself, which has no `.tables`. Go through `get_output_files()`, whose `"tables"` entry is a plain dict of table objects:

```python
from biopipelines import Load, list_tables, get_table_path, table_exists

prev = Load("/path/to/job/003_Boltz2")
tables = prev.get_output_files()["tables"]

list_tables(tables)                                  # ['structures', 'confidence', ...]
table_exists(tables, "affinity")                     # True / False
df = pd.read_csv(get_table_path(tables, "confidence"))
```

### The Five Table Helpers

`get_table`, `get_table_path`, `get_indexed_table`, `list_tables` and `table_exists` are exported from `biopipelines`, and all five accept **any container you might be holding**, not just the dict above:

| You have | Where it comes from |
|---|---|
| a plain `{name: table}` dict | `Load(folder).get_output_files()["tables"]` |
| a `TableContainer` | any tool's `.tables` |
| a `StandardizedOutput` | what a tool constructor returns inside a `Pipeline` block |
| a `ToolOutput` | the framework's own wrapper around a registered tool (`ToolOutput(config)`); accepted so code holding one need not reach for its `.output` |

All five route through one resolver, so a container that one helper accepts is accepted by all of them — there is no need to unwrap `.tables` yourself, and no helper is fussier than its neighbors:

```python
get_table_path(boltz, "confidence")           # the StandardizedOutput a tool returns
get_table_path(boltz.tables, "confidence")    # its TableContainer
get_table_path(tables, "confidence")          # the plain dict from get_output_files()
```

**Anything else raises, and the message says what you passed.** Hand a helper a list, a string, or an unrelated object and you get a `TypeError` naming the type and listing the four accepted shapes. That includes `table_exists`, which raises rather than answering `False` — `False` here would read as "this run has no such table", which is a wrong answer rather than an error, and wrong answers are the ones that cost you an afternoon.

**Ordinary tables and indexed tables have separate accessors.** Most tables are one CSV, so `get_table_path(source, name)` returns its path. A tool that writes *one table per input ID* declares the collection as a single named entry (an `IndexedTableContainer`), and there is no single path to return — you have to say which entry you want:

```python
list_tables(prev_output)                             # ['rmsf', ...]
get_indexed_table(prev_output, "rmsf", "1gfl")       # -> the TableInfo for that ID
get_indexed_table(prev_output, "rmsf", "1gfl").info.path

get_table(prev_output, "rmsf")                       # -> the IndexedTableContainer itself
```

Call the wrong one and the error points at the right one, in both directions:

- `get_table_path(source, "rmsf")` on an indexed collection → `TypeError: Table 'rmsf' is an IndexedTableContainer with 2 entries. Use get_indexed_table(source, 'rmsf', entry_id) to get one entry's path.`
- `get_indexed_table(source, "scores", "1gfl")` on an ordinary table → `TypeError: Table 'scores' is a TableInfo, not an IndexedTableContainer. Use get_table_path(source, 'scores') for its single path.`

`get_table` is the one helper that is indifferent to the distinction: it returns whatever is stored under that name, a `TableInfo` or an `IndexedTableContainer`, so use it when you want to inspect (`.info.columns`, `.ids`) rather than to read a file. A name that is not there raises `KeyError` listing the names that are.

This is the whole loop: run the pipeline, `Load` the step you care about, read its tables with `pandas`, and check `missing.csv` for anything that is not there.

---

## Filesystem Structure

Each tool's output folder follows a predictable sub-layout so you always
know where to look for a given artefact. Stream map_tables live inside
their own stream's folder; standalone tables live under `tables/`.

```
<biopipelines_output>/<project>/<job>_<NNN>/
├── RunTime/                    # Execution scripts and the run page
│   ├── pipeline.sh
│   ├── pipeline.html           # self-contained record of the run — open it in a browser
│   ├── pipeline_graph.json     # the same content as data, for regenerating the page
│   ├── 001_<tool>.sh
│   ├── 002_<tool>.sh
│   └── .internal/              # scripts for auto-generated internal tools
│       └── 001_<tool>.sh
├── Logs/                       # Execution logs
│   ├── 001_<tool>.log
│   ├── 002_<tool>.log
│   └── .internal/
│       └── 001_<tool>.log
├── ToolOutputs/                # Tool output predictions (JSON manifests)
│   ├── 001_<tool>.json
│   ├── 002_<tool>.json
│   └── .internal/
│       └── 001_<tool>.json
├── 001_<Tool>/                 # Tool outputs — canonical sub-layout
│   ├── .expected_outputs.json  # status/manifest at the folder root
│   ├── _configuration/         # config-time inputs: JSONs, YAMLs, CSVs
│   ├── _execution/             # raw model dumps (CIFs, intermediate files)
│   ├── <stream_name>/          # one folder per declared output stream
│   │   ├── <stream>_map.csv    # map_table lives WITH its stream
│   │   └── <id>.pdb / .fasta / ...
│   ├── tables/                 # standalone TableInfo CSVs
│   └── _extras/                # catch-all (plots, logs, sessions)
├── 002_<Tool>/
│   └── ...
├── <folder>/                   # Folder("<folder>") groups public tool outputs
│   └── 003_<Tool>/
└── .internal/                  # framework-owned internal tools (hidden)
    └── 001_<Tool>/
```

**Where to find each kind of output:**

| You want | Look under |
|---|---|
| Per-ID files (PDBs, FASTAs, images) | `<stream>/` (e.g. `structures/`, `sequences/`) |
| A stream's map_table | `<stream>/<stream>_map.csv` (or `<stream>/<stream>.csv` for content-bearing streams like Sequence) |
| Standalone metric/analysis tables | `tables/<name>.csv` |
| Input JSONs / YAMLs sent to the tool's CLI | `_configuration/` |
| Raw model dumps (Boltz's `boltz_results_*`, ColabFold's Folding dump, etc.) | `_execution/` |
| The tool's quick-glance status manifest | `.expected_outputs.json` (root) |
| The tool's own console log | `_log` (root) — the same stream as `Logs/<NNN>_<Tool>.log`, kept beside the outputs so you don't have to leave the folder |
| Outputs grouped under `Folder("x")` | `x/<NNN>_<Tool>/` |
| Auto-generated internal tools (e.g. a `Ligand` from `ligand="LIG"`) | `.internal/<NNN>_<Tool>/` |

### Visualizing one step: `bp-visualize`

`pipeline.html` covers the whole run. `bp-visualize` covers one step, and takes selection and ordering arguments the run page has no way to express:

```bash
bp-visualize /path/to/job/003_Boltz2
bp-visualize /path/to/job/003_Boltz2 --descending confidence.confidence_score --max-items 5
bp-visualize /path/to/job/002_ESMFold --ascending confidence.plddt --streams structures
```

It writes `<tool folder>/_extras/<step>_view.html` by default (`-o` to choose), and the page is self-contained exactly as `pipeline.html` is — the vendored py3Dmol copy is inlined, so the 3D viewers work from a `file://` URL with no network. Copy the one file off the cluster and open it.

Two properties make it usable mid-run: it reads the step's exported `ToolOutputs/<step>.json` plus whatever is on disk, so it needs neither the scheduler nor the pipeline script, and it does not care whether the rest of the job has finished. Run it on a login node against a job that is still queued.

| Option | Meaning |
|---|---|
| `--descending TABLE.COLUMN` / `--ascending TABLE.COLUMN` | Order items by a column of one of the step's tables |
| `--max-items N` | Render at most N items; **default 10**. `0` means no cap, leaving each renderer its own sampling |
| `--ids id1,id2` | Render only these ids |
| `--streams a,b` / `--tables a,b` | Render only these streams / tables |
| `-o PATH` | Where to write the page |
| `--open` | Open the page in a browser (for local runs, not a login node) |
| `--allow-external` | Keep a renderer's CDN `<script>` when no vendored copy exists |

**The sort key must be qualified.** `--descending confidence.plddt`, not `--descending plddt`: one tool's tables can carry the same column name twice, and picking one silently is a wrong answer rather than an error. A bare name is refused and the message names the qualified spellings that would work. The available columns are in `docs/tool_reference.md`, and in the step's own `.expected_outputs.json`.

Ordering and capping apply to the streams *and* the tables on the page, so a top-5-by-pLDDT view shows the same five items in the viewer and in every table beside it.

The default is **10**, not the run page's 5. A single-step page carries one step's structures rather than every node's, and 5 hides exactly one item of a 6-item stream — the least useful place to stop. Pass `--max-items 0` for no cap.

**A cap alongside an ordering means "the top N of what was ranked."** It applies only to the streams the sort table actually keys — a Boltz2 step ordered by `confidence.confidence_score` caps its structures and leaves its `sequences`, `compounds` and `msas` whole, since the confidence table keys structures and cutting the others would drop items on no criterion at all. The page says which streams were spared and why. With no ordering, `--max-items` is purely a page-size limit and applies to every stream.

### The run page

`RunTime/pipeline.html` is a record of the run you can open in a browser. It lays the pipeline out as one lane per batch with that batch's resources, one collapsible card per step carrying its paths, tool version, environment and completion marker, and arrows for the dataflow between steps. A donut summarises completed against failed, and a Machine panel records the scheduler and environment manager the run actually used.

It is **self-contained**: no network, no sibling files. Copy that one file off a cluster and it still works, which is what it is for.

Two things worth knowing:

- **It is written twice.** `save()` writes it at configuration time, when nothing has run yet — useful for checking the wiring, but every step reads "pending". Refresh it afterwards with `regenerate_pipeline_page("<job>/RunTime")` to fill in the results. On-the-fly runs refresh it after each step automatically.
- **The 3D viewers work with no network.** The structure and grid renderers need the py3Dmol library, and a copy is vendored in the repo at `biopipelines/renderers/vendor/3Dmol-min.js`, inlined once per page however many structure nodes it has. A page with viewers is around 695 KB against 28 KB without — still one file you can copy off a cluster and open anywhere. If that vendored copy is missing, the page falls back to the CDN when you pass `allow_external=True`, and to a metadata table otherwise, naming what it dropped either way.

`allow_external` therefore does **not** decide whether structures render. It means "allow renderer output that fetches from the network", which now only matters for a renderer that fetches *data*.

An arrow appears only where the consuming tool holds an object stamped with its producer. A stream derived inside `DataStream` itself (`stream.chunks(...)`, iterating a bare stream) is a fresh object with no back-reference, so its edge is left out rather than guessed — the graph shows what the framework could prove, not everything that is true.

Each public tool carries the same step number across its output folder, its `RunTime/`/`Logs/` scripts, and its `ToolOutputs/` manifest, so `002_LigandMPNN/` pairs with `RunTime/002_LigandMPNN.sh`. Internal tools (auto-generated, e.g. a `Ligand` from `ligand="LIG"`) are numbered separately under their own `.internal/` subdirectories and don't consume a public step number.

Configure paths and environments in the repository root's `config.<variant>.yaml` — there is no single `config.yaml`. One file ships per site variant (`config.cluster.yaml`, `config.local.yaml`, `config.container.yaml`, `config.daint.yaml`, `config.colab.yaml`); a gitignored `.config.<variant>.yaml` overlay next to it holds your machine-specific paths and wins over the committed defaults. Edit them with the `bp-config` command rather than by hand. The active variant is detected automatically and can be forced with `BIOPIPELINES_CONFIG_VARIANT=<variant>`.

---

## Reproducibility

BioPipelines distinguishes two layers of dependency information.

**Install-time floors.** `pyproject.toml` and `environments/biopipelines.yaml`
declare `>=` minimum versions for the framework's Python dependencies. They
guarantee a clean install lands on a known-good baseline but leave room for
forward upgrades; the per-tool `environments/<tool>.<variant>.yaml` files stay
loose for the same reason — site-specific CUDA drivers and ML-stack
constraints make a single locked version unportable across clusters.

**Per-run records.** A pipeline constructed with `debug=True` writes a complete
runtime snapshot under `<output>/_debug_capture/` when the job runs:

- `_debug_capture/environments/<env>.yaml` — `mamba env export --no-builds`
  per environment used in the pipeline.
- `_debug_capture/environments/<env>.pip.txt` — `pip freeze` per environment.
- `_debug_capture/system/` — `uname`, `nvidia-smi`, scheduler version,
  container-runtime version.

```python
with Pipeline(project="Project", job="Job", debug=True):
    ...
# the pipeline.sh now exports BIOPIPELINES_DEBUG=1 and snapshots the runtime
```

This pair (loose `>=` for install, exact freeze per run) lets a third party
reproduce a specific result by recreating the captured environment, while not
forcing every clean install to use bleeding-edge versions. A reference
artefact captured on UZH S3IT lives at
`docs/reviewer_evidence/a6_debug_capture/environments/`.

---

## Troubleshooting

**Path errors**: Run from BioPipelines root directory.

**Tool installation issues**: The variety of HPC configurations makes it very difficult to define a tool installer that is portable accross all. Please refer to the official documentation (references are in the README file). If installation fails on Colab, please open a Git Issue.

**Job killed/00MM**: Most likely the jobs requires more resources. Adjust GPU/memory in `Resources()`.

**Missing files**: Check `Logs/<NNN>_<tool>.log`

**Local output**: Write results to the current directory instead of the config-defined path:

```python
with Pipeline("Test", "Debug", "Testing", local_output=True):
    ...
```

Note: `local_output` defaults to `True` automatically when `on_the_fly` is enabled (i.e., in Jupyter notebooks).

**Load previous outputs**:

```python
from biopipelines import Load, LoadMultiple

# Single output (pass the tool's output folder)
prev = Load("/path/to/job/001_Boltz2")

# Multiple outputs by tool name (pass the job folder)
all_boltz = LoadMultiple("/path/to/job/", tool="Boltz2")
```
