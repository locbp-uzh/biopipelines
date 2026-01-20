# BioPipelines User Manual

## Index

### BioPipelines

- [Architecture Overview](#architecture-overview)
  - [What is BioPipelines?](#what-is-biopipelines)
  - [Installation](#installation)
  - [Repository Structure](#repository-structure)
  - [Core Concepts](#core-concepts)
  - [Resources](#resources)
  - [Job Submission](#job-submission)
  - [Filesystem Structure](#filesystem-structure)
  - [Environment Management](#environment-management)
- [Tool I/O Reference Guide](#tool-io-reference-guide)
  - [Overview](#overview)
  - [Table Organization](#table-organization)
  - [Table Column References](#table-column-references)
- [Troubleshooting](#troubleshooting)
  - [Common Issues](#common-issues)
  - [Debug Mode](#debug-mode)

## Architecture Overview

### What is BioPipelines?

BioPipelines is a Python framework that generates bash scripts for bioinformatics workflows. It does not execute computations directly - instead, it predicts the filesystem structure and creates scripts that will be executed on SLURM clusters.

### Installation

Login to your cluster via terminal or the website ([S3IT Apps](https://apps.s3it.uzh.ch) > Clusters > Shell access). In your home directory, clone the biopipelines repository:
```bash
git clone https://gitlab.uzh.ch/locbp/public/biopipelines
```

### Repository Structure

The repository is organized as follows:

```
/biopipelines/
├── Docs/                       # Markdown files with usage information
│   ├── UserManual.md           # General biopipelines manual
│   ├── ToolReference.md        # Reference for installation, input, and expected output from each tool
├── Environments/               # Definition of most important and hard-to-replicate conda environments 
├── ExamplePipelines/           # Usage cases
├── HelpScripts/                # Scripts that run at slurm runtime.
│   ├── pipe_<tool&function>.py # Generic pipe scripts. 
├── Ligands/                    # Folder for ligand structural files. The tool Ligand will download here.
├── MyPipelines/                # Folder for user-written pipelines.
├── PDBs/                       # Folder for protein structural files. The tool PDB will download here.
├── PipelineScripts/            # Scripts that run at pipeline runtime, thus generating and coordinating other scripts.
│   ├── pipeline.py             # Contains fundamental elements like Pipeline, Resources, ...
│   ├── base_config.py          # Base class for tools and tool outputs.
│   ├── <tool>.py               # Definition and declaration of a tool.
│   ├── ...                     # Other utilities like folders.py, config_manager.py, ...
├── config.yaml                 # Here the user can choose which environments to use and define where the models are.
├── [re]submit                  # Bash scripts for [re]submission of pipeline scripts

```

### Core Concepts

**Pipeline**: The main coordinator that orchestrates tools, manages the folder structure, and switches environments. Pipelines can be used with context manager syntax (preferred) or traditional syntax. Both generate a unique folder /shares/*USER*/MyProject/JobName_*NNN* where all the output will be. **WARNING**: Don't include blank spaces in the project and job names.

```python
# Context manager syntax (preferred)
from PipelineScripts.pipeline import *

with Pipeline("MyProject", "JobName", "Description"):
    Resources(gpu="V100", time="24:00:00", memory="16GB")
    tool1 = Tool1(param1=value1,
                  param2=value2)
    tool2 = Tool2(input=tool1)
# Auto generates SLURM script on exit

# Traditional syntax (equivalent, also supported)
from PipelineScripts.pipeline import Pipeline

pipeline = Pipeline("MyProject", "JobName", "Description")
pipeline.resources(gpu="V100", time="24:00:00", memory="16GB")
tool1 = pipeline.add(Tool1(param1=value1, param2=value2))
tool2 = pipeline.add(Tool2(input=tool1))
pipeline.slurm()
```

**Tools**: Individual bioinformatics operations like running models (RFdiffusion, LigandMPNN, Boltz2, ...) or analyzing results (Filter, SelectBest, ...) that generate bash scripts and predict their outputs. Tools return an object containing predictions of the filesystem after SLURM execution. The prediction can be used as input in subsequent tools. One can access the prediction with the default `print(<prediction>)` method.

```python
from PipelineScripts.pipeline import *
from PipelineScripts.rfdiffusion import RFdiffusion
with Pipeline("TestProject","Test","Some test"):
  rfd = RFdiffusion(contigs="50-100", num_designs=5)
  print(rfd)
"""
structures:
    – '<output_folder>/Test_1.pdb'
    – ...
    – '<output_folder>/Test_5.pdb'
structure_ids:
    – Test_1
    – ...
    – Test_5
tables:
    $structures (id, source_id, pdb, fixed, designed, contigs, time, status):
        – '<output_folder>/rfdiffusion_results.csv'
output_folder:
    – 'path/to/001_RFdiffusion'
main:
    – '<output_folder>/rfdiffusion_results.csv'
#Aliases
structures=pdbs
"""
```
The predicted output can then be used by other tools as input, thus enabling further prediction prior to the actual generation of the files:
```python
from PipelineScripts.pipeline import *
from PipelineScripts.rfdiffusion import RFdiffusion
from PipelineScripts.protein_mpnn import ProteinMPNN
with Pipeline("TestProject","Test","Some test"):
  rfd = RFdiffusion(contigs="50-100", num_designs=5)
  pmd = ProteinMPNN(structures=rfd,num_sequences=2)
  print(pmd)
"""
sequences:
    – '<output_folder>/Test_queries.csv'
sequence_ids:
    – Test_1_1
    – Test_1_2
    – ...
    – Test_5_1
    – Test_5_2
tables:
    $sequences (id, source_id, source_pdb, sequence, score, seq_recovery, rmsd):
        – '<output_folder>/Test_queries.csv'
output_folder:
    – 'path/to/002_ProteinMPNN'
fa_files:
    – '<output_folder>/seqs/Test_sequences.fa'
main:
    – '<output_folder>/proteinmpnn_results.csv'
queries_fasta:
    – '<output_folder>/Test_queries.fasta'
seqs_folder:
    – '<output_folder>/seqs'
#Aliases
sequences=queries_csv
"""
```

All the outputs are standardized such that the identity of the previous tool is not needed for the prediction to be used as input, and in general tools have to be developed agnostic of previous tool identities. Furthermore, this systems allows:

1. Verification of the success or failure of a given tool. At slurm runtime, after running a tool, the pipeline coordinator check for the presence of the predicted output and creates a file `<NNN>_<Tool>_COMPLETED` or `<NNN>_<Tool>_FAILED`. One can resubmit the same slurm bash script(for example, if the time ran out) and the completed steps will be skipped.
2. Standardized saving and loading of tool outputs. All the predictions are saved in the folder `ToolOutputs` within the job folder, and can be loaded with `LoadOutput` (single file) or `LoadOutputs` (multiple files). In the following example, rfd behaves identically to the one in the previous snippet:
```python
from PipelineScripts.pipeline import *
from PipelineScripts.load_output import LoadOutput, LoadOutputs
from PipelineScripts.protein_mpnn import ProteinMPNN

with Pipeline("TestProject","Test","Some test"):
  # Load single output
  rfd = LoadOutput('/path/to/ToolOutputs/001_RFdiffusion.json')
  pmd = ProteinMPNN(structures=rfd, num_sequences=2)

  # Or load multiple outputs at once
  all_merge_tables = LoadOutputs('/path/to/job', tool="MergeTables")
```

### Resources

Resources are set before using tools as follows:

```python
with Pipeline("Project", "Job", "Description"):
    Resources(gpu="A100", memory="32GB", time="24:00:00")              # Specific model
    Resources(gpu="32GB|80GB|96GB", memory="32GB", time="24:00:00")   # V100, A100, H100, or H200
    Resources(memory="128GB", time="24:00:00", cpus=32)                # CPU-only
```

**GPU parameter options:**
- `"T4"`, `"L4"`, `"V100"`, `"A100"`, `"H100"`, `"H200"` - Specific GPU models
- `"24GB"`, `"32GB"`, `"80GB"`, `"96GB"`, `"32GB|80GB|96GB"` - Memory-based selection
- `"gpu"` or `"any"` - Any available GPU
- `"high-memory"` - Equivalent to `"32GB|80GB|96GB"`
- Omit parameter for CPU-only jobs

### Job submission

Login to your cluster via terminal or the website ([S3IT Apps](https://apps.s3it.uzh.ch) > Clusters > Shell access). Move to the biopipelines repository, then submit using the submit script:

```bash
cd biopipelines

./submit /path/to/<pipeline.py>
```

The script will activate the biopipelines environment for proper execution of the pipeline. If not available or not up-to-date, it will install it. If installation fails due to memory errors (e.g. std::bad_alloc) you can first run a console with more memory, and then submit:

```bash
srun --mem=8G --time=1:00:00 --pty bash

./submit /path/to/<pipeline.py>
```

If for any reason you want to resubmit a pipeline another time, use the resubmit script with the path to the slurm shell script:
```bash
cd biopipelines

./resubmit /shares/<group>/<user>/BioPipelines/<project>/<job>/RunTime/slurm.sh
```

Alternatively, running the pipeline script with python will result in the following output:
```bash
"""
==============================Job============================== # Suggested job name
Pipeline: Test (002)
==============================Slurm Script============================== # Slurm script for manual submission
#!/usr/bin/bash
#SBATCH --mem=16GB
#SBATCH --time=24:00:00
#SBATCH --output=job.out
#SBATCH --begin=now+0hour
# Make all files group-writable by default
umask 002
module load mamba singularityce
# Execute pipeline
path/to/RunTime/pipeline.sh    # Do not run this: it is part of the slurm script, it will most likely fail from console unless you have sufficient resources (e.g. GPU)
==============================SBATCH============================== # Simple submission using sbatch command
sbatch --job-name=Test --output path/to/RunTime/slurm.out path/to/RunTime/slurm.sh
"""
```

**Multiple job submission**
When submitting jobs with `./submit /path/to/<pipeline.py>`, all pipelines created within the script will be submitted. For example, one can test the influence of a given parameter on the pipeline as follows:

```python
#MyPipelines/beta_test.py
#imports
for beta in [0, 1, 10, 100]:
  with Pipeline("Diffusion",
                f"Rhodamine_Beta{beta}",
                f"Run diffusion simulation with beta value {beta}"):
      #tools
```
Running `./submit MyPipelines/beta_test.py` will then result in the submission of four jobs, *Rhodamine_Beta0*, *Rhodamine_Beta1*, *Rhodamine_Beta10*, *Rhodamine_Beta100*.

**Batch dependencies**
Calling `Resources()` multiple times within a pipeline creates sequential batches with automatic dependency management. Each batch waits for the previous one to complete successfully before starting:

```python
with Pipeline("Project", "Job", "Description"):
    Resources(gpu="V100", time="4:00:00")    # Batch 1: GPU tools
    tool1 = RFdiffusion(...)
    tool2 = ProteinMPNN(...)

    Resources(time="2:00:00")                # Batch 2: CPU-only analysis (waits for Batch 1)
    tool3 = Filter(...)
```

The `./submit` script automatically chains batch jobs using `--dependency=afterok:<jobid>`, ensuring proper execution order and resource optimization.

### Filesystem Structure

After the execution, the filesystem will look somewhat like this:

```
/shares/<user>/BioPipelines/<project>/<job>_<NNN>/
├── RunTime/                    # Execution scripts only
│   ├── pipeline.sh             # Main orchestration script
│   ├── 001_<tool 1>.sh         # Tool-specific script
│   └── 002_<tool 2>.sh
│   ├── ...
├── Logs/                       # Execution logs
│   ├── 001_<tool 1>.log
│   └── 002_<tool 2>.log
│   ├── ...
├── ToolOutputs/                # Tool output predictions
│   ├── 001_<tool 1>.json
│   └── 002_<tool 2>.json
│   ├── ...
├── 001_<Tool 1>/               # Tool outputs only
│   ├── <Some temporary folder>
│   ├── structure_1.pdb
│   ├── structure_2.pdb
│   ├── ...
│   └── <results>.csv
├── 002_LigandMPNN/
│   ├── <Some temporary folder>
│   ├── <sequence queries>.csv
│   └── <sequence queries>.fasta
│   ├── ...
```

Base folder paths (like `/shares/<group>`, tool data folders, etc.) are configured in `config.yaml` at the repository root. Edit this file to match your cluster's filesystem layout.

### Environment Management

Most tools require a conda environment to run. Default environments are defined in `config.yaml` at the repository root. To change default environments edit the `environments` section in `config.yaml`. Overriding a specific tool in a given pipeline is not recommended, but can be done using the `env` parameter in `pipeline.add(...)`:

```python
# Use default environment from config.yaml
tool1 =Tool1(...)
# Override with custom environment
tool2 = pipeline.add(Tool2(...), env="CustomEnv") #pipeline must be defined
```

Edit `config.yaml` in the repository root to change which conda environment each tool uses by default. Tools with `null` environments don't require conda activation.

## Tool I/O Reference Guide

### Overview
All tools return outputs through a unified `StandardizedOutput` object with general structure:
```python
structures        # List[str] - File paths e.g. /path/to/structure_1.pdb, ...
structure_ids     # List[str] - Structure identifiers e.g. structure_1, ...
compounds         # [str]     - Path to csv file with columns id, format, smiles, ccd
compound_ids      # List[str] - Compound identifiers
sequences         # [str]     - Path to csv file with columns id, sequence, ...
sequence_ids      # List[str] - Sequence identifiers
msas              # [str]     - Path to csv file with columns id, sequence, msa, ...
msa_ids           # List[str] - MSA identifiers
tables            # List[DI]  – Contains objects of class TableInfo
tool_folder       # [str]     – Path to output directory
```

### Table Organization
Table names are tool-dependent. Accessing a table with dot notation will return the path to the csv file. For accessing the metadata, use the info function.

```python
# Simple path
path = <tool output>.tables.analysis  # path/to/analysis.csv

# Rich metadata access
info = <tool output>.tables.info("analysis")
print(info.path)         # path/to/analysis.csv
print(info.columns)      # ["id", ...]
print(info.description)  # "Results from ..."
print(info.count)        # Number of entries
```

### Table Column References
Tools can explicitly reference specific columns from upstream tables using dot notation <tool output>.tables.<table name>.<column name>, which returns the tuple (TableInfo,str). 

```python
lmpnn = LigandMPNN(
    structures=rfdaa, #output from RFdiffusionAllAtom
    ligand='LIG',
    redesigned=rfdaa.tables.structures.designed     #structures has columns id, fixed, designed, ... 
)
```

## Troubleshooting

### Common Issues

**Path errors**: Ensure running from BioPipelines root directory.

**Environment issues**: Check conda environment availability and tool compatibility.

**Resource limits**: Adjust GPU/memory requirements in `pipeline.resources()`.

**Missing files**: Check logs in `Logs/<NNN>_<tool>.log` files.


### Debug Mode

You can test the filesystem prediction locally by instantiating the pipeline object with the debug flag set as True:

```python
pipeline = Pipeline(..., debug=True)
```
