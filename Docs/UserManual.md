# BioPipelines User Manual

## Index

### BioPipelines

- [Architecture Overview](#architecture-overview)
  - [What is BioPipelines?](#what-is-biopipelines)
  - [Core Concepts](#core-concepts)
  - [Job Submission](#job-submission)
  - [Filesystem Structure](#filesystem-structure)
  - [Environment Management](#environment-management)
- [Tool I/O Reference Guide](#tool-io-reference-guide)
  - [Overview](#overview)
  - [Datasheet Organization](#datasheet-organization)
  - [Datasheet Column References](#datasheet-column-references)
- [Troubleshooting](#troubleshooting)
  - [Common Issues](#common-issues)
  - [Debug Mode](#debug-mode)

### Tool Reference

- [BioPipelines User Manual](#biopipelines-user-manual)
  - [Index](#index)
    - [BioPipelines](#biopipelines)
    - [Tool Reference](#tool-reference)
  - [Architecture Overview](#architecture-overview)
    - [What is BioPipelines?](#what-is-biopipelines)
    - [Core Concepts](#core-concepts)
    - [Job submission](#job-submission)
    - [Filesystem Structure](#filesystem-structure)
    - [Environment Management](#environment-management)
  - [Tool I/O Reference Guide](#tool-io-reference-guide)
    - [Overview](#overview)
    - [Datasheet Organization](#datasheet-organization)
    - [Datasheet Column References](#datasheet-column-references)
  - [Troubleshooting](#troubleshooting)
    - [Common Issues](#common-issues)
    - [Debug Mode](#debug-mode)
- [Complete BioPipelines Tool Reference](#complete-biopipelines-tool-reference)
  - [Structure Generation](#structure-generation)
    - [RFdiffusion](#rfdiffusion)
    - [RFdiffusionAllAtom](#rfdiffusionallatom)
    - [BoltzGen](#boltzgen)
  - [Sequence Design](#sequence-design)
    - [ProteinMPNN](#proteinmpnn)
    - [LigandMPNN](#ligandmpnn)
    - [MutationComposer](#mutationcomposer)
    - [SDM (SiteDirectedMutagenesis)](#sdm-sitedirectedmutagenesis)
    - [Fuse](#fuse)
    - [StitchSequences](#stitchsequences)
  - [Structure Prediction](#structure-prediction)
    - [AlphaFold](#alphafold)
    - [ESMFold](#esmfold)
    - [Boltz2](#boltz2)
    - [RF3](#rf3)
    - [OnionNet](#onionnet)
    - [OnionNet2](#onionnet2)
  - [Analysis](#analysis)
    - [DynamicBind](#dynamicbind)
    - [ResidueAtomDistance](#residueatomdistance)
    - [PLIP (Protein-Ligand Interaction Profiler)](#plip-protein-ligand-interaction-profiler)
    - [DistanceSelector](#distanceselector)
    - [SelectionEditor](#selectioneditor)
    - [ConformationalChange](#conformationalchange)
    - [MutationProfiler](#mutationprofiler)
    - [ProteinLigandContacts](#proteinligandcontacts)
  - [Data Management](#data-management)
    - [Filter](#filter)
    - [SelectBest](#selectbest)
    - [RemoveDuplicates](#removeduplicates)
    - [MergeDatasheets](#mergedatasheets)
    - [ConcatenateDatasheets](#concatenatedatasheets)
    - [SliceDatasheet](#slicedatasheet)
    - [ExtractMetrics](#extractmetrics)
    - [AverageByDatasheet](#averagebydatasheet)
  - [Utilities](#utilities)
    - [LoadOutput](#loadoutput)
    - [MMseqs2Server](#mmseqs2server)
    - [CompoundLibrary](#compoundlibrary)
    - [PDB](#pdb)
    - [PyMOL](#pymol)

## Architecture Overview

### What is BioPipelines?

BioPipelines is a Python framework that generates bash scripts for bioinformatics workflows. It does not execute computations directly - instead, it predicts the filesystem structure and creates scripts that will be executed on SLURM clusters.

### Core Concepts

**Pipeline**: The main coordinator that orchestrates tools, manages the folder structure, and switches environments. The following code generates a unique folder /shares/*USER*/MyPipeline/JobName_*NNN* where all the output will be. **WARNING**: Don't include blank spaces in the pipeline and job names. 

```python
from PipelineScripts.pipeline import Pipeline
pipeline = Pipeline("MyPipeline", "JobName", "Description")
```

**Tools**: Individual bioinformatics operations like running models (RFdiffusion, LigandMPNN, Boltz2, ...) or analyzing results (Filter, SelectBest, ...) that generate bash scripts and predict their outputs. They are added to the pipeline by calling `pipeline.add(<Tool>(<Parameters>))`, which returns an object containing predictions of the filesystem after SLURM execution. The prediction can be used as input in subsequent tools. One can access the prediction with the default `print(<prediction>)` method.

```python
from PipelineScripts.pipeline import Pipeline
from PipelineScripts.rfdiffusion import RFdiffusion
pipeline = Pipeline("Pipeline","Test","Some test")
rfd = pipeline.add(RFdiffusion(contigs="50-100", num_designs=5))
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
datasheets:
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
from PipelineScripts.pipeline import Pipeline
from PipelineScripts.rfdiffusion import RFdiffusion
from PipelineScripts.protein_mpnn import ProteinMPNN
pipeline = Pipeline("Pipeline","Test","Some test")
rfd = pipeline.add(RFdiffusion(contigs="50-100", num_designs=5))
pmd = pipeline.add(ProteinMPNN(structures=rfd,num_sequences=2))
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
datasheets:
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
2. Standardized saving and loading of tool outputs. All the predictions are saved in the folder `ToolOutputs` within the job folder, and can be loaded with the tool LoadOutput. In the following example, rfd behaves identically to the one in the previous snippet:
```python
from PipelineScripts.pipeline import Pipeline
from PipelineScripts.load_output import LoadOutput
from PipelineScripts.protein_mpnn import ProteinMPNN
pipeline = Pipeline("Pipeline","Test","Some test")
rfd = pipeline.add(LoadOutput('/path/to/ToolOutputs/001_RFdiffusion.json'))
pmd = pipeline.add(ProteinMPNN(structures=rfd,num_sequences=2))
```

### Job submission

Configure computational resources before submitting:

```python
# Examples
pipeline.resources(gpu="A100", memory="32GB", time="24:00:00")      # Specific model
pipeline.resources(gpu="32GB|80GB", memory="32GB", time="24:00:00") # V100 or A100/H100
pipeline.resources(gpu="!T4", memory="32GB", time="24:00:00")       # Exclude T4
pipeline.resources(memory="128GB", time="24:00:00", cpus=32)        # CPU-only
```

**GPU parameter options:**
- `"T4"`, `"V100"`, `"A100"`, `"H100"` - Specific GPU models
- `"32GB"`, `"80GB"`, `"32GB|80GB"` - Memory-based selection
- `"!T4"`, `"!V100"` - Exclude specific models
- `"gpu"` - Any available GPU
- `"high-memory"` - Equivalent to `"32GB|80GB"`
- Omit parameter for CPU-only jobs

Call `pipeline.slurm()` to validate the pipeline, generate the full scripts, and generate a slurm file that can be executed on the cluster. Importantly, this will not result in the submission to slurm. For this you have to either:
1. Run from console the pipeline with python (requires packages pandas and yaml to be executable) and then run the sbatch command in output;
2. Run from console the pipeline with python and then copy-paste the job name and slurm content in the job composer form (https://apps.s3it.uzh.ch/pun/sys/myjobs);
3. Instead of running the pipeline with `python /path/to/<pipeline>.py`, use `./submit /path/to/<pipeline.py>` to both run and submit. If not available, this will install and activate a biopipelines environment with pandas and yaml.

```python
pipeline.slurm()
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

### Filesystem Structure

After the execution, the filesystem will look somewhat like this:

```
/shares/<user>/BioPipelines/<pipeline>/<job>_<NNN>/
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

Most tools require a conda environment to run. Default environments are defined in `config.yaml` at the repository root. To change default environments system-wide, edit the `environments` section in `config.yaml`. To override for a specific tool, use the `env` parameter in `pipeline.add(...)`:

```python
# Use default environment from config.yaml
tool1 = pipeline.add(Tool1(...))
# Override with custom environment
tool2 = pipeline.add(Tool2(...), env="CustomEnv")
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
datasheets        # List[DI]  – Contains objects of class DatasheetInfo
tool_folder       # [str]     – Path to output directory
```

### Datasheet Organization
Datasheet names are tool-dependent. Accessing a datasheet with dot notation will return the path to the csv file. For accessing the metadata, use the info function.

```python
# Simple path
path = <tool output>.datasheets.analysis  # path/to/analysis.csv

# Rich metadata access
info = <tool output>.datasheets.info("analysis")
print(info.path)         # path/to/analysis.csv
print(info.columns)      # ["id", ...]
print(info.description)  # "Results from ..."
print(info.count)        # Number of entries
```

### Datasheet Column References
Tools can explicitly reference specific columns from upstream datasheets using dot notation <tool output>.datasheets.<datasheet name>.<column name>, which returns the tuple (DatasheetInfo,str). 

```python
lmpnn = LigandMPNN(
    structures=rfdaa, #output from RFdiffusionAllAtom
    ligand='LIG',
    redesigned=rfdaa.datasheets.structures.designed     #structures has columns id, fixed, designed, ... 
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

# Complete BioPipelines Tool Reference

This document contains the complete, verified tool reference with all parameters, default environments, and output specifications.

## Structure Generation

### RFdiffusion

Generates novel protein backbone structures using diffusion models. Designs de novo proteins or scaffolds functional motifs into new contexts.

**References**: https://www.nature.com/articles/s41586-023-06415-8

**Installation**: Go to the data directory, then clone the RFdiffusion git repository and download the weights as indicated in https://github.com/RosettaCommons/RFdiffusion. The environment they provide won't work on our cluster, so instead go to the biopipelines directory and run this:
```bash
conda env create -f Environments/ProteinEnv.yaml
conda activate ProteinEnv
pip install -r Environments/ProteinEnv_pip_requirements.txt
```
Followed by the original instructions:
```bash
cd env/SE3Transformer
pip install --no-cache-dir -r requirements.txt
python setup.py install
cd ../.. # change into the root directory of the repository
pip install -e . # install the rfdiffusion module from the root of the repository
```

**Parameters**:
- `pdb`: Union[str, ToolOutput, StandardizedOutput] = "" - Input PDB template (optional)
- `contigs`: str (required) - Contig specification (e.g., "A1-100", "10-20,A6-140")
- `inpaint`: str = "" - Inpainting specification
- `num_designs`: int = 1 - Number of backbone designs to generate
- `active_site`: bool = False - Use active site model for small motifs (<30 residues)
- `steps`: int = 50 - Number of diffusion steps
- `partial_steps`: int = 0 - Partial diffusion steps
- `reproducible`: bool = False - Use deterministic sampling
- `design_startnum`: int = 1 - Starting number for design IDs

**Outputs**:
- `structures`: List of generated PDB files
- `datasheets.structures`:

  | id | source_id | pdb | fixed | designed | contigs | time | status |
  |----|-----------|-----|-------|----------|---------|------|--------|

**Example**:
```python
rfd = pipeline.add(RFdiffusion(
    contigs="50-100",
    num_designs=10
))
```

---

### RFdiffusionAllAtom

Generates protein structures with explicit modeling of ligands and small molecules. Diffusion model that handles all-atom representation including non-protein entities.

**References**: https://www.science.org/doi/10.1126/science.adl2528.

**Installation**: Works in ProteinEnv environment as installed in RFdiffusion.

**Parameters**:
- `ligand`: str (required) - Ligand identifier in PDB (e.g., 'LIG', 'ATP')
- `pdb`: Union[str, ToolOutput, StandardizedOutput] = "" - Input PDB template
- `contigs`: str (required) - Contig specification
- `inpaint`: str = "" - Inpainting specification
- `num_designs`: int = 1 - Number of designs
- `active_site`: bool = False - Use active site model
- `steps`: int = 200 - Diffusion steps (default higher for AllAtom)
- `partial_steps`: int = 0 - Partial diffusion steps
- `reproducible`: bool = False - Deterministic sampling
- `design_startnum`: int = 1 - Starting design number
- `ppi_design`: bool = False - Enable protein-protein interface design
- `ppi_hotspot_residues`: List[str] = None - Hotspot residues for PPI
- `ppi_binder_length`: int = None - Length of PPI binder
- `autogenerate_contigs`: bool = False - Auto-infer fixed contig segments
- `model_only_neighbors`: bool = False - Only remodel residues near ligand
- `num_recycles`: int = 1 - Number of diffusion recycles
- `scaffold_guided`: bool = False - Enforce strict adherence to scaffold
- `align_motif`: bool = True - Pre-align functional motif
- `deterministic`: bool = False - Fixed RNG seeds
- `inpaint_str`: str = None - Secondary structure pattern for inpainting
- `inpaint_seq`: str = None - Sequence pattern for inpainting
- `inpaint_length`: int = None - Target length for inpainted regions
- `guiding_potentials`: str = None - Custom external potentials

**Outputs**:
- `structures`: List of generated PDB files
- `datasheets.structures`:

  | id | source_id | pdb | fixed | designed | contigs | time | status |
  |----|-----------|-----|-------|----------|---------|------|--------|

**Example**:
```python
rfdaa = pipeline.add(RFdiffusionAllAtom(
    pdb=template,
    ligand='LIG',
    contigs='10-20,A6-140',
    num_designs=5
))
```

---

### BoltzGen

Generates protein binders (proteins, peptides, or nanobodies) targeting specified molecules using an end-to-end pipeline that combines diffusion-based backbone generation, inverse folding, structure prediction, and multi-metric filtering.

**References**: https://github.com/HannesStark/boltzgen

**Installation**: Install BoltzGen via pip with Python ≥3.11:
```bash
conda create -n boltzgen python=3.11
conda activate boltzgen
pip install boltzgen
```

Models (~6GB) download automatically to `~/.cache` on first run, or specify custom cache with `cache_dir` parameter.

**Parameters**:
- `design_spec`: Union[str, Dict] (required) - YAML configuration string/dict or path to YAML file defining:
  - Target entities (proteins, ligands from .cif/.pdb files)
  - Binder specification (sequence ranges like "80..140")
  - Binding site constraints (binding_types, not_binding)
  - Secondary structure specifications (helix/sheet/loop)
  - Structural constraints (disulfide bonds, covalent connections)
- `protocol`: str = "protein-anything" - Design protocol:
  - "protein-anything": General protein binder design
  - "peptide-anything": Peptide binder design
  - "protein-small_molecule": Protein-small molecule binder design
  - "nanobody-anything": Nanobody design
- `num_designs`: int = 10000 - Total designs to generate
- `budget`: int = 100 - Final number after diversity filtering
- `design_checkpoints`: Optional[List[str]] = None - Model checkpoint paths
- `step_scale`: Optional[float] = None - Diffusion step scale
- `noise_scale`: Optional[float] = None - Diffusion noise scale
- `diffusion_batch_size`: Optional[int] = None - Samples per trunk run (auto if None)
- `inverse_fold_num_sequences`: int = 1 - Sequences per backbone
- `skip_inverse_folding`: bool = False - Skip sequence redesign
- `alpha`: float = 0.5 - Quality/diversity trade-off (0.0=quality only, 1.0=diversity only)
- `filter_biased`: bool = True - Remove amino acid composition outliers
- `additional_filters`: Optional[List[str]] = None - Hard threshold expressions (e.g., ["design_ALA>0.3"])
- `metrics_override`: Optional[Dict[str, float]] = None - Per-metric ranking weights
- `devices`: Optional[int] = None - Number of GPUs (auto-detect if None)
- `reuse`: bool = False - Resume interrupted runs
- `steps`: Optional[List[str]] = None - Run only specific pipeline steps
- `cache_dir`: Optional[str] = None - Model download location

**Outputs**:
- `structures`: Final designed structure files (in `final_ranked_designs/final_<budget>_designs/`)
- `datasheets.all_metrics`:

  | id | RMSD | hydrogen_bonds | packing_quality | interface_contacts | binding_energy | design_plddt |
  |----|------|----------------|-----------------|-------------------|----------------|--------------|

- `datasheets.final_metrics`:

  | id | RMSD | hydrogen_bonds | packing_quality | interface_contacts | binding_energy | design_plddt | rank |
  |----|------|----------------|-----------------|-------------------|----------------|--------------|------|

- `datasheets.aggregate_metrics`:

  | metric | mean | std | min | max |
  |--------|------|-----|-----|-----|

- `datasheets.per_target_metrics`:

  | target_id | num_designs | avg_rmsd | avg_plddt |
  |-----------|-------------|----------|-----------|

**Important Notes**:
- **Residue indexing**: All residue indices start at 1 and use canonical mmcif `label_asym_id`, not `auth_asym_id`
- **Sequence specification**: Use ranges like "80..140" for random length, "15..20AAAA" for random prefix + fixed tail
- **Output structure**: Pipeline generates intermediate designs, inverse-folded structures, and final ranked designs with comprehensive metrics

**Example**:
```python
# YAML design specification
design_yaml = """
entities:
  - protein:
      id: C
      sequence: 80..140
  - file:
      path: target.cif
      include:
        - chain:
            id: A
"""

boltzgen = pipeline.add(BoltzGen(
    design_spec=design_yaml,
    protocol="protein-anything",
    num_designs=10000,
    budget=100,
    alpha=0.5
))

# Access final designs
best_designs = boltzgen.datasheets.final_metrics
```

**Example with File**:
```python
# Using external YAML file
boltzgen = pipeline.add(BoltzGen(
    design_spec="designs/my_binder.yaml",
    protocol="peptide-anything",
    num_designs=5000,
    budget=50,
    additional_filters=["design_ALA<0.3", "filter_rmsd_design<2.5"]
))
```

---

## Sequence Design

### ProteinMPNN

Designs protein sequences for given backbone structures. Uses graph neural networks to optimize sequences for structure stability while respecting fixed/designed region constraints.

**References**: https://www.science.org/doi/10.1126/science.add2187.

**Installation**: Go to your data folder and clone the official repository (https://github.com/dauparas/ProteinMPNN). The model will then work in the same environment as RFdiffusion.
```bash
git clone https://github.com/dauparas/ProteinMPN
```

**Parameters**:
- `structures`: Union[str, List[str], ToolOutput] (required) - Input structures
- `datasheets`: Optional[List[str]] = None - Input datasheet files
- `num_sequences`: int = 1 - Number of sequences per structure
- `fixed`: str = "" - Fixed positions (PyMOL selection or datasheet reference)
- `redesigned`: str = "" - Redesigned positions (PyMOL selection or datasheet reference)
- `fixed_chain`: str = "A" - Chain to apply fixed positions
- `plddt_threshold`: float = 100.0 - pLDDT threshold for automatic fixing (residues above threshold are fixed)
- `sampling_temp`: float = 0.1 - Sampling temperature
- `model_name`: str = "v_48_020" - ProteinMPNN model variant
- `soluble_model`: bool = True - Use soluble protein model

**Outputs**:
- `sequences`: CSV file with generated sequences
- `datasheets.sequences`:

  | id | source_id | source_pdb | sequence | score | seq_recovery | rmsd |
  |----|-----------|------------|----------|-------|--------------|------|

**Note**: Sample 0 is the original/template sequence, samples 1+ are designs.

**Example**:
```python
pmpnn = pipeline.add(ProteinMPNN(
    structures=rfd,
    num_sequences=10,
    fixed="1-10+50-60",
    redesigned="20-40"
))
```

---

### LigandMPNN

Designs protein sequences optimized for ligand binding. Specialized version of ProteinMPNN that considers protein-ligand interactions during sequence design.

**References**: https://www.nature.com/articles/s41592-025-02626-1.

**Installation**: As from the official repository (https://github.com/dauparas/LigandMPNN), go to your data folder then run:
```bash
git clone https://github.com/dauparas/LigandMPNN.git
cd LigandMPNN
bash get_model_params.sh "./model_params"
conda create -n ligandmpnn_env python=3.11
pip3 install -r requirements.txt
```

**Parameters**:
- `structures`: Union[str, List[str], ToolOutput] (required) - Input structures
- `ligand`: str (required) - Ligand identifier for binding site focus
- `datasheets`: Optional[List[str]] = None - Input datasheet files
- `name`: str = "" - Job name for output files
- `num_sequences`: int = 1 - Number of sequences to generate per structure
- `fixed`: str = "" - Fixed positions (LigandMPNN format "A3 A4 A5" or datasheet reference)
- `redesigned`: str = "" - Designed positions (LigandMPNN format or datasheet reference)
- `design_within`: float = 5.0 - Distance in Angstroms from ligand for post-generation analysis only (does not control design). For actually designing residues within a distance, use [DistanceSelector](#distanceselector) to select positions first.
- `model`: str = "v_32_010" - LigandMPNN model version (v_32_005, v_32_010, v_32_020, v_32_025)
- `batch_size`: int = 1 - Batch size for processing

**Outputs**:
- `sequences`: CSV file with generated sequences
- `datasheets.sequences`:

  | id | sequence | sample | T | seed | overall_confidence | ligand_confidence | seq_rec |
  |----|----------|--------|---|------|-------------------|-------------------|---------|

**Example**:
```python
lmpnn = pipeline.add(LigandMPNN(
    structures=rfdaa,
    ligand="LIG",
    num_sequences=5,
    redesigned=rfdaa.datasheets.structures.designed
))
```

---

### MutationComposer

Generates new protein sequences by composing mutations based on frequency analysis. Creates combinatorial mutants from mutation profiles with different sampling strategies.

**Installation**: Same environment as MutationProfiler.

**Parameters**:
- `frequencies`: Union[List, ToolOutput, StandardizedOutput, DatasheetInfo, str] (required) - Mutation frequency datasheet(s) from MutationProfiler
- `num_sequences`: int = 10 - Number of sequences to generate
- `mode`: str = "single_point" - Generation strategy:
  - "single_point": One mutation per sequence
  - "weighted_random": Random mutations weighted by frequency
  - "hotspot_focused": Focus on high-frequency positions
  - "top_mutations": Use only top N mutations
- `min_frequency`: float = 0.01 - Minimum frequency threshold for mutations
- `max_mutations`: int = None - Maximum mutations per sequence
- `random_seed`: int = None - Random seed for reproducibility
- `prefix`: str = "" - Prefix for sequence IDs
- `hotspot_count`: int = 10 - Number of top hotspot positions (for hotspot_focused mode)
- `combination_strategy`: str = "average" - Strategy for combining multiple datasheets (average, maximum, stack, round_robin)

**Outputs**:
- `sequences`: CSV file with composed sequences
- `datasheets.sequences`:

  | id | sequence | mutations | mutation_positions |
  |----|----------|-----------|-------------------|

**Example**:
```python
profiler = pipeline.add(MutationProfiler(original=ref, mutants=variants))
composer = pipeline.add(MutationComposer(
    frequencies=profiler.datasheets.relative_frequencies,
    num_sequences=50,
    mode="weighted_random",
    max_mutations=5
))
```

---

### SDM (SiteDirectedMutagenesis)

Performs site-directed mutagenesis at specified positions. Generates systematic amino acid substitutions for experimental library design or computational scanning.

**Environment**: `ProteinEnv`

**Parameters**:
- `original`: Union[str, ToolOutput, StandardizedOutput] (required) - Input structure/sequence
- `position`: int (required) - Target position for mutagenesis (1-indexed)
- `mode`: str = "saturation" - Mutagenesis strategy:
  - "saturation": All 20 amino acids
  - "hydrophobic": Hydrophobic residues only
  - "hydrophilic": Hydrophilic residues only
  - "charged": Charged residues only
  - "polar": Polar residues only
  - "nonpolar": Nonpolar residues only
  - "aromatic": Aromatic residues only
  - "aliphatic": Aliphatic residues only
  - "positive": Positively charged residues only
  - "negative": Negatively charged residues only
- `include_original`: bool = False - Include original amino acid in output
- `exclude`: str = "" - Amino acids to exclude (single letter codes as string, e.g., "CP")
- `prefix`: str = "" - Prefix for sequence IDs

**Outputs**:
- `sequences`: CSV file with mutant sequences
- `datasheets.sequences`:

  | id | sequence | mutation | position | original_aa | new_aa |
  |----|----------|----------|----------|-------------|--------|

- `datasheets.missing_sequences`:

  | id | sequence | reason |
  |----|----------|--------|

**Example**:
```python
sdm = pipeline.add(SDM(
    original=template,
    position=42,
    mode="saturation",
    exclude="CP"
))
```

---

### Fuse

Concatenates multiple protein sequences with flexible linkers. Creates fusion proteins with customizable linker lengths for domain engineering.

**Environment**: `ProteinEnv`

**Parameters**:
- `proteins`: Union[List[str], str] (required) - List of protein sequences or PDB file paths
- `sequences`: Union[List[str], str] = None - Alias for proteins
- `name`: str = "" - Job name for output files
- `linker`: str = "GGGGSGGGGSGGGGSGGGGS" - Linker sequence that will be cut based on `linker_lengths` if specified
- `linker_lengths`: List[str] = None - List of length ranges for each junction to generate multiple variants by cutting the linker (e.g., ["1-6", "1-6"])

**Outputs**:
- `sequences`: CSV file with fused sequences
- `datasheets.sequences`:

  | id | sequence | lengths |
  |----|----------|---------|

**Example**:
```python
fused = pipeline.add(Fuse(
    proteins=[domain1, domain2, domain3],
    linker="GGGGS"
))
```

---

### StitchSequences

Combines different regions from multiple sequence sources into chimeric proteins. Useful for creating hybrid sequences from different design outputs.

**Environment**: `ProteinEnv`

**Parameters**:
- `sequences`: List[Union[ToolOutput, StandardizedOutput]] (required) - List of sequence outputs to stitch
- `selections`: Union[List[Union[str, ToolOutput]], str] = None - Position specifications for each sequence (e.g., ["1-50", "51-100"])
- `id_map`: Dict[str, str] = None - ID mapping pattern (default: {"*": "*_N"})

**Outputs**:
- `sequences`: CSV file with stitched sequences
- `datasheets.sequences`:

  | id | sequence |
  |----|----------|

**Example**:
```python
stitched = pipeline.add(StitchSequences(
    sequences=[lmpnn1, lmpnn2],
    selections=["1-50", "51-100"]
))
```

---

## Structure Prediction

### AlphaFold

Predicts protein structures from amino acid sequences using AlphaFold2. Generates high-confidence 3D models with optional relaxation.

**Environment**: `ProteinEnv`

**Parameters**:
- `sequences`: Union[str, List[str], ToolOutput, Dict[str, Any]] (required) - Input sequences or dict with sequences
- `datasheets`: Optional[List[str]] = None - Input datasheet files
- `name`: str = "" - Job name
- `num_relax`: int = 0 - Number of best models to relax with AMBER
- `num_recycle`: int = 3 - Number of recycling iterations
- `rand_seed`: int = 0 - Random seed (0 = random)

**Outputs**:
- `structures`: List of predicted PDB files
- `datasheets.structures`:

  | id | source_id | sequence |
  |----|-----------|----------|

**Example**:
```python
af = pipeline.add(AlphaFold(
    sequences=lmpnn,
    num_relax=1,
    num_recycle=5
))
```

---

### ESMFold

Predicts protein structures using Meta's ESM-2 with ESMFold. Fast single-sequence prediction without requiring MSAs. Models are cached in shared folder for reuse.

**References**: https://github.com/facebookresearch/esm

**Installation**: ESMFold works in the ProteinEnv environment. Install requirements:
```bash
conda activate ProteinEnv
pip install "fair-esm[esmfold]"
pip install 'dllogger @ git+https://github.com/NVIDIA/dllogger.git'
pip install 'openfold @ git+https://github.com/aqlaboratory/openfold.git@4b41059694619831a7db195b7e0988fc4ff3a307'
```

**Environment**: `ProteinEnv`

**Parameters**:
- `sequences`: Union[str, List[str], ToolOutput, StandardizedOutput] (required) - Input sequences
- `datasheets`: Optional[List[str]] = None - Input datasheet files
- `name`: str = "" - Job name
- `num_recycle`: int = 4 - Number of recycling iterations
- `chunk_size`: Optional[int] = None - Chunk size for long sequences (auto if None)

**Outputs**:
- `structures`: List of predicted PDB files
- `datasheets.structures`:

  | id | sequence |
  |----|----------|

**Example**:
```python
esm = pipeline.add(ESMFold(
    sequences=lmpnn,
    num_recycle=4
))
```

---

### Boltz2

Predicts biomolecular complexes including proteins, nucleic acids, and small molecules. State-of-the-art model for protein-ligand and protein-protein complex prediction.

**Environment**: `Boltz2Env`

**Parameters**:
- `config`: Optional[str] = None - Direct YAML configuration string
- `proteins`: Union[str, List[str], ToolOutput] (required) - Protein sequences
- `ligands`: Union[str, ToolOutput, StandardizedOutput, None] = None - Ligand SMILES string or compound library ToolOutput
- `msas`: Optional[Union[str, ToolOutput]] = None - Pre-computed MSA files for recycling (pass entire ToolOutput, not .msas)
- `sequences`: Union[str, List[str], ToolOutput] = None - Legacy parameter (use proteins instead)
- `ligand_smiles`: Optional[str] = None - Legacy parameter (use ligands instead)
- `ligand_library`: Optional[str] = None - Path to CSV file with ligand library (deprecated, use CompoundLibrary tool)
- `primary_key`: Optional[str] = None - Key column in library to filter by
- `library_repr`: str = "SMILES" - Ligand representation (SMILES, CCD)
- `library_type`: str = "noncovalent" - Binding type (noncovalent, covalent)
- `affinity`: bool = True - Calculate binding affinity predictions
- `output_format`: str = "pdb" - Output format (pdb, mmcif)
- `msa_server`: str = "public" - MSA generation (public, local)
- `global_msas_cache`: bool = False - Enable global MSA caching across jobs
- `recycling_steps`: Optional[int] = None - Number of recycling steps (default: model-specific)
- `diffusion_samples`: Optional[int] = None - Number of diffusion samples (default: model-specific)
- `use_potentials`: bool = False - Enable external potentials

**Outputs**:
- `structures`: List of predicted complex PDB files
- `datasheets.confidence`:

  | id | input_file | confidence_score | ptm | iptm | complex_plddt | complex_iplddt |
  |----|------------|------------------|-----|------|---------------|----------------|

- `datasheets.affinity`:

  | id | input_file | affinity_pred_value | affinity_probability_binary |
  |----|------------|---------------------|----------------------------|

**Example**:
```python
boltz_apo = pipeline.add(Boltz2(proteins=lmpnn))
boltz_holo = pipeline.add(Boltz2(
    proteins=lmpnn,
    ligands="CC(=O)OC1=CC=CC=C1C(=O)O",  # Aspirin SMILES
    msas=boltz_apo,  # Pass entire ToolOutput
    affinity=True
), env="Boltz2Env")
```

---

### RF3

⚠️ **UNDER DEVELOPMENT** - This tool is still being tested and refined.

Predicts biomolecular structures using RoseTTAFold3. Supports protein-only and protein-ligand complex prediction with batch processing capabilities.

**Important**: RF3 requires MSAs to be provided. You can obtain MSAs either by:
1. Running Boltz2 first on your proteins (which uses public MMseqs2 server automatically)
2. Using our MMseqs2 implementation to generate MSAs separately

**Environment**: `modelforge`

**Installation**:
The official RF3 repository uses `uv` for installation, but for consistency with BioPipelines we use conda. Run the following in your data folder:
```bash
cd /data/$USER
git clone https://github.com/RosettaCommons/modelforge.git
cd modelforge
conda create -n modelforge python=3.12
conda activate modelforge
pip install -e .
```

**Parameters**:
- `proteins`: Union[str, List[str], ToolOutput, StandardizedOutput] (required) - Protein sequences
- `ligands`: Union[str, ToolOutput, StandardizedOutput, None] = None - Ligand SMILES string or compound library ToolOutput
- `msas`: Optional[Union[str, ToolOutput]] = None - Pre-computed MSA files (optional)
- `output_format`: str = "pdb" - Output format ("pdb" or "cif")
- `checkpoint_path`: Optional[str] = None - Path to RF3 checkpoint file
- `early_stopping_plddt`: Optional[float] = None - pLDDT threshold for early stopping
- `use_templates`: bool = False - Enable template-based prediction

**Note**: The `num_models` parameter is not currently supported due to configuration override limitations in RF3.

**Outputs**:
- `structures`: List of predicted structure files
- `datasheets.structures`:

  | id | model_id | file_path | plddt_score |
  |----|----------|-----------|-------------|

- `datasheets.confidence`:

  | id | model_id | plddt_score | ptm_score |
  |----|----------|-------------|-----------|

**Example**:
```python
# Apo prediction
rf3_apo = pipeline.add(RF3(
    proteins=sequences
))

# Protein-ligand complex prediction
rf3_holo = pipeline.add(RF3(
    proteins=sequences,
    ligands="CC(=O)OC1=CC=CC=C1C(=O)O",  # Aspirin SMILES
    early_stopping_plddt=85.0
))

# Batch prediction with compound library
compounds = pipeline.add(CompoundLibrary({
    "AspA": "CC(=O)OC1=CC=CC=C1C(=O)O",
    "AspB": "CC(=O)OC1=CC=CC=C1C(=O)O"
}))
rf3_batch = pipeline.add(RF3(
    proteins=protein_sequences,
    ligands=compounds
))
```

---

### OnionNet

⚠️ **UNDER DEVELOPMENT** - This tool is still being tested and refined.

Predicts protein-ligand binding affinities from complex structures using OnionNet CNN-based model with rotation-free element-pair-specific contacts.

**Environment**: `OnionNetEnv`

**Installation**:
```bash
cd /data/$USER
git lfs clone https://github.com/zhenglz/onionnet.git
cd onionnet
conda env create -f onet_env.yaml
```

**Parameters**:
- `structures`: Union[str, ToolOutput, StandardizedOutput] (required) - Protein-ligand complex structures
- `model_weights`: Optional[str] = None - Path to model weights file (.h5)
- `scaler_model`: Optional[str] = None - Path to scaler model file
- `output_format`: str = "csv" - Output format ("csv" or "json")

**Outputs**:
- `datasheets.affinities`:

  | id | structure_path | predicted_affinity_pKa |
  |----|----------------|------------------------|

**Example**:
```python
# Predict affinities from Boltz2 structures
affinity = pipeline.add(OnionNet(
    structures=boltz_output,
    model_weights="/path/to/weights.h5",
    scaler_model="/path/to/scaler.model"
))
```

---

### OnionNet2

⚠️ **UNDER DEVELOPMENT** - This tool is still being tested and refined.

Predicts protein-ligand binding affinities using OnionNet-2, an improved version with higher accuracy and lower computational cost. Uses residue-atom contacting shells in CNN architecture.

**Environment**: `OnionNet2Env`

**Installation**:
```bash
cd /data/$USER
git clone https://github.com/zchwang/OnionNet-2.git
cd OnionNet-2
conda create -n OnionNet2Env python=3.8
conda activate OnionNet2Env
pip install tensorflow==2.3 pandas==1.3.4 scikit-learn==0.22.1 numpy==1.18.5 scipy==1.4.1
```

**Parameters**:
- `structures`: Union[str, ToolOutput, StandardizedOutput] (required) - Protein-ligand complex structures
- `model_path`: Optional[str] = None - Path to trained model file
- `scaler_path`: Optional[str] = None - Path to scaler file
- `shells`: int = 62 - Number of contacting shells
- `output_format`: str = "csv" - Output format ("csv" or "json")

**Outputs**:
- `datasheets.affinities`:

  | id | structure_path | predicted_affinity_pKa |
  |----|----------------|------------------------|

**Example**:
```python
# Predict affinities with OnionNet-2
affinity = pipeline.add(OnionNet2(
    structures=boltz_output,
    model_path="/path/to/model.h5",
    scaler_path="/path/to/scaler.pkl",
    shells=62
))
```

---

## Analysis

### DynamicBind

⚠️ **UNDER DEVELOPMENT** - Installation failed.

Predicts ligand-specific protein-ligand complex structures using equivariant diffusion models. Recovers ligand-specific conformations from unbound protein structures.

**Key Features:**
- Predicts protein-ligand binding conformations
- Handles multiple ligands per protein
- Generates multiple poses with affinity predictions
- Outputs relaxed structures with confidence scores

**Installation:**

Requires two conda environments. 

```bash
cd /data/$USER
git clone https://github.com/luwei0917/DynamicBind.git
cd DynamicBind
wget https://zenodo.org/records/10183369/files/workdir.zip
unzip workdir.zip

# Request additional resources during installation e.g.:
srun --mem=32G --cpus-per-task=8 --time=04:00:00 --pty bash

conda create -n dynamicbind python=3.9 -y
conda activate dynamicbind
pip install torch==1.13.1+cu117 torchvision==0.14.1+cu117 torchaudio==0.13.1 \
    --extra-index-url https://download.pytorch.org/whl/cu117

# Setup LD_LIBRARY_PATH inside the environment
mkdir -p $CONDA_PREFIX/etc/conda/activate.d
echo 'export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH' > \
     $CONDA_PREFIX/etc/conda/activate.d/env_vars.sh
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH

# Continue installation
conda install -c conda-forge rdkit biopython scipy pandas pyyaml networkx -y
pip install numpy==1.26.4 e3nn==0.4.4 spyrmsd fair-esm tqdm
pip install \
    torch-scatter==2.1.1+pt113cu117 \
    torch-sparse==0.6.17+pt113cu117 \
    torch-cluster==1.6.1+pt113cu117 \
    torch-spline-conv \
    pyg_lib \
    -f https://data.pyg.org/whl/torch-1.13.1+cu117.html \
    --no-cache-dir
pip install torch-geometric

# Create relax environment
conda create --name relax python=3.8
conda activate relax
conda install -c conda-forge openmm pdbfixer biopython openmmforcefields openff-toolkit ambertools=22 compilers -y
```

Config: `DynamicBind: null` (tool manages environments)

Model weights location: `/data/{username}/DynamicBind/workdir`

**Parameters:**
- `proteins`: Input protein(s) - PDB filename, list of PDBs, or tool output (all structures used)
- `ligands`: SMILES string or tool output with compounds datasheet (must have 'smiles' column)
- `samples_per_complex`: Number of samples per complex (default: 10)
- `inference_steps`: Diffusion steps (default: 20)
- `no_relax`: Skip OpenMM relaxation (default: False)
- `hts`: High-throughput screening mode (default: False)
- `seed`: Random seed (default: 42)

**Example:**

```python
# Basic usage - single protein, SMILES string
db = pipeline.add(DynamicBind(
    proteins="target.pdb",
    ligands="CCO",  # SMILES string
    samples_per_complex=10
))

# Multiple proteins from PDBs folder
db = pipeline.add(DynamicBind(
    proteins=["protein1.pdb", "protein2.pdb"],
    ligands="CC(=O)O"
))

# Chain with other tools - all structures from tool output
rfdaa = pipeline.add(RFdiffusionAllAtom(ligand="ZIT", contigs="A1-100", num_designs=5))
compounds = pipeline.add(CompoundLibrary(smiles_list=["CCO", "CC(=O)O"]))
db = pipeline.add(DynamicBind(proteins=rfdaa, ligands=compounds))
```

**Outputs:**
- Structures: SDF files with pattern `rank{N}_ligand_lddt{score}_affinity{score}_relaxed.sdf`
- Datasheets:
  - `dynamicbind_results.csv`: columns `id`, `ligand_id`, `structure`, `affinity`, `lddt`, `rank`
  - `compounds.csv`: compounds used (available as `output.compounds`)

---

### ResidueAtomDistance

Measures distances between specific atoms and residues in structures. Useful for tracking ligand-protein interactions or structural features.

**Environment**: `ProteinEnv`

**Parameters**:
- `structures`: Union[ToolOutput, StandardizedOutput] (required) - Input structures
- `atom`: str (required) - Atom selection (e.g., 'LIG.Cl', 'name CA', 'A10.CA')
- `residue`: str (required) - Residue selection (e.g., 'D in IGDWG', '145', 'resn ALA')
- `method`: str = "min" - Distance calculation method (min, max, mean, closest)
- `metric_name`: str = None - Custom name for distance column in output (default: "distance")

**Outputs**:
- `datasheets.analysis`:

  | id | source_structure | {metric_name} |
  |----|------------------|---------------|

**Example**:
```python
distances = pipeline.add(ResidueAtomDistance(
    structures=boltz,
    atom="LIG.Cl",
    residue="D in IGDWG",
    method="min",
    metric_name="chlorine_distance"
))
```

---

### PLIP (Protein-Ligand Interaction Profiler)

⚠️ **UNDER DEVELOPMENT** - This tool is still being tested and refined.

Analyzes protein-ligand interactions using PLIP. Identifies hydrogen bonds, hydrophobic contacts, salt bridges, and other interaction types.

**Environment**: `ProteinEnv`

**Note**: This tool is not fully debugged yet and may require adjustments.

**Parameters**:
- `structures`: Union[str, List[str], ToolOutput, StandardizedOutput] (required) - Input structures
- `ligand`: str = "" - Specific ligand identifier (if empty, analyzes all ligands)
- `output_format`: List[str] = None - Output formats (default: ['xml', 'txt', 'pymol'])
- `create_pymol`: bool = True - Generate PyMOL session files
- `create_images`: bool = False - Generate ray-traced images
- `analyze_peptides`: bool = False - Include protein-peptide interactions
- `analyze_intra`: bool = False - Include intra-chain interactions
- `analyze_dna`: bool = False - Include DNA/RNA interactions
- `max_threads`: int = 4 - Maximum threads for parallel processing
- `verbose`: bool = True - Enable verbose output

**Outputs**:
- `datasheets.interactions`:

  | id | ligand_id | interaction_type | residue | distance | angle | energy |
  |----|-----------|------------------|---------|----------|-------|--------|

**Example**:
```python
plip = pipeline.add(PLIP(
    structures=boltz,
    ligand="LIG",
    create_pymol=True
))
```

---

### DistanceSelector

Selects protein residues based on proximity to ligands or other reference points. Generates position specifications for downstream design tools. 

**Installation**: It requires an environment containing pandas (e.g. biopipelines).

**Parameters**:
- `structures`: Union[str, List[str], ToolOutput] (required) - Input structures
- `ligand`: str (required) - Ligand identifier for distance reference
- `distance`: float = 5.0 - Distance cutoff in Angstroms
- `reference_type`: str = "ligand" - Type of reference (ligand, atoms, residues)
- `reference_selection`: str = "" - Specific PyMOL selection if not using ligand

**Outputs**:
- `datasheets.selections`:

  | id | pdb | within | beyond | distance_cutoff | reference_ligand |
  |----|-----|--------|--------|-----------------|------------------|

**Example**:
```python
selector = pipeline.add(DistanceSelector(
    structures=boltz,
    ligand="ATP",
    distance=8.0
))
```

---

### SelectionEditor

Modifies PyMOL-formatted selection strings (e.g., "3-45+58-60") with structure-aware operations. Validates all operations against actual PDB residue numbering and automatically merges overlapping/adjacent ranges.

**Installation**: Requires an environment containing pandas (e.g. biopipelines).

**Parameters**:
- `selection`: tuple (required) - Datasheet column reference (e.g., `tool.datasheets.structures.designed`)
- `structures`: Optional[Union[ToolOutput, List[str]]] = None - Input structures (auto-detected from selection source if not provided)
- `expand`: int = 0 - Number of residues to add on each side of intervals
- `shrink`: int = 0 - Number of residues to remove from each side of intervals
- `shift`: int = 0 - Number of residues to shift all intervals (+/-)
- `invert`: bool = False - Whether to invert the selection (select complement)

**Outputs**:
- `datasheets.selections`:

  | id | pdb | {column_name} | original_{column_name} |
  |----|-----|---------------|------------------------|

**Example**:
```python
# Expand binding site selection by 2 residues
distances = pipeline.add(DistanceSelector(
    structures=rfdaa,
    ligand="LIG",
    distance=5
))

expanded = pipeline.add(SelectionEditor(
    selection=distances.datasheets.selections.within,
    expand=2
))

# Use expanded selection with LigandMPNN
lmpnn = pipeline.add(LigandMPNN(
    structures=rfdaa,
    ligand="LIG",
    redesigned=expanded.datasheets.selections.within
))

# Invert selection to get everything except binding site
fixed_region = pipeline.add(SelectionEditor(
    selection=distances.datasheets.selections.within,
    invert=True
))

# Shrink selection
tighter = pipeline.add(SelectionEditor(
    selection=distances.datasheets.selections.within,
    shrink=1
))
```

**Key Features**:
- **Structure-aware**: Validates against actual PDB residue numbers (handles gaps, non-sequential numbering)
- **Range merging**: Automatically merges overlapping/adjacent intervals (e.g., "1-4+6-10" with expand=1 becomes "1-11")
- **Auto-detection**: Structures parameter is optional; can auto-detect from selection source

---

### ConformationalChange

⚠️ **UNDER DEVELOPMENT** - This tool is still being tested and refined.

Quantifies structural changes between reference and target structures. Calculates RMSD and distance metrics for specified regions.

**Environment**: `ProteinEnv`

**Note**: This tool is not fully debugged yet and may require adjustments.

**Parameters**:
- `reference_structures`: Union[str, ToolOutput, StandardizedOutput] (required) - Reference structures
- `target_structures`: Union[ToolOutput, StandardizedOutput] (required) - Target structures to compare
- `selection`: Union[str, ToolOutput] (required) - Region specification (PyMOL selection or datasheet reference)
- `alignment`: str = "align" - Alignment method (align, super, cealign)

**Outputs**:
- `datasheets.conformational_analysis`:

  | id | reference_structure | target_structure | selection | num_residues | RMSD | max_distance | mean_distance | sum_over_square_root |
  |----|---------------------|------------------|-----------|--------------|------|--------------|---------------|---------------------|

**Example**:
```python
conf_change = pipeline.add(ConformationalChange(
    reference_structures=apo_structures,
    target_structures=holo_structures,
    selection="resi 10-50",  # PyMOL selection
    alignment="super"
))
```

---

### MutationProfiler

Analyzes mutation patterns across sequence sets. Calculates position-specific amino acid frequencies for understanding sequence diversity.

**Installation**: This tool only needs the a small environment:
```bash
conda create -n MutationEnv seaborn matplotlib pandas logomaker scipy
```

**Parameters**:
- `original`: Union[ToolOutput, StandardizedOutput] (required) - Original/reference sequences
- `mutants`: Union[ToolOutput, StandardizedOutput] (required) - Mutant sequences to analyze
- `include_original`: bool = True - Include original sequence in frequency analysis

**Outputs**:
- `datasheets.profile`:

  | position | original | count | frequency |
  |----------|----------|-------|-----------|

- `datasheets.mutations`:

  | position | original | A | C | D | E | F | G | H | I | K | L | M | N | P | Q | R | S | T | V | W | Y |
  |----------|----------|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|

- `datasheets.absolute_frequencies`:

  | position | original | A | C | ... | Y |
  |----------|----------|---|---|-----|---|

- `datasheets.relative_frequencies`:

  | position | original | A | C | ... | Y |
  |----------|----------|---|---|-----|---|

**Example**:
```python
profiler = pipeline.add(MutationProfiler(
    original=template,
    mutants=lmpnn
))
```

---

### ProteinLigandContacts
**Environment**: `ProteinEnv`

Analyzes contacts between selected protein regions and ligands. For each selected residue, calculates the minimum distance to any ligand atom. Returns contact counts and distance statistics.

**Parameters**:
- `structures`: Union[ToolOutput, StandardizedOutput] (required) - Input structures
- `selections`: Union[str, ToolOutput] = None - Protein region selections (string format: '10-20+30-40', datasheet reference, or None for all protein)
- `ligand`: str (required) - Ligand residue name (3-letter code, e.g., 'LIG', 'ATP', 'GDP')
- `contact_threshold`: float = 5.0 - Distance threshold for counting contacts (Angstroms)
- `contact_metric_name`: str = None - Custom name for contact count column (default: "contacts")

**Outputs**:
- `datasheets.contact_analysis`:

  | id | source_structure | selections | ligand | contacts | min_distance | max_distance | mean_distance | sum_distances_sqrt_normalized |
  |----|------------------|------------|--------|----------|--------------|--------------|---------------|-------------------------------|

**Output Columns**:
- `id`: Structure identifier
- `source_structure`: Path to input structure file
- `selections`: Protein residues analyzed (e.g., '10-20+30-40' or 'all_protein')
- `ligand`: Ligand residue name
- `contacts`: Number of residues within contact_threshold distance
- `min_distance`: Minimum distance from any selected residue to ligand (Å)
- `max_distance`: Maximum distance from any selected residue to ligand (Å)
- `mean_distance`: Mean distance from selected residues to ligand (Å)
- `sum_distances_sqrt_normalized`: Sum of distances divided by √(number of residues)

**Example**:
```python
# Analyze contacts with specific protein regions
contacts = pipeline.add(ProteinLigandContacts(
    structures=rfdaa,
    selections=rfdaa.datasheets.structures.designed,
    ligand="LIG",
    contact_threshold=5.0
))

# Use fixed selection for all structures
contacts = pipeline.add(ProteinLigandContacts(
    structures=boltz,
    selections='10-20+30-40',
    ligand="ATP",
    contact_threshold=4.0
))

# Analyze all protein residues
contacts = pipeline.add(ProteinLigandContacts(
    structures=boltz,
    ligand="GDP"
))
```

---

### PoseDistance

Measures ligand pose distance between reference holo structure and sample structures. Calculates RMSD and geometric metrics to quantify how well designed structures reproduce known binding poses.

**Environment**: `ProteinEnv`

**Parameters**:
- `reference_structure`: Union[str, ToolOutput, StandardizedOutput] (required) - Reference holo structure (e.g., XRC structure)
- `sample_structures`: Union[str, List[str], ToolOutput, StandardizedOutput] (required) - Designed/predicted structures to compare
- `reference_ligand`: str (required) - Ligand residue name in reference structure (e.g., 'LIG', 'ATP')
- `sample_ligand`: Optional[str] = None - Ligand residue name in sample structures (default: same as reference_ligand)
- `alignment_selection`: str = "protein" - PyMOL selection for protein alignment (e.g., "chain A", "backbone")
- `calculate_centroid`: bool = True - Calculate ligand centroid distance
- `calculate_orientation`: bool = False - Calculate orientation angle difference

**Outputs**:
- `datasheets.analysis`:

  | id | target_structure | reference_structure | ligand_rmsd | centroid_distance | alignment_rmsd | num_ligand_atoms | alignment_selection |
  |----|------------------|---------------------|-------------|-------------------|----------------|------------------|---------------------|

**Output Columns**:
- `ligand_rmsd`: RMSD between ligand poses after protein alignment (Å)
- `centroid_distance`: Distance between ligand centroids (Å)
- `alignment_rmsd`: RMSD of protein alignment (Å)
- `num_ligand_atoms`: Number of atoms in ligand

**Example**:
```python
# Compare designed structures to XRC reference
xrc = pipeline.add(PDB(pdbs="4ufc", ids="reference"))
designed = pipeline.add(Boltz2(proteins=sequences, ligands="CCO"))

pose_analysis = pipeline.add(PoseDistance(
    reference_structure=xrc,
    sample_structures=designed,
    reference_ligand="ATP",
    sample_ligand="LIG",
    alignment_selection="chain A and backbone"
))

# Filter structures with RMSD < 2.0 Å
good_poses = pipeline.add(Filter(
    data=pose_analysis.datasheets.analysis,
    pool=designed,
    expression="ligand_rmsd < 2.0"
))
```

---

## Data Management

### Filter

Filters structures or sequences based on metric criteria. Uses pandas query expressions to select items meeting specified thresholds.

**Environment**: `ProteinEnv`

**Parameters**:
- `data`: Union[ToolOutput, StandardizedOutput] (required) - Datasheet input to filter
- `pool`: Union[ToolOutput, StandardizedOutput] = None - Structure/sequence pool for copying filtered items
- `expression`: str (required) - Pandas query-style filter expression (e.g., "distance < 3.5 and confidence > 0.8")
- `max_items`: Optional[int] = None - Maximum items to keep after filtering
- `sort_by`: Optional[str] = None - Column name to sort by before applying max_items
- `sort_ascending`: bool = True - Sort order (True = ascending, False = descending)

**Outputs**:
- Filtered pool with same structure as input
- All the datasheets of the upstream tool given as input will be copied and fitered based on the expression, and maintain the same datasheet name. For example, after filtering an output of MergeDatasheets, you can access filtered.datasheets.merged.
- `datasheets.missing`:

  | id | structure | msa |
  |----|-----------|-----|

**Example**:
```python
filtered = pipeline.add(Filter(
    data=distances.datasheets.analysis,
    pool=boltz,
    expression="distance < 3.5 and confidence_score > 0.85",
    max_items=10,
    sort_by="distance"
))
```

---

### SelectBest

Selects the single best structure or sequence based on optimization criteria. Supports single or multi-objective optimization with configurable weights.

**Environment**: `ProteinEnv`

**Parameters**:
- `pool`: Union[ToolOutput, StandardizedOutput, List[Union[ToolOutput, StandardizedOutput]]] (required) - Single or list of tool outputs to select from
- `datasheets`: Union[List[Union[ToolOutput, StandardizedOutput, DatasheetInfo, str]], List[str]] (required) - Datasheets to evaluate for selection
- `metric`: str (required) - Primary metric to optimize
- `mode`: str = "max" - Optimization direction ("max" or "min")
- `weights`: Optional[Dict[str, float]] = None - Dictionary of {metric_name: weight} for multi-metric selection
- `tie_breaker`: str = "first" - How to break ties ("first", "random", or metric name)
- `composite_function`: str = "weighted_sum" - How to combine metrics (weighted_sum, product, min, max)
- `name`: str = "best" - Name for output structure file

**Outputs**:
- Single best structure/sequence with same format as input pool

**Example**:
```python
best = pipeline.add(SelectBest(
    pool=boltz,
    datasheets=[distances.datasheets.analysis],
    metric="distance",
    mode="min"
))

# Multi-objective selection
best_multi = pipeline.add(SelectBest(
    pool=boltz,
    datasheets=[analysis.datasheets.merged],
    metric="composite_score",
    weights={"binding_affinity": 0.6, "pLDDT": 0.4},
    mode="max"
))
```

---

### RemoveDuplicates

Removes duplicate structures or sequences from a pool. Supports deduplication by sequence, structure similarity, or ID matching.

**Environment**: `ProteinEnv`

**Parameters**:
- `pool`: Union[ToolOutput, StandardizedOutput] (required) - Items to deduplicate
- `history`: Optional[Union[ToolOutput, StandardizedOutput, List]] = None - Previous datasheets for cross-cycle deduplication
- `compare`: str = "sequence" - Comparison method (sequence, structure, id)
- `similarity_threshold`: float = 1.0 - Similarity threshold for structure comparison (1.0 = exact match)

**Outputs**:
- Deduplicated pool with same structure as input
- `datasheets.removed`:

  | id | reason |
  |----|--------|

**Example**:
```python
unique = pipeline.add(RemoveDuplicates(
    pool=lmpnn,
    compare="sequence"
))
```

---

### MergeDatasheets

Combines multiple datasheets by joining on a common key column. Enables integration of metrics from different analysis tools.

**Environment**: `ProteinEnv`

**Parameters**:
- `datasheets`: List[Union[ToolOutput, StandardizedOutput, DatasheetInfo, str]] (required) - List of datasheets to merge
- `key`: str = "id" - Join column name
- `prefixes`: Optional[List[str]] = None - Prefixes for columns from each datasheet
- `suffixes`: Optional[List[str]] = None - Suffixes for columns from each datasheet
- `how`: str = "inner" - Join type (inner, outer, left, right)
- `calculate`: Optional[Dict[str, str]] = None - Derived column expressions {new_col: expression}

**Outputs**:
- `datasheets.merged`: Combined datasheet with columns from all inputs

**Example**:
```python
merged = pipeline.add(MergeDatasheets(
    datasheets=[distances.datasheets.analysis, plip.datasheets.interactions],
    prefixes=["dist_", "plip_"],
    key="id",
    calculate={"score": "dist_distance + plip_energy"}
))
```

---

### ConcatenateDatasheets

Stacks multiple datasheets vertically (row-wise). Useful for combining results from multiple cycles or parallel runs.

**Environment**: `ProteinEnv`

**Parameters**:
- `datasheets`: List[Union[ToolOutput, StandardizedOutput, DatasheetInfo, str]] (required) - List of datasheets to concatenate
- `fill`: str = "N/A" - Value for missing columns
- `ignore_index`: bool = True - Reset index in concatenated output

**Outputs**:
- `datasheets.concatenated`: Row-wise concatenation of all input datasheets

**Example**:
```python
concat = pipeline.add(ConcatenateDatasheets(
    datasheets=[cycle1_results, cycle2_results, cycle3_results],
    fill="N/A"
))
```

---

### SliceDatasheet

Extracts a subset of rows and/or columns from a datasheet. Enables data sampling and column selection.

**Environment**: `ProteinEnv`

**Parameters**:
- `datasheet`: Union[ToolOutput, StandardizedOutput, DatasheetInfo, str] (required) - Input datasheet to slice
- `start`: int = 0 - Starting row index
- `end`: Optional[int] = None - Ending row index (None = to end)
- `step`: int = 1 - Step size for slicing
- `columns`: Optional[List[str]] = None - Specific columns to keep (None = all columns)

**Outputs**:
- `datasheets.sliced`: Sliced datasheet

**Example**:
```python
sliced = pipeline.add(SliceDatasheet(
    datasheet=results.datasheets.analysis,
    start=0,
    end=100,
    columns=["id", "distance", "confidence"]
))
```

---

### ExtractMetrics

Extracts and aggregates specific metrics from datasheets. Supports grouping and various aggregation functions for data summarization.

**Environment**: `ProteinEnv`

**Parameters**:
- `datasheets`: List[Union[ToolOutput, StandardizedOutput, DatasheetInfo, str]] (required) - Input datasheets
- `metrics`: List[str] (required) - Metric column names to extract
- `group_by`: Optional[str] = None - Column to group by for aggregation
- `aggregation`: str = "mean" - Aggregation function (mean, median, min, max, sum, std)
- `pivot`: bool = False - Pivot metrics to columns

**Outputs**:
- `datasheets.extracted`: Extracted metrics datasheet

**Example**:
```python
metrics = pipeline.add(ExtractMetrics(
    datasheets=[boltz.datasheets.confidence],
    metrics=["complex_plddt", "ptm"],
    group_by="input_file",
    aggregation="mean"
))
```

---

### AverageByDatasheet

Computes averages of metrics grouped by a specified column. Useful for summarizing results across multiple structures or cycles.

**Environment**: `ProteinEnv`

**Parameters**:
- `datasheets`: List[Union[ToolOutput, StandardizedOutput, DatasheetInfo, str]] (required) - Input datasheets
- `group_by`: str (required) - Column to group by
- `metrics`: List[str] (required) - Metric columns to average
- `weights`: Optional[Dict[str, float]] = None - Weights for each metric

**Outputs**:
- `datasheets.averaged`: Averaged metrics by group

**Example**:
```python
averaged = pipeline.add(AverageByDatasheet(
    datasheets=[cycle1.datasheets.analysis, cycle2.datasheets.analysis],
    group_by="structure_id",
    metrics=["distance", "confidence"]
))
```

---

## Utilities

### LoadOutput

Loads previously saved tool outputs for reuse in new pipelines. Enables incremental development and filtering of existing results at pipeline runtime.

**Environment**: `ProteinEnv`

**Parameters**:
- `output_json`: str (required) - Path to tool output JSON file (in ToolOutputs folder)
- `filter`: Optional[str] = None - Pandas query-style filter expression
- `validate_files`: bool = True - Check file existence when loading

**Outputs**:
- Same structure as original tool that created the output

**Example**:
```python
previous_boltz = pipeline.add(LoadOutput(
    output_json="/path/to/job/ToolOutputs/003_Boltz2.json",
    filter="confidence_score > 0.8"
))
```

**PracticalTips**:
- Instead of copy-pasting ligand SMILES across pipelines, you can create a compound library, and then load the smiles passing an id filter:
```python
# Pipeline 1 to store the library
# imports, pipeline instantiation, ...
pipeline.add(CompoundLibrary({
  "Compound1": "CCNCNNCC(=O)C",
  "Compound2": "CCNCNNCC(=O)C",
  ...
}))
# submit

# Pipeline 2 running calculatios with one of the compounds
# imports, pipeline instantiation, ...
compound1 = pipeline.add(LoadOutput(
    output_json="/path/to/job/ToolOutputs/001_CompoundLibrary.json",
    filter='id == "Compound1"' #quotes are important for proper pandas query here: x is a column name; "x" is a string.
))
boltz = pipeline.add(Boltz2(proteins=HaloTag,
                            ligands=compound1))
#submit

---

### MMseqs2

Generates multiple sequence alignments (MSAs) for protein sequences. Used for improving structure prediction quality by providing evolutionary information.

**Environment**: `ProteinEnv`

**Parameters**:
- `sequences`: Union[str, List[str], ToolOutput, StandardizedOutput] (required) - Input sequences
- `output_format`: str = "csv" - Output format (csv, a3m)
- `timeout`: int = 3600 - Timeout in seconds for server response

**Outputs**:
- `datasheets.msas`:

  | id | sequence_id | sequence | msa_file |
  |----|-------------|----------|----------|

**Example**:
```python
msas = pipeline.add(MMseqs2(
    sequences=lmpnn,
    timeout=7200
))
```

---

### MMseqs2Server

Runs an MMseqs2 server for local MSA generation. Automatically started by MMseqs2 client when needed; manual setup typically not required.

**Environment**: None (doesn't require conda)

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
- `datasheets.compounds`:

  | id | format | smiles | ccd | {branching_keys} |
  |----|--------|--------|-----|------------------|

**Examples**:
```python
# With primary_key for combinatorial expansion
library = pipeline.add(CompoundLibrary(
    library={
        "scaffold": "<linker><fluorophore>",
        "linker": ["CCOCC", "CCOCCOCC"],
        "fluorophore": ["c1ccc(N)cc1"]
    },
    primary_key="scaffold"
))

# Without primary_key - direct SMILES list
library = pipeline.add(CompoundLibrary(
    library={
        "compounds": ["CCO", "CCCO", "CCCCO", "c1ccccc1"]
    }
))
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
- `datasheets.structures`:

  | id | pdb_id | file_path | format | file_size | source |
  |----|--------|-----------|--------|-----------|--------|

- `datasheets.sequences`:

  | id | sequence |
  |----|----------|

- `datasheets.failed`:

  | pdb_id | error_message | source | attempted_path |
  |--------|---------------|--------|----------------|

**Examples**:
```python
# Automatic fallback: checks PDBs/, then downloads
pdb = pipeline.add(PDB(
    pdbs=["4ufc", "1abc"],
    ids=["POI1", "POI2"]
))

# Check custom folder first, then PDBs/, then download
pdb = pipeline.add(PDB(
    pdbs=["my_structure"],
    local_folder="/path/to/my/pdbs"
))

# Download biological assembly
pdb = pipeline.add(PDB(
    pdbs="4ufc",
    ids="MyProtein",
    biological_assembly=True,
    remove_waters=False
))
```

---

### PyMOL

Creates PyMOL visualizations and session files for structural analysis. Supports alignment, coloring schemes, and ray-traced image generation.

**Environment**: `ProteinEnv`

**Note**: This tool is not fully debugged yet and may require adjustments.

**Parameters**:
- `structures`: Union[str, List[str], ToolOutput, StandardizedOutput] (required) - Structures to visualize
- `reference_structure`: Optional[Union[str, ToolOutput]] = None - Reference for alignment
- `color_by`: str = "chain" - Coloring scheme (chain, b_factor, element, ss, spectrum)
- `alignment`: str = "align" - Alignment method (align, super, cealign)
- `ray_trace`: bool = False - Generate ray-traced images
- `image_size`: Tuple[int, int] = (1920, 1080) - Image dimensions if ray tracing

**Outputs**:
- PyMOL session files (.pse)
- Images (if ray_trace=True)

**Example**:
```python
pymol = pipeline.add(PyMOL(
    structures=best,
    reference_structure=template,
    color_by="b_factor",
    alignment="super"
))
```

---
