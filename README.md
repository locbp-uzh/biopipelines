# BioPipelines

## Overview

A Python framework for automated computational protein design workflows that generates SLURM-executable bash scripts for high-performance computing clusters. BioPipelines creates modular, reproducible workflows by connecting bioinformatics tools through standardized interfaces. 

## Installation

```bash
git clone https://gitlab.uzh.ch/locbp/public/biopipelines
```
Specific environments for using the tools have to be installed separately. Good luck.

## Available Tools

### Structure Generation
- **RFdiffusion** - Backbone structure generation
- **RFdiffusionAllAtom** - All-atom structure generation with ligands

### Sequence Generation
- **LigandMPNN** - Ligand-aware sequence design
- **ProteinMPNN** - Protein sequence design (0-based indexing: sample 0 = template)
- **MutationComposer** - Mutation-guided sequence generation
- **Fuse** - Protein-protein fusion design
- **StitchSequences** - Combine sequences from multiple tools
- **SiteDirectedMutagenesis** - Introduce specific mutations

### Folding/Cofolding
- **AlphaFold** - Structure prediction (ColabFold backend)
- **Boltz** - Advanced structure prediction with noncovalent interactions

### Analysis of MPNN Output
- **MutationProfiler** - Statistical analysis of sequence variations

### Structure Analysis
- **ResidueAtomDistance** - Distance measurements and contacts
- **ProteinLigandContacts** - Protein-ligand contact analysis
- **PoseDistance** - Measure ligand pose RMSD between reference and designed structures
- **ConformationalChange** - Quantify structural changes with multiple metrics (RMSD, distances)
- **DistanceSelector** - Distance-based residue selection
- **PLIP** - Protein-ligand interaction profiling

### Data Management & Filtering
- **Filter** - Expression-based result filtering
- **SelectBest** - Select top-performing structures/sequences
- **MergeDatasheets** - Combine multiple analysis results
- **ConcatenateDatasheets** - Merge datasets across cycles
- **RemoveDuplicates** - Sequence deduplication
- **SliceDatasheet** - Extract subset of rows
- **AverageByDatasheet** - Calculate average metrics
- **ExtractMetrics** - Extract specific metrics from datasheets

### Visualization
- **PyMOL** - Automated molecular visualization sessions

### Utilities
- **LoadOutput** - Import results from previous pipelines
- **MMseqs2** - Multiple sequence alignment and homology search
- **CompoundLibrary** - Generate compound libraries and properties
- **PDB** - Fetch structures from local folders or RCSB PDB

## Documentation

- **[User Manual](UserManual.md)** - Complete usage guide with tool reference
- **[Examples](ExamplePipelines/)** - Pipeline examples
- **[Developer Manual](DeveloperManual.md)** - Complete usage guide with tool reference

## Key Features

- **SLURM Integration**: Automatic job submission and resource management
- **Modular Design**: Mix and match tools for custom workflows  
- **Standardized Interfaces**: CSV datasheets for seamless tool chaining
- **Environment Management**: Automatic conda environment switching
- **Filesystem Prediction**: Pre-determines output structure before execution
- **Cluster Optimization**: Designed for university HPC environments


