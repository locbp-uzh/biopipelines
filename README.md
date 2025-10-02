# BioPipelines

A Python framework for automated computational protein design workflows that generates SLURM-executable bash scripts for high-performance computing clusters.

## Overview

BioPipelines creates modular, reproducible workflows by connecting bioinformatics tools through standardized interfaces. The system **generates bash scripts** for cluster execution rather than running computations directly.

## Installation

```bash
git clone https://gitlab.uzh.ch/locbp/public/biopipelines
```

**⚠️ Important**: Always run pipelines from the BioPipelines root directory for proper path resolution.

## Quick Start

```python
from PipelineScripts.pipeline import Pipeline
from PipelineScripts.ligand_mpnn import LigandMPNN
from PipelineScripts.boltz2 import Boltz2

# Create pipeline
pipeline = Pipeline("MyPipeline", "TestJob")
pipeline.resources(gpu="V100", time="24:00:00", memory="16GB")

# Add tools
lmpnn = pipeline.add(LigandMPNN(structures="protein.pdb", num_sequences=10))
boltz = pipeline.add(Boltz2(proteins=lmpnn.output, ligands="SMILES_STRING"))

# Generate and submit to cluster
pipeline.save()
pipeline.slurm(email="user@domain.com")
```

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
- **FetchStructure** - Download structures from databases

## Documentation

- **[User Manual](UserManual.md)** - Complete usage guide with tool reference
- **[Architecture](Docs/Architecture.md)** - System design and concepts
- **[Tool Development](Docs/ToolDevelopment.md)** - Guide for developing new tools
- **[Tool Data](Docs/ToolData.md)** - Data format specifications
- **[Examples](ExamplePipelines/)** - Pipeline examples

## Key Features

- **SLURM Integration**: Automatic job submission and resource management
- **Modular Design**: Mix and match tools for custom workflows  
- **Standardized Interfaces**: CSV datasheets for seamless tool chaining
- **Environment Management**: Automatic conda environment switching
- **Filesystem Prediction**: Pre-determines output structure before execution
- **Cluster Optimization**: Designed for university HPC environments


## Support

- Issues and feature requests: [GitHub Issues](https://github.com/your-repo/issues)
- Documentation: [User Manual](UserManual.md)
- Examples: [ExamplePipelines/](ExamplePipelines/)