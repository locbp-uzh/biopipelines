# BioPipelines

A Python framework for automated computational protein design workflows that generates SLURM-executable bash scripts for high-performance computing clusters.

## Overview

BioPipelines creates modular, reproducible workflows by connecting bioinformatics tools through standardized interfaces. The system **generates bash scripts** for cluster execution rather than running computations directly.

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
- **ProteinMPNN** - Protein sequence design
- **MutationComposer** - Mutation-guided sequence generation
- **Fuse** - Protein-protein fusion design

### Folding/Cofolding
- **AlphaFold** - Structure prediction (ColabFold backend)
- **Boltz** - Advanced structure prediction with noncovalent interactions

### Analysis of MPNN Output
- **MutationProfiler** - Statistical analysis of sequence variations

### Structure Analysis
- **ResidueAtomDistance** - Distance measurements and contacts

### Filtering
- **Filter** - Expression-based result filtering
- **MergeDatasheets** - Combine multiple analysis results
- **ConcatenateDatasheets** - Merge datasets across cycles
- **RemoveDuplicates** - Sequence deduplication

### Visualization
- **PyMOL** - Automated molecular visualization sessions

### Miscellaneous
- **LoadOutput** - Import results from previous pipelines

## Documentation

- **[User Manual](UserManual.md)** - Complete usage guide
- **[Tool Documentation](Docs/)** - Detailed tool-specific guides
- **[Architecture](Docs/Architecture.md)** - System design and concepts
- **[Examples](ExamplePipelines/)** - Pipeline examples

## Key Features

- **SLURM Integration**: Automatic job submission and resource management
- **Modular Design**: Mix and match tools for custom workflows  
- **Standardized Interfaces**: CSV datasheets for seamless tool chaining
- **Environment Management**: Automatic conda environment switching
- **Filesystem Prediction**: Pre-determines output structure before execution
- **Cluster Optimization**: Designed for university HPC environments

## Installation

```bash
git clone <repository-url>
cd biopipelines
# Run notebooks directly from this directory for optimal functionality
jupyter notebook
```

**⚠️ Important**: Always run pipelines from the BioPipelines root directory for proper path resolution.

## Support

- Issues and feature requests: [GitHub Issues](https://github.com/your-repo/issues)
- Documentation: [User Manual](UserManual.md)
- Examples: [ExamplePipelines/](ExamplePipelines/)