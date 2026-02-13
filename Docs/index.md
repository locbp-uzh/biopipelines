# BioPipelines

A Python framework for automated computational protein design workflows on SLURM clusters and Jupyter notebooks.

---

## What is BioPipelines?

BioPipelines provides standardized interfaces to connect bioinformatics tools into reproducible workflows. It does not execute computations directly -- instead, it generates bash scripts and predicts output file paths, which are then executed on SLURM clusters or interactively in notebooks.

```python
from biopipelines.pipeline import *
from biopipelines.entities import *
from biopipelines.rfdiffusion import RFdiffusion
from biopipelines.protein_mpnn import ProteinMPNN
from biopipelines.alphafold import AlphaFold

with Pipeline(project="Examples",
              job="RFD-ProteinMPNN-AlphaFold2",
              description="Redesign of N terminus domain of lysozyme"):
    Resources(gpu="A100", time="4:00:00", memory="16GB")
    lysozyme = PDB("168L")
    rfd = RFdiffusion(pdb=lysozyme,
                      contigs='50-70/A81-140',
                      num_designs=3)
    pmpnn = ProteinMPNN(structures=rfd, num_sequences=2)
    af = AlphaFold(proteins=pmpnn)
```

---

## Key Features

**30+ integrated tools** -- Structure generation (RFdiffusion, BoltzGen), sequence design (ProteinMPNN, LigandMPNN), structure prediction (AlphaFold, Boltz2), analysis, and more.

**Two execution modes** -- Submit to SLURM clusters with `biopipelines-submit`, or run interactively in Jupyter/Colab notebooks with on-the-fly execution.

**Combinatorics** -- Cartesian products (`Each`) and grouping (`Bundle`) to systematically explore protein-ligand combinations.

**Data management** -- DataStreams for file tracking, Tables for metrics, and Panda for pandas-style transformations (filter, sort, merge, concat).

**Visualization** -- Automated PyMOL sessions and plots.

---

## Quick Start

### Installation

```bash
git clone https://github.com/locbp-uzh/biopipelines
cd biopipelines
mamba env create -f Environments/biopipelines.yaml
mamba activate biopipelines
pip install -e .
ipython kernel install --user --name biopipelines
```

Some clusters are configured to give low memory to the default bash shell, which might result in failure of the above procedure (std_alloc). You can avoid this by running the following prior to the installation:

```bash
srun --mem=16GB --time=1:00:00 --pty bash
```

Edit `config.yaml` to match your cluster configuration.

Individual models have to be installed separately. We provide a pipeline (ExamplePipelines/install_tools.py) to install all the tools used in the repository at once, but please refer to the respective official documentation in case your particular cluster configuration requires adjustments:

```bash
cd ExamplePipelines
biopipelines-submit install_tools.py
```

### Run an Example

```bash
biopipelines-submit rfd_pmpnn_af2.py
```

---

## Documentation

- **[User Guide](UserManual.md)** -- Core concepts, installation, and usage
- **[Tool Reference](ToolReference.md)** -- Complete reference for all tools
- **[Examples](examples.md)** -- Example pipelines with explanations
- **[Developer Guide](DeveloperManual.md)** -- Architecture and tool development
