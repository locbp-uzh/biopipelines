# BioPipelines

[![Documentation](https://img.shields.io/badge/docs-readthedocs-blue)](https://biopipelines.readthedocs.io/en/latest/)
[![Preprint](https://img.shields.io/badge/preprint-bioRxiv-B31B1B)](https://www.biorxiv.org/content/10.64898/2026.03.11.711024v1)

## Overview

A Python framework for automated computational protein design workflows that can run in Jupyter/Colab notebooks as well as on SLURM-based computing clusters. BioPipelines provides standardized interfaces to connect bioinformatics tools.

## Google Colab Notebooks for example pipelines

| Notebook | Link | Description | Tools |
|----------|------|-------------|-------|
| **Inverse Folding** | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/locbp-uzh/biopipelines/blob/main/ExamplePipelines/notebooks/ubiquitin_results.ipynb) | Inverse folding of ubiquitin, AlphaFold2 refolding, RMSD/pLDDT filter, codon optimisation for *E. coli* | ProteinMPNN · AlphaFold · DNAEncoder |
| **Kinase LID Redesign** | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/locbp-uzh/biopipelines/blob/main/ExamplePipelines/notebooks/kinase_LID_redesign_results.ipynb) | De novo backbone design of the adenylate kinase LID domain, filtered by RMSD on the fixed scaffold | RFdiffusion · ProteinMPNN · AlphaFold · ConformationalChange |
| **FRET Biosensor Design** | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/locbp-uzh/biopipelines/blob/main/ExamplePipelines/notebooks/FRET_results.ipynb) | Linker length optimisation for a Ca²⁺-responsive EBFP–CaM–EYFP FRET sensor | Fuse · Boltz2 · Distance · Panda · Plot |
| **Compound Library Screening** | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/locbp-uzh/biopipelines/blob/main/ExamplePipelines/notebooks/compound_screening_results.ipynb) | SAR screening of a compound library against a protein–DNA complex; affinities plotted by substituent | Boltz2 · CompoundLibrary · Panda · Plot |
| **Iterative Binding Optimization** | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/locbp-uzh/biopipelines/blob/main/ExamplePipelines/notebooks/iterative_binding_optimization_results.ipynb) | 5-cycle directed evolution loop to improve ligand binding affinity via LigandMPNN + Boltz2 | Boltz2 · LigandMPNN · DistanceSelector · MutationProfiler · MutationComposer |
| **Boltz2 Showcase** | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/locbp-uzh/biopipelines/blob/main/ExamplePipelines/notebooks/boltz2_results.ipynb) | All Boltz2 input modes: sequences, PDB, ligands, DNA, glycosylation, covalent linkages, Bundle/Each combinatorics | Boltz2 · CompoundLibrary · Bundle · Each |

## Example: Inverse Folding Pipeline

<p align="center">
  <img src="Docs/images/figure1_ubiquitin_pipeline.png" alt="Ubiquitin inverse folding pipeline: PDB → ProteinMPNN → AlphaFold + DNAEncoder" width="800">
</p>

## Example: Compound Library Screening with Boltz2

<p align="center">
  <img src="Docs/images/figure3_boltz2_compound_screening.png" alt="Boltz2 compound library screening with protein-DNA-ligand co-folding" width="800">
</p>

## Documentation

Full documentation is available at **[biopipelines.readthedocs.io](https://biopipelines.readthedocs.io/en/latest/)**.

- **[User Manual](Docs/UserManual.md)**
- **[Tool Reference](Docs/ToolReference.md)**
- **[Examples](ExamplePipelines/)**
- **[Developer Manual](Docs/DeveloperManual.md)**

## Tool Compatibility

✅ Working · 🟡 Partial / caveats · ⏳ In progress / testing · ❌ Not supported

| Category | Tool | Linux (SLURM / Jupyter) | Google Colab |
|----------|------|:-----------------------:|:------------:|
| **Utilities** | PDB | ✅ | ✅ |
| | RCSB | ✅ | ✅ |
| | Sequence | ✅ | ✅ |
| | Ligand | ✅ | ✅ |
| | CompoundLibrary | ✅ | ✅ |
| | MMseqs2 | ⏳ Partial uniref30; will move to colabfold_search | ❌ Databases too large |
| **Visualization** | PyMOL | ✅ | 🟡 Works, no inline display |
| | Plot | ✅ | ✅ |
| **Structure Generation** | RFdiffusion | ✅ | ✅ |
| | RFdiffusionAllAtom | ✅ | ✅ |
| | RFdiffusion3 | ✅ | ⏳ |
| | BoltzGen | ✅ | 🟡 High compute time |
| **Sequence Design** | ProteinMPNN | ✅ | ✅ |
| | LigandMPNN | ✅ | ✅ |
| | MutationComposer | ✅ | ✅ |
| | Mutagenesis | ✅ | ✅ |
| | Fuse | ✅ | ✅ |
| | StitchSequences | ✅ | ✅ |
| | DNAEncoder | ✅ | ✅ |
| | RBSDesigner | ✅ | ✅ |
| **Structure Prediction** | AlphaFold | ✅ | ✅ |
| | Boltz2 | ✅ | ✅ |
| | Gnina | ✅ | ⏳ |
| **Analysis** | Distance | ✅ | ✅ |
| | Angle | ✅ | ✅ |
| | DistanceSelector | ✅ | ✅ |
| | ConformationalChange | ✅ | ✅ |
| | Contacts | ✅ | ✅ |
| | PoseBusters | ✅ | ⏳ |
| | PoseChange | ✅ | ✅ |
| | CABSflex | ✅ | ⏳ |
| **Statistics** | MutationProfiler | ✅ | ✅ |
| | SequenceMetricCorrelation | ✅ | ⏳ |
| | BayesianAdjuster | ✅ | ⏳ |
| **Data Management** | Panda | ✅ | ✅ |