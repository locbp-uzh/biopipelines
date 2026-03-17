# BioPipelines

[![Documentation](https://img.shields.io/badge/docs-readthedocs-blue)](https://biopipelines.readthedocs.io/en/latest/)
[![Preprint](https://img.shields.io/badge/preprint-bioRxiv-B31B1B)](https://www.biorxiv.org/content/10.64898/2026.03.11.711024v1)

## Overview

A Python framework for automated computational protein design workflows that can run in Jupyter/Colab notebooks as well as on SLURM-based computing clusters. BioPipelines provides standardized interfaces to connect bioinformatics tools.

## Example: Inverse Folding Pipeline

<p align="center">
  <img src="Docs/images/figure1_ubiquitin_pipeline.png" alt="Ubiquitin inverse folding pipeline: PDB → ProteinMPNN → AlphaFold + DNAEncoder" width="800">
</p>

## Example: Compound Library Screening with Boltz2

<p align="center">
  <img src="Docs/images/figure3_boltz2_compound_screening.png" alt="Boltz2 compound library screening with protein-DNA-ligand co-folding" width="800">
</p>

## Example Pipelines

| Notebook | Description | Tools |
|----------|-------------|-------|
| [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/locbp-uzh/biopipelines/blob/main/ExamplePipelines/notebooks/ubiquitin.ipynb) **Inverse Folding** | Inverse folding of ubiquitin, AlphaFold2 refolding, RMSD/pLDDT filter, codon optimisation for *E. coli* | ProteinMPNN · AlphaFold · DNAEncoder |
| [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/locbp-uzh/biopipelines/blob/main/ExamplePipelines/notebooks/kinase_LID_redesign.ipynb) **Kinase LID Redesign** | De novo backbone design of the adenylate kinase LID domain, filtered by RMSD on the fixed scaffold | RFdiffusion · ProteinMPNN · AlphaFold · ConformationalChange |
| [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/locbp-uzh/biopipelines/blob/main/ExamplePipelines/notebooks/FRET.ipynb) **FRET Biosensor Design** | Linker length optimisation for a Ca²⁺-responsive EBFP–CaM–EYFP FRET sensor | Fuse · Boltz2 · Distance · Panda · Plot · PyMOL |
| [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/locbp-uzh/biopipelines/blob/main/ExamplePipelines/notebooks/compound_screening.ipynb) **Compound Library Screening** | SAR screening of a compound library against a protein–DNA complex; affinities plotted by substituent | Boltz2 · CompoundLibrary · Panda · Plot |
| [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/locbp-uzh/biopipelines/blob/main/ExamplePipelines/notebooks/iterative_binding_optimization.ipynb) **Iterative Binding Optimization** | 5-cycle directed evolution loop to improve ligand binding affinity via LigandMPNN + Boltz2 | Boltz2 · LigandMPNN · DistanceSelector · MutationProfiler · MutationComposer |
| [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/locbp-uzh/biopipelines/blob/main/ExamplePipelines/notebooks/boltz2.ipynb) **Boltz2 Showcase** | All Boltz2 input modes: sequences, PDB, ligands, DNA, glycosylation, covalent linkages, Bundle/Each combinatorics | Boltz2 · CompoundLibrary · Bundle · Each |

## Documentation

Full documentation is available at **[biopipelines.readthedocs.io](https://biopipelines.readthedocs.io/en/latest/)**.

- **[User Manual](Docs/UserManual.md)**
- **[Tool Reference](Docs/ToolReference.md)**
- **[Examples](ExamplePipelines/)**
- **[Developer Manual](Docs/DeveloperManual.md)**
