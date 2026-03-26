# BioPipelines

[![Documentation](https://img.shields.io/badge/docs-readthedocs-blue)](https://biopipelines.readthedocs.io/en/latest/)
[![Preprint](https://img.shields.io/badge/preprint-bioRxiv-B31B1B)](https://www.biorxiv.org/content/10.64898/2026.03.11.711024v1)

## Repository Structure

This repository is hosted on the lab's internal GitLab instance and is **not** the same as the public GitHub mirror.

### Branches

| Branch | Purpose |
|--------|---------|
| `main` | Internal lab branch — shared within the group; may contain work-in-progress tools and pipelines not yet ready for public release |
| `public` | Stable, release-ready code — automatically mirrored to [GitHub](https://github.com/locbp-uzh/biopipelines) |
| `<name>-main` | Personal development branch for a lab member (e.g. `gianluca-main`) |
| `<name>-tool-<toolname>` | Isolated feature branch for a specific tool or component (e.g. `gianluca-tool-esmfold`) |
| `release/<version>` | Short-lived branch for batching cherry-picks before a public release (e.g. `release/1.2.0`) — deleted after the merge request to `public` is accepted |

### Workflow

```
<name>-tool-<feature>          # develop a new tool or larger feature in isolation
        │
        └──► <name>-main       # merge feature branch back when ready
                │
                └──► main      # squash and push when the work is stable enough to share internally
                        │
                        └──► public   # cherry-pick specific commits; open a merge request for review
                                │
                                └──► GitHub (automatic mirror)
```

### Guidelines

- **Personal `<name>-main` branches** are the default place for day-to-day development. Commits can be messy here.
- **Feature branches** (`<name>-tool-<toolname>`) should be used for larger, self-contained additions. Merge back into your personal branch when done.
- **Merging to `main`**: squash your commits into one clean commit before pushing to `main`. Open a merge request, so that it will be properly evaluated for potential transfer to `public`. — `main` is for internal sharing, not gated releases.
- **Promoting to `public`**: cherry-pick only the commits that are stable and ready for external users. Open a merge request so at least one other person reviews before the code lands on GitHub.
  - If you have several features to release at once, create a `release/<version>` branch, cherry-pick onto it, and open a single merge request to `public`. Delete the release branch afterwards.
- **Bug fixes** follow the same path: fix on `<name>-main`, squash to `main`, cherry-pick to `public` with a merge request.
- `main` and `public` are intentionally kept separate — **not everything on `main` belongs in `public`**.

---

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