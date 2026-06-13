# BioPipelines

[![Documentation](https://img.shields.io/badge/docs-readthedocs-blue)](https://biopipelines.readthedocs.io/en/latest/)
[![Preprint](https://img.shields.io/badge/preprint-bioRxiv-B31B1B)](https://www.biorxiv.org/content/10.64898/2026.03.11.711024v1)

## Overview

A Python framework for automated computational protein design workflows that can run in Jupyter/Colab notebooks as well as on SLURM-based computing clusters. BioPipelines provides standardized interfaces to connect bioinformatics tools.

<table>
  <tr>
    <td align="center" width="50%"><b>Inverse Folding Pipeline</b></td>
    <td align="center" width="50%"><b>Compound Library Screening with Boltz2</b></td>
  </tr>
  <tr>
    <td align="center"><img src="docs/images/figure1_ubiquitin_pipeline.png" alt="Ubiquitin inverse folding pipeline: PDB → ProteinMPNN → AlphaFold + DNAEncoder" width="100%"></td>
    <td align="center"><img src="docs/images/figure3_boltz2_compound_screening.png" alt="Boltz2 compound library screening with protein-DNA-ligand co-folding" width="100%"></td>
  </tr>
</table>

## Quick Start (LLM)

**1. Install [Git](https://git-scm.com/downloads).**

**2. Clone the repository.**

```bash
git clone https://github.com/locbp-uzh/biopipelines
cd biopipelines
```

**3. Open an AI coding assistant inside the repo.** Install a coding agent like [Claude Code](https://claude.com/claude-code) or [Codex](https://openai.com/codex/). Start it from the repository root. The framework provides work contracts under the `llm/` folder: `pipelines.md` (writing workflows), `development.md` (extending the framework), `cluster.md` (automated runs on HPCs), `colab.md` (automated runs on Colab).

```bash
claude #or: codex
> Read and follow `llm/pipelines.md`.
```

For framework development please fork the repository. For manual editing we recommend a visual IDE with autocompletion like [VS Code](http://code.visualstudio.com/).

## Google Colab Notebooks for example pipelines

The **Open in Colab** badge launches the clean notebook to run yourself; the **preview** badge opens a pre-executed copy (outputs and plots already rendered).

| Notebook | Link | Description | Tools |
|----------|------|-------------|-------|
| **Inverse Folding** | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/locbp-uzh/biopipelines/blob/main/example_pipelines/notebooks/ubiquitin.ipynb) [![preview](https://img.shields.io/badge/%F0%9F%91%81%EF%B8%8F%20preview-white?style=flat-square&labelColor=white&color=white)](https://colab.research.google.com/github/locbp-uzh/biopipelines/blob/main/example_pipelines/notebooks/ubiquitin_results.ipynb) | Inverse folding of ubiquitin, AlphaFold2 refolding, RMSD/pLDDT filter, codon optimisation for *E. coli* | ProteinMPNN · AlphaFold · ConformationalChange · DNAEncoder |
| **Kinase LID Redesign** | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/locbp-uzh/biopipelines/blob/main/example_pipelines/notebooks/kinase_LID_redesign.ipynb) [![preview](https://img.shields.io/badge/%F0%9F%91%81%EF%B8%8F%20preview-white?style=flat-square&labelColor=white&color=white)](https://colab.research.google.com/github/locbp-uzh/biopipelines/blob/main/example_pipelines/notebooks/kinase_LID_redesign_results.ipynb) | De novo backbone design of the adenylate kinase LID domain, filtered by RMSD on the fixed scaffold | RFdiffusion · ProteinMPNN · AlphaFold · ConformationalChange |
| **FRET Biosensor Design** | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/locbp-uzh/biopipelines/blob/main/example_pipelines/notebooks/FRET.ipynb) [![preview](https://img.shields.io/badge/%F0%9F%91%81%EF%B8%8F%20preview-white?style=flat-square&labelColor=white&color=white)](https://colab.research.google.com/github/locbp-uzh/biopipelines/blob/main/example_pipelines/notebooks/FRET_results.ipynb) | Linker length optimisation for a Ca²⁺-responsive EBFP–CaM–EYFP FRET sensor | Fuse · Boltz2 · Distance · Panda |
| **Compound Library Screening** | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/locbp-uzh/biopipelines/blob/main/example_pipelines/notebooks/compound_screening.ipynb) [![preview](https://img.shields.io/badge/%F0%9F%91%81%EF%B8%8F%20preview-white?style=flat-square&labelColor=white&color=white)](https://colab.research.google.com/github/locbp-uzh/biopipelines/blob/main/example_pipelines/notebooks/compound_screening_results.ipynb) | SAR screening of a compound library against a protein–DNA complex; affinities plotted by substituent | Boltz2 · CompoundLibrary · Panda |
| **Iterative Binding Optimization** | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/locbp-uzh/biopipelines/blob/main/example_pipelines/notebooks/iterative_binding_optimization.ipynb) [![preview](https://img.shields.io/badge/%F0%9F%91%81%EF%B8%8F%20preview-white?style=flat-square&labelColor=white&color=white)](https://colab.research.google.com/github/locbp-uzh/biopipelines/blob/main/example_pipelines/notebooks/iterative_binding_optimization_results.ipynb) | Multi-cycle directed-evolution loop to improve ligand binding affinity via LigandMPNN + Boltz2 | Boltz2 · LigandMPNN · DistanceSelector · MutationProfiler · MutationComposer |
| **Bayesian Directed Evolution** | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/locbp-uzh/biopipelines/blob/main/example_pipelines/notebooks/bayesian_directed_evolution.ipynb) [![preview](https://img.shields.io/badge/%F0%9F%91%81%EF%B8%8F%20preview-white?style=flat-square&labelColor=white&color=white)](https://colab.research.google.com/github/locbp-uzh/biopipelines/blob/main/example_pipelines/notebooks/bayesian_directed_evolution_results.ipynb) | Data-driven directed evolution of an SH3–peptide interface: correlate ΔG with mutations, then Bayesian-reweight the mutation distribution | Boltz2 · Prodigy · Mutagenesis · MutationProfiler · SequenceMetricCorrelation · BayesianAdjuster · MutationComposer |
| **Compound ADMET Triage** | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/locbp-uzh/biopipelines/blob/main/example_pipelines/notebooks/compound_admet_triage.ipynb) [![preview](https://img.shields.io/badge/%F0%9F%91%81%EF%B8%8F%20preview-white?style=flat-square&labelColor=white&color=white)](https://colab.research.google.com/github/locbp-uzh/biopipelines/blob/main/example_pipelines/notebooks/compound_admet_triage_results.ipynb) | Cheap-to-expensive funnel: RDKit + ADMET filter a kinase-inhibitor library, then Boltz2 co-fold and xTB QM interaction energy on the survivors | RDKit · ADMETAI · Boltz2 · DistanceSelector · XTB |
| **Conformational Flexibility** | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/locbp-uzh/biopipelines/blob/main/example_pipelines/notebooks/conformational_flexibility.ipynb) [![preview](https://img.shields.io/badge/%F0%9F%91%81%EF%B8%8F%20preview-white?style=flat-square&labelColor=white&color=white)](https://colab.research.google.com/github/locbp-uzh/biopipelines/blob/main/example_pipelines/notebooks/conformational_flexibility_results.ipynb) | Two flexibility samplers head-to-head (CABSflex Monte Carlo vs BioEmu) compared by per-residue RMSF on T4 lysozyme | OpenMM · CABSflex · BioEmu · EnsembleAnalysis |
| **Interaction Profiling** | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/locbp-uzh/biopipelines/blob/main/example_pipelines/notebooks/interaction_profiling.ipynb) [![preview](https://img.shields.io/badge/%F0%9F%91%81%EF%B8%8F%20preview-white?style=flat-square&labelColor=white&color=white)](https://colab.research.google.com/github/locbp-uzh/biopipelines/blob/main/example_pipelines/notebooks/interaction_profiling_results.ipynb) | Dissect a co-folded protein–ligand binding mode four ways: ProLIF, PLIP, per-residue contacts and buried SASA | Boltz2 · Reduce · ProLIF · PLIP · Contacts · SASA |
| **MSA-Free Folding** | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/locbp-uzh/biopipelines/blob/main/example_pipelines/notebooks/msa_free_folding.ipynb) [![preview](https://img.shields.io/badge/%F0%9F%91%81%EF%B8%8F%20preview-white?style=flat-square&labelColor=white&color=white)](https://colab.research.google.com/github/locbp-uzh/biopipelines/blob/main/example_pipelines/notebooks/msa_free_folding_results.ipynb) | MSA-free fold → inverse-fold → re-fold round trip comparing Frame2Seq vs ProteinMPNN self-consistency | ESMFold · Frame2Seq · ProteinMPNN · ConformationalChange |
| **Pocket Validation** | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/locbp-uzh/biopipelines/blob/main/example_pipelines/notebooks/pocket_validation.ipynb) [![preview](https://img.shields.io/badge/%F0%9F%91%81%EF%B8%8F%20preview-white?style=flat-square&labelColor=white&color=white)](https://colab.research.google.com/github/locbp-uzh/biopipelines/blob/main/example_pipelines/notebooks/pocket_validation_results.ipynb) | Does a co-folded ligand land in a detectable pocket? P2Rank and AF2BIND scored by recall/precision against a consensus ground-truth pocket | Boltz2 · DistanceSelector · Consensus · Selection · P2Rank · AF2BIND |
| **PPI Interface Energetics** | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/locbp-uzh/biopipelines/blob/main/example_pipelines/notebooks/ppi_interface_energetics.ipynb) [![preview](https://img.shields.io/badge/%F0%9F%91%81%EF%B8%8F%20preview-white?style=flat-square&labelColor=white&color=white)](https://colab.research.google.com/github/locbp-uzh/biopipelines/blob/main/example_pipelines/notebooks/ppi_interface_energetics_results.ipynb) | Characterise the barnase–barstar interface: Prodigy ΔG/Kᴅ, APBS per-chain electrostatics, PLIP interface contacts | UniProt · PDB · Prodigy · APBS · PLIP |
| **Urolithin A Target Search** | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/locbp-uzh/biopipelines/blob/main/example_pipelines/notebooks/urolithin_target_search.ipynb) [![preview](https://img.shields.io/badge/%F0%9F%91%81%EF%B8%8F%20preview-white?style=flat-square&labelColor=white&color=white)](https://colab.research.google.com/github/locbp-uzh/biopipelines/blob/main/example_pipelines/notebooks/urolithin_target_search_results.ipynb) | Find candidate targets for a metabolite: RCSB chemical-similarity search, Boltz2 co-fold, rank by binder probability | RCSB · Boltz2 · PoseBusters · Panda |






## Tools

Each tool lists its references, the compute resources it uses, and the platforms it has been tested on. *References:* <img src="https://img.shields.io/badge/-181717?style=flat-square&logo=github&logoColor=white" alt="repo"> repository, <img src="https://img.shields.io/badge/-paper-B31B1B?style=flat-square&logo=readthedocs&logoColor=white" alt="paper"> paper, <img src="https://img.shields.io/badge/BP-1ABC9C?style=flat-square" alt="BP"> native BioPipelines tool. *Resources:* <img src="https://img.shields.io/badge/CPU-4C6EF5?style=flat-square&logo=data%3Aimage%2Fsvg%2Bxml%3Bbase64%2CPHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNCAyNCIgZmlsbD0id2hpdGUiPjxwYXRoIGQ9Ik05IDJ2Mkg3LjVBMS41IDEuNSAwIDAgMCA2IDUuNVY3SDR2MmgydjJINHYyaDJ2Mkg0djJoMnYxLjVBMS41IDEuNSAwIDAgMCA3LjUgMjBIOXYyaDJ2LTJoMnYyaDJ2LTJoMS41YTEuNSAxLjUgMCAwIDAgMS41LTEuNVYxOGgydi0yaC0ydi0yaDJ2LTJoLTJWOWgyVjdoLTJWNS41QTEuNSAxLjUgMCAwIDAgMTYuNSA0SDE1VjJoLTJ2MmgtMlYySDl6bS0xIDZoOHY4SDhWOHoiLz48L3N2Zz4%3D&logoColor=white" alt="CPU"> runs on CPU, <img src="https://img.shields.io/badge/GPU-76B900?style=flat-square&logo=nvidia&logoColor=white" alt="GPU"> GPU-accelarated. *Platforms:* <img src="https://img.shields.io/badge/HPC-2C3E50?style=flat-square&logo=linux&logoColor=white" alt="HPC"> = Linux HPC SLURM cluster / Jupyter, <img src="https://img.shields.io/badge/Colab-F9AB00?style=flat-square&logo=googlecolab&logoColor=white" alt="Colab"> = Google Colab.

<table>
<colgroup>
<col>
<col width="100">
<col width="80">
<col width="115">
</colgroup>
<tr><td colspan="4"><h3>🧬 Structure Generation</h3></td></tr>

<tr>
  <td><sub><b>BoltzGen</b><sup><i>a</i></sup><br>Design protein, peptide, or nanobody binders against a chosen target.</sub></td>
  <td width="100" align="center"><sub><a href="https://github.com/jwohlwend/boltz"><img src="https://img.shields.io/badge/-181717?style=flat-square&logo=github&logoColor=white" alt="repo"></a> <a href="https://www.biorxiv.org/content/10.1101/2025.06.14.659707v1"><img src="https://img.shields.io/badge/-paper-B31B1B?style=flat-square&logo=readthedocs&logoColor=white" alt="paper"></a></sub></td>
  <td width="80" align="center"><sub><img src="https://img.shields.io/badge/GPU-76B900?style=flat-square&logo=nvidia&logoColor=white" alt="GPU"></sub></td>
  <td width="115" align="center"><sub><span style="white-space:nowrap"><img src="https://img.shields.io/badge/HPC-2C3E50?style=flat-square&logo=linux&logoColor=white" alt="HPC ok"></span>&nbsp;<span style="white-space:nowrap"><img src="https://img.shields.io/badge/Colab-F9AB00?style=flat-square&logo=googlecolab&logoColor=white" alt="Colab ok"></span></sub></td>
</tr>
<tr>
  <td><sub><b>HBDesigner</b><br>Design hydrogen-bond networks onto a protein backbone.</sub></td>
  <td width="100" align="center"><sub><a href="https://github.com/Kuhlman-Lab/HBDesigner"><img src="https://img.shields.io/badge/-181717?style=flat-square&logo=github&logoColor=white" alt="repo"></a> <a href="https://rosettacommons.github.io/HBDesigner/"><img src="https://img.shields.io/badge/-docs-B31B1B?style=flat-square&logo=readthedocs&logoColor=white" alt="docs"></a></sub></td>
  <td width="80" align="center"><sub><img src="https://img.shields.io/badge/GPU-76B900?style=flat-square&logo=nvidia&logoColor=white" alt="GPU"> <img src="https://img.shields.io/badge/CPU-4C6EF5?style=flat-square&logo=data%3Aimage%2Fsvg%2Bxml%3Bbase64%2CPHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNCAyNCIgZmlsbD0id2hpdGUiPjxwYXRoIGQ9Ik05IDJ2Mkg3LjVBMS41IDEuNSAwIDAgMCA2IDUuNVY3SDR2MmgydjJINHYyaDJ2Mkg0djJoMnYxLjVBMS41IDEuNSAwIDAgMCA3LjUgMjBIOXYyaDJ2LTJoMnYyaDJ2LTJoMS41YTEuNSAxLjUgMCAwIDAgMS41LTEuNVYxOGgydi0yaC0ydi0yaDJ2LTJoLTJWOWgyVjdoLTJWNS41QTEuNSAxLjUgMCAwIDAgMTYuNSA0SDE1VjJoLTJ2MmgtMlYySDl6bS0xIDZoOHY4SDhWOHoiLz48L3N2Zz4%3D&logoColor=white" alt="CPU"><sup><i>c</i></sup></sub></td>
  <td width="115" align="center"><sub><span style="white-space:nowrap"><img src="https://img.shields.io/badge/HPC-2C3E50?style=flat-square&logo=linux&logoColor=white" alt="HPC ok"></span>&nbsp;<span style="white-space:nowrap"><img src="https://img.shields.io/badge/Colab-F9AB00?style=flat-square&logo=googlecolab&logoColor=white" alt="Colab ok"></span></sub></td>
</tr>
<tr>
  <td><sub><b>PocketGen</b><br>Design a protein binding pocket tailored to a chosen target ligand.</sub></td>
  <td width="100" align="center"><sub><a href="https://github.com/zaixizhang/PocketGen"><img src="https://img.shields.io/badge/-181717?style=flat-square&logo=github&logoColor=white" alt="repo"></a> <a href="https://www.nature.com/articles/s42256-024-00920-9"><img src="https://img.shields.io/badge/-paper-B31B1B?style=flat-square&logo=readthedocs&logoColor=white" alt="paper"></a></sub></td>
  <td width="80" align="center"><sub><img src="https://img.shields.io/badge/GPU-76B900?style=flat-square&logo=nvidia&logoColor=white" alt="GPU"></sub></td>
  <td width="115" align="center"><sub><span style="white-space:nowrap"><img src="https://img.shields.io/badge/Colab-F9AB00?style=flat-square&logo=googlecolab&logoColor=white" alt="Colab ok"></span></sub></td>
</tr>
<tr>
  <td><sub><b>RFdiffusion</b><br>Generate novel protein backbones for scaffolds, motifs, and target binders.</sub></td>
  <td width="100" align="center"><sub><a href="https://github.com/RosettaCommons/RFdiffusion"><img src="https://img.shields.io/badge/-181717?style=flat-square&logo=github&logoColor=white" alt="repo"></a> <a href="https://www.nature.com/articles/s41586-023-06415-8"><img src="https://img.shields.io/badge/-paper-B31B1B?style=flat-square&logo=readthedocs&logoColor=white" alt="paper"></a></sub></td>
  <td width="80" align="center"><sub><img src="https://img.shields.io/badge/GPU-76B900?style=flat-square&logo=nvidia&logoColor=white" alt="GPU"></sub></td>
  <td width="115" align="center"><sub><span style="white-space:nowrap"><img src="https://img.shields.io/badge/HPC-2C3E50?style=flat-square&logo=linux&logoColor=white" alt="HPC ok"></span>&nbsp;<span style="white-space:nowrap"><img src="https://img.shields.io/badge/Colab-F9AB00?style=flat-square&logo=googlecolab&logoColor=white" alt="Colab ok"></span></sub></td>
</tr>
<tr>
  <td><sub><b>RFdiffusion2</b><br>Atom-level diffusion for enzyme active-site scaffolding around catalytic residues and ligands.</sub></td>
  <td width="100" align="center"><sub><a href="https://github.com/RosettaCommons/RFdiffusion2"><img src="https://img.shields.io/badge/-181717?style=flat-square&logo=github&logoColor=white" alt="repo"></a> <a href="https://www.biorxiv.org/content/10.1101/2025.11.20.689494v1"><img src="https://img.shields.io/badge/-paper-B31B1B?style=flat-square&logo=readthedocs&logoColor=white" alt="paper"></a></sub></td>
  <td width="80" align="center"><sub><img src="https://img.shields.io/badge/GPU-76B900?style=flat-square&logo=nvidia&logoColor=white" alt="GPU"></sub></td>
  <td width="115" align="center"><sub><span style="white-space:nowrap"><img src="https://img.shields.io/badge/HPC-2C3E50?style=flat-square&logo=linux&logoColor=white" alt="HPC ok"></span></sub></td>
</tr>
<tr>
  <td><sub><b>RFdiffusion3</b><br>Latest RFdiffusion generation of protein backbones, scaffolds, and binders.</sub></td>
  <td width="100" align="center"><sub><a href="https://github.com/RosettaCommons/foundry"><img src="https://img.shields.io/badge/-181717?style=flat-square&logo=github&logoColor=white" alt="repo"></a> <a href="https://www.biorxiv.org/content/10.1101/2024.11.13.623358v1"><img src="https://img.shields.io/badge/-paper-B31B1B?style=flat-square&logo=readthedocs&logoColor=white" alt="paper"></a></sub></td>
  <td width="80" align="center"><sub><img src="https://img.shields.io/badge/GPU-76B900?style=flat-square&logo=nvidia&logoColor=white" alt="GPU"></sub></td>
  <td width="115" align="center"><sub><span style="white-space:nowrap"><img src="https://img.shields.io/badge/HPC-2C3E50?style=flat-square&logo=linux&logoColor=white" alt="HPC ok"></span>&nbsp;<span style="white-space:nowrap"><img src="https://img.shields.io/badge/Colab-F9AB00?style=flat-square&logo=googlecolab&logoColor=white" alt="Colab ok"></span></sub></td>
</tr>
<tr>
  <td><sub><b>RFdiffusionAllAtom</b><br>Generate protein backbones around small-molecule, metal, or cofactor ligands.</sub></td>
  <td width="100" align="center"><sub><a href="https://github.com/baker-laboratory/rf_diffusion_all_atom"><img src="https://img.shields.io/badge/-181717?style=flat-square&logo=github&logoColor=white" alt="repo"></a> <a href="https://www.science.org/doi/10.1126/science.adl2528"><img src="https://img.shields.io/badge/-paper-B31B1B?style=flat-square&logo=readthedocs&logoColor=white" alt="paper"></a></sub></td>
  <td width="80" align="center"><sub><img src="https://img.shields.io/badge/GPU-76B900?style=flat-square&logo=nvidia&logoColor=white" alt="GPU"></sub></td>
  <td width="115" align="center"><sub><span style="white-space:nowrap"><img src="https://img.shields.io/badge/HPC-2C3E50?style=flat-square&logo=linux&logoColor=white" alt="HPC ok"></span>&nbsp;<span style="white-space:nowrap"><img src="https://img.shields.io/badge/Colab-F9AB00?style=flat-square&logo=googlecolab&logoColor=white" alt="Colab ok"></span></sub></td>
</tr>

<tr><td colspan="4"><h3>✏️ Sequence Design</h3></td></tr>

<tr>
  <td><sub><b>DNAEncoder</b><br>Reverse-translate a protein into codon-optimised DNA for your host.</sub></td>
  <td width="100" align="center"><sub><img src="https://img.shields.io/badge/BP-1ABC9C?style=flat-square" alt="BioPipelines native tool"></sub></td>
  <td width="80" align="center"><sub><img src="https://img.shields.io/badge/CPU-4C6EF5?style=flat-square&logo=data%3Aimage%2Fsvg%2Bxml%3Bbase64%2CPHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNCAyNCIgZmlsbD0id2hpdGUiPjxwYXRoIGQ9Ik05IDJ2Mkg3LjVBMS41IDEuNSAwIDAgMCA2IDUuNVY3SDR2MmgydjJINHYyaDJ2Mkg0djJoMnYxLjVBMS41IDEuNSAwIDAgMCA3LjUgMjBIOXYyaDJ2LTJoMnYyaDJ2LTJoMS41YTEuNSAxLjUgMCAwIDAgMS41LTEuNVYxOGgydi0yaC0ydi0yaDJ2LTJoLTJWOWgyVjdoLTJWNS41QTEuNSAxLjUgMCAwIDAgMTYuNSA0SDE1VjJoLTJ2MmgtMlYySDl6bS0xIDZoOHY4SDhWOHoiLz48L3N2Zz4%3D&logoColor=white" alt="CPU"></sub></td>
  <td width="115" align="center"><sub><span style="white-space:nowrap"><img src="https://img.shields.io/badge/HPC-2C3E50?style=flat-square&logo=linux&logoColor=white" alt="HPC ok"></span>&nbsp;<span style="white-space:nowrap"><img src="https://img.shields.io/badge/Colab-F9AB00?style=flat-square&logo=googlecolab&logoColor=white" alt="Colab ok"></span></sub></td>
</tr>
<tr>
  <td><sub><b>Frame2Seq</b><br>Inverse folding that quickly generates many sequences per backbone.</sub></td>
  <td width="100" align="center"><sub><a href="https://github.com/dakpinaroglu/Frame2seq"><img src="https://img.shields.io/badge/-181717?style=flat-square&logo=github&logoColor=white" alt="repo"></a> <a href="https://arxiv.org/abs/2312.02447"><img src="https://img.shields.io/badge/-paper-B31B1B?style=flat-square&logo=readthedocs&logoColor=white" alt="paper"></a></sub></td>
  <td width="80" align="center"><sub><img src="https://img.shields.io/badge/CPU-4C6EF5?style=flat-square&logo=data%3Aimage%2Fsvg%2Bxml%3Bbase64%2CPHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNCAyNCIgZmlsbD0id2hpdGUiPjxwYXRoIGQ9Ik05IDJ2Mkg3LjVBMS41IDEuNSAwIDAgMCA2IDUuNVY3SDR2MmgydjJINHYyaDJ2Mkg0djJoMnYxLjVBMS41IDEuNSAwIDAgMCA3LjUgMjBIOXYyaDJ2LTJoMnYyaDJ2LTJoMS41YTEuNSAxLjUgMCAwIDAgMS41LTEuNVYxOGgydi0yaC0ydi0yaDJ2LTJoLTJWOWgyVjdoLTJWNS41QTEuNSAxLjUgMCAwIDAgMTYuNSA0SDE1VjJoLTJ2MmgtMlYySDl6bS0xIDZoOHY4SDhWOHoiLz48L3N2Zz4%3D&logoColor=white" alt="CPU"> <img src="https://img.shields.io/badge/GPU-76B900?style=flat-square&logo=nvidia&logoColor=white" alt="GPU"></sub></td>
  <td width="115" align="center"><sub><span style="white-space:nowrap"><img src="https://img.shields.io/badge/Colab-F9AB00?style=flat-square&logo=googlecolab&logoColor=white" alt="Colab ok"></span></sub></td>
</tr>
<tr>
  <td><sub><b>Fuse</b><br>Join multiple sequences into one fusion protein with chosen linkers.</sub></td>
  <td width="100" align="center"><sub><img src="https://img.shields.io/badge/BP-1ABC9C?style=flat-square" alt="BioPipelines native tool"></sub></td>
  <td width="80" align="center"><sub><img src="https://img.shields.io/badge/CPU-4C6EF5?style=flat-square&logo=data%3Aimage%2Fsvg%2Bxml%3Bbase64%2CPHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNCAyNCIgZmlsbD0id2hpdGUiPjxwYXRoIGQ9Ik05IDJ2Mkg3LjVBMS41IDEuNSAwIDAgMCA2IDUuNVY3SDR2MmgydjJINHYyaDJ2Mkg0djJoMnYxLjVBMS41IDEuNSAwIDAgMCA3LjUgMjBIOXYyaDJ2LTJoMnYyaDJ2LTJoMS41YTEuNSAxLjUgMCAwIDAgMS41LTEuNVYxOGgydi0yaC0ydi0yaDJ2LTJoLTJWOWgyVjdoLTJWNS41QTEuNSAxLjUgMCAwIDAgMTYuNSA0SDE1VjJoLTJ2MmgtMlYySDl6bS0xIDZoOHY4SDhWOHoiLz48L3N2Zz4%3D&logoColor=white" alt="CPU"></sub></td>
  <td width="115" align="center"><sub><span style="white-space:nowrap"><img src="https://img.shields.io/badge/HPC-2C3E50?style=flat-square&logo=linux&logoColor=white" alt="HPC ok"></span>&nbsp;<span style="white-space:nowrap"><img src="https://img.shields.io/badge/Colab-F9AB00?style=flat-square&logo=googlecolab&logoColor=white" alt="Colab ok"></span></sub></td>
</tr>
<tr>
  <td><sub><b>LigandMPNN</b><br>Design sequences around bound ligands, nucleic acids, and metal ions.</sub></td>
  <td width="100" align="center"><sub><a href="https://github.com/dauparas/LigandMPNN"><img src="https://img.shields.io/badge/-181717?style=flat-square&logo=github&logoColor=white" alt="repo"></a> <a href="https://www.biorxiv.org/content/10.1101/2023.12.22.573103v1"><img src="https://img.shields.io/badge/-paper-B31B1B?style=flat-square&logo=readthedocs&logoColor=white" alt="paper"></a></sub></td>
  <td width="80" align="center"><sub><img src="https://img.shields.io/badge/GPU-76B900?style=flat-square&logo=nvidia&logoColor=white" alt="GPU"></sub></td>
  <td width="115" align="center"><sub><span style="white-space:nowrap"><img src="https://img.shields.io/badge/HPC-2C3E50?style=flat-square&logo=linux&logoColor=white" alt="HPC ok"></span>&nbsp;<span style="white-space:nowrap"><img src="https://img.shields.io/badge/Colab-F9AB00?style=flat-square&logo=googlecolab&logoColor=white" alt="Colab ok"></span></sub></td>
</tr>
<tr>
  <td><sub><b>Mutagenesis</b><br>Generate systematic single- or multi-point mutant variants of a sequence.</sub></td>
  <td width="100" align="center"><sub><img src="https://img.shields.io/badge/BP-1ABC9C?style=flat-square" alt="BioPipelines native tool"></sub></td>
  <td width="80" align="center"><sub><img src="https://img.shields.io/badge/CPU-4C6EF5?style=flat-square&logo=data%3Aimage%2Fsvg%2Bxml%3Bbase64%2CPHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNCAyNCIgZmlsbD0id2hpdGUiPjxwYXRoIGQ9Ik05IDJ2Mkg3LjVBMS41IDEuNSAwIDAgMCA2IDUuNVY3SDR2MmgydjJINHYyaDJ2Mkg0djJoMnYxLjVBMS41IDEuNSAwIDAgMCA3LjUgMjBIOXYyaDJ2LTJoMnYyaDJ2LTJoMS41YTEuNSAxLjUgMCAwIDAgMS41LTEuNVYxOGgydi0yaC0ydi0yaDJ2LTJoLTJWOWgyVjdoLTJWNS41QTEuNSAxLjUgMCAwIDAgMTYuNSA0SDE1VjJoLTJ2MmgtMlYySDl6bS0xIDZoOHY4SDhWOHoiLz48L3N2Zz4%3D&logoColor=white" alt="CPU"></sub></td>
  <td width="115" align="center"><sub><span style="white-space:nowrap"><img src="https://img.shields.io/badge/HPC-2C3E50?style=flat-square&logo=linux&logoColor=white" alt="HPC ok"></span>&nbsp;<span style="white-space:nowrap"><img src="https://img.shields.io/badge/Colab-F9AB00?style=flat-square&logo=googlecolab&logoColor=white" alt="Colab ok"></span></sub></td>
</tr>
<tr>
  <td><sub><b>MutationComposer</b><br>Combine candidate point mutations into all desired multi-mutant variants.</sub></td>
  <td width="100" align="center"><sub><img src="https://img.shields.io/badge/BP-1ABC9C?style=flat-square" alt="BioPipelines native tool"></sub></td>
  <td width="80" align="center"><sub><img src="https://img.shields.io/badge/CPU-4C6EF5?style=flat-square&logo=data%3Aimage%2Fsvg%2Bxml%3Bbase64%2CPHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNCAyNCIgZmlsbD0id2hpdGUiPjxwYXRoIGQ9Ik05IDJ2Mkg3LjVBMS41IDEuNSAwIDAgMCA2IDUuNVY3SDR2MmgydjJINHYyaDJ2Mkg0djJoMnYxLjVBMS41IDEuNSAwIDAgMCA3LjUgMjBIOXYyaDJ2LTJoMnYyaDJ2LTJoMS41YTEuNSAxLjUgMCAwIDAgMS41LTEuNVYxOGgydi0yaC0ydi0yaDJ2LTJoLTJWOWgyVjdoLTJWNS41QTEuNSAxLjUgMCAwIDAgMTYuNSA0SDE1VjJoLTJ2MmgtMlYySDl6bS0xIDZoOHY4SDhWOHoiLz48L3N2Zz4%3D&logoColor=white" alt="CPU"></sub></td>
  <td width="115" align="center"><sub><span style="white-space:nowrap"><img src="https://img.shields.io/badge/HPC-2C3E50?style=flat-square&logo=linux&logoColor=white" alt="HPC ok"></span>&nbsp;<span style="white-space:nowrap"><img src="https://img.shields.io/badge/Colab-F9AB00?style=flat-square&logo=googlecolab&logoColor=white" alt="Colab ok"></span></sub></td>
</tr>
<tr>
  <td><sub><b>ProteinMPNN</b><br>Design amino-acid sequences predicted to fold to your input backbone.</sub></td>
  <td width="100" align="center"><sub><a href="https://github.com/dauparas/ProteinMPNN"><img src="https://img.shields.io/badge/-181717?style=flat-square&logo=github&logoColor=white" alt="repo"></a> <a href="https://www.science.org/doi/10.1126/science.add2187"><img src="https://img.shields.io/badge/-paper-B31B1B?style=flat-square&logo=readthedocs&logoColor=white" alt="paper"></a></sub></td>
  <td width="80" align="center"><sub><img src="https://img.shields.io/badge/CPU-4C6EF5?style=flat-square&logo=data%3Aimage%2Fsvg%2Bxml%3Bbase64%2CPHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNCAyNCIgZmlsbD0id2hpdGUiPjxwYXRoIGQ9Ik05IDJ2Mkg3LjVBMS41IDEuNSAwIDAgMCA2IDUuNVY3SDR2MmgydjJINHYyaDJ2Mkg0djJoMnYxLjVBMS41IDEuNSAwIDAgMCA3LjUgMjBIOXYyaDJ2LTJoMnYyaDJ2LTJoMS41YTEuNSAxLjUgMCAwIDAgMS41LTEuNVYxOGgydi0yaC0ydi0yaDJ2LTJoLTJWOWgyVjdoLTJWNS41QTEuNSAxLjUgMCAwIDAgMTYuNSA0SDE1VjJoLTJ2MmgtMlYySDl6bS0xIDZoOHY4SDhWOHoiLz48L3N2Zz4%3D&logoColor=white" alt="CPU"> <img src="https://img.shields.io/badge/GPU-76B900?style=flat-square&logo=nvidia&logoColor=white" alt="GPU"></sub></td>
  <td width="115" align="center"><sub><span style="white-space:nowrap"><img src="https://img.shields.io/badge/HPC-2C3E50?style=flat-square&logo=linux&logoColor=white" alt="HPC ok"></span>&nbsp;<span style="white-space:nowrap"><img src="https://img.shields.io/badge/Colab-F9AB00?style=flat-square&logo=googlecolab&logoColor=white" alt="Colab ok"></span></sub></td>
</tr>
<tr>
  <td><sub><b>RBSDesigner</b><sup><i>b</i></sup><br>Design ribosome-binding sites for a chosen bacterial expression level.</sub></td>
  <td width="100" align="center"><sub><a href="https://doi.org/10.1038/nbt.1568"><img src="https://img.shields.io/badge/-paper-B31B1B?style=flat-square&logo=readthedocs&logoColor=white" alt="paper"></a></sub></td>
  <td width="80" align="center"><sub><img src="https://img.shields.io/badge/CPU-4C6EF5?style=flat-square&logo=data%3Aimage%2Fsvg%2Bxml%3Bbase64%2CPHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNCAyNCIgZmlsbD0id2hpdGUiPjxwYXRoIGQ9Ik05IDJ2Mkg3LjVBMS41IDEuNSAwIDAgMCA2IDUuNVY3SDR2MmgydjJINHYyaDJ2Mkg0djJoMnYxLjVBMS41IDEuNSAwIDAgMCA3LjUgMjBIOXYyaDJ2LTJoMnYyaDJ2LTJoMS41YTEuNSAxLjUgMCAwIDAgMS41LTEuNVYxOGgydi0yaC0ydi0yaDJ2LTJoLTJWOWgyVjdoLTJWNS41QTEuNSAxLjUgMCAwIDAgMTYuNSA0SDE1VjJoLTJ2MmgtMlYySDl6bS0xIDZoOHY4SDhWOHoiLz48L3N2Zz4%3D&logoColor=white" alt="CPU"></sub></td>
  <td width="115" align="center"><sub><span style="white-space:nowrap"><img src="https://img.shields.io/badge/HPC-2C3E50?style=flat-square&logo=linux&logoColor=white" alt="HPC ok"></span>&nbsp;<span style="white-space:nowrap"><img src="https://img.shields.io/badge/Colab-F9AB00?style=flat-square&logo=googlecolab&logoColor=white" alt="Colab ok"></span></sub></td>
</tr>
<tr>
  <td><sub><b>StitchSequences</b><br>Concatenate sequence fragments across streams into full-length constructs.</sub></td>
  <td width="100" align="center"><sub><img src="https://img.shields.io/badge/BP-1ABC9C?style=flat-square" alt="BioPipelines native tool"></sub></td>
  <td width="80" align="center"><sub><img src="https://img.shields.io/badge/CPU-4C6EF5?style=flat-square&logo=data%3Aimage%2Fsvg%2Bxml%3Bbase64%2CPHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNCAyNCIgZmlsbD0id2hpdGUiPjxwYXRoIGQ9Ik05IDJ2Mkg3LjVBMS41IDEuNSAwIDAgMCA2IDUuNVY3SDR2MmgydjJINHYyaDJ2Mkg0djJoMnYxLjVBMS41IDEuNSAwIDAgMCA3LjUgMjBIOXYyaDJ2LTJoMnYyaDJ2LTJoMS41YTEuNSAxLjUgMCAwIDAgMS41LTEuNVYxOGgydi0yaC0ydi0yaDJ2LTJoLTJWOWgyVjdoLTJWNS41QTEuNSAxLjUgMCAwIDAgMTYuNSA0SDE1VjJoLTJ2MmgtMlYySDl6bS0xIDZoOHY4SDhWOHoiLz48L3N2Zz4%3D&logoColor=white" alt="CPU"></sub></td>
  <td width="115" align="center"><sub><span style="white-space:nowrap"><img src="https://img.shields.io/badge/HPC-2C3E50?style=flat-square&logo=linux&logoColor=white" alt="HPC ok"></span>&nbsp;<span style="white-space:nowrap"><img src="https://img.shields.io/badge/Colab-F9AB00?style=flat-square&logo=googlecolab&logoColor=white" alt="Colab ok"></span></sub></td>
</tr>

<tr><td colspan="4"><h3>🔮 Structure Prediction & Docking</h3></td></tr>

<tr>
  <td><sub><b>AlphaFold</b><br>Predict monomer or multimer protein structures from sequence.</sub></td>
  <td width="100" align="center"><sub><a href="https://github.com/sokrypton/ColabFold"><img src="https://img.shields.io/badge/-181717?style=flat-square&logo=github&logoColor=white" alt="repo"></a> <a href="https://www.nature.com/articles/s41586-021-03819-2"><img src="https://img.shields.io/badge/-paper-B31B1B?style=flat-square&logo=readthedocs&logoColor=white" alt="paper"></a></sub></td>
  <td width="80" align="center"><sub><img src="https://img.shields.io/badge/GPU-76B900?style=flat-square&logo=nvidia&logoColor=white" alt="GPU"></sub></td>
  <td width="115" align="center"><sub><span style="white-space:nowrap"><img src="https://img.shields.io/badge/HPC-2C3E50?style=flat-square&logo=linux&logoColor=white" alt="HPC ok"></span>&nbsp;<span style="white-space:nowrap"><img src="https://img.shields.io/badge/Colab-F9AB00?style=flat-square&logo=googlecolab&logoColor=white" alt="Colab ok"></span></sub></td>
</tr>
<tr>
  <td><sub><b>Boltz2</b><br>Co-fold proteins, nucleic acids, and ligands, with binding-affinity output.</sub></td>
  <td width="100" align="center"><sub><a href="https://github.com/jwohlwend/boltz"><img src="https://img.shields.io/badge/-181717?style=flat-square&logo=github&logoColor=white" alt="repo"></a> <a href="https://www.biorxiv.org/content/10.1101/2025.06.14.659707v1"><img src="https://img.shields.io/badge/-paper-B31B1B?style=flat-square&logo=readthedocs&logoColor=white" alt="paper"></a></sub></td>
  <td width="80" align="center"><sub><img src="https://img.shields.io/badge/GPU-76B900?style=flat-square&logo=nvidia&logoColor=white" alt="GPU"></sub></td>
  <td width="115" align="center"><sub><span style="white-space:nowrap"><img src="https://img.shields.io/badge/HPC-2C3E50?style=flat-square&logo=linux&logoColor=white" alt="HPC ok"></span>&nbsp;<span style="white-space:nowrap"><img src="https://img.shields.io/badge/Colab-F9AB00?style=flat-square&logo=googlecolab&logoColor=white" alt="Colab ok"></span></sub></td>
</tr>
<tr>
  <td><sub><b>DiffDock</b><br>Blind-dock a ligand to a protein without specifying the pocket.</sub></td>
  <td width="100" align="center"><sub><a href="https://github.com/gcorso/DiffDock"><img src="https://img.shields.io/badge/-181717?style=flat-square&logo=github&logoColor=white" alt="repo"></a> <a href="https://arxiv.org/abs/2210.01776"><img src="https://img.shields.io/badge/-paper-B31B1B?style=flat-square&logo=readthedocs&logoColor=white" alt="paper"></a></sub></td>
  <td width="80" align="center"><sub><img src="https://img.shields.io/badge/GPU-76B900?style=flat-square&logo=nvidia&logoColor=white" alt="GPU"></sub></td>
  <td width="115" align="center"><sub><span style="white-space:nowrap"><img src="https://img.shields.io/badge/Colab-F9AB00?style=flat-square&logo=googlecolab&logoColor=white" alt="Colab ok"></span></sub></td>
</tr>
<tr>
  <td><sub><b>DynamicBind</b><br>Dock a ligand while flexing the apo protein toward its bound state.</sub></td>
  <td width="100" align="center"><sub><a href="https://github.com/luwei0917/DynamicBind"><img src="https://img.shields.io/badge/-181717?style=flat-square&logo=github&logoColor=white" alt="repo"></a> <a href="https://www.nature.com/articles/s41467-024-45461-2"><img src="https://img.shields.io/badge/-paper-B31B1B?style=flat-square&logo=readthedocs&logoColor=white" alt="paper"></a></sub></td>
  <td width="80" align="center"><sub><img src="https://img.shields.io/badge/GPU-76B900?style=flat-square&logo=nvidia&logoColor=white" alt="GPU"></sub></td>
  <td width="115" align="center"><sub><span style="white-space:nowrap"><img src="https://img.shields.io/badge/Colab-F9AB00?style=flat-square&logo=googlecolab&logoColor=white" alt="Colab ok"></span></sub></td>
</tr>
<tr>
  <td><sub><b>ESMFold</b><br>Predict a protein structure from a single sequence, no MSA needed.</sub></td>
  <td width="100" align="center"><sub><a href="https://github.com/facebookresearch/esm"><img src="https://img.shields.io/badge/-181717?style=flat-square&logo=github&logoColor=white" alt="repo"></a> <a href="https://www.science.org/doi/10.1126/science.ade2574"><img src="https://img.shields.io/badge/-paper-B31B1B?style=flat-square&logo=readthedocs&logoColor=white" alt="paper"></a></sub></td>
  <td width="80" align="center"><sub><img src="https://img.shields.io/badge/GPU-76B900?style=flat-square&logo=nvidia&logoColor=white" alt="GPU"></sub></td>
  <td width="115" align="center"><sub><span style="white-space:nowrap"><img src="https://img.shields.io/badge/HPC-2C3E50?style=flat-square&logo=linux&logoColor=white" alt="HPC ok"></span>&nbsp;<span style="white-space:nowrap"><img src="https://img.shields.io/badge/Colab-F9AB00?style=flat-square&logo=googlecolab&logoColor=white" alt="Colab ok"></span></sub></td>
</tr>
<tr>
  <td><sub><b>Gnina</b><br>Dock ligands into a pocket and score poses with deep learning.</sub></td>
  <td width="100" align="center"><sub><a href="https://github.com/gnina/gnina"><img src="https://img.shields.io/badge/-181717?style=flat-square&logo=github&logoColor=white" alt="repo"></a> <a href="https://doi.org/10.1186/s13321-021-00522-2"><img src="https://img.shields.io/badge/-paper-B31B1B?style=flat-square&logo=readthedocs&logoColor=white" alt="paper"></a></sub></td>
  <td width="80" align="center"><sub><img src="https://img.shields.io/badge/CPU-4C6EF5?style=flat-square&logo=data%3Aimage%2Fsvg%2Bxml%3Bbase64%2CPHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNCAyNCIgZmlsbD0id2hpdGUiPjxwYXRoIGQ9Ik05IDJ2Mkg3LjVBMS41IDEuNSAwIDAgMCA2IDUuNVY3SDR2MmgydjJINHYyaDJ2Mkg0djJoMnYxLjVBMS41IDEuNSAwIDAgMCA3LjUgMjBIOXYyaDJ2LTJoMnYyaDJ2LTJoMS41YTEuNSAxLjUgMCAwIDAgMS41LTEuNVYxOGgydi0yaC0ydi0yaDJ2LTJoLTJWOWgyVjdoLTJWNS41QTEuNSAxLjUgMCAwIDAgMTYuNSA0SDE1VjJoLTJ2MmgtMlYySDl6bS0xIDZoOHY4SDhWOHoiLz48L3N2Zz4%3D&logoColor=white" alt="CPU"> <img src="https://img.shields.io/badge/GPU-76B900?style=flat-square&logo=nvidia&logoColor=white" alt="GPU"></sub></td>
  <td width="115" align="center"><sub><span style="white-space:nowrap"><img src="https://img.shields.io/badge/HPC-2C3E50?style=flat-square&logo=linux&logoColor=white" alt="HPC ok"></span>&nbsp;<span style="white-space:nowrap"><img src="https://img.shields.io/badge/Colab-F9AB00?style=flat-square&logo=googlecolab&logoColor=white" alt="Colab ok"></span></sub></td>
</tr>
<tr>
  <td><sub><b>NeuralPLexer</b></span><sup><i>n</i></sup><br>Predict protein-ligand complex structures from sequence and ligand graph.</sub></td>
  <td width="100" align="center"><sub><a href="https://github.com/zrqiao/NeuralPLexer"><img src="https://img.shields.io/badge/-181717?style=flat-square&logo=github&logoColor=white" alt="repo"></a> <a href="https://www.nature.com/articles/s42256-024-00792-z"><img src="https://img.shields.io/badge/-paper-B31B1B?style=flat-square&logo=readthedocs&logoColor=white" alt="paper"></a></sub></td>
  <td width="80" align="center"><sub><img src="https://img.shields.io/badge/GPU-76B900?style=flat-square&logo=nvidia&logoColor=white" alt="GPU"></sub></td>
  <td width="115" align="center"><sub><span style="white-space:nowrap"><img src="https://img.shields.io/badge/Colab-F9AB00?style=flat-square&logo=googlecolab&logoColor=white" alt="Colab ok"></sub></td>
</tr>
<tr>
  <td><sub><b>PLACER</b><br>Generate scored ligand-pose or sidechain ensembles in a protein pocket.</sub></td>
  <td width="100" align="center"><sub><a href="https://github.com/baker-laboratory/PLACER"><img src="https://img.shields.io/badge/-181717?style=flat-square&logo=github&logoColor=white" alt="repo"></a> <a href="https://www.pnas.org/doi/10.1073/pnas.2427161122"><img src="https://img.shields.io/badge/-paper-B31B1B?style=flat-square&logo=readthedocs&logoColor=white" alt="paper"></a></sub></td>
  <td width="80" align="center"><sub><img src="https://img.shields.io/badge/GPU-76B900?style=flat-square&logo=nvidia&logoColor=white" alt="GPU"></sub></td>
  <td width="115" align="center"><sub><span style="white-space:nowrap"><img src="https://img.shields.io/badge/HPC-2C3E50?style=flat-square&logo=linux&logoColor=white" alt="HPC ok"></span>&nbsp;<span style="white-space:nowrap"><img src="https://img.shields.io/badge/Colab-F9AB00?style=flat-square&logo=googlecolab&logoColor=white" alt="Colab ok"></span></sub></td>
</tr>

<tr><td colspan="4"><h3>📐 Analysis</h3></td></tr>

<tr>
  <td><sub><b>ADMETAI</b><br>Predict ADMET and physicochemical properties for your compounds.</sub></td>
  <td width="100" align="center"><sub><a href="https://github.com/swansonk14/admet_ai"><img src="https://img.shields.io/badge/-181717?style=flat-square&logo=github&logoColor=white" alt="repo"></a> <a href="https://doi.org/10.1093/bioinformatics/btae416"><img src="https://img.shields.io/badge/-paper-B31B1B?style=flat-square&logo=readthedocs&logoColor=white" alt="paper"></a></sub></td>
  <td width="80" align="center"><sub><img src="https://img.shields.io/badge/CPU-4C6EF5?style=flat-square&logo=data%3Aimage%2Fsvg%2Bxml%3Bbase64%2CPHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNCAyNCIgZmlsbD0id2hpdGUiPjxwYXRoIGQ9Ik05IDJ2Mkg3LjVBMS41IDEuNSAwIDAgMCA2IDUuNVY3SDR2MmgydjJINHYyaDJ2Mkg0djJoMnYxLjVBMS41IDEuNSAwIDAgMCA3LjUgMjBIOXYyaDJ2LTJoMnYyaDJ2LTJoMS41YTEuNSAxLjUgMCAwIDAgMS41LTEuNVYxOGgydi0yaC0ydi0yaDJ2LTJoLTJWOWgyVjdoLTJWNS41QTEuNSAxLjUgMCAwIDAgMTYuNSA0SDE1VjJoLTJ2MmgtMlYySDl6bS0xIDZoOHY4SDhWOHoiLz48L3N2Zz4%3D&logoColor=white" alt="CPU"></sub></td>
  <td width="115" align="center"><sub><span style="white-space:nowrap"><img src="https://img.shields.io/badge/HPC-2C3E50?style=flat-square&logo=linux&logoColor=white" alt="HPC ok"></span>&nbsp;<span style="white-space:nowrap"><img src="https://img.shields.io/badge/Colab-F9AB00?style=flat-square&logo=googlecolab&logoColor=white" alt="Colab ok"></span></sub></td>
</tr>
<tr>
  <td><sub><b>AF2BIND</b><sup><i>q</i></sup><br>Predict small-molecule binding-site residues on a protein structure.</sub></td>
  <td width="100" align="center"><sub><a href="https://github.com/sokrypton/af2bind"><img src="https://img.shields.io/badge/-181717?style=flat-square&logo=github&logoColor=white" alt="repo"></a> <a href="https://www.biorxiv.org/content/10.1101/2023.10.15.562410v1"><img src="https://img.shields.io/badge/-paper-B31B1B?style=flat-square&logo=readthedocs&logoColor=white" alt="paper"></a></sub></td>
  <td width="80" align="center"><sub><img src="https://img.shields.io/badge/GPU-76B900?style=flat-square&logo=nvidia&logoColor=white" alt="GPU"></sub></td>
  <td width="115" align="center"><sub><span style="white-space:nowrap"><img src="https://img.shields.io/badge/HPC-2C3E50?style=flat-square&logo=linux&logoColor=white" alt="HPC ok"></span>&nbsp;<span style="white-space:nowrap"><img src="https://img.shields.io/badge/Colab-F9AB00?style=flat-square&logo=googlecolab&logoColor=white" alt="Colab ok"></span></sub></td>
</tr>
<tr>
  <td><sub><b>Aggrescan3D</b><br>Score per-residue and global protein aggregation propensity from structure.</sub></td>
  <td width="100" align="center"><sub><a href="https://bitbucket.org/lcbio/aggrescan3d"><img src="https://img.shields.io/badge/-181717?style=flat-square&logo=github&logoColor=white" alt="repo"></a> <a href="https://academic.oup.com/nar/article/47/W1/W300/5485072"><img src="https://img.shields.io/badge/-paper-B31B1B?style=flat-square&logo=readthedocs&logoColor=white" alt="paper"></a></sub></td>
  <td width="80" align="center"><sub><img src="https://img.shields.io/badge/CPU-4C6EF5?style=flat-square&logo=data%3Aimage%2Fsvg%2Bxml%3Bbase64%2CPHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNCAyNCIgZmlsbD0id2hpdGUiPjxwYXRoIGQ9Ik05IDJ2Mkg3LjVBMS41IDEuNSAwIDAgMCA2IDUuNVY3SDR2MmgydjJINHYyaDJ2Mkg0djJoMnYxLjVBMS41IDEuNSAwIDAgMCA3LjUgMjBIOXYyaDJ2LTJoMnYyaDJ2LTJoMS41YTEuNSAxLjUgMCAwIDAgMS41LTEuNVYxOGgydi0yaC0ydi0yaDJ2LTJoLTJWOWgyVjdoLTJWNS41QTEuNSAxLjUgMCAwIDAgMTYuNSA0SDE1VjJoLTJ2MmgtMlYySDl6bS0xIDZoOHY4SDhWOHoiLz48L3N2Zz4%3D&logoColor=white" alt="CPU"></sub></td>
  <td width="115" align="center"><sub><span style="white-space:nowrap"><img src="https://img.shields.io/badge/HPC-2C3E50?style=flat-square&logo=linux&logoColor=white" alt="HPC ok"></span>&nbsp;<span style="white-space:nowrap"><img src="https://img.shields.io/badge/Colab-F9AB00?style=flat-square&logo=googlecolab&logoColor=white" alt="Colab ok"></span></sub></td>
</tr>
<tr>
  <td><sub><b>Angle</b><br>Measure bond, torsion, and vector angles within structures.</sub></td>
  <td width="100" align="center"><sub><img src="https://img.shields.io/badge/BP-1ABC9C?style=flat-square" alt="BioPipelines native tool"></sub></td>
  <td width="80" align="center"><sub><img src="https://img.shields.io/badge/CPU-4C6EF5?style=flat-square&logo=data%3Aimage%2Fsvg%2Bxml%3Bbase64%2CPHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNCAyNCIgZmlsbD0id2hpdGUiPjxwYXRoIGQ9Ik05IDJ2Mkg3LjVBMS41IDEuNSAwIDAgMCA2IDUuNVY3SDR2MmgydjJINHYyaDJ2Mkg0djJoMnYxLjVBMS41IDEuNSAwIDAgMCA3LjUgMjBIOXYyaDJ2LTJoMnYyaDJ2LTJoMS41YTEuNSAxLjUgMCAwIDAgMS41LTEuNVYxOGgydi0yaC0ydi0yaDJ2LTJoLTJWOWgyVjdoLTJWNS41QTEuNSAxLjUgMCAwIDAgMTYuNSA0SDE1VjJoLTJ2MmgtMlYySDl6bS0xIDZoOHY4SDhWOHoiLz48L3N2Zz4%3D&logoColor=white" alt="CPU"></sub></td>
  <td width="115" align="center"><sub><span style="white-space:nowrap"><img src="https://img.shields.io/badge/HPC-2C3E50?style=flat-square&logo=linux&logoColor=white" alt="HPC ok"></span>&nbsp;<span style="white-space:nowrap"><img src="https://img.shields.io/badge/Colab-F9AB00?style=flat-square&logo=googlecolab&logoColor=white" alt="Colab ok"></span></sub></td>
</tr>
<tr>
  <td><sub><b>APBS</b><br>Compute electrostatic potential and surface charge of a structure.</sub></td>
  <td width="100" align="center"><sub><a href="https://github.com/Electrostatics/apbs"><img src="https://img.shields.io/badge/-181717?style=flat-square&logo=github&logoColor=white" alt="repo"></a> <a href="https://doi.org/10.1002/pro.3280"><img src="https://img.shields.io/badge/-paper-B31B1B?style=flat-square&logo=readthedocs&logoColor=white" alt="paper"></a></sub></td>
  <td width="80" align="center"><sub><img src="https://img.shields.io/badge/CPU-4C6EF5?style=flat-square&logo=data%3Aimage%2Fsvg%2Bxml%3Bbase64%2CPHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNCAyNCIgZmlsbD0id2hpdGUiPjxwYXRoIGQ9Ik05IDJ2Mkg3LjVBMS41IDEuNSAwIDAgMCA2IDUuNVY3SDR2MmgydjJINHYyaDJ2Mkg0djJoMnYxLjVBMS41IDEuNSAwIDAgMCA3LjUgMjBIOXYyaDJ2LTJoMnYyaDJ2LTJoMS41YTEuNSAxLjUgMCAwIDAgMS41LTEuNVYxOGgydi0yaC0ydi0yaDJ2LTJoLTJWOWgyVjdoLTJWNS41QTEuNSAxLjUgMCAwIDAgMTYuNSA0SDE1VjJoLTJ2MmgtMlYySDl6bS0xIDZoOHY4SDhWOHoiLz48L3N2Zz4%3D&logoColor=white" alt="CPU"></sub></td>
  <td width="115" align="center"><sub><span style="white-space:nowrap"><img src="https://img.shields.io/badge/HPC-2C3E50?style=flat-square&logo=linux&logoColor=white" alt="HPC ok"></span>&nbsp;<span style="white-space:nowrap"><img src="https://img.shields.io/badge/Colab-F9AB00?style=flat-square&logo=googlecolab&logoColor=white" alt="Colab ok"></span></sub></td>
</tr>
<tr>
  <td><sub><b>BioEmu</b><br>Sample a protein's equilibrium conformational ensemble from its sequence.</sub></td>
  <td width="100" align="center"><sub><a href="https://github.com/microsoft/bioemu"><img src="https://img.shields.io/badge/-181717?style=flat-square&logo=github&logoColor=white" alt="repo"></a> <a href="https://www.science.org/doi/10.1126/science.adv9817"><img src="https://img.shields.io/badge/-paper-B31B1B?style=flat-square&logo=readthedocs&logoColor=white" alt="paper"></a></sub></td>
  <td width="80" align="center"><sub><img src="https://img.shields.io/badge/GPU-76B900?style=flat-square&logo=nvidia&logoColor=white" alt="GPU"></sub></td>
  <td width="115" align="center"><sub><span style="white-space:nowrap"><img src="https://img.shields.io/badge/HPC-2C3E50?style=flat-square&logo=linux&logoColor=white" alt="HPC ok"></span>&nbsp;<span style="white-space:nowrap"><img src="https://img.shields.io/badge/Colab-F9AB00?style=flat-square&logo=googlecolab&logoColor=white" alt="Colab ok"></span></sub></td>
</tr>
<tr>
  <td><sub><b>CABSflex</b><br>Estimate protein backbone flexibility and per-residue fluctuation profiles.</sub></td>
  <td width="100" align="center"><sub><a href="https://github.com/lcbio/CABSflex"><img src="https://img.shields.io/badge/-181717?style=flat-square&logo=github&logoColor=white" alt="repo"></a> <a href="https://doi.org/10.1093/nar/gky356"><img src="https://img.shields.io/badge/-paper-B31B1B?style=flat-square&logo=readthedocs&logoColor=white" alt="paper"></a></sub></td>
  <td width="80" align="center"><sub><img src="https://img.shields.io/badge/CPU-4C6EF5?style=flat-square&logo=data%3Aimage%2Fsvg%2Bxml%3Bbase64%2CPHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNCAyNCIgZmlsbD0id2hpdGUiPjxwYXRoIGQ9Ik05IDJ2Mkg3LjVBMS41IDEuNSAwIDAgMCA2IDUuNVY3SDR2MmgydjJINHYyaDJ2Mkg0djJoMnYxLjVBMS41IDEuNSAwIDAgMCA3LjUgMjBIOXYyaDJ2LTJoMnYyaDJ2LTJoMS41YTEuNSAxLjUgMCAwIDAgMS41LTEuNVYxOGgydi0yaC0ydi0yaDJ2LTJoLTJWOWgyVjdoLTJWNS41QTEuNSAxLjUgMCAwIDAgMTYuNSA0SDE1VjJoLTJ2MmgtMlYySDl6bS0xIDZoOHY4SDhWOHoiLz48L3N2Zz4%3D&logoColor=white" alt="CPU"></sub></td>
  <td width="115" align="center"><sub><span style="white-space:nowrap"><img src="https://img.shields.io/badge/HPC-2C3E50?style=flat-square&logo=linux&logoColor=white" alt="HPC ok"></span>&nbsp;<span style="white-space:nowrap"><img src="https://img.shields.io/badge/Colab-F9AB00?style=flat-square&logo=googlecolab&logoColor=white" alt="Colab ok"></span></sub></td>
</tr>
<tr>
  <td><sub><b>ConformationalChange</b><br>Quantify backbone differences (RMSD, displacement) between paired structures.</sub></td>
  <td width="100" align="center"><sub><img src="https://img.shields.io/badge/BP-1ABC9C?style=flat-square" alt="BioPipelines native tool"></sub></td>
  <td width="80" align="center"><sub><img src="https://img.shields.io/badge/GPU-76B900?style=flat-square&logo=nvidia&logoColor=white" alt="GPU"></sub></td>
  <td width="115" align="center"><sub><span style="white-space:nowrap"><img src="https://img.shields.io/badge/HPC-2C3E50?style=flat-square&logo=linux&logoColor=white" alt="HPC ok"></span>&nbsp;<span style="white-space:nowrap"><img src="https://img.shields.io/badge/Colab-F9AB00?style=flat-square&logo=googlecolab&logoColor=white" alt="Colab ok"></span></sub></td>
</tr>
<tr>
  <td><sub><b>Contacts</b><br>List atom and residue contacts within and between chains.</sub></td>
  <td width="100" align="center"><sub><img src="https://img.shields.io/badge/BP-1ABC9C?style=flat-square" alt="BioPipelines native tool"></sub></td>
  <td width="80" align="center"><sub><img src="https://img.shields.io/badge/CPU-4C6EF5?style=flat-square&logo=data%3Aimage%2Fsvg%2Bxml%3Bbase64%2CPHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNCAyNCIgZmlsbD0id2hpdGUiPjxwYXRoIGQ9Ik05IDJ2Mkg3LjVBMS41IDEuNSAwIDAgMCA2IDUuNVY3SDR2MmgydjJINHYyaDJ2Mkg0djJoMnYxLjVBMS41IDEuNSAwIDAgMCA3LjUgMjBIOXYyaDJ2LTJoMnYyaDJ2LTJoMS41YTEuNSAxLjUgMCAwIDAgMS41LTEuNVYxOGgydi0yaC0ydi0yaDJ2LTJoLTJWOWgyVjdoLTJWNS41QTEuNSAxLjUgMCAwIDAgMTYuNSA0SDE1VjJoLTJ2MmgtMlYySDl6bS0xIDZoOHY4SDhWOHoiLz48L3N2Zz4%3D&logoColor=white" alt="CPU"></sub></td>
  <td width="115" align="center"><sub><span style="white-space:nowrap"><img src="https://img.shields.io/badge/HPC-2C3E50?style=flat-square&logo=linux&logoColor=white" alt="HPC ok"></span>&nbsp;<span style="white-space:nowrap"><img src="https://img.shields.io/badge/Colab-F9AB00?style=flat-square&logo=googlecolab&logoColor=white" alt="Colab ok"></span></sub></td>
</tr>
<tr>
  <td><sub><b>Distance</b><br>Measure distances between selected atoms or residues across structures.</sub></td>
  <td width="100" align="center"><sub><img src="https://img.shields.io/badge/BP-1ABC9C?style=flat-square" alt="BioPipelines native tool"></sub></td>
  <td width="80" align="center"><sub><img src="https://img.shields.io/badge/CPU-4C6EF5?style=flat-square&logo=data%3Aimage%2Fsvg%2Bxml%3Bbase64%2CPHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNCAyNCIgZmlsbD0id2hpdGUiPjxwYXRoIGQ9Ik05IDJ2Mkg3LjVBMS41IDEuNSAwIDAgMCA2IDUuNVY3SDR2MmgydjJINHYyaDJ2Mkg0djJoMnYxLjVBMS41IDEuNSAwIDAgMCA3LjUgMjBIOXYyaDJ2LTJoMnYyaDJ2LTJoMS41YTEuNSAxLjUgMCAwIDAgMS41LTEuNVYxOGgydi0yaC0ydi0yaDJ2LTJoLTJWOWgyVjdoLTJWNS41QTEuNSAxLjUgMCAwIDAgMTYuNSA0SDE1VjJoLTJ2MmgtMlYySDl6bS0xIDZoOHY4SDhWOHoiLz48L3N2Zz4%3D&logoColor=white" alt="CPU"></sub></td>
  <td width="115" align="center"><sub><span style="white-space:nowrap"><img src="https://img.shields.io/badge/HPC-2C3E50?style=flat-square&logo=linux&logoColor=white" alt="HPC ok"></span>&nbsp;<span style="white-space:nowrap"><img src="https://img.shields.io/badge/Colab-F9AB00?style=flat-square&logo=googlecolab&logoColor=white" alt="Colab ok"></span></sub></td>
</tr>
<tr>
  <td><sub><b>DistanceSelector</b><br>Select residues within a cutoff of an atom, residue, or ligand.</sub></td>
  <td width="100" align="center"><sub><img src="https://img.shields.io/badge/BP-1ABC9C?style=flat-square" alt="BioPipelines native tool"></sub></td>
  <td width="80" align="center"><sub><img src="https://img.shields.io/badge/CPU-4C6EF5?style=flat-square&logo=data%3Aimage%2Fsvg%2Bxml%3Bbase64%2CPHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNCAyNCIgZmlsbD0id2hpdGUiPjxwYXRoIGQ9Ik05IDJ2Mkg3LjVBMS41IDEuNSAwIDAgMCA2IDUuNVY3SDR2MmgydjJINHYyaDJ2Mkg0djJoMnYxLjVBMS41IDEuNSAwIDAgMCA3LjUgMjBIOXYyaDJ2LTJoMnYyaDJ2LTJoMS41YTEuNSAxLjUgMCAwIDAgMS41LTEuNVYxOGgydi0yaC0ydi0yaDJ2LTJoLTJWOWgyVjdoLTJWNS41QTEuNSAxLjUgMCAwIDAgMTYuNSA0SDE1VjJoLTJ2MmgtMlYySDl6bS0xIDZoOHY4SDhWOHoiLz48L3N2Zz4%3D&logoColor=white" alt="CPU"></sub></td>
  <td width="115" align="center"><sub><span style="white-space:nowrap"><img src="https://img.shields.io/badge/HPC-2C3E50?style=flat-square&logo=linux&logoColor=white" alt="HPC ok"></span>&nbsp;<span style="white-space:nowrap"><img src="https://img.shields.io/badge/Colab-F9AB00?style=flat-square&logo=googlecolab&logoColor=white" alt="Colab ok"></span></sub></td>
</tr>
<tr>
  <td><sub><b>DSSP</b><br>Assign per-residue secondary structure and overall helix, sheet, and coil fractions.</sub></td>
  <td width="100" align="center"><sub><a href="https://github.com/PDB-REDO/dssp"><img src="https://img.shields.io/badge/-181717?style=flat-square&logo=github&logoColor=white" alt="repo"></a> <a href="https://doi.org/10.1002/bip.360221211"><img src="https://img.shields.io/badge/-paper-B31B1B?style=flat-square&logo=readthedocs&logoColor=white" alt="paper"></a></sub></td>
  <td width="80" align="center"><sub><img src="https://img.shields.io/badge/CPU-4C6EF5?style=flat-square&logo=data%3Aimage%2Fsvg%2Bxml%3Bbase64%2CPHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNCAyNCIgZmlsbD0id2hpdGUiPjxwYXRoIGQ9Ik05IDJ2Mkg3LjVBMS41IDEuNSAwIDAgMCA2IDUuNVY3SDR2MmgydjJINHYyaDJ2Mkg0djJoMnYxLjVBMS41IDEuNSAwIDAgMCA3LjUgMjBIOXYyaDJ2LTJoMnYyaDJ2LTJoMS41YTEuNSAxLjUgMCAwIDAgMS41LTEuNVYxOGgydi0yaC0ydi0yaDJ2LTJoLTJWOWgyVjdoLTJWNS41QTEuNSAxLjUgMCAwIDAgMTYuNSA0SDE1VjJoLTJ2MmgtMlYySDl6bS0xIDZoOHY4SDhWOHoiLz48L3N2Zz4%3D&logoColor=white" alt="CPU"></sub></td>
  <td width="115" align="center"><sub><span style="white-space:nowrap"><img src="https://img.shields.io/badge/HPC-2C3E50?style=flat-square&logo=linux&logoColor=white" alt="HPC ok"></span>&nbsp;<span style="white-space:nowrap"><img src="https://img.shields.io/badge/Colab-F9AB00?style=flat-square&logo=googlecolab&logoColor=white" alt="Colab ok"></span></sub></td>
</tr>
<tr>
  <td><sub><b>EnsembleAnalysis</b><br>Compute per-residue RMSF and ensemble metrics from a conformer set.</sub></td>
  <td width="100" align="center"><sub><img src="https://img.shields.io/badge/BP-1ABC9C?style=flat-square" alt="BioPipelines native tool"></sub></td>
  <td width="80" align="center"><sub><img src="https://img.shields.io/badge/CPU-4C6EF5?style=flat-square&logo=data%3Aimage%2Fsvg%2Bxml%3Bbase64%2CPHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNCAyNCIgZmlsbD0id2hpdGUiPjxwYXRoIGQ9Ik05IDJ2Mkg3LjVBMS41IDEuNSAwIDAgMCA2IDUuNVY3SDR2MmgydjJINHYyaDJ2Mkg0djJoMnYxLjVBMS41IDEuNSAwIDAgMCA3LjUgMjBIOXYyaDJ2LTJoMnYyaDJ2LTJoMS41YTEuNSAxLjUgMCAwIDAgMS41LTEuNVYxOGgydi0yaC0ydi0yaDJ2LTJoLTJWOWgyVjdoLTJWNS41QTEuNSAxLjUgMCAwIDAgMTYuNSA0SDE1VjJoLTJ2MmgtMlYySDl6bS0xIDZoOHY4SDhWOHoiLz48L3N2Zz4%3D&logoColor=white" alt="CPU"></sub></td>
  <td width="115" align="center"><sub><span style="white-space:nowrap"><img src="https://img.shields.io/badge/HPC-2C3E50?style=flat-square&logo=linux&logoColor=white" alt="HPC ok"></span>&nbsp;<span style="white-space:nowrap"><img src="https://img.shields.io/badge/Colab-F9AB00?style=flat-square&logo=googlecolab&logoColor=white" alt="Colab ok"></span></sub></td>
</tr>
<tr>
  <td><sub><b>FPocket</b><br>Detect binding pockets and score their volume and druggability.</sub></td>
  <td width="100" align="center"><sub><a href="https://github.com/Discngine/fpocket"><img src="https://img.shields.io/badge/-181717?style=flat-square&logo=github&logoColor=white" alt="repo"></a> <a href="https://doi.org/10.1186/1471-2105-10-168"><img src="https://img.shields.io/badge/-paper-B31B1B?style=flat-square&logo=readthedocs&logoColor=white" alt="paper"></a></sub></td>
  <td width="80" align="center"><sub><img src="https://img.shields.io/badge/CPU-4C6EF5?style=flat-square&logo=data%3Aimage%2Fsvg%2Bxml%3Bbase64%2CPHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNCAyNCIgZmlsbD0id2hpdGUiPjxwYXRoIGQ9Ik05IDJ2Mkg3LjVBMS41IDEuNSAwIDAgMCA2IDUuNVY3SDR2MmgydjJINHYyaDJ2Mkg0djJoMnYxLjVBMS41IDEuNSAwIDAgMCA3LjUgMjBIOXYyaDJ2LTJoMnYyaDJ2LTJoMS41YTEuNSAxLjUgMCAwIDAgMS41LTEuNVYxOGgydi0yaC0ydi0yaDJ2LTJoLTJWOWgyVjdoLTJWNS41QTEuNSAxLjUgMCAwIDAgMTYuNSA0SDE1VjJoLTJ2MmgtMlYySDl6bS0xIDZoOHY4SDhWOHoiLz48L3N2Zz4%3D&logoColor=white" alt="CPU"></sub></td>
  <td width="115" align="center"><sub><span style="white-space:nowrap"><img src="https://img.shields.io/badge/HPC-2C3E50?style=flat-square&logo=linux&logoColor=white" alt="HPC ok"></span>&nbsp;<span style="white-space:nowrap"><img src="https://img.shields.io/badge/Colab-F9AB00?style=flat-square&logo=googlecolab&logoColor=white" alt="Colab ok"></span></sub></td>
</tr>
<tr>
  <td><sub><b>GEMS</b><br>Predict protein-ligand binding affinity from a complex structure.</sub></td>
  <td width="100" align="center"><sub><a href="https://github.com/camlab-ethz/GEMS"><img src="https://img.shields.io/badge/-181717?style=flat-square&logo=github&logoColor=white" alt="repo"></a> <a href="https://doi.org/10.1021/acs.jcim.4c01979"><img src="https://img.shields.io/badge/-paper-B31B1B?style=flat-square&logo=readthedocs&logoColor=white" alt="paper"></a></sub></td>
  <td width="80" align="center"><sub><img src="https://img.shields.io/badge/GPU-76B900?style=flat-square&logo=nvidia&logoColor=white" alt="GPU"></sub></td>
  <td width="115" align="center"><sub><span style="white-space:nowrap"><img src="https://img.shields.io/badge/Colab-F9AB00?style=flat-square&logo=googlecolab&logoColor=white" alt="Colab ok"></span></sub></td>
</tr>
<tr>
  <td><sub><b>OpenMM</b><br>Energy-minimise structures to relieve clashes before downstream metrics.</sub></td>
  <td width="100" align="center"><sub><a href="https://github.com/openmm/openmm"><img src="https://img.shields.io/badge/-181717?style=flat-square&logo=github&logoColor=white" alt="repo"></a> <a href="https://doi.org/10.1371/journal.pcbi.1005659"><img src="https://img.shields.io/badge/-paper-B31B1B?style=flat-square&logo=readthedocs&logoColor=white" alt="paper"></a></sub></td>
  <td width="80" align="center"><sub><img src="https://img.shields.io/badge/CPU-4C6EF5?style=flat-square&logo=data%3Aimage%2Fsvg%2Bxml%3Bbase64%2CPHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNCAyNCIgZmlsbD0id2hpdGUiPjxwYXRoIGQ9Ik05IDJ2Mkg3LjVBMS41IDEuNSAwIDAgMCA2IDUuNVY3SDR2MmgydjJINHYyaDJ2Mkg0djJoMnYxLjVBMS41IDEuNSAwIDAgMCA3LjUgMjBIOXYyaDJ2LTJoMnYyaDJ2LTJoMS41YTEuNSAxLjUgMCAwIDAgMS41LTEuNVYxOGgydi0yaC0ydi0yaDJ2LTJoLTJWOWgyVjdoLTJWNS41QTEuNSAxLjUgMCAwIDAgMTYuNSA0SDE1VjJoLTJ2MmgtMlYySDl6bS0xIDZoOHY4SDhWOHoiLz48L3N2Zz4%3D&logoColor=white" alt="CPU"> <img src="https://img.shields.io/badge/GPU-76B900?style=flat-square&logo=nvidia&logoColor=white" alt="GPU"></sub></td>
  <td width="115" align="center"><sub><span style="white-space:nowrap"><img src="https://img.shields.io/badge/HPC-2C3E50?style=flat-square&logo=linux&logoColor=white" alt="HPC ok"></span>&nbsp;<span style="white-space:nowrap"><img src="https://img.shields.io/badge/Colab-F9AB00?style=flat-square&logo=googlecolab&logoColor=white" alt="Colab ok"></span></sub></td>
</tr>
<tr>
  <td><sub><b>P2Rank</b><sup><i>p</i></sup><br>Predict ligand binding sites on a protein surface from structure.</sub></td>
  <td width="100" align="center"><sub><a href="https://github.com/rdk/p2rank"><img src="https://img.shields.io/badge/-181717?style=flat-square&logo=github&logoColor=white" alt="repo"></a> <a href="https://doi.org/10.1186/s13321-018-0285-8"><img src="https://img.shields.io/badge/-paper-B31B1B?style=flat-square&logo=readthedocs&logoColor=white" alt="paper"></a></sub></td>
  <td width="80" align="center"><sub><img src="https://img.shields.io/badge/CPU-4C6EF5?style=flat-square&logo=data%3Aimage%2Fsvg%2Bxml%3Bbase64%2CPHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNCAyNCIgZmlsbD0id2hpdGUiPjxwYXRoIGQ9Ik05IDJ2Mkg3LjVBMS41IDEuNSAwIDAgMCA2IDUuNVY3SDR2MmgydjJINHYyaDJ2Mkg0djJoMnYxLjVBMS41IDEuNSAwIDAgMCA3LjUgMjBIOXYyaDJ2LTJoMnYyaDJ2LTJoMS41YTEuNSAxLjUgMCAwIDAgMS41LTEuNVYxOGgydi0yaC0ydi0yaDJ2LTJoLTJWOWgyVjdoLTJWNS41QTEuNSAxLjUgMCAwIDAgMTYuNSA0SDE1VjJoLTJ2MmgtMlYySDl6bS0xIDZoOHY4SDhWOHoiLz48L3N2Zz4%3D&logoColor=white" alt="CPU"></sub></td>
  <td width="115" align="center"><sub><span style="white-space:nowrap"><img src="https://img.shields.io/badge/HPC-2C3E50?style=flat-square&logo=linux&logoColor=white" alt="HPC ok"></span>&nbsp;<span style="white-space:nowrap"><img src="https://img.shields.io/badge/Colab-F9AB00?style=flat-square&logo=googlecolab&logoColor=white" alt="Colab ok"></span></sub></td>
</tr>
<tr>
  <td><sub><b>PLIP</b><sup><i>k</i></sup><br>Profile non-covalent protein-ligand interactions directly from a complex.</sub></td>
  <td width="100" align="center"><sub><a href="https://github.com/pharmai/plip"><img src="https://img.shields.io/badge/-181717?style=flat-square&logo=github&logoColor=white" alt="repo"></a> <a href="https://doi.org/10.1093/nar/gkab294"><img src="https://img.shields.io/badge/-paper-B31B1B?style=flat-square&logo=readthedocs&logoColor=white" alt="paper"></a></sub></td>
  <td width="80" align="center"><sub><img src="https://img.shields.io/badge/CPU-4C6EF5?style=flat-square&logo=data%3Aimage%2Fsvg%2Bxml%3Bbase64%2CPHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNCAyNCIgZmlsbD0id2hpdGUiPjxwYXRoIGQ9Ik05IDJ2Mkg3LjVBMS41IDEuNSAwIDAgMCA2IDUuNVY3SDR2MmgydjJINHYyaDJ2Mkg0djJoMnYxLjVBMS41IDEuNSAwIDAgMCA3LjUgMjBIOXYyaDJ2LTJoMnYyaDJ2LTJoMS41YTEuNSAxLjUgMCAwIDAgMS41LTEuNVYxOGgydi0yaC0ydi0yaDJ2LTJoLTJWOWgyVjdoLTJWNS41QTEuNSAxLjUgMCAwIDAgMTYuNSA0SDE1VjJoLTJ2MmgtMlYySDl6bS0xIDZoOHY4SDhWOHoiLz48L3N2Zz4%3D&logoColor=white" alt="CPU"></sub></td>
  <td width="115" align="center"><sub><span style="white-space:nowrap"><img src="https://img.shields.io/badge/HPC-2C3E50?style=flat-square&logo=linux&logoColor=white" alt="HPC ok"></span>&nbsp;<span style="white-space:nowrap"><img src="https://img.shields.io/badge/Colab-F9AB00?style=flat-square&logo=googlecolab&logoColor=white" alt="Colab ok"></span></sub></td>
</tr>
<tr>
  <td><sub><b>PLM_Sol</b><br>Predict protein solubility from sequence alone via protein-language-model embeddings.</sub></td>
  <td width="100" align="center"><sub><a href="https://github.com/Violet969/PLM_Sol"><img src="https://img.shields.io/badge/-181717?style=flat-square&logo=github&logoColor=white" alt="repo"></a> <a href="https://academic.oup.com/bib/article/25/5/bbae404/7739950"><img src="https://img.shields.io/badge/-paper-B31B1B?style=flat-square&logo=readthedocs&logoColor=white" alt="paper"></a></sub></td>
  <td width="80" align="center"><sub><img src="https://img.shields.io/badge/GPU-76B900?style=flat-square&logo=nvidia&logoColor=white" alt="GPU"></sub></td>
  <td width="115" align="center"><sub><span style="white-space:nowrap"><img src="https://img.shields.io/badge/HPC-2C3E50?style=flat-square&logo=linux&logoColor=white" alt="HPC ok"></span>&nbsp;<span style="white-space:nowrap"><img src="https://img.shields.io/badge/Colab-F9AB00?style=flat-square&logo=googlecolab&logoColor=white" alt="Colab ok"></span></sub></td>
</tr>
<tr>
  <td><sub><b>PoseBusters</b><br>Flag physically implausible poses: clashes, bad geometry, wrong stereochemistry.</sub></td>
  <td width="100" align="center"><sub><a href="https://github.com/maabuu/posebusters"><img src="https://img.shields.io/badge/-181717?style=flat-square&logo=github&logoColor=white" alt="repo"></a> <a href="https://doi.org/10.1039/D3SC04185A"><img src="https://img.shields.io/badge/-paper-B31B1B?style=flat-square&logo=readthedocs&logoColor=white" alt="paper"></a></sub></td>
  <td width="80" align="center"><sub><img src="https://img.shields.io/badge/CPU-4C6EF5?style=flat-square&logo=data%3Aimage%2Fsvg%2Bxml%3Bbase64%2CPHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNCAyNCIgZmlsbD0id2hpdGUiPjxwYXRoIGQ9Ik05IDJ2Mkg3LjVBMS41IDEuNSAwIDAgMCA2IDUuNVY3SDR2MmgydjJINHYyaDJ2Mkg0djJoMnYxLjVBMS41IDEuNSAwIDAgMCA3LjUgMjBIOXYyaDJ2LTJoMnYyaDJ2LTJoMS41YTEuNSAxLjUgMCAwIDAgMS41LTEuNVYxOGgydi0yaC0ydi0yaDJ2LTJoLTJWOWgyVjdoLTJWNS41QTEuNSAxLjUgMCAwIDAgMTYuNSA0SDE1VjJoLTJ2MmgtMlYySDl6bS0xIDZoOHY4SDhWOHoiLz48L3N2Zz4%3D&logoColor=white" alt="CPU"></sub></td>
  <td width="115" align="center"><sub><span style="white-space:nowrap"><img src="https://img.shields.io/badge/HPC-2C3E50?style=flat-square&logo=linux&logoColor=white" alt="HPC ok"></span>&nbsp;<span style="white-space:nowrap"><img src="https://img.shields.io/badge/Colab-F9AB00?style=flat-square&logo=googlecolab&logoColor=white" alt="Colab ok"></span></sub></td>
</tr>
<tr>
  <td><sub><b>PoseChange</b><br>Compare a ligand's pose between two structures of one complex.</sub></td>
  <td width="100" align="center"><sub><img src="https://img.shields.io/badge/BP-1ABC9C?style=flat-square" alt="BioPipelines native tool"></sub></td>
  <td width="80" align="center"><sub><img src="https://img.shields.io/badge/GPU-76B900?style=flat-square&logo=nvidia&logoColor=white" alt="GPU"></sub></td>
  <td width="115" align="center"><sub><span style="white-space:nowrap"><img src="https://img.shields.io/badge/HPC-2C3E50?style=flat-square&logo=linux&logoColor=white" alt="HPC ok"></span>&nbsp;<span style="white-space:nowrap"><img src="https://img.shields.io/badge/Colab-F9AB00?style=flat-square&logo=googlecolab&logoColor=white" alt="Colab ok"></span></sub></td>
</tr>
<tr>
  <td><sub><b>Prodigy</b><br>Predict protein-protein binding affinity (Kd, dG) from contacts.</sub></td>
  <td width="100" align="center"><sub><a href="https://github.com/haddocking/prodigy"><img src="https://img.shields.io/badge/-181717?style=flat-square&logo=github&logoColor=white" alt="repo"></a> <a href="https://doi.org/10.7554/eLife.07454"><img src="https://img.shields.io/badge/-paper-B31B1B?style=flat-square&logo=readthedocs&logoColor=white" alt="paper"></a></sub></td>
  <td width="80" align="center"><sub><img src="https://img.shields.io/badge/CPU-4C6EF5?style=flat-square&logo=data%3Aimage%2Fsvg%2Bxml%3Bbase64%2CPHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNCAyNCIgZmlsbD0id2hpdGUiPjxwYXRoIGQ9Ik05IDJ2Mkg3LjVBMS41IDEuNSAwIDAgMCA2IDUuNVY3SDR2MmgydjJINHYyaDJ2Mkg0djJoMnYxLjVBMS41IDEuNSAwIDAgMCA3LjUgMjBIOXYyaDJ2LTJoMnYyaDJ2LTJoMS41YTEuNSAxLjUgMCAwIDAgMS41LTEuNVYxOGgydi0yaC0ydi0yaDJ2LTJoLTJWOWgyVjdoLTJWNS41QTEuNSAxLjUgMCAwIDAgMTYuNSA0SDE1VjJoLTJ2MmgtMlYySDl6bS0xIDZoOHY4SDhWOHoiLz48L3N2Zz4%3D&logoColor=white" alt="CPU"></sub></td>
  <td width="115" align="center"><sub><span style="white-space:nowrap"><img src="https://img.shields.io/badge/HPC-2C3E50?style=flat-square&logo=linux&logoColor=white" alt="HPC ok"></span>&nbsp;<span style="white-space:nowrap"><img src="https://img.shields.io/badge/Colab-F9AB00?style=flat-square&logo=googlecolab&logoColor=white" alt="Colab ok"></span></sub></td>
</tr>
<tr>
  <td><sub><b>ProLIF</b><br>Build interaction fingerprints (H-bonds, hydrophobic, pi-stacking) for complexes.</sub></td>
  <td width="100" align="center"><sub><a href="https://github.com/chemosim-lab/ProLIF"><img src="https://img.shields.io/badge/-181717?style=flat-square&logo=github&logoColor=white" alt="repo"></a> <a href="https://doi.org/10.1186/s13321-021-00548-6"><img src="https://img.shields.io/badge/-paper-B31B1B?style=flat-square&logo=readthedocs&logoColor=white" alt="paper"></a></sub></td>
  <td width="80" align="center"><sub><img src="https://img.shields.io/badge/CPU-4C6EF5?style=flat-square&logo=data%3Aimage%2Fsvg%2Bxml%3Bbase64%2CPHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNCAyNCIgZmlsbD0id2hpdGUiPjxwYXRoIGQ9Ik05IDJ2Mkg3LjVBMS41IDEuNSAwIDAgMCA2IDUuNVY3SDR2MmgydjJINHYyaDJ2Mkg0djJoMnYxLjVBMS41IDEuNSAwIDAgMCA3LjUgMjBIOXYyaDJ2LTJoMnYyaDJ2LTJoMS41YTEuNSAxLjUgMCAwIDAgMS41LTEuNVYxOGgydi0yaC0ydi0yaDJ2LTJoLTJWOWgyVjdoLTJWNS41QTEuNSAxLjUgMCAwIDAgMTYuNSA0SDE1VjJoLTJ2MmgtMlYySDl6bS0xIDZoOHY4SDhWOHoiLz48L3N2Zz4%3D&logoColor=white" alt="CPU"></sub></td>
  <td width="115" align="center"><sub><span style="white-space:nowrap"><img src="https://img.shields.io/badge/HPC-2C3E50?style=flat-square&logo=linux&logoColor=white" alt="HPC ok"></span>&nbsp;<span style="white-space:nowrap"><img src="https://img.shields.io/badge/Colab-F9AB00?style=flat-square&logo=googlecolab&logoColor=white" alt="Colab ok"></span></sub></td>
</tr>
<tr>
  <td><sub><b>Reduce</b><br>Add and optimise explicit hydrogens on protein and ligand atoms.</sub></td>
  <td width="100" align="center"><sub><a href="https://github.com/rlabduke/reduce"><img src="https://img.shields.io/badge/-181717?style=flat-square&logo=github&logoColor=white" alt="repo"></a> <a href="https://doi.org/10.1006/jmbi.1998.2401"><img src="https://img.shields.io/badge/-paper-B31B1B?style=flat-square&logo=readthedocs&logoColor=white" alt="paper"></a></sub></td>
  <td width="80" align="center"><sub><img src="https://img.shields.io/badge/CPU-4C6EF5?style=flat-square&logo=data%3Aimage%2Fsvg%2Bxml%3Bbase64%2CPHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNCAyNCIgZmlsbD0id2hpdGUiPjxwYXRoIGQ9Ik05IDJ2Mkg3LjVBMS41IDEuNSAwIDAgMCA2IDUuNVY3SDR2MmgydjJINHYyaDJ2Mkg0djJoMnYxLjVBMS41IDEuNSAwIDAgMCA3LjUgMjBIOXYyaDJ2LTJoMnYyaDJ2LTJoMS41YTEuNSAxLjUgMCAwIDAgMS41LTEuNVYxOGgydi0yaC0ydi0yaDJ2LTJoLTJWOWgyVjdoLTJWNS41QTEuNSAxLjUgMCAwIDAgMTYuNSA0SDE1VjJoLTJ2MmgtMlYySDl6bS0xIDZoOHY4SDhWOHoiLz48L3N2Zz4%3D&logoColor=white" alt="CPU"></sub></td>
  <td width="115" align="center"><sub><span style="white-space:nowrap"><img src="https://img.shields.io/badge/HPC-2C3E50?style=flat-square&logo=linux&logoColor=white" alt="HPC ok"></span>&nbsp;<span style="white-space:nowrap"><img src="https://img.shields.io/badge/Colab-F9AB00?style=flat-square&logo=googlecolab&logoColor=white" alt="Colab ok"></span></sub></td>
</tr>
<tr>
  <td><sub><b>RTMScore</b><br>Score and rank docking poses by predicted binding likelihood.</sub></td>
  <td width="100" align="center"><sub><a href="https://github.com/sc8668/RTMScore"><img src="https://img.shields.io/badge/-181717?style=flat-square&logo=github&logoColor=white" alt="repo"></a> <a href="https://doi.org/10.1021/acs.jmedchem.2c00991"><img src="https://img.shields.io/badge/-paper-B31B1B?style=flat-square&logo=readthedocs&logoColor=white" alt="paper"></a></sub></td>
  <td width="80" align="center"><sub><img src="https://img.shields.io/badge/GPU-76B900?style=flat-square&logo=nvidia&logoColor=white" alt="GPU"></sub></td>
  <td width="115" align="center"><sub><span style="white-space:nowrap"><img src="https://img.shields.io/badge/Colab-F9AB00?style=flat-square&logo=googlecolab&logoColor=white" alt="Colab ok"></span></sub></td>
</tr>
<tr>
  <td><sub><b>SASA</b><br>Compute per-residue and total solvent-accessible surface area of a structure.</sub></td>
  <td width="100" align="center"><sub><img src="https://img.shields.io/badge/BP-1ABC9C?style=flat-square" alt="BioPipelines native tool"></sub></td>
  <td width="80" align="center"><sub><img src="https://img.shields.io/badge/GPU-76B900?style=flat-square&logo=nvidia&logoColor=white" alt="GPU"></sub></td>
  <td width="115" align="center"><sub><span style="white-space:nowrap"><img src="https://img.shields.io/badge/HPC-2C3E50?style=flat-square&logo=linux&logoColor=white" alt="HPC ok"></span>&nbsp;<span style="white-space:nowrap"><img src="https://img.shields.io/badge/Colab-F9AB00?style=flat-square&logo=googlecolab&logoColor=white" alt="Colab ok"></span></sub></td>
</tr>
<tr>
  <td><sub><b>ThermoMPNN</b><br>Predict fold-stability change (ddG) for each point mutation.</sub></td>
  <td width="100" align="center"><sub><a href="https://github.com/Kuhlman-Lab/ThermoMPNN"><img src="https://img.shields.io/badge/-181717?style=flat-square&logo=github&logoColor=white" alt="repo"></a> <a href="https://www.pnas.org/doi/10.1073/pnas.2314853121"><img src="https://img.shields.io/badge/-paper-B31B1B?style=flat-square&logo=readthedocs&logoColor=white" alt="paper"></a></sub></td>
  <td width="80" align="center"><sub><img src="https://img.shields.io/badge/CPU-4C6EF5?style=flat-square&logo=data%3Aimage%2Fsvg%2Bxml%3Bbase64%2CPHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNCAyNCIgZmlsbD0id2hpdGUiPjxwYXRoIGQ9Ik05IDJ2Mkg3LjVBMS41IDEuNSAwIDAgMCA2IDUuNVY3SDR2MmgydjJINHYyaDJ2Mkg0djJoMnYxLjVBMS41IDEuNSAwIDAgMCA3LjUgMjBIOXYyaDJ2LTJoMnYyaDJ2LTJoMS41YTEuNSAxLjUgMCAwIDAgMS41LTEuNVYxOGgydi0yaC0ydi0yaDJ2LTJoLTJWOWgyVjdoLTJWNS41QTEuNSAxLjUgMCAwIDAgMTYuNSA0SDE1VjJoLTJ2MmgtMlYySDl6bS0xIDZoOHY4SDhWOHoiLz48L3N2Zz4%3D&logoColor=white" alt="CPU"> <img src="https://img.shields.io/badge/GPU-76B900?style=flat-square&logo=nvidia&logoColor=white" alt="GPU"></sub></td>
  <td width="115" align="center"><sub><span style="white-space:nowrap"><img src="https://img.shields.io/badge/HPC-2C3E50?style=flat-square&logo=linux&logoColor=white" alt="HPC ok"></span>&nbsp;<span style="white-space:nowrap"><img src="https://img.shields.io/badge/Colab-F9AB00?style=flat-square&logo=googlecolab&logoColor=white" alt="Colab ok"></span></sub></td>
</tr>
<tr>
  <td><sub><b>VespaG</b><sup><i>v</i></sup><br>Predict single-substitution fitness effects from sequence, GPU-fast and MSA-free.</sub></td>
  <td width="100" align="center"><sub><a href="https://github.com/JSchlensok/VespaG"><img src="https://img.shields.io/badge/-181717?style=flat-square&logo=github&logoColor=white" alt="repo"></a> <a href="https://academic.oup.com/bioinformatics/article/40/11/btae621/7907184"><img src="https://img.shields.io/badge/-paper-B31B1B?style=flat-square&logo=readthedocs&logoColor=white" alt="paper"></a></sub></td>
  <td width="80" align="center"><sub><img src="https://img.shields.io/badge/GPU-76B900?style=flat-square&logo=nvidia&logoColor=white" alt="GPU"></sub></td>
  <td width="115" align="center"><sub><span style="white-space:nowrap"><img src="https://img.shields.io/badge/HPC-2C3E50?style=flat-square&logo=linux&logoColor=white" alt="HPC ok"></span>&nbsp;<span style="white-space:nowrap"><img src="https://img.shields.io/badge/Colab-F9AB00?style=flat-square&logo=googlecolab&logoColor=white" alt="Colab ok"></span></sub></td>
</tr>
<tr>
  <td><sub><b>XTB</b><br>Score protein-ligand interaction energy with semi-empirical quantum chemistry.</sub></td>
  <td width="100" align="center"><sub><a href="https://github.com/grimme-lab/xtb"><img src="https://img.shields.io/badge/-181717?style=flat-square&logo=github&logoColor=white" alt="repo"></a> <a href="https://doi.org/10.1002/wcms.1493"><img src="https://img.shields.io/badge/-paper-B31B1B?style=flat-square&logo=readthedocs&logoColor=white" alt="paper"></a></sub></td>
  <td width="80" align="center"><sub><img src="https://img.shields.io/badge/CPU-4C6EF5?style=flat-square&logo=data%3Aimage%2Fsvg%2Bxml%3Bbase64%2CPHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNCAyNCIgZmlsbD0id2hpdGUiPjxwYXRoIGQ9Ik05IDJ2Mkg3LjVBMS41IDEuNSAwIDAgMCA2IDUuNVY3SDR2MmgydjJINHYyaDJ2Mkg0djJoMnYxLjVBMS41IDEuNSAwIDAgMCA3LjUgMjBIOXYyaDJ2LTJoMnYyaDJ2LTJoMS41YTEuNSAxLjUgMCAwIDAgMS41LTEuNVYxOGgydi0yaC0ydi0yaDJ2LTJoLTJWOWgyVjdoLTJWNS41QTEuNSAxLjUgMCAwIDAgMTYuNSA0SDE1VjJoLTJ2MmgtMlYySDl6bS0xIDZoOHY4SDhWOHoiLz48L3N2Zz4%3D&logoColor=white" alt="CPU"></sub></td>
  <td width="115" align="center"><sub><span style="white-space:nowrap"><img src="https://img.shields.io/badge/HPC-2C3E50?style=flat-square&logo=linux&logoColor=white" alt="HPC ok"></span>&nbsp;<span style="white-space:nowrap"><img src="https://img.shields.io/badge/Colab-F9AB00?style=flat-square&logo=googlecolab&logoColor=white" alt="Colab ok"></span></sub></td>
</tr>

<tr><td colspan="4"><h3>🧪 Cheminformatics</h3></td></tr>

<tr>
  <td><sub><b>OpenBabel</b><br>Convert molecules between formats; add pH-aware hydrogens and 3D coordinates.</sub></td>
  <td width="100" align="center"><sub><a href="https://github.com/openbabel/openbabel"><img src="https://img.shields.io/badge/-181717?style=flat-square&logo=github&logoColor=white" alt="repo"></a> <a href="https://doi.org/10.1186/1758-2946-3-33"><img src="https://img.shields.io/badge/-paper-B31B1B?style=flat-square&logo=readthedocs&logoColor=white" alt="paper"></a></sub></td>
  <td width="80" align="center"><sub><img src="https://img.shields.io/badge/CPU-4C6EF5?style=flat-square&logo=data%3Aimage%2Fsvg%2Bxml%3Bbase64%2CPHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNCAyNCIgZmlsbD0id2hpdGUiPjxwYXRoIGQ9Ik05IDJ2Mkg3LjVBMS41IDEuNSAwIDAgMCA2IDUuNVY3SDR2MmgydjJINHYyaDJ2Mkg0djJoMnYxLjVBMS41IDEuNSAwIDAgMCA3LjUgMjBIOXYyaDJ2LTJoMnYyaDJ2LTJoMS41YTEuNSAxLjUgMCAwIDAgMS41LTEuNVYxOGgydi0yaC0ydi0yaDJ2LTJoLTJWOWgyVjdoLTJWNS41QTEuNSAxLjUgMCAwIDAgMTYuNSA0SDE1VjJoLTJ2MmgtMlYySDl6bS0xIDZoOHY4SDhWOHoiLz48L3N2Zz4%3D&logoColor=white" alt="CPU"></sub></td>
  <td width="115" align="center"><sub><span style="white-space:nowrap"><img src="https://img.shields.io/badge/HPC-2C3E50?style=flat-square&logo=linux&logoColor=white" alt="HPC ok"></span>&nbsp;<span style="white-space:nowrap"><img src="https://img.shields.io/badge/Colab-F9AB00?style=flat-square&logo=googlecolab&logoColor=white" alt="Colab ok"></span></sub></td>
</tr>
<tr>
  <td><sub><b>RDKit</b><br>Compute molecular descriptors (MW, logP, TPSA, QED) from SMILES.</sub></td>
  <td width="100" align="center"><sub><a href="https://github.com/rdkit/rdkit"><img src="https://img.shields.io/badge/-181717?style=flat-square&logo=github&logoColor=white" alt="repo"></a></sub></td>
  <td width="80" align="center"><sub><img src="https://img.shields.io/badge/CPU-4C6EF5?style=flat-square&logo=data%3Aimage%2Fsvg%2Bxml%3Bbase64%2CPHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNCAyNCIgZmlsbD0id2hpdGUiPjxwYXRoIGQ9Ik05IDJ2Mkg3LjVBMS41IDEuNSAwIDAgMCA2IDUuNVY3SDR2MmgydjJINHYyaDJ2Mkg0djJoMnYxLjVBMS41IDEuNSAwIDAgMCA3LjUgMjBIOXYyaDJ2LTJoMnYyaDJ2LTJoMS41YTEuNSAxLjUgMCAwIDAgMS41LTEuNVYxOGgydi0yaC0ydi0yaDJ2LTJoLTJWOWgyVjdoLTJWNS41QTEuNSAxLjUgMCAwIDAgMTYuNSA0SDE1VjJoLTJ2MmgtMlYySDl6bS0xIDZoOHY4SDhWOHoiLz48L3N2Zz4%3D&logoColor=white" alt="CPU"></sub></td>
  <td width="115" align="center"><sub><span style="white-space:nowrap"><img src="https://img.shields.io/badge/HPC-2C3E50?style=flat-square&logo=linux&logoColor=white" alt="HPC ok"></span>&nbsp;<span style="white-space:nowrap"><img src="https://img.shields.io/badge/Colab-F9AB00?style=flat-square&logo=googlecolab&logoColor=white" alt="Colab ok"></span></sub></td>
</tr>

<tr><td colspan="4"><h3>📊 Sequence Statistics</h3></td></tr>

<tr>
  <td><sub><b>BayesianAdjuster</b><br>Refine observed mutation frequencies using correlation-informed Bayesian updates across positions.</sub></td>
  <td width="100" align="center"><sub><img src="https://img.shields.io/badge/BP-1ABC9C?style=flat-square" alt="BioPipelines native tool"></sub></td>
  <td width="80" align="center"><sub><img src="https://img.shields.io/badge/CPU-4C6EF5?style=flat-square&logo=data%3Aimage%2Fsvg%2Bxml%3Bbase64%2CPHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNCAyNCIgZmlsbD0id2hpdGUiPjxwYXRoIGQ9Ik05IDJ2Mkg3LjVBMS41IDEuNSAwIDAgMCA2IDUuNVY3SDR2MmgydjJINHYyaDJ2Mkg0djJoMnYxLjVBMS41IDEuNSAwIDAgMCA3LjUgMjBIOXYyaDJ2LTJoMnYyaDJ2LTJoMS41YTEuNSAxLjUgMCAwIDAgMS41LTEuNVYxOGgydi0yaC0ydi0yaDJ2LTJoLTJWOWgyVjdoLTJWNS41QTEuNSAxLjUgMCAwIDAgMTYuNSA0SDE1VjJoLTJ2MmgtMlYySDl6bS0xIDZoOHY4SDhWOHoiLz48L3N2Zz4%3D&logoColor=white" alt="CPU"></sub></td>
  <td width="115" align="center"><sub><span style="white-space:nowrap"><img src="https://img.shields.io/badge/HPC-2C3E50?style=flat-square&logo=linux&logoColor=white" alt="HPC ok"></span>&nbsp;<span style="white-space:nowrap"><img src="https://img.shields.io/badge/Colab-F9AB00?style=flat-square&logo=googlecolab&logoColor=white" alt="Colab ok"></span></sub></td>
</tr>
<tr>
  <td><sub><b>MutationProfiler</b><br>Profile mutation patterns and per-position frequencies across a sequence set.</sub></td>
  <td width="100" align="center"><sub><img src="https://img.shields.io/badge/BP-1ABC9C?style=flat-square" alt="BioPipelines native tool"></sub></td>
  <td width="80" align="center"><sub><img src="https://img.shields.io/badge/CPU-4C6EF5?style=flat-square&logo=data%3Aimage%2Fsvg%2Bxml%3Bbase64%2CPHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNCAyNCIgZmlsbD0id2hpdGUiPjxwYXRoIGQ9Ik05IDJ2Mkg3LjVBMS41IDEuNSAwIDAgMCA2IDUuNVY3SDR2MmgydjJINHYyaDJ2Mkg0djJoMnYxLjVBMS41IDEuNSAwIDAgMCA3LjUgMjBIOXYyaDJ2LTJoMnYyaDJ2LTJoMS41YTEuNSAxLjUgMCAwIDAgMS41LTEuNVYxOGgydi0yaC0ydi0yaDJ2LTJoLTJWOWgyVjdoLTJWNS41QTEuNSAxLjUgMCAwIDAgMTYuNSA0SDE1VjJoLTJ2MmgtMlYySDl6bS0xIDZoOHY4SDhWOHoiLz48L3N2Zz4%3D&logoColor=white" alt="CPU"></sub></td>
  <td width="115" align="center"><sub><span style="white-space:nowrap"><img src="https://img.shields.io/badge/HPC-2C3E50?style=flat-square&logo=linux&logoColor=white" alt="HPC ok"></span>&nbsp;<span style="white-space:nowrap"><img src="https://img.shields.io/badge/Colab-F9AB00?style=flat-square&logo=googlecolab&logoColor=white" alt="Colab ok"></span></sub></td>
</tr>
<tr>
  <td><sub><b>SequenceMetricCorrelation</b><br>Correlate individual mutations with downstream pipeline metrics to find drivers.</sub></td>
  <td width="100" align="center"><sub><img src="https://img.shields.io/badge/BP-1ABC9C?style=flat-square" alt="BioPipelines native tool"></sub></td>
  <td width="80" align="center"><sub><img src="https://img.shields.io/badge/CPU-4C6EF5?style=flat-square&logo=data%3Aimage%2Fsvg%2Bxml%3Bbase64%2CPHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNCAyNCIgZmlsbD0id2hpdGUiPjxwYXRoIGQ9Ik05IDJ2Mkg3LjVBMS41IDEuNSAwIDAgMCA2IDUuNVY3SDR2MmgydjJINHYyaDJ2Mkg0djJoMnYxLjVBMS41IDEuNSAwIDAgMCA3LjUgMjBIOXYyaDJ2LTJoMnYyaDJ2LTJoMS41YTEuNSAxLjUgMCAwIDAgMS41LTEuNVYxOGgydi0yaC0ydi0yaDJ2LTJoLTJWOWgyVjdoLTJWNS41QTEuNSAxLjUgMCAwIDAgMTYuNSA0SDE1VjJoLTJ2MmgtMlYySDl6bS0xIDZoOHY4SDhWOHoiLz48L3N2Zz4%3D&logoColor=white" alt="CPU"></sub></td>
  <td width="115" align="center"><sub><span style="white-space:nowrap"><img src="https://img.shields.io/badge/HPC-2C3E50?style=flat-square&logo=linux&logoColor=white" alt="HPC ok"></span>&nbsp;<span style="white-space:nowrap"><img src="https://img.shields.io/badge/Colab-F9AB00?style=flat-square&logo=googlecolab&logoColor=white" alt="Colab ok"></span></sub></td>
</tr>

<tr><td colspan="4"><h3>🧷 MSAs</h3></td></tr>

<tr>
  <td><sub><b>MMseqs2 / MMseqs2Server</b><br>Build multiple-sequence alignments for folding by fast homology search.</sub></td>
  <td width="100" align="center"><sub><a href="https://github.com/soedinglab/MMseqs2"><img src="https://img.shields.io/badge/-181717?style=flat-square&logo=github&logoColor=white" alt="repo"></a> <a href="https://www.nature.com/articles/nbt.3988"><img src="https://img.shields.io/badge/-paper-B31B1B?style=flat-square&logo=readthedocs&logoColor=white" alt="paper"></a></sub></td>
  <td width="80" align="center"><sub><img src="https://img.shields.io/badge/CPU-4C6EF5?style=flat-square&logo=data%3Aimage%2Fsvg%2Bxml%3Bbase64%2CPHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNCAyNCIgZmlsbD0id2hpdGUiPjxwYXRoIGQ9Ik05IDJ2Mkg3LjVBMS41IDEuNSAwIDAgMCA2IDUuNVY3SDR2MmgydjJINHYyaDJ2Mkg0djJoMnYxLjVBMS41IDEuNSAwIDAgMCA3LjUgMjBIOXYyaDJ2LTJoMnYyaDJ2LTJoMS41YTEuNSAxLjUgMCAwIDAgMS41LTEuNVYxOGgydi0yaC0ydi0yaDJ2LTJoLTJWOWgyVjdoLTJWNS41QTEuNSAxLjUgMCAwIDAgMTYuNSA0SDE1VjJoLTJ2MmgtMlYySDl6bS0xIDZoOHY4SDhWOHoiLz48L3N2Zz4%3D&logoColor=white" alt="CPU"> </sub></td>
  <td width="115" align="center"><sub><span style="white-space:nowrap"><img src="https://img.shields.io/badge/HPC-2C3E50?style=flat-square&logo=linux&logoColor=white" alt="HPC ok"></span></sub></td>
</tr>
<tr>
  <td><sub><b>MSA</b><br>Convert and manipulate multiple-sequence alignments between CSV and A3M.</sub></td>
  <td width="100" align="center"><sub><img src="https://img.shields.io/badge/BP-1ABC9C?style=flat-square" alt="BioPipelines native tool"></sub></td>
  <td width="80" align="center"><sub><img src="https://img.shields.io/badge/CPU-4C6EF5?style=flat-square&logo=data%3Aimage%2Fsvg%2Bxml%3Bbase64%2CPHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNCAyNCIgZmlsbD0id2hpdGUiPjxwYXRoIGQ9Ik05IDJ2Mkg3LjVBMS41IDEuNSAwIDAgMCA2IDUuNVY3SDR2MmgydjJINHYyaDJ2Mkg0djJoMnYxLjVBMS41IDEuNSAwIDAgMCA3LjUgMjBIOXYyaDJ2LTJoMnYyaDJ2LTJoMS41YTEuNSAxLjUgMCAwIDAgMS41LTEuNVYxOGgydi0yaC0ydi0yaDJ2LTJoLTJWOWgyVjdoLTJWNS41QTEuNSAxLjUgMCAwIDAgMTYuNSA0SDE1VjJoLTJ2MmgtMlYySDl6bS0xIDZoOHY4SDhWOHoiLz48L3N2Zz4%3D&logoColor=white" alt="CPU"></sub></td>
  <td width="115" align="center"><sub><span style="white-space:nowrap"><img src="https://img.shields.io/badge/HPC-2C3E50?style=flat-square&logo=linux&logoColor=white" alt="HPC ok"></span>&nbsp;<span style="white-space:nowrap"><img src="https://img.shields.io/badge/Colab-F9AB00?style=flat-square&logo=googlecolab&logoColor=white" alt="Colab ok"></span></sub></td>
</tr>

<tr><td colspan="4"><h3>🗂️ Data Management</h3></td></tr>
<tr>
  <td><sub><b>Consensus</b><br>Aggregate a per-residue stream across grouped structures into consensus columns.</sub></td>
  <td width="100" align="center"><sub><img src="https://img.shields.io/badge/BP-1ABC9C?style=flat-square" alt="BioPipelines native tool"></sub></td>
  <td width="80" align="center"><sub><img src="https://img.shields.io/badge/CPU-4C6EF5?style=flat-square&logo=data%3Aimage%2Fsvg%2Bxml%3Bbase64%2CPHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNCAyNCIgZmlsbD0id2hpdGUiPjxwYXRoIGQ9Ik05IDJ2Mkg3LjVBMS41IDEuNSAwIDAgMCA2IDUuNVY3SDR2MmgydjJINHYyaDJ2Mkg0djJoMnYxLjVBMS41IDEuNSAwIDAgMCA3LjUgMjBIOXYyaDJ2LTJoMnYyaDJ2LTJoMS41YTEuNSAxLjUgMCAwIDAgMS41LTEuNVYxOGgydi0yaC0ydi0yaDJ2LTJoLTJWOWgyVjdoLTJWNS41QTEuNSAxLjUgMCAwIDAgMTYuNSA0SDE1VjJoLTJ2MmgtMlYySDl6bS0xIDZoOHY4SDhWOHoiLz48L3N2Zz4%3D&logoColor=white" alt="CPU"></sub></td>
  <td width="115" align="center"><sub><span style="white-space:nowrap"><img src="https://img.shields.io/badge/HPC-2C3E50?style=flat-square&logo=linux&logoColor=white" alt="HPC ok"></span>&nbsp;<span style="white-space:nowrap"><img src="https://img.shields.io/badge/Colab-F9AB00?style=flat-square&logo=googlecolab&logoColor=white" alt="Colab ok"></span></sub></td>
</tr>
<tr>
  <td><sub><b>ExtractMetrics</b><br>Extract and reshape pipeline metrics into a Prism-ready layout.</sub></td>
  <td width="100" align="center"><sub><img src="https://img.shields.io/badge/BP-1ABC9C?style=flat-square" alt="BioPipelines native tool"></sub></td>
  <td width="80" align="center"><sub><img src="https://img.shields.io/badge/CPU-4C6EF5?style=flat-square&logo=data%3Aimage%2Fsvg%2Bxml%3Bbase64%2CPHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNCAyNCIgZmlsbD0id2hpdGUiPjxwYXRoIGQ9Ik05IDJ2Mkg3LjVBMS41IDEuNSAwIDAgMCA2IDUuNVY3SDR2MmgydjJINHYyaDJ2Mkg0djJoMnYxLjVBMS41IDEuNSAwIDAgMCA3LjUgMjBIOXYyaDJ2LTJoMnYyaDJ2LTJoMS41YTEuNSAxLjUgMCAwIDAgMS41LTEuNVYxOGgydi0yaC0ydi0yaDJ2LTJoLTJWOWgyVjdoLTJWNS41QTEuNSAxLjUgMCAwIDAgMTYuNSA0SDE1VjJoLTJ2MmgtMlYySDl6bS0xIDZoOHY4SDhWOHoiLz48L3N2Zz4%3D&logoColor=white" alt="CPU"></sub></td>
  <td width="115" align="center"><sub><span style="white-space:nowrap"><img src="https://img.shields.io/badge/HPC-2C3E50?style=flat-square&logo=linux&logoColor=white" alt="HPC ok"></span>&nbsp;<span style="white-space:nowrap"><img src="https://img.shields.io/badge/Colab-F9AB00?style=flat-square&logo=googlecolab&logoColor=white" alt="Colab ok"></span></sub></td>
</tr>
<tr>
  <td><sub><b>Panda</b><br>Filter, sort, merge, concatenate, and rank pipeline tables.</sub></td>
  <td width="100" align="center"><sub><img src="https://img.shields.io/badge/BP-1ABC9C?style=flat-square" alt="BioPipelines native tool"></sub></td>
  <td width="80" align="center"><sub><img src="https://img.shields.io/badge/CPU-4C6EF5?style=flat-square&logo=data%3Aimage%2Fsvg%2Bxml%3Bbase64%2CPHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNCAyNCIgZmlsbD0id2hpdGUiPjxwYXRoIGQ9Ik05IDJ2Mkg3LjVBMS41IDEuNSAwIDAgMCA2IDUuNVY3SDR2MmgydjJINHYyaDJ2Mkg0djJoMnYxLjVBMS41IDEuNSAwIDAgMCA3LjUgMjBIOXYyaDJ2LTJoMnYyaDJ2LTJoMS41YTEuNSAxLjUgMCAwIDAgMS41LTEuNVYxOGgydi0yaC0ydi0yaDJ2LTJoLTJWOWgyVjdoLTJWNS41QTEuNSAxLjUgMCAwIDAgMTYuNSA0SDE1VjJoLTJ2MmgtMlYySDl6bS0xIDZoOHY4SDhWOHoiLz48L3N2Zz4%3D&logoColor=white" alt="CPU"></sub></td>
  <td width="115" align="center"><sub><span style="white-space:nowrap"><img src="https://img.shields.io/badge/HPC-2C3E50?style=flat-square&logo=linux&logoColor=white" alt="HPC ok"></span>&nbsp;<span style="white-space:nowrap"><img src="https://img.shields.io/badge/Colab-F9AB00?style=flat-square&logo=googlecolab&logoColor=white" alt="Colab ok"></span></sub></td>
</tr>
<tr>
  <td><sub><b>Pool</b><br>Gather outputs of parallel runs into one standardised output.</sub></td>
  <td width="100" align="center"><sub><img src="https://img.shields.io/badge/BP-1ABC9C?style=flat-square" alt="BioPipelines native tool"></sub></td>
  <td width="80" align="center"><sub><img src="https://img.shields.io/badge/CPU-4C6EF5?style=flat-square&logo=data%3Aimage%2Fsvg%2Bxml%3Bbase64%2CPHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNCAyNCIgZmlsbD0id2hpdGUiPjxwYXRoIGQ9Ik05IDJ2Mkg3LjVBMS41IDEuNSAwIDAgMCA2IDUuNVY3SDR2MmgydjJINHYyaDJ2Mkg0djJoMnYxLjVBMS41IDEuNSAwIDAgMCA3LjUgMjBIOXYyaDJ2LTJoMnYyaDJ2LTJoMS41YTEuNSAxLjUgMCAwIDAgMS41LTEuNVYxOGgydi0yaC0ydi0yaDJ2LTJoLTJWOWgyVjdoLTJWNS41QTEuNSAxLjUgMCAwIDAgMTYuNSA0SDE1VjJoLTJ2MmgtMlYySDl6bS0xIDZoOHY4SDhWOHoiLz48L3N2Zz4%3D&logoColor=white" alt="CPU"></sub></td>
  <td width="115" align="center"><sub><span style="white-space:nowrap"><img src="https://img.shields.io/badge/HPC-2C3E50?style=flat-square&logo=linux&logoColor=white" alt="HPC ok"></span>&nbsp;<span style="white-space:nowrap"><img src="https://img.shields.io/badge/Colab-F9AB00?style=flat-square&logo=googlecolab&logoColor=white" alt="Colab ok"></span></sub></td>
</tr>
<tr>
  <td><sub><b>ReMap</b><br>Rename identifiers consistently across multiple streams and tables in a run.</sub></td>
  <td width="100" align="center"><sub><img src="https://img.shields.io/badge/BP-1ABC9C?style=flat-square" alt="BioPipelines native tool"></sub></td>
  <td width="80" align="center"><sub><img src="https://img.shields.io/badge/CPU-4C6EF5?style=flat-square&logo=data%3Aimage%2Fsvg%2Bxml%3Bbase64%2CPHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNCAyNCIgZmlsbD0id2hpdGUiPjxwYXRoIGQ9Ik05IDJ2Mkg3LjVBMS41IDEuNSAwIDAgMCA2IDUuNVY3SDR2MmgydjJINHYyaDJ2Mkg0djJoMnYxLjVBMS41IDEuNSAwIDAgMCA3LjUgMjBIOXYyaDJ2LTJoMnYyaDJ2LTJoMS41YTEuNSAxLjUgMCAwIDAgMS41LTEuNVYxOGgydi0yaC0ydi0yaDJ2LTJoLTJWOWgyVjdoLTJWNS41QTEuNSAxLjUgMCAwIDAgMTYuNSA0SDE1VjJoLTJ2MmgtMlYySDl6bS0xIDZoOHY4SDhWOHoiLz48L3N2Zz4%3D&logoColor=white" alt="CPU"></sub></td>
  <td width="115" align="center"><sub><span style="white-space:nowrap"><img src="https://img.shields.io/badge/HPC-2C3E50?style=flat-square&logo=linux&logoColor=white" alt="HPC ok"></span>&nbsp;<span style="white-space:nowrap"><img src="https://img.shields.io/badge/Colab-F9AB00?style=flat-square&logo=googlecolab&logoColor=white" alt="Colab ok"></span></sub></td>
</tr>
<tr>
  <td><sub><b>Selection</b><br>Build and manipulate residue and atom selection strings.</sub></td>
  <td width="100" align="center"><sub><img src="https://img.shields.io/badge/BP-1ABC9C?style=flat-square" alt="BioPipelines native tool"></sub></td>
  <td width="80" align="center"><sub><img src="https://img.shields.io/badge/CPU-4C6EF5?style=flat-square&logo=data%3Aimage%2Fsvg%2Bxml%3Bbase64%2CPHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNCAyNCIgZmlsbD0id2hpdGUiPjxwYXRoIGQ9Ik05IDJ2Mkg3LjVBMS41IDEuNSAwIDAgMCA2IDUuNVY3SDR2MmgydjJINHYyaDJ2Mkg0djJoMnYxLjVBMS41IDEuNSAwIDAgMCA3LjUgMjBIOXYyaDJ2LTJoMnYyaDJ2LTJoMS41YTEuNSAxLjUgMCAwIDAgMS41LTEuNVYxOGgydi0yaC0ydi0yaDJ2LTJoLTJWOWgyVjdoLTJWNS41QTEuNSAxLjUgMCAwIDAgMTYuNSA0SDE1VjJoLTJ2MmgtMlYySDl6bS0xIDZoOHY4SDhWOHoiLz48L3N2Zz4%3D&logoColor=white" alt="CPU"></sub></td>
  <td width="115" align="center"><sub><span style="white-space:nowrap"><img src="https://img.shields.io/badge/HPC-2C3E50?style=flat-square&logo=linux&logoColor=white" alt="HPC ok"></span>&nbsp;<span style="white-space:nowrap"><img src="https://img.shields.io/badge/Colab-F9AB00?style=flat-square&logo=googlecolab&logoColor=white" alt="Colab ok"></span></sub></td>
</tr>

<tr><td colspan="4"><h3>🧰 Inputs & I/O</h3></td></tr>

<tr>
  <td><sub><b>CompoundLibrary</b><br>Build a compound collection from a list, file, or SMILES set.</sub></td>
  <td width="100" align="center"><sub><img src="https://img.shields.io/badge/BP-1ABC9C?style=flat-square" alt="BioPipelines native tool"></sub></td>
  <td width="80" align="center"><sub><img src="https://img.shields.io/badge/CPU-4C6EF5?style=flat-square&logo=data%3Aimage%2Fsvg%2Bxml%3Bbase64%2CPHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNCAyNCIgZmlsbD0id2hpdGUiPjxwYXRoIGQ9Ik05IDJ2Mkg3LjVBMS41IDEuNSAwIDAgMCA2IDUuNVY3SDR2MmgydjJINHYyaDJ2Mkg0djJoMnYxLjVBMS41IDEuNSAwIDAgMCA3LjUgMjBIOXYyaDJ2LTJoMnYyaDJ2LTJoMS41YTEuNSAxLjUgMCAwIDAgMS41LTEuNVYxOGgydi0yaC0ydi0yaDJ2LTJoLTJWOWgyVjdoLTJWNS41QTEuNSAxLjUgMCAwIDAgMTYuNSA0SDE1VjJoLTJ2MmgtMlYySDl6bS0xIDZoOHY4SDhWOHoiLz48L3N2Zz4%3D&logoColor=white" alt="CPU"></sub></td>
  <td width="115" align="center"><sub><span style="white-space:nowrap"><img src="https://img.shields.io/badge/HPC-2C3E50?style=flat-square&logo=linux&logoColor=white" alt="HPC ok"></span>&nbsp;<span style="white-space:nowrap"><img src="https://img.shields.io/badge/Colab-F9AB00?style=flat-square&logo=googlecolab&logoColor=white" alt="Colab ok"></span></sub></td>
</tr>
<tr>
  <td><sub><b>Ligand</b><br>Fetch small molecules from identifiers or SMILES into a stream.</sub></td>
  <td width="100" align="center"><sub><img src="https://img.shields.io/badge/BP-1ABC9C?style=flat-square" alt="BioPipelines native tool"></sub></td>
  <td width="80" align="center"><sub><img src="https://img.shields.io/badge/CPU-4C6EF5?style=flat-square&logo=data%3Aimage%2Fsvg%2Bxml%3Bbase64%2CPHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNCAyNCIgZmlsbD0id2hpdGUiPjxwYXRoIGQ9Ik05IDJ2Mkg3LjVBMS41IDEuNSAwIDAgMCA2IDUuNVY3SDR2MmgydjJINHYyaDJ2Mkg0djJoMnYxLjVBMS41IDEuNSAwIDAgMCA3LjUgMjBIOXYyaDJ2LTJoMnYyaDJ2LTJoMS41YTEuNSAxLjUgMCAwIDAgMS41LTEuNVYxOGgydi0yaC0ydi0yaDJ2LTJoLTJWOWgyVjdoLTJWNS41QTEuNSAxLjUgMCAwIDAgMTYuNSA0SDE1VjJoLTJ2MmgtMlYySDl6bS0xIDZoOHY4SDhWOHoiLz48L3N2Zz4%3D&logoColor=white" alt="CPU"></sub></td>
  <td width="115" align="center"><sub><span style="white-space:nowrap"><img src="https://img.shields.io/badge/HPC-2C3E50?style=flat-square&logo=linux&logoColor=white" alt="HPC ok"></span>&nbsp;<span style="white-space:nowrap"><img src="https://img.shields.io/badge/Colab-F9AB00?style=flat-square&logo=googlecolab&logoColor=white" alt="Colab ok"></span></sub></td>
</tr>
<tr>
  <td><sub><b>Load / LoadMultiple</b><br>Reload previously produced pipeline outputs back into a run.</sub></td>
  <td width="100" align="center"><sub><img src="https://img.shields.io/badge/BP-1ABC9C?style=flat-square" alt="BioPipelines native tool"></sub></td>
  <td width="80" align="center"><sub><img src="https://img.shields.io/badge/CPU-4C6EF5?style=flat-square&logo=data%3Aimage%2Fsvg%2Bxml%3Bbase64%2CPHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNCAyNCIgZmlsbD0id2hpdGUiPjxwYXRoIGQ9Ik05IDJ2Mkg3LjVBMS41IDEuNSAwIDAgMCA2IDUuNVY3SDR2MmgydjJINHYyaDJ2Mkg0djJoMnYxLjVBMS41IDEuNSAwIDAgMCA3LjUgMjBIOXYyaDJ2LTJoMnYyaDJ2LTJoMS41YTEuNSAxLjUgMCAwIDAgMS41LTEuNVYxOGgydi0yaC0ydi0yaDJ2LTJoLTJWOWgyVjdoLTJWNS41QTEuNSAxLjUgMCAwIDAgMTYuNSA0SDE1VjJoLTJ2MmgtMlYySDl6bS0xIDZoOHY4SDhWOHoiLz48L3N2Zz4%3D&logoColor=white" alt="CPU"></sub></td>
  <td width="115" align="center"><sub><span style="white-space:nowrap"><img src="https://img.shields.io/badge/HPC-2C3E50?style=flat-square&logo=linux&logoColor=white" alt="HPC ok"></span>&nbsp;<span style="white-space:nowrap"><img src="https://img.shields.io/badge/Colab-F9AB00?style=flat-square&logo=googlecolab&logoColor=white" alt="Colab ok"></span></sub></td>
</tr>
<tr>
  <td><sub><b>PDB</b><br>Load protein structures from local files or by PDB code.</sub></td>
  <td width="100" align="center"><sub><img src="https://img.shields.io/badge/BP-1ABC9C?style=flat-square" alt="BioPipelines native tool"></sub></td>
  <td width="80" align="center"><sub><img src="https://img.shields.io/badge/CPU-4C6EF5?style=flat-square&logo=data%3Aimage%2Fsvg%2Bxml%3Bbase64%2CPHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNCAyNCIgZmlsbD0id2hpdGUiPjxwYXRoIGQ9Ik05IDJ2Mkg3LjVBMS41IDEuNSAwIDAgMCA2IDUuNVY3SDR2MmgydjJINHYyaDJ2Mkg0djJoMnYxLjVBMS41IDEuNSAwIDAgMCA3LjUgMjBIOXYyaDJ2LTJoMnYyaDJ2LTJoMS41YTEuNSAxLjUgMCAwIDAgMS41LTEuNVYxOGgydi0yaC0ydi0yaDJ2LTJoLTJWOWgyVjdoLTJWNS41QTEuNSAxLjUgMCAwIDAgMTYuNSA0SDE1VjJoLTJ2MmgtMlYySDl6bS0xIDZoOHY4SDhWOHoiLz48L3N2Zz4%3D&logoColor=white" alt="CPU"></sub></td>
  <td width="115" align="center"><sub><span style="white-space:nowrap"><img src="https://img.shields.io/badge/HPC-2C3E50?style=flat-square&logo=linux&logoColor=white" alt="HPC ok"></span>&nbsp;<span style="white-space:nowrap"><img src="https://img.shields.io/badge/Colab-F9AB00?style=flat-square&logo=googlecolab&logoColor=white" alt="Colab ok"></span></sub></td>
</tr>
<tr>
  <td><sub><b>Plot</b><br>Generate publication-style plots directly from your pipeline tables.</sub></td>
  <td width="100" align="center"><sub><img src="https://img.shields.io/badge/BP-1ABC9C?style=flat-square" alt="BioPipelines native tool"></sub></td>
  <td width="80" align="center"><sub><img src="https://img.shields.io/badge/CPU-4C6EF5?style=flat-square&logo=data%3Aimage%2Fsvg%2Bxml%3Bbase64%2CPHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNCAyNCIgZmlsbD0id2hpdGUiPjxwYXRoIGQ9Ik05IDJ2Mkg3LjVBMS41IDEuNSAwIDAgMCA2IDUuNVY3SDR2MmgydjJINHYyaDJ2Mkg0djJoMnYxLjVBMS41IDEuNSAwIDAgMCA3LjUgMjBIOXYyaDJ2LTJoMnYyaDJ2LTJoMS41YTEuNSAxLjUgMCAwIDAgMS41LTEuNVYxOGgydi0yaC0ydi0yaDJ2LTJoLTJWOWgyVjdoLTJWNS41QTEuNSAxLjUgMCAwIDAgMTYuNSA0SDE1VjJoLTJ2MmgtMlYySDl6bS0xIDZoOHY4SDhWOHoiLz48L3N2Zz4%3D&logoColor=white" alt="CPU"></sub></td>
  <td width="115" align="center"><sub><span style="white-space:nowrap"><img src="https://img.shields.io/badge/HPC-2C3E50?style=flat-square&logo=linux&logoColor=white" alt="HPC ok"></span>&nbsp;<span style="white-space:nowrap"><img src="https://img.shields.io/badge/Colab-F9AB00?style=flat-square&logo=googlecolab&logoColor=white" alt="Colab ok"></span></sub></td>
</tr>
<tr>
  <td><sub><b>PyMOL</b><br>Build PyMOL sessions and render publication-quality structure figures.</sub></td>
  <td width="100" align="center"><sub><a href="https://github.com/schrodinger/pymol-open-source"><img src="https://img.shields.io/badge/-181717?style=flat-square&logo=github&logoColor=white" alt="repo"></a></sub></td>
  <td width="80" align="center"><sub><img src="https://img.shields.io/badge/GPU-76B900?style=flat-square&logo=nvidia&logoColor=white" alt="GPU"></sub></td>
  <td width="115" align="center"><sub><span style="white-space:nowrap"><img src="https://img.shields.io/badge/HPC-2C3E50?style=flat-square&logo=linux&logoColor=white" alt="HPC ok"></span>&nbsp;<span style="white-space:nowrap"><img src="https://img.shields.io/badge/Colab-F9AB00?style=flat-square&logo=googlecolab&logoColor=white" alt="Colab ok"></span></sub></td>
</tr>
<tr>
  <td><sub><b>RCSB</b><br>Search the RCSB PDB and download matching protein structures.</sub></td>
  <td width="100" align="center"><sub><a href="https://doi.org/10.1093/nar/gkac1077"><img src="https://img.shields.io/badge/-paper-B31B1B?style=flat-square&logo=readthedocs&logoColor=white" alt="paper"></a></sub></td>
  <td width="80" align="center"><sub><img src="https://img.shields.io/badge/CPU-4C6EF5?style=flat-square&logo=data%3Aimage%2Fsvg%2Bxml%3Bbase64%2CPHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNCAyNCIgZmlsbD0id2hpdGUiPjxwYXRoIGQ9Ik05IDJ2Mkg3LjVBMS41IDEuNSAwIDAgMCA2IDUuNVY3SDR2MmgydjJINHYyaDJ2Mkg0djJoMnYxLjVBMS41IDEuNSAwIDAgMCA3LjUgMjBIOXYyaDJ2LTJoMnYyaDJ2LTJoMS41YTEuNSAxLjUgMCAwIDAgMS41LTEuNVYxOGgydi0yaC0ydi0yaDJ2LTJoLTJWOWgyVjdoLTJWNS41QTEuNSAxLjUgMCAwIDAgMTYuNSA0SDE1VjJoLTJ2MmgtMlYySDl6bS0xIDZoOHY4SDhWOHoiLz48L3N2Zz4%3D&logoColor=white" alt="CPU"></sub></td>
  <td width="115" align="center"><sub><span style="white-space:nowrap"><img src="https://img.shields.io/badge/HPC-2C3E50?style=flat-square&logo=linux&logoColor=white" alt="HPC ok"></span>&nbsp;<span style="white-space:nowrap"><img src="https://img.shields.io/badge/Colab-F9AB00?style=flat-square&logo=googlecolab&logoColor=white" alt="Colab ok"></span></sub></td>
</tr>
<tr>
  <td><sub><b>Scripting</b><br>Run a custom two-phase script as a typed step, with its own declared streams and tables.</sub></td>
  <td width="100" align="center"><sub><img src="https://img.shields.io/badge/BP-1ABC9C?style=flat-square" alt="BioPipelines native tool"></sub></td>
  <td width="80" align="center"><sub><img src="https://img.shields.io/badge/CPU-4C6EF5?style=flat-square&logo=data%3Aimage%2Fsvg%2Bxml%3Bbase64%2CPHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNCAyNCIgZmlsbD0id2hpdGUiPjxwYXRoIGQ9Ik05IDJ2Mkg3LjVBMS41IDEuNSAwIDAgMCA2IDUuNVY3SDR2MmgydjJINHYyaDJ2Mkg0djJoMnYxLjVBMS41IDEuNSAwIDAgMCA3LjUgMjBIOXYyaDJ2LTJoMnYyaDJ2LTJoMS41YTEuNSAxLjUgMCAwIDAgMS41LTEuNVYxOGgydi0yaC0ydi0yaDJ2LTJoLTJWOWgyVjdoLTJWNS41QTEuNSAxLjUgMCAwIDAgMTYuNSA0SDE1VjJoLTJ2MmgtMlYySDl6bS0xIDZoOHY4SDhWOHoiLz48L3N2Zz4%3D&logoColor=white" alt="CPU"></sub></td>
  <td width="115" align="center"><sub><span style="white-space:nowrap"><img src="https://img.shields.io/badge/HPC-2C3E50?style=flat-square&logo=linux&logoColor=white" alt="HPC ok"></span>&nbsp;<span style="white-space:nowrap"><img src="https://img.shields.io/badge/Colab-F9AB00?style=flat-square&logo=googlecolab&logoColor=white" alt="Colab ok"></span></sub></td>
</tr>
<tr>
  <td><sub><b>Sequence</b><br>Create a sequences stream directly from one or more strings.</sub></td>
  <td width="100" align="center"><sub><img src="https://img.shields.io/badge/BP-1ABC9C?style=flat-square" alt="BioPipelines native tool"></sub></td>
  <td width="80" align="center"><sub><img src="https://img.shields.io/badge/CPU-4C6EF5?style=flat-square&logo=data%3Aimage%2Fsvg%2Bxml%3Bbase64%2CPHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNCAyNCIgZmlsbD0id2hpdGUiPjxwYXRoIGQ9Ik05IDJ2Mkg3LjVBMS41IDEuNSAwIDAgMCA2IDUuNVY3SDR2MmgydjJINHYyaDJ2Mkg0djJoMnYxLjVBMS41IDEuNSAwIDAgMCA3LjUgMjBIOXYyaDJ2LTJoMnYyaDJ2LTJoMS41YTEuNSAxLjUgMCAwIDAgMS41LTEuNVYxOGgydi0yaC0ydi0yaDJ2LTJoLTJWOWgyVjdoLTJWNS41QTEuNSAxLjUgMCAwIDAgMTYuNSA0SDE1VjJoLTJ2MmgtMlYySDl6bS0xIDZoOHY4SDhWOHoiLz48L3N2Zz4%3D&logoColor=white" alt="CPU"></sub></td>
  <td width="115" align="center"><sub><span style="white-space:nowrap"><img src="https://img.shields.io/badge/HPC-2C3E50?style=flat-square&logo=linux&logoColor=white" alt="HPC ok"></span>&nbsp;<span style="white-space:nowrap"><img src="https://img.shields.io/badge/Colab-F9AB00?style=flat-square&logo=googlecolab&logoColor=white" alt="Colab ok"></span></sub></td>
</tr>
<tr>
  <td><sub><b>Table</b><br>Construct a table directly from inline values or a file.</sub></td>
  <td width="100" align="center"><sub><img src="https://img.shields.io/badge/BP-1ABC9C?style=flat-square" alt="BioPipelines native tool"></sub></td>
  <td width="80" align="center"><sub><img src="https://img.shields.io/badge/CPU-4C6EF5?style=flat-square&logo=data%3Aimage%2Fsvg%2Bxml%3Bbase64%2CPHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNCAyNCIgZmlsbD0id2hpdGUiPjxwYXRoIGQ9Ik05IDJ2Mkg3LjVBMS41IDEuNSAwIDAgMCA2IDUuNVY3SDR2MmgydjJINHYyaDJ2Mkg0djJoMnYxLjVBMS41IDEuNSAwIDAgMCA3LjUgMjBIOXYyaDJ2LTJoMnYyaDJ2LTJoMS41YTEuNSAxLjUgMCAwIDAgMS41LTEuNVYxOGgydi0yaC0ydi0yaDJ2LTJoLTJWOWgyVjdoLTJWNS41QTEuNSAxLjUgMCAwIDAgMTYuNSA0SDE1VjJoLTJ2MmgtMlYySDl6bS0xIDZoOHY4SDhWOHoiLz48L3N2Zz4%3D&logoColor=white" alt="CPU"></sub></td>
  <td width="115" align="center"><sub><span style="white-space:nowrap"><img src="https://img.shields.io/badge/HPC-2C3E50?style=flat-square&logo=linux&logoColor=white" alt="HPC ok"></span>&nbsp;<span style="white-space:nowrap"><img src="https://img.shields.io/badge/Colab-F9AB00?style=flat-square&logo=googlecolab&logoColor=white" alt="Colab ok"></span></sub></td>
</tr>
<tr>
  <td><sub><b>UniProt</b><br>Fetch sequences and annotations (GO, domains, organism) from UniProt.</sub></td>
  <td width="100" align="center"><sub><a href="https://www.uniprot.org/"><img src="https://img.shields.io/badge/-181717?style=flat-square&logo=github&logoColor=white" alt="repo"></a> <a href="https://doi.org/10.1093/nar/gkac1052"><img src="https://img.shields.io/badge/-paper-B31B1B?style=flat-square&logo=readthedocs&logoColor=white" alt="paper"></a></sub></td>
  <td width="80" align="center"><sub><img src="https://img.shields.io/badge/CPU-4C6EF5?style=flat-square&logo=data%3Aimage%2Fsvg%2Bxml%3Bbase64%2CPHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNCAyNCIgZmlsbD0id2hpdGUiPjxwYXRoIGQ9Ik05IDJ2Mkg3LjVBMS41IDEuNSAwIDAgMCA2IDUuNVY3SDR2MmgydjJINHYyaDJ2Mkg0djJoMnYxLjVBMS41IDEuNSAwIDAgMCA3LjUgMjBIOXYyaDJ2LTJoMnYyaDJ2LTJoMS41YTEuNSAxLjUgMCAwIDAgMS41LTEuNVYxOGgydi0yaC0ydi0yaDJ2LTJoLTJWOWgyVjdoLTJWNS41QTEuNSAxLjUgMCAwIDAgMTYuNSA0SDE1VjJoLTJ2MmgtMlYySDl6bS0xIDZoOHY4SDhWOHoiLz48L3N2Zz4%3D&logoColor=white" alt="CPU"></sub></td>
  <td width="115" align="center"><sub><span style="white-space:nowrap"><img src="https://img.shields.io/badge/HPC-2C3E50?style=flat-square&logo=linux&logoColor=white" alt="HPC ok"></span>&nbsp;<span style="white-space:nowrap"><img src="https://img.shields.io/badge/Colab-F9AB00?style=flat-square&logo=googlecolab&logoColor=white" alt="Colab ok"></span></sub></td>
</tr>

</table>

<sub><sup><i><b>a</b></i></sup> Only protein-small molecule protocol implemented <sup><i><b>b</b></i></sup> RBSDesigner2 outperforms RBSDesigner, but is not publicly available; shown for illustrative purposes. <sup><i><b>c</b></i></sup> A CPU environment (`HBDesigner.install(device="cpu")`) is provided, but only the GPU path has been tested; CPU inference is upstream-marked "not recommended" and is slow. <sup><i><b>k</b></i></sup> On HPC, runs via the apptainer container. &nbsp; <sup><i><b>p</b></i></sup> On Colab requires a high-RAM runtime: the P2Rank JVM misreads the default runtime's cgroup memory limit and crashes the kernel. &nbsp; <sup><i><b>q</b></i></sup> On Colab runs on CPU; slower for large/many inputs. &nbsp; <sup><i><b>v</b></i></sup> Generates ESM-2 3B embeddings internally (~11 GB weights, downloaded lazily on first run); needs a GPU runtime — impractically slow on CPU. &nbsp; <sup><i><b>n</b></i></sup> On Colab needs a high-RAM runtime: the openfold attention + complex sampling OOM the free T4's 12 GB system RAM; verified end-to-end on an A100 high-RAM runtime.</sub>


## Documentation

Full documentation is available at **[biopipelines.readthedocs.io](https://biopipelines.readthedocs.io/en/latest/)**.

- **[User Manual](docs/user_manual.md)**
- **[Tool Reference](docs/tool_reference.md)**
- **[Examples](example_pipelines/)**
- **[Developer Manual](docs/developer_manual.md)**