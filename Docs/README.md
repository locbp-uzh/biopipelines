# ProteinNotebooks Documentation

A comprehensive guide to computational protein design workflows using integrated modeling tools.

## Quick Navigation

- [User Manual](user_manual.md) - For general usage
- [Developer Guide](developer_guide.md) - For coding and pipeline development  
- [Tool I/O Reference](tool_io_reference.md) - Detailed input/output formats
- [Pipeline Architecture](pipeline_architecture.md) - New modular system design

## Overview

This repository contains Jupyter notebooks and Python scripts for automated protein modeling workflows using:

- **RFdiffusion/RFdiffusion-AllAtom**: Structure generation and design
- **ProteinMPNN/LigandMPNN**: Sequence design and optimization
- **AlphaFold2/ColabFold**: Structure prediction and validation
- **Boltz1/Boltz2**: Advanced structure prediction with ligand support
- **OmegaFold**: Alternative folding method

## Architecture

The system is transitioning from Jupyter notebooks to a modular Python pipeline architecture:

```
ProteinNotebooks/
├── ConfigScripts/     # Tool configuration modules
├── IOScripts/         # Input/output conversion utilities  
├── Utilities/         # Core pipeline management
├── Docs/              # Documentation (this folder)
├── HelpScripts/       # Legacy utility scripts
└── [notebooks]/       # Legacy Jupyter notebooks
```

## Getting Started

1. **For general usage**: Start with the [User Manual](user_manual.md)
2. **For development**: See the [Developer Guide](developer_guide.md)
3. **For pipeline creation**: Check [Pipeline Architecture](pipeline_architecture.md)

## Environment Requirements

- **ProteinEnv**: Python 3.9 (RFdiffusion, ProteinMPNN, AlphaFold2, OmegaFold)
- **Boltz2Env**: Python 3.11 (Boltz2 predictions)
- **ligandmpnn_env**: Python 3.11 (LigandMPNN, incompatible with Boltz2Env)
- **GPU**: T4/V100/A100 depending on workload complexity