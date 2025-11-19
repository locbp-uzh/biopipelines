# BioPipelines

## Overview

A Python framework for automated computational protein design workflows that generates SLURM-executable bash scripts for high-performance computing clusters. BioPipelines creates modular, reproducible workflows by connecting bioinformatics tools through standardized interfaces. 

## Installation

```bash
git clone https://gitlab.uzh.ch/locbp/public/biopipelines
```
Specific environments for using the tools have to be installed separately. Good luck.

## More info

Check out the following resources for a detailed explanation of it.

- **[User Manual](UserManual.md)** - Complete usage guide with tool reference
- **[Examples](ExamplePipelines/)** - Pipeline examples
- **[Developer Manual](DeveloperManual.md)** - Complete usage guide with tool reference

## Key Features

- **SLURM Integration**: Automatic job submission and resource management
- **Modular Design**: Mix and match tools for custom workflows  
- **Standardized Interfaces**: CSV tables for seamless tool chaining
- **Environment Management**: Automatic conda environment switching
- **Filesystem Prediction**: Pre-determines output structure before execution
- **Cluster Optimization**: Designed for university HPC environments


