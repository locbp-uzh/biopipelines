# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Important code guideline

**Never** implement fallback methods or values. If the intended function is not working, the code has to crash. If the intended value is not found, the code has to crash. Do not guess variables values or file paths, rather specify them only once and pass them appropriately.

## Project Overview

BioPipelines is a system that generates bash scripts for bioinformatics workflows. The pipeline itself does not execute computations directly - instead, it predicts and prepares the filesystem structure (output directories and files) that will exist at SLURM runtime, then generates bash scripts to be executed when the SLURM job runs.

### Pipeline Execution Model
- **Pipeline Role**: Generates bash scripts and predicts filesystem structure
- **Execution Time**: Scripts are executed only at SLURM runtime
- **HelpScripts**: Utility scripts designed to be executed during SLURM job execution
- **Filesystem Prediction**: Pipeline pre-determines output paths and directory structure before job submission

## Environment Setup

Each tool can depend on a specific conda environment. It is up to the user to specify one that works.

## Architecture & Workflow Structure

### Key Workflow

At "pipeline runtime" a Tool generates bash scripts and predict the filesystem structure of the output these bash scripts will produce based on its input. Tools are written agnostic of other tools, as they predict based on tool specific input parameters and output of other tools, which is given in a standard tool-aspecific format. At "slurm runtime" (on cluster supercomputer) those bash scripts are executed. The coordination of tools and managing of folders is done via a Pipeline class. Any python code that is executed at slurm time is present as a script in HelpScripts with name pipe_<purpose>.py. 

Pipeline | Slurm | Filesystem
Tool1 -> <tool1>.sh -> <job folder>/<Execution order>_<Tool name>/<output files>
Tool2 -> <tool2>.sh -> <job folder>/<Execution order>_<Tool name>/<output files>
...

#### Core Directories

#### File Organization Patterns

### Execution Environments
The biopipelines are designed for University of Zürich's cluster computing environment with:
- GPU allocation: L4, V100, A100, H100, H200 options
- SLURM job scheduling integration
- Shared storage at `/shares/locbp.chem.uzh`
- User data directories at `/home/{username}/data` and `/home/{username}`

## Common Development Commands

## Development Guidelines

### Error Handling
- All pipelines log to both console and log files using `tee`
- Failed jobs should be debugged using log files in job output directories
- Check GPU availability with `nvidia-smi` before execution

### File Management
- The system auto-generates unique names to avoid overwriting outputs
- Clean up intermediate files periodically to manage storage
- Archive completed analyses to free up working space

# Communication Guidelines

## Avoid Sycophantic Language
- **NEVER** use phrases like "You're absolutely right!", "You're absolutely correct!", "Excellent point!", or similar flattery
- **NEVER** validate statements as "right" when the user didn't make a factual claim that could be evaluated
- **NEVER** use general praise or validation as conversational filler