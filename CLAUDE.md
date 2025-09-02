# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Important code guideline

**Never** implement fallback methods. If the intended function is not working, the code has to crash.

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

### Key Components

#### Core Directories

#### File Organization Patterns

### Execution Environments
The notebooks are designed for University of Zürich's cluster computing environment with:
- GPU allocation: T4, V100, A100 options
- SLURM job scheduling integration
- Shared storage at `/shares/locbp.chem.uzh`
- User data directories at `/data/{username}` and `/home/{username}`

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