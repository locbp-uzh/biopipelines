# BioPipelines

A Python framework for automated computational protein design workflows on HPC clusters (SLURM, LSF, PBS) and Jupyter notebooks.

---

## What is BioPipelines?

BioPipelines provides standardized interfaces to connect bioinformatics tools into reproducible workflows. It does not execute computations directly -- instead, it generates bash scripts and predicts output file paths, which are then executed on HPC clusters (SLURM, LSF, PBS) or interactively in notebooks.

<figure markdown>
  ![Ubiquitin inverse folding pipeline](images/figure1_ubiquitin_pipeline.png){ width="800" }
</figure>

<figure markdown>
  ![Boltz2 compound library screening](images/figure3_boltz2_compound_screening.png){ width="800" }
</figure>

---

## Key Features

**70+ integrated tools** -- Structure generation (RFdiffusion, BoltzGen), sequence design (ProteinMPNN, LigandMPNN), structure prediction (AlphaFold, Boltz2), analysis, and more.

**Three ways to run** -- Let an AI coding assistant author and run pipelines for you, submit to HPC clusters (SLURM, LSF, PBS) with `biopipelines-submit`, or run interactively in Jupyter/Colab notebooks with on-the-fly execution.

**Combinatorics** -- Cartesian products (`Each`) and grouping (`Bundle`) to systematically explore protein-ligand combinations.

**Data management** -- DataStreams for file tracking, Tables for metrics, and Panda for pandas-style transformations (filter, sort, merge, concat).

**Visualization** -- Automated PyMOL sessions and plots.

---

## Quick Start

=== "AI coding assistant"

    The easiest way to use BioPipelines, especially if you don't write Python. Clone the repo on your computer and run an AI coding assistant (Claude Code, Codex, …) inside it; describe your protocol in plain language and let it author and run the pipeline. We recommend forking the repository prior to start working.

    ```bash
    #git clone https://github.com/<your-project>/biopipelines
    git clone https://github.com/locbp-uzh/biopipelines
    cd biopipelines
    claude   # or: codex
    ```

    Then open the session by pointing the assistant at the pipeline prompt:

    > Read and follow `llm/pipelines.md`. <then describe the protocol you want>

    The assistant reads the framework's contract from `llm/`, interviews you about any open choices, and writes and runs the pipeline. Four contracts available: `pipelines.md` (workflows), `development.md` (tool implementation), `cluster.md` (automated debugging on HPCs), `colab.md` (automated debugging on Google Colab).

=== "Cluster (conda/mamba)"

    ```bash
    git clone https://github.com/locbp-uzh/biopipelines
    cd biopipelines
    mamba env create -f environments/biopipelines.yaml
    mamba activate biopipelines
    pip install -e .
    ipython kernel install --user --name biopipelines
    ```

    Some clusters are configured to give low memory to the default bash shell, which might result in failure of the above procedure (std_alloc). You can avoid this by running the following prior to the installation:

    ```bash
    srun --mem=16GB --time=1:00:00 --pty bash
    ```

    Edit `config.yaml` to match your cluster configuration.

    Individual models have to be installed separately. We provide a pipeline (example_pipelines/install_tools.py) to install all the tools used in the repository at once, but please refer to the respective official documentation in case your particular cluster configuration requires adjustments:

    ```bash
    cd example_pipelines
    biopipelines-submit install_tools.py
    ```

    ### Run an Example

    ```bash
    biopipelines-submit ubiquitin.py
    ```

=== "Google Colab"

    Run these two cells at the top of your Colab notebook:

    ```python
    # Cell 1: Install BioPipelines and micromamba
    !git clone https://github.com/locbp-uzh/biopipelines
    %cd biopipelines
    !pip install -e ".[all]"
    !wget -q https://github.com/mamba-org/micromamba-releases/releases/latest/download/micromamba-linux-64 -O /usr/local/bin/micromamba && chmod +x /usr/local/bin/micromamba
    !micromamba create -f environments/biopipelines.yaml -y
    ```

    ```python
    # Cell 2: Install the tools you need (one-time per session)
    from biopipelines.pipeline import *
    from biopipelines.rfdiffusion import RFdiffusion

    with Pipeline("Setup", "install", description="Install tools"):
        RFdiffusion.install()
    ```

    The Colab configuration (`colab.yaml`) is detected automatically — no manual config needed. Tools are installed via `micromamba` into isolated environments, matching the cluster behavior. See [Google Colab](user_manual.md#installation-google-colab) in the User Guide for details.

---

## Documentation

- **[User Guide](user_manual.md)** -- Core concepts, installation, and usage
- **[Tool Reference](tool_reference.md)** -- Complete reference for all tools
- **[Developer Guide](developer_manual.md)** -- Architecture and tool development
