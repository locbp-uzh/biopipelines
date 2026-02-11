# Copyright (c) 2026 Gianluca Quargnali @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Installation pipeline for BioPipelines tools.

Each tool.install() call adds an installation step to the pipeline.
When submitted, each step runs the tool's installation script (git clone,
environment creation, model downloads, etc.) as a SLURM job.

Uncomment/comment tools as needed for your setup.
"""

from PipelineScripts.pipeline import *
from PipelineScripts.rfdiffusion import RFdiffusion
from PipelineScripts.rfdiffusion_allatom import RFdiffusionAllAtom
from PipelineScripts.rfdiffusion3 import RFdiffusion3
from PipelineScripts.protein_mpnn import ProteinMPNN
from PipelineScripts.ligand_mpnn import LigandMPNN
from PipelineScripts.alphafold import AlphaFold
from PipelineScripts.boltz2 import Boltz2
from PipelineScripts.boltzgen import BoltzGen
from PipelineScripts.mutation_profiler import MutationProfiler

with Pipeline(project="Setup",
              job="InstallTools",
              description="Install external tools and environments"):

    Resources(time="8:00:00",
              memory="32GB")

    # Structure generation
    RFdiffusion.install()         # Creates ProteinEnv, clones repo, downloads weights
    RFdiffusionAllAtom.install()  # Clones repo (uses ProteinEnv from RFdiffusion)
    RFdiffusion3.install()        # Creates foundry env, downloads checkpoints
    BoltzGen.install()            # Creates boltzgen env

    # Sequence design
    ProteinMPNN.install()         # Clones repo (uses ProteinEnv from RFdiffusion)
    LigandMPNN.install()          # Creates ligandmpnn_env, clones repo

    # Structure prediction
    AlphaFold.install()           # Downloads LocalColabFold installer
    Boltz2.install()              # Creates Boltz2Env

    # Analysis
    MutationProfiler.install()    # Creates MutationEnv
