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

from biopipelines.pipeline import *
from biopipelines.rfdiffusion import RFdiffusion
from biopipelines.rfdiffusion_allatom import RFdiffusionAllAtom
from biopipelines.rfdiffusion3 import RFdiffusion3
from biopipelines.protein_mpnn import ProteinMPNN
from biopipelines.ligand_mpnn import LigandMPNN
from biopipelines.alphafold import AlphaFold
from biopipelines.boltz2 import Boltz2
from biopipelines.boltzgen import BoltzGen
from biopipelines.mutation_profiler import MutationProfiler
from biopipelines.pymol import PyMOL

with Pipeline(project="Setup",
              job="InstallTools",
              description="Install external tools and environments"):

    Resources(time="8:00:00",
              memory="32GB")

    # Structure generation
    RFdiffusion.install()         # Creates SE3nv, clones repo, downloads weights
    RFdiffusionAllAtom.install()  # Clones repo (uses SE3nv from RFdiffusion)
    RFdiffusion3.install()        # Creates foundry env, downloads checkpoints
    BoltzGen.install()            # Creates boltzgen env

    # Sequence design
    ProteinMPNN.install()         # Clones repo (uses SE3nv from RFdiffusion)
    LigandMPNN.install()          # Creates ligandmpnn_env, clones repo

    # Structure prediction
    AlphaFold.install()           # LocalColabFold installation
    Boltz2.install()              # Creates Boltz2Env

    # Visualization & Analysis
    PyMOL.install()               # Creates ProteinEnv (PyMOL, pandas, biopython, rdkit)
    MutationProfiler.install()    # Creates MutationEnv

