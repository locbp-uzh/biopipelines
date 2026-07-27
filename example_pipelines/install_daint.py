# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

# tested:

"""
Installation pipeline for the tools verified on CSCS Alps/Daint (aarch64/GH200).

Submit with the daint variant so every path, the billing account and the venv
root come from your config rather than this file:

    BIOPIPELINES_CONFIG_VARIANT=daint biopipelines-submit install_daint.py

Environments land under machine.env_manager.venv_root. Where that points is a
config decision — on scratch it survives until the purge window and needs
re-running afterwards; on project storage it persists but costs inodes.

Daint runs a subset of the catalogue. Tools are absent here when they cannot be
built for aarch64 at all: RFdiffusion needs DGL, whose aarch64 wheel ships no
C++ library matching the installed torch, and the RosettaCommons images are
x86-64 only. See llm/daint.md for the evidence behind each.
"""

from biopipelines.pipeline import *
from biopipelines import (
    ADMETAI,
    AlphaFold,
    Boltz2,
    BoltzGen,
    DSSP,
    ESMFold2,
    LigandMPNN,
    MMseqs2Server,
    PoseBusters,
    ProteinMPNN,
    PyMOL,
    RFdiffusion3,
)

with Pipeline(project="Setup", job="InstallDaint",
              description="Install the Daint-verified tool set"):
    Resources(gpu="gh", time="8:00:00", memory="32GB")

    # Structure generation and inverse folding
    RFdiffusion3.install()   # foundry: rc-foundry[all] over the NGC image
    ProteinMPNN.install()    # torch + numpy only; no container needed
    LigandMPNN.install()     # compiles ProDy, patches the vendored openfold

    # Folding and co-folding
    Boltz2.install()         # boltz[cuda] layered on NGC torch
    BoltzGen.install()
    ESMFold2.install()
    AlphaFold.install()      # LocalColabFold; also what MMseqs2Server calls

    # Analysis
    PyMOL.install()          # ProteinEnv: also ConformationalChange, Contacts, SASA
    DSSP.install()
    PoseBusters.install()
    ADMETAI.install()        # CPU-only; the install downloads the bundled weights

    # ColabFold MSA databases: ~216 GB of downloads, no GPU and no compute, so
    # they go to the transfer partition rather than holding a GH200 node idle
    # for the ~80 minutes this takes. A new Resources() opens a second batch.
    Resources(gpu="none", partition="xfer", time="24:00:00", memory="32GB")
    MMseqs2Server.install(step="databases")

    # Index build is deliberately absent: it is CPU-heavy and its cost has not
    # been measured yet. Run MMseqs2Server.install(step="build", mode="cpu")
    # on a compute node once the downloads are in.
