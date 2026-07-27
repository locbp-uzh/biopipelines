# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

# tested:

"""
Design and fold on one exclusive node, with the MSA server warming in parallel.

    BIOPIPELINES_CONFIG_VARIANT=daint biopipelines-submit daint_packed_design.py

Daint allocates whole nodes, so a job that folds on one GPU is billed for four
and the other three idle. This packs the node instead: five concurrent tasks in
one allocation, one of them the MMseqs2 server.

The server needs ~2.6 h to make its 719 GB index resident (see llm/daint.md),
which is dead time if nothing else runs. Here it overlaps with backbone
generation and inverse folding, so by the time there are sequences to align the
server is warm. Its clients block until it advertises readiness, so no explicit
coordination is needed.

Tasks cannot consume each other's outputs — each is a separate subshell — so
this pipeline only produces per-task outputs. Gather them afterwards with a
second pipeline that Loads the run folder and Pools the results.
"""

from biopipelines.pipeline import *
from biopipelines import (
    MMseqs2Server,
    MMseqs2,
    RFdiffusion3,
    ProteinMPNN,
    Boltz2,
    PDB,
)

with Pipeline(project="Daint", job="PackedDesign",
              description="MSA server warms while designs are generated, then fold"):

    scaffold = PDB("1AKI", ids="LYZ")

    # pack=5 on a 4-GPU node: --ntasks-per-node=5 with the node's GPUs requested
    # once, so the CPU-only server does not overcommit them. cpus= splits the
    # node's 288 cores deliberately — a step with no --cpus-per-task gets ONE.
    with Parallel(pack=5):
        Resources(gpu="gh", time="12:00:00", memory="800GB")

        # The server: no GPU, most of the cores, and an idle timeout long enough
        # that it outlives the design work rather than exiting mid-pipeline.
        with Run(cpus=160, gpus=0):
            MMseqs2Server(mode="cpu", idle_timeout=86400)

        # Four design tasks, one GPU each. Tools inside a task run in sequence,
        # so each backbone set is inverse-folded and folded within its own task.
        for i in range(4):
            with Run(cpus=32, gpus=1):
                backbones = RFdiffusion3(pdb=scaffold,
                                         contigs="50-70/A81-140",
                                         num_designs=2)
                seqs = ProteinMPNN(structures=backbones, 
                                   num_sequences=2)
                msas = MMseqs2(sequences=seqs)
                Boltz2(proteins=seqs, 
                       msas=msas)
