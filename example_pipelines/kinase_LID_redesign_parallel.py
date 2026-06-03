# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

# tested:

"""
Parallelized version of kinase_LID_redesign.

Splits the 10 RFdiffusion designs across N parallel SLURM jobs (one per
GPU), pools the resulting backbones into a single DataStream, and runs
the downstream sequence-design / refold / filter chain once on the
gathered set. Compared with the single-job version the wall-clock for
the diffusion step drops by ~N (each branch generates num_designs/N
structures concurrently), at the cost of N independent sbatch
allocations instead of one.
"""

from biopipelines.pipeline import *
from biopipelines import (
    RFdiffusion,
    ProteinMPNN,
    AlphaFold,
    ConformationalChange,
    Panda,
    Pool,
    PyMOL,
)

N_BRANCHES = 5            # number of parallel RFdiffusion jobs
DESIGNS_PER_BRANCH = 2    # → 10 backbones total, same as the serial version

with Pipeline(project="AdenylateKinase", job="LID_Redesign_Parallel"):
    Resources(gpu="A100", time="0:30:00", memory="8GB")
    kinase = PDB("4AKE")

    # Fan out: N independent RFdiffusion runs, each on its own A100.
    # Resources() inside the block opens a sibling SLURM job per branch.
    backbones_runs = []
    with Parallel():
        for i in range(N_BRANCHES):
            Resources(gpu="A100", time="4:00:00", memory="16GB")
            backbones_runs.append(RFdiffusion(pdb=kinase,
                                  contigs='A1-117/50-70/A161-214',
                                  num_designs=DESIGNS_PER_BRANCH))

    # Fan in: the next Resources() call waits for every sibling to finish,
    # and Pool gathers the per-branch structures into one DataStream with
    # collision-free ids.
    Resources(gpu="A100", time="4:00:00", memory="16GB")
    backbones = Pool(runs=backbones_runs,
                     recount_prefix="kinase_design")

    sequences = ProteinMPNN(structures=backbones,
                            num_sequences=2,
                            redesigned=backbones.tables.structures.designed)
    refolded = AlphaFold(proteins=sequences)
    conf_change = ConformationalChange(reference_structures=kinase,
                                       target_structures=refolded,
                                       selection=backbones.tables.structures.fixed)
    top3 = Panda(tables=[refolded.tables.confidence,
                         conf_change.tables.changes],
                 operations=[Panda.merge(),
                             Panda.filter("RMSD < 1.5"),
                             Panda.sort("plddt", ascending=False),
                             Panda.head(3)],
                 pool=refolded)
    PyMOL(PyMOL.Load(top3),
          PyMOL.ColorAF(),
          PyMOL.Color("white", selection=backbones.tables.structures.fixed))
