# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

# tested:


"""
De novo antibiotic binding protein design using BoltzGen.

Runs 10x1000 designs in parallel SLURM batches via ``Parallel()``, then
merges and filters the combined output in a fan-in batch. The whole
sweep (design + analysis + filtering) is one pipeline submission, so
no manual job-id bookkeeping is needed.
"""

from biopipelines.pipeline import *
from biopipelines.ligand import Ligand
from biopipelines.boltzgen import BoltzGen, BoltzGenMerge

NUM_BATCHES = 4
NUM_DESIGNS_PER_BATCH = 10
NUM_FINAL = 10

METRICS_OVERRIDE = {
    "design_iiptm": 1.1,
    "design_ptm": 1.1,
    "neg_min_design_to_target_pae": 1.1,
    "plip_hbonds_refolded": 2.0,
    "plip_saltbridge_refolded": 2.0,
    "delta_sasa_refolded": 2.0,
}

with Pipeline(project="Examples",
              job="Dopamine_BoltzGen_10x1000designs",
              description="BoltzGen de novo binder design against dopamine — parallel sweep + merge + filter"):
    # Anchor batch: fetch the ligand once on a cheap CPU job; the parallel
    # BoltzGen siblings below all reuse this single output.
    Resources(time="1:00:00", memory="4GB")
    ligand = Ligand("dopamine")

    # 10 parallel BoltzGen runs of 1000 designs each. Recommended is 10-20k total
    runs = []
    with Parallel():
        for batch in range(NUM_BATCHES):
            Resources(gpu="80GB|96GB", time="24:00:00", memory="16GB")
            Suffix(f"batch{batch+1:02d}") # Important for followup merge, to have different IDs per batch. BoltzGen assigns IDs based on the folder name.
            runs.append(BoltzGen(ligand=ligand,
                                 binder_spec="140-180",
                                 protocol="protein-small_molecule",
                                 num_designs=NUM_DESIGNS_PER_BATCH,
                                 steps=["design", 
                                        "inverse_folding", 
                                        "folding",
                                        "design_folding", 
                                        "affinity"]))

    # Run analysis, merge, then filter. 
    Resources(gpu="any", time="24:00:00", memory="64GB")
    Suffix()
    analyzed = []
    for run in runs:
        analyzed.append(BoltzGen(reuse=run,
                                 ligand=ligand,
                                 binder_spec="140-180",
                                 protocol="protein-small_molecule",
                                 num_designs=NUM_DESIGNS_PER_BATCH,
                                 steps=["analysis"]))
        
    merged = BoltzGenMerge(sources=analyzed)

    BoltzGen(reuse=merged,
             steps=["filtering"],
             binder_spec="140-180",
             protocol="protein-small_molecule",
             num_designs=NUM_BATCHES*NUM_DESIGNS_PER_BATCH,
             budget=NUM_FINAL,       
             refolding_rmsd_threshold=2.5,
             metrics_override=METRICS_OVERRIDE)
