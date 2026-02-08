"""
De novo antibiotic binding protein design using BoltzGen.
"""

from PipelineScripts.pipeline import *
from PipelineScripts.ligand import Ligand
from PipelineScripts.boltzgen import BoltzGen


#We run in parallel 10x 1000 designs -> 10 000 designs. Recommended is 10 to 20 thousands.
for _ in range(10):
    with Pipeline(project="Examples",
                job="Gentamicin_BoltzGen_1000designs",
                description="BoltzGen-based de novo protein binder design against gentamicin - analysis and filtering"):
        # GPU-heavy steps: design through affinity
        Resources(gpu="80GB|96GB", time="24:00:00", memory="16GB")
        designs = BoltzGen(
            ligand=Ligand(lookup="gentamicin", ids="gentamicin", codes="LIG"), # here we give a pdb in input
            binder_spec="140-180", # Preprint uses this range of lengths
            protocol="protein-small_molecule",
            num_designs=1000,
            steps=["design", "inverse_folding", "folding", "design_folding", "affinity"]) # We stop at affinity as subsequent steps require more memory and less GPU

# By running the pipeline with ./submit /path/to/<pipeline>.py, 10 jobs will be submitted.
# Note down the job ids e.g. 517802, 517803, ..., 517811, 517812
# These can be used to make sure part 2 only starts after all jobs for part 1 are completed