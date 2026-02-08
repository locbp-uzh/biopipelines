"""
This pipeline shows how to use boltz2
"""

from PipelineScripts.pipeline import *
from PipelineScripts.entities import * # import PDB, Sequence, Ligand
from PipelineScripts.boltz2 import Boltz2

"""
when running via submit, all the sbatch commands generated in the python script will be submitted
./submit /ExamplePipelines/multiple_pipelines_submission.py
will result in 3 jobs
"""

Linkers = {
    "original": r"ClCCCCCCOCCOCC",
    "peg2": r"ClCCOCCOCC",
    "peg3": r"ClCCOCCOCCOCC"
}

for linker in Linkers.keys():
    with Pipeline(project="Examples",
                  job=f"HT_Linker_{linker}",
                  description="Folding of HaloTag7 with different chloroalkane linkers"):
        Resources(gpu="any",
                  time="24:00:00",
                  memory="16GB")
        HaloTag = PDB("6U32",ids="HT")
        Boltz2(proteins=HaloTag,
        ligands=Ligand(smiles=Linkers[linker])) # we never pass smiles directly

        
