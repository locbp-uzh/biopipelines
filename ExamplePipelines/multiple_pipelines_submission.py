"""
This pipeline shows how to use boltz2
"""

import os, sys
sys.path.insert(0, os.getcwd()) #to see scripts in current folder

from PipelineScripts.pipeline import *
from PipelineScripts.boltz2 import Boltz2
from PipelineScripts.entities import * # import PDB, Sequence, Ligand

"""
when running via submit, all the sbatch commands generated in the python script will be submitted
./submit /ExamplePipelines/multiple_pipelines_submission.py
will result in 4 jobs
"""

Cy5s = {
    "Cy5_R":r"CC/1(C)C2=C(C=CC=C2)N(C)\C1=C\C=C\C=C\C3=[N+](C)C4=C(C=CC=C4)[C@@]3(CC5=CN(CCOCCOCCCCCCCl)N=N5)CC(=O)NC",
    "Cy5_S":r"CC/1(C)C2=C(C=CC=C2)N(C)\C1=C\C=C\C=C\C3=[N+](C)C4=C(C=CC=C4)[C@]3(CC5=CN(CCOCCOCCCCCCCl)N=N5)CC(=O)NC",
    "Cy5_RR":r"CC/1(C)C2=C(C=CC=C2)N(C)\C1=C\C=C\C=C\[C@]34[C@](CC5=CN(CCOCCOCCCCCCCl)N=N5)(CC(=O)N3C)C6=C(C=CC=C6)N4C",
    "Cy5_SS":r"CC/1(C)C2=C(C=CC=C2)N(C)\C1=C\C=C\C=C\[C@@]34[C@@](CC5=CN(CCOCCOCCCCCCCl)N=N5)(CC(=O)N3C)C6=C(C=CC=C6)N4C"
}

for cy5 in Cy5s.keys():
    with Pipeline(project="Examples",
                  job=f"HT_{cy5}",
                  description="Folding of HaloTag7 with Cy5 methyl amide, several enantiomers"):
        Resources(gpu="any",
                  time="24:00:00",
                  memory="16GB")
        HaloTag = PDB("6U32",ids="HT")
        Boltz2(proteins=HaloTag,
        ligands=Ligand(smiles=Cy5s[cy5])) # we never pass smiles directly

        
