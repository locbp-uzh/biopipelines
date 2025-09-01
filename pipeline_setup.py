#!/usr/bin/env python3
import os, sys
sys.path.insert(0, os.getcwd()) #to see scripts in current folder

from PipelineScripts.pipeline import Pipeline
from PipelineScripts.ligand_mpnn import LigandMPNN
from PipelineScripts.boltz2 import Boltz2
from PipelineScripts.filter import Filter
from PipelineScripts.residue_atom_distance import ResidueAtomDistance
from PipelineScripts.confidence import Confidence

pipeline = Pipeline(
    pipeline_name="LigandMPNN-Boltz", #Will create a folder in /shares/USER/...
    job_name="HT7_Cy7_ChlorineFilter", #Unique job folder in /shares/USER/.../job_name_NNN
    job_description="Test on filter based on distance between chlorine and aspartate")

pipeline.resources(
    gpu="V100",
    time="24:00:00",
    memory="16GB"
)

lmpnn = pipeline.add(
    LigandMPNN(structures="HT7_Cy7_CHF2_S_noncov.pdb",
        num_sequences=3, 
        design="145-180", 
        design_within=4),
    env="ProteinEnv")
print("="*30+"LigandMPNN"+"="*30)
print(lmpnn.output)

boltz = pipeline.add(
    Boltz2(
        proteins=lmpnn.output,
        ligand="CN(C(/C=C/C=C/C=C/C=C(C1(C)C)/N(C)C2=C1C=CC=C2)(N3CC(F)F)C4(CC3=O)CC5=CN(CCOCCOCCCCCCCl)N=N5)C6=C4C=CC=C6",
        affinity=True
    )
)
print("="*30+"Boltz2"+"="*30)
print(boltz.output)

filter = pipeline.add(
    Filter(
        criteria=[
            ResidueAtomDistance(
                atom = "ligand.Cl",
                residue = "protein.D in IGDWG",
                expression = "distance<=5"
            ),
            Confidence(
                expression="pLDDT>80"
            )
        ]
    )
)
print("="*30+"Filter"+"="*30)
print(filter.output)

pipeline.save()
pipeline.slurm() #prints the instructions for slurm submission


