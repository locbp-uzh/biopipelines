"""
This pipeline shows how predict the structure of the protein HaloTag7 bound to a ligand
"""

#!/usr/bin/env python3
import os, sys
sys.path.insert(0, os.getcwd()) #to see scripts in current folder

from PipelineScripts.pipeline import Pipeline
from PipelineScripts.load_output import LoadOutput
from PipelineScripts.boltz2 import Boltz2
from PipelineScripts.site_directed_mutagenesis import SDM

pp = Pipeline(
    pipeline_name="SaturationMutagenesis", #Will create a folder in /shares/USER/<pipeline_name>
    job_name="HT7_Cy7_C_R", #Unique job folder in /shares/USER/<pipeline_name>/job_name_NNN
    job_description="Saturation mutagenesis of HaloTag7 at position 175 with Cy7 methyl amide open enantiomer R")

pp.resources(
    gpu="32GB",
    time="24:00:00",
    memory="16GB"
)

#we load this
best_R = pp.add(LoadOutput('/shares/locbp.chem.uzh/public/BioPipelines/Boltz/HT7_Cy7_C_R_001/ToolOutputs/1_Boltz2_output.json'))


single_point_mutants = pp.add(SDM(original=best_R,
                                  position=175,
                                  mode="saturation"))

boltz2 = pp.add(Boltz2(proteins=single_point_mutants,
                        ligands="CCC(=O)N(CC1=C2C=CC=C(CC(=O)N(O)C(=O)N[C@@H]3CC4=CC=C(C=C4)C5=CC=C(C[C@H]6CC[C@H](O)CC6)C=C5)C=C2S1)N1C(=O)CCC1=O"))

pp.save()
pp.slurm(email="")