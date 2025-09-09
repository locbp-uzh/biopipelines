"""
This pipeline shows how predict the structure of the protein HaloTag7 bound to a ligand
"""

#!/usr/bin/env python3
import os, sys
sys.path.insert(0, os.getcwd()) #to see scripts in current folder

from PipelineScripts.pipeline import Pipeline
from PipelineScripts.boltz2 import Boltz2

pipeline = Pipeline(
    pipeline_name="Boltz", #Will create a folder in /shares/USER/<pipeline_name>
    job_name="HT7_Cy7_CHF2_S", #Unique job folder in /shares/USER/<pipeline_name>/job_name_NNN
    job_description="Folding of HaloTag7 with Cy7-CHF2 open for enantiomer S")

pipeline.resources(
    gpu="V100",
    time="24:00:00",
    memory="16GB"
)

"""
Importantly, when writing the ligand smiles as a string, use r prefix (raw string) to avoid misinterpreting any backslash \ character. 
We set global_msas_cache to True to save the protein msa for any later use of it – or to retrieve it if it was already computed.
We add pipeline=pipeline only so the tool can access the job name. This is not done when folding sequences from previous tools.
"""
boltz2 = pipeline.add(Boltz2(proteins="MAEIGTGFPFDPHYVEVLGERMHYVDVGPRDGTPVLFLHGNPTSSYVWRNIIPHVAPTHRCIAPDLIGMGKSDKPDLGYFFDDHVRFMDAFIEALGLEEVVLVIHDWGSALGFHWAKRNPERVKGIAFMEFIRPIPTWDEWPEFARETFQAFRTTDVGRKLIIDQNVFIEGTLPMGVVRPLTEVEMDHYREPFLNPVDREPLWRFPNELPIAGEPANIVALVEEYMDWLHQSPVPKLLFWGTPGVLIPPAEAARLAKSLPNCKAVDIGPGLNLLQEDNPDLIGSEIARWLSTLEISG",
ligands=r"C[N+](C1=C2C=CC=C1)=C([C@@]2(CC3=CN(N=N3)CCOCCOCCCCCCCl)CC(NCC(F)F)=O)/C=C/C=C/C=C/C=C(N(C4=C5C=CC=C4)C)\C5(C)C",
affinity=True,
global_msas_cache=True,
pipeline=pipeline))

pipeline.save()
pipeline.slurm(email="")