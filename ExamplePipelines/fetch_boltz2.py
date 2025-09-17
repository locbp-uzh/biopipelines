"""
This pipeline shows how predict the structure of the protein HaloTag7 bound to a ligand
"""

#!/usr/bin/env python3
import os, sys
sys.path.insert(0, os.getcwd()) #to see scripts in current folder

from PipelineScripts.pipeline import Pipeline
from PipelineScripts.fetch_structure import FetchStructure
from PipelineScripts.boltz2 import Boltz2

pipeline = Pipeline(
    pipeline_name="Boltz", #Will create a folder in /shares/USER/<pipeline_name>
    job_name="HT7_Cy7_C_S", #Unique job folder in /shares/USER/<pipeline_name>/job_name_NNN
    job_description="Folding of HaloTag7 with Cy7 methyl amide open enantiomer S")

pipeline.resources(
    gpu="V100",
    time="24:00:00",
    memory="16GB"
)

HaloTag = pipeline.add(FetchStructure("6U32","HT"))
"""
Importantly, when writing the ligand smiles as a string, use r prefix (raw string) to avoid misinterpreting any backslash \ character. 
We set global_msas_cache to True to save the protein msa for any later use of it – or to retrieve it if it was already computed.
We add pipeline=pipeline only so the tool can access the job name. This is not done when folding sequences from previous tools.
"""
boltz2 = pipeline.add(Boltz2(proteins=HaloTag.output,
ligands=r"CC/1(C)C2=C(C=CC=C2)N(C)\C1=C\C=C\C=C\C=C\C3=[N+](C)C4=C(C=CC=C4)[C@]3(CC5=CN(CCOCCOCCCCCCCl)N=N5)CC(=O)NC",
affinity=True,
global_msas_cache=True,
pipeline=pipeline))

pipeline.save()
pipeline.slurm(email="")