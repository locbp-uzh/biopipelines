"""
This pipeline shows how predict the structure of the protein HaloTag7 bound to several ligands
"""

#!/usr/bin/env python3
import os, sys
sys.path.insert(0, os.getcwd()) #to see scripts in current folder

from PipelineScripts.pipeline import Pipeline
from PipelineScripts.boltz2 import Boltz2
from PipelineScripts.compound_library import CompoundLibrary

pipeline = Pipeline(
    pipeline_name="Boltz", #Will create a folder in /shares/USER/<pipeline_name>
    job_name="HT7_Cy7_CHF2_RESONANCE_EFFECT", #Unique job folder in /shares/USER/<pipeline_name>/job_name_NNN
    job_description="Folding of HaloTag7 with Cy7-CHF2, several enantiomers, with the positive charge being either on close to the ring closing group (RCG+) or the other cap (CAP+)")

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
compounds = pipeline.add(CompoundLibrary(
    {
        '(S)Cy7-CH2F-RCG+':r'C[N+]1=C(/C=C/C=C/C=C/C=C(C2(C)C)/N(C)C3=C2C=CC=C3)[C@@](CC4=CN(CCOCCOCCCCCCCl)N=N4)(CC(NCCF)=O)C5=C1C=CC=C5',
        '(S)Cy7-CH2F-CAP+':r'CN(/C([C@]1(CC(NCCF)=O)CC2=CN(CCOCCOCCCCCCCl)N=N2)=C\C=C\C=C\C=C\C(C3(C)C)=[N+](C)C4=C3C=CC=C4)C5=C1C=CC=C5',
        '(R)Cy7-CH2F-RCG+':r'C[N+]1=C(/C=C/C=C/C=C/C=C(C2(C)C)/N(C)C3=C2C=CC=C3)[C@](CC4=CN(CCOCCOCCCCCCCl)N=N4)(CC(NCCF)=O)C5=C1C=CC=C5',
        '(R)Cy7-CH2F-CAP+':r'CN(/C([C@@]1(CC(NCCF)=O)CC2=CN(CCOCCOCCCCCCCl)N=N2)=C\C=C\C=C\C=C\C(C3(C)C)=[N+](C)C4=C3C=CC=C4)C5=C1C=CC=C5'
    }
))

boltz2 = pipeline.add(Boltz2(proteins="MAEIGTGFPFDPHYVEVLGERMHYVDVGPRDGTPVLFLHGNPTSSYVWRNIIPHVAPTHRCIAPDLIGMGKSDKPDLGYFFDDHVRFMDAFIEALGLEEVVLVIHDWGSALGFHWAKRNPERVKGIAFMEFIRPIPTWDEWPEFARETFQAFRTTDVGRKLIIDQNVFIEGTLPMGVVRPLTEVEMDHYREPFLNPVDREPLWRFPNELPIAGEPANIVALVEEYMDWLHQSPVPKLLFWGTPGVLIPPAEAARLAKSLPNCKAVDIGPGLNLLQEDNPDLIGSEIARWLSTLEISG",
ligands=compounds.output,
affinity=True,
global_msas_cache=True,
pipeline=pipeline))

pipeline.save()
pipeline.slurm(email="")