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
    job_name="HT7_Cy7_C_RESONANCE_CANONICAL", #Unique job folder in /shares/USER/<pipeline_name>/job_name_NNN
    job_description="Folding of HaloTag7 with Cy7, several enantiomers, with the positive charge being either on close to the ring closing group (RCG+) or the other cap (CAP+)")

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
        '(S)Cy7-C-RCG+':r'CC/1(C)C2=C(C=CC=C2)N(C)\C1=C\C=C\C=C\C=C\C3=[N+](C)C4=C(C=CC=C4)[C@]3(CC5=CN(CCOCCOCCCCCCCl)N=N5)CC(=O)NC',
        '(S)Cy7-C-CAP+':r'CC1(C)C2=C(C=CC=C2)[N+](=C1/C=C/C=C/C=C/C=C\3/[C@@](CC4=CN(CCOCCOCCCCCCCl)N=N4)(CC(=O)NC)C5=C(C=CC=C5)N3C)C',
        '(R)Cy7-C-RCG+':r'CC/1(C)C2=C(C=CC=C2)N(C)\C1=C\C=C\C=C\C=C\C3=[N+](C)C4=C(C=CC=C4)[C@@]3(CC5=CN(CCOCCOCCCCCCCl)N=N5)CC(=O)NC',
        '(R)Cy7-C-CAP+':r'CC1(C)C2=C(C=CC=C2)[N+](=C1/C=C/C=C/C=C/C=C\3/[C@](CC4=CN(CCOCCOCCCCCCCl)N=N4)(CC(=O)NC)C5=C(C=CC=C5)N3C)C'
    }
))

boltz2 = pipeline.add(Boltz2(proteins="MAEIGTGFPFDPHYVEVLGERMHYVDVGPRDGTPVLFLHGNPTSSYVWRNIIPHVAPTHRCIAPDLIGMGKSDKPDLGYFFDDHVRFMDAFIEALGLEEVVLVIHDWGSALGFHWAKRNPERVKGIAFMEFIRPIPTWDEWPEFARETFQAFRTTDVGRKLIIDQNVFIEGTLPMGVVRPLTEVEMDHYREPFLNPVDREPLWRFPNELPIAGEPANIVALVEEYMDWLHQSPVPKLLFWGTPGVLIPPAEAARLAKSLPNCKAVDIGPGLNLLQEDNPDLIGSEIARWLSTLEISG",
ligands=compounds,
affinity=True,
global_msas_cache=True,
pipeline=pipeline))

pipeline.save()
pipeline.slurm()