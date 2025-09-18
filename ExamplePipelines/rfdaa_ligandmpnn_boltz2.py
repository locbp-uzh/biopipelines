"""
This pipeline shows how to run RFdiffusion-AllAtom, generate sequences for the designed part with LigandMPNN, and fold with Boltz2
"""

#!/usr/bin/env python3
import os, sys
sys.path.insert(0, os.getcwd()) #to see scripts in current folder

from PipelineScripts.pipeline import Pipeline
from PipelineScripts.rfdiffusion_allatom import RFdiffusionAllAtom
from PipelineScripts.ligand_mpnn import LigandMPNN
from PipelineScripts.boltz2 import Boltz2
from PipelineScripts.residue_atom_distance import ResidueAtomDistance
from PipelineScripts.merge_datasheets import MergeDatasheets
from PipelineScripts.filter import Filter

pipeline = Pipeline(
    pipeline_name="RFDAA-LigandMPNN-Boltz", #Will create a folder in /shares/USER/<pipeline_name>
    job_name="RifampicinN", #Unique job folder in /shares/USER/<pipeline_name>/job_name_NNN
    job_description="N terminus redesign of rifampicin binder")

pipeline.resources(
    gpu="V100",
    time="24:00:00",
    memory="16GB"
)

"""
Design a couple of backbones with RFD-AA. The output structures will still contain the ligand with the same name as in the original pdb.
"""
rfdaa = pipeline.add(RFdiffusionAllAtom(
    pdb='rifampicin_2hw2_clean.pdb', #copied within the folder PDBs in the biopipelines folder. 
    ligand='RFP', #in rfdaa always specify the ligand name
    contigs='10-20,A6-140',
    num_designs=2
))

"""
Diversify with LigandMPNN
"""
lmpnn = pipeline.add(LigandMPNN(structures=rfdaa, 
                                ligand="RFP", #in ligand mpnn you should always specify the ligand name. 
                                num_sequences=2, 
                                design_within=4))

"""
We run the Apo version first. One can also extract confidence parameters from it, and in general here is where we calculate the MSAs, which will be recycled later on with the msas input parameter.
"""
boltz_apo = pipeline.add(Boltz2(proteins=lmpnn))

"""
MSAs are passed with <tool>, not with <tool>.msas
"""
boltz_holo = pipeline.add(Boltz2(proteins=lmpnn,
                                ligands=r"C[C@H]1/C=C/C=C(C)\C(NC2=C(/C=N/N3CCN(C)CC3)C(O)=C4C(C(O)=C(C)C5=C4C([C@](O/C=C/[C@H](OC)[C@@H](C)[C@@H](OC(C)=O)[C@H](C)[C@H](O)[C@H](C)[C@H]1O)(C)O5)=O)=C2O)=O",
                                msas=boltz_apo,
                                affinity=True))

#Prints
pipeline.save()
pipeline.slurm(email="") 
