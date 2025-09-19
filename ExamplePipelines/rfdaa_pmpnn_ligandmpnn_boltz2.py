###IN PROGRESS
"""
This pipeline shows how to run RFdiffusion-AllAtom, generate sequences for the backbone with ProteinMPNN, add the top mutations from LigandMPNN, fold the sequences with Boltz2.
"""

#!/usr/bin/env python3
import os, sys
sys.path.insert(0, os.getcwd()) #to see scripts in current folder

from PipelineScripts.pipeline import Pipeline
from PipelineScripts.load_output import LoadOutput
from PipelineScripts.rfdiffusion_allatom import RFdiffusionAllAtom
from PipelineScripts.ligand_mpnn import LigandMPNN
from PipelineScripts.boltz2 import Boltz2
from PipelineScripts.residue_atom_distance import ResidueAtomDistance
from PipelineScripts.merge_datasheets import MergeDatasheets
from PipelineScripts.filter import Filter

pipeline = Pipeline(
    pipeline_name="RFDAA-ProteinMPNN-LigandMPNN-Boltz", #Will create a folder in /shares/USER/<pipeline_name>
    job_name="rifampicin", #Unique job folder in /shares/USER/<pipeline_name>/job_name_NNN
    job_description="redesign of N terminus of rifampicin binding protein")

pipeline.resources(
    gpu="80GB", #ask for A100-80GB or H100-80GB
    time="24:00:00",
    memory="16GB"
)

rifampicin = pipeline.add(LoadOutput("/shares/locbp.chem.uzh/public/BioPipelines/Boltz/rifampicin_001/ToolOutputs/1_Boltz2_output.json"))


rfdaa = pipeline.add(RFdiffusionAllAtom(pdb=rifampicin, #can also be a path, preferentially to PDBs folder inside biopipelines folder
                                        ligand='LIG', #in rfdaa always specify the ligand name
                                        contigs='10-20,A6-140',
                                        num_designs=2))


lmpnn = pipeline.add(LigandMPNN(structures=rfdaa, 
                                ligand="RFP", #in ligand mpnn you should always specify the ligand name. 
                                num_sequences=2, 
                                design_within=4))

boltz_apo = pipeline.add(Boltz2(proteins=lmpnn))
boltz_holo = pipeline.add(Boltz2(proteins=lmpnn,
                                ligands=rifampicin, #ligand smiles taken from original
                                msas=boltz_apo)) #MSAs are passed with <tool output>, not with <tool output>.msas

#Prints
pipeline.save()
pipeline.slurm() 
