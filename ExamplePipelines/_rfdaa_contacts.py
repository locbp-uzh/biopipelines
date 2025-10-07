###
###
###
###
###
###
####IN PROGRESS
###
###
###
###
###
###
"""
This pipeline shows how to run RFdiffusionAllAtom to redesign the N-terminus of a rifampicin-binding protein and then filter structures in based on the proximity of the newly generated structure towards the ligand.
"""

import os, sys
sys.path.insert(0, os.getcwd()) #to see scripts in current folder

from PipelineScripts.pipeline import Pipeline
from PipelineScripts.load_output import LoadOutput
from PipelineScripts.rfdiffusion_allatom import RFdiffusionAllAtom
from PipelineScripts.protein_ligand_contacts import ProteinLigandContacts

pipeline = Pipeline(
    pipeline_name="Examples",
    job_name="RFDAA-Contacts",
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
                                        num_designs=2,
                                        steps=20))

contacts = pipeline.add(ProteinLigandContacts(structures=rfdaa,
                                              selections=rfdaa.datasheets.structures.designed,
                                              ligand="LIG"))


pipeline.save()
pipeline.slurm() 
