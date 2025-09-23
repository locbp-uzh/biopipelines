"""
This pipeline shows how to run RFdiffusion-AllAtom, generate sequences from LigandMPNN, fold the sequences with Boltz2.
"""

#!/usr/bin/env python3
import os, sys
sys.path.insert(0, os.getcwd()) #to see scripts in current folder

from PipelineScripts.pipeline import Pipeline
from PipelineScripts.load_output import LoadOutput
from PipelineScripts.rfdiffusion_allatom import RFdiffusionAllAtom
from PipelineScripts.distance_selector import DistanceSelector
from PipelineScripts.mmseqs2 import MMseqs2
from PipelineScripts.ligand_mpnn import LigandMPNN
from PipelineScripts.boltz2 import Boltz2

pipeline = Pipeline(
    pipeline_name="RFDAA-LigandMPNN-MMseqs2-Boltz", #Will create a folder in /shares/USER/<pipeline_name>
    job_name="rifampicin", #Unique job folder in /shares/USER/<pipeline_name>/job_name_NNN
    job_description="redesign of N terminus of rifampicin binding protein")

pipeline.resources(
    gpu="80GB", #ask for A100-80GB or H100-80GB
    time="24:00:00",
    memory="16GB"
)

#contains both the protein-ligand complex and the ligand smiles
rifampicin = pipeline.add(LoadOutput("/shares/locbp.chem.uzh/public/BioPipelines/Boltz/rifampicin_001/ToolOutputs/1_Boltz2_output.json"))


rfdaa = pipeline.add(RFdiffusionAllAtom(pdb=rifampicin, #can also be a path, preferentially to PDBs folder inside biopipelines folder
                                        ligand='LIG', #in rfdaa always specify the ligand name
                                        contigs='10-20,A6-140',
                                        num_designs=2,
                                        steps=200))


lmpnn = pipeline.add(LigandMPNN(structures=rfdaa,
                                ligand="LIG",
                                num_sequences=2,
                                redesigned=rfdaa.datasheets.structures.designed))

msas = pipeline.add(MMseqs2(sequences=lmpnn))

boltz_holo = pipeline.add(Boltz2(proteins=lmpnn,
                                ligands=rifampicin,
                                msas=msas))

#Prints
pipeline.save()
pipeline.slurm() 
