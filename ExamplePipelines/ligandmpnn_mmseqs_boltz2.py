"""
This pipeline shows how to improve the difference in predicted binding affinity between open and close form of a carbocyanine 7 chloride and halotag7 starting from a Boltz model of the open form.
"""

import os, sys
sys.path.insert(0, os.getcwd()) #to see scripts in current folder

from PipelineScripts.pipeline import *
from PipelineScripts.load_output import LoadOutput
from PipelineScripts.mmseqs2 import MMseqs2
from PipelineScripts.ligand_mpnn import LigandMPNN
from PipelineScripts.boltz2 import Boltz2
from PipelineScripts.residue_atom_distance import ResidueAtomDistance
from PipelineScripts.merge_tables import MergeTables
from PipelineScripts.filter import Filter

with Pipeline(project="Examples",
              job="LigandMPNN-MSA-Boltz",
              description="Test of MMseqs2 local"):

    Resources(gpu="V100",
              time="24:00:00",
              memory="16GB")

    original = LoadOutput('/shares/locbp.chem.uzh/public/BioPipelines/Boltz/HT7_Cy7_C_R_001/ToolOutputs/1_Boltz2_output.json')

    lmpnn = LigandMPNN(structures=original, #this is equivalent to boltz2
                       ligand="LIG", #in ligand mpnn you should always specify the ligand name, which is LIG if from Boltz
                       num_sequences=3)

    msas = MMseqs2(lmpnn)

    boltz_holo_open = Boltz2(proteins=lmpnn,
                             ligands=original,
                             msas=msas) 
