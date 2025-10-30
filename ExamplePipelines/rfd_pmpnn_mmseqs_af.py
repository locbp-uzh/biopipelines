"""
This pipeline shows how to run RFdiffusion-AllAtom, generate sequences far from the ligand with ProteinMPNN and close with LigandMPNN, obtain alignments with MMseqs2 and fold the sequences with Boltz2.
"""

import os, sys
sys.path.insert(0, os.getcwd()) #to see scripts in current folder

from PipelineScripts.pipeline import *
from PipelineScripts.load_output import LoadOutput
from PipelineScripts.rfdiffusion import RFdiffusion
from PipelineScripts.protein_mpnn import ProteinMPNN
from PipelineScripts.mmseqs2 import MMseqs2
from PipelineScripts.alphafold import AlphaFold

with Pipeline(project="Examples",
              job="RFD-ProteinMPNN-MMseqs-AlphaFold",
              description="redesign of N terminus of rifampicin binding protein"):

    Resources(gpu="80GB", #ask for A100-80GB or H100-80GB
              time="24:00:00",
              memory="16GB")

    rifampicin = LoadOutput("/shares/locbp.chem.uzh/public/BioPipelines/Boltz/rifampicin_001/ToolOutputs/1_Boltz2_output.json")


    rfd = RFdiffusion(pdb=rifampicin,
                      contigs='10-20,A6-140',
                      num_designs=2,
                      steps=20)

    pmpnn = ProteinMPNN(structures=rfd,
                        num_sequences=2)

    msas = MMseqs2(sequences=pmpnn)
    
    af = AlphaFold(sequences=pmpnn,
                   msas=msas)


