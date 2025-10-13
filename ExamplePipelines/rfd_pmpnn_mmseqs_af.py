"""
This pipeline shows how to run RFdiffusion-AllAtom, generate sequences far from the ligand with ProteinMPNN and close with LigandMPNN, obtain alignments with MMseqs2 and fold the sequences with Boltz2.
"""

import os, sys
sys.path.insert(0, os.getcwd()) #to see scripts in current folder

from PipelineScripts.pipeline import Pipeline
from PipelineScripts.load_output import LoadOutput
from PipelineScripts.rfdiffusion import RFdiffusion
from PipelineScripts.distance_selector import DistanceSelector
from PipelineScripts.protein_mpnn import ProteinMPNN
from PipelineScripts.mmseqs2 import MMseqs2
from PipelineScripts.alphafold import AlphaFold

pipeline = Pipeline(
    pipeline_name="Examples",
    job_name="RFD-ProteinMPNN-MMseqs-AlphaFold",
    job_description="redesign of N terminus of rifampicin binding protein")

pipeline.resources(
    gpu="80GB", #ask for A100-80GB or H100-80GB
    time="24:00:00",
    memory="16GB"
)

rifampicin = pipeline.add(LoadOutput("/shares/locbp.chem.uzh/public/BioPipelines/Boltz/rifampicin_001/ToolOutputs/1_Boltz2_output.json"))


rfd = pipeline.add(RFdiffusion(pdb=rifampicin, 
                                contigs='10-20,A6-140',
                                num_designs=2,
                                steps=20))

pmpnn = pipeline.add(ProteinMPNN(structures=rfd,
                                 num_sequences=2))

msas = pipeline.add(MMseqs2(sequences=pmpnn))
af = pipeline.add(AlphaFold(proteins=pmpnn,
                            msas=msas))


pipeline.slurm() 
