"""
This pipeline shows how to run RFdiffusion-AllAtom, generate sequences far from the ligand with ProteinMPNN and close with LigandMPNN, obtain alignments with MMseqs2 and fold the sequences with Boltz2.
"""

import os, sys
sys.path.insert(0, os.getcwd()) #to see scripts in current folder

from PipelineScripts.pipeline import Pipeline
from PipelineScripts.load_output import LoadOutput
from PipelineScripts.rfdiffusion_allatom import RFdiffusionAllAtom
from PipelineScripts.distance_selector import DistanceSelector
from PipelineScripts.protein_mpnn import ProteinMPNN
from PipelineScripts.ligand_mpnn import LigandMPNN
from PipelineScripts.stitch_sequences import StitchSequences
from PipelineScripts.mmseqs2 import MMseqs2
from PipelineScripts.boltz2 import Boltz2

pipeline = Pipeline(
    pipeline_name="Examples",
    job_name="RFDAA-ProteinMPNN-LigandMPNN-MMseqs-Boltz",
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
#this generates a table showing for each structure id a pymol selection for residues within and beyond the distance from the ligand
distances = pipeline.add(DistanceSelector(structures=rfdaa,
                                          ligand="LIG",
                                          distance=5))

pmpnn = pipeline.add(ProteinMPNN(structures=rfdaa,
                                 num_sequences=2,
                                 redesigned=distances.datasheets.selections.beyond))

lmpnn = pipeline.add(LigandMPNN(structures=rfdaa,
                                ligand="LIG", #in ligand mpnn you should always specify the ligand name.
                                num_sequences=2,
                                redesigned=distances.datasheets.selections.within))

sequences = pipeline.add(StitchSequences(sequences=[pmpnn,lmpnn],
                                         selections=["",distances.datasheets.selections.within]))

msas = pipeline.add(MMseqs2(sequences=sequences))
boltz_holo = pipeline.add(Boltz2(proteins=sequences,
                                ligands=rifampicin, #ligand smiles taken from original
                                msas=msas)) #MSAs are passed with <tool output>, not with <tool output>.msas


pipeline.save()
pipeline.slurm() 
