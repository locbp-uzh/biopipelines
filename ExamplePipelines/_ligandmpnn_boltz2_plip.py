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
This pipeline shows how to integrate PLIP analysis in a design workflow.

Generates sequences with LigandMPNN, folds them with Boltz2, then analyzes
protein-ligand interactions with PLIP to profile binding characteristics.
"""

import os, sys
sys.path.insert(0, os.getcwd()) #to see scripts in current folder

from PipelineScripts.pipeline import Pipeline
from PipelineScripts.load_output import LoadOutput
from PipelineScripts.ligand_mpnn import LigandMPNN
from PipelineScripts.boltz2 import Boltz2
from PipelineScripts.plip import PLIP

pipeline = Pipeline(
    pipeline_name="Examples",
    job_name="LigandMPNN-Boltz2-PLIP",
    job_description="LigandMPNN sequence design followed by Boltz2 folding and PLIP interaction analysis")

pipeline.resources(
    gpu="80GB", #ask for A100-80GB or H100-80GB
    time="24:00:00",
    memory="16GB"
)

# Load original protein-ligand complex
original = pipeline.add(LoadOutput("/shares/locbp.chem.uzh/public/BioPipelines/Boltz/rifampicin_001/ToolOutputs/1_Boltz2_output.json"))

# Generate diverse sequences with LigandMPNN
lmpnn = pipeline.add(LigandMPNN(structures=original,
                                ligand="RFP",  # Rifampicin ligand name
                                num_sequences=10,
                                design_within=4))

# Fold sequences with Boltz2
boltz_holo = pipeline.add(Boltz2(proteins=lmpnn,
                                 ligands=original)) # Use ligands from original structure

# Analyze protein-ligand interactions with PLIP
plip_analysis = pipeline.add(PLIP(
    structures=boltz_holo,
    ligand="RFP",  # Focus on rifampicin interactions
    output_format=['xml', 'txt', 'pymol'],  # Comprehensive output
    create_pymol=True,  # Generate PyMOL sessions for visualization
    create_images=True,  # Generate publication-quality images
    analyze_peptides=False,
    analyze_intra=False,
    analyze_dna=False,
    max_threads=4,
    verbose=True
))

pipeline.slurm()