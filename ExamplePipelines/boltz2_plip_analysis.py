#!/usr/bin/env python3
"""
This pipeline shows how to use PLIP for protein-ligand interaction profiling.

Takes protein-ligand complexes (e.g., from Boltz2), analyzes interactions with PLIP,
and generates interaction profiles, visualizations, and reports.
"""

import os, sys
sys.path.insert(0, os.getcwd()) #to see scripts in current folder

from PipelineScripts.pipeline import Pipeline
from PipelineScripts.load_output import LoadOutput
from PipelineScripts.plip import PLIP

pipeline = Pipeline(
    pipeline_name="Boltz2-PLIP-Analysis", #Will create a folder in /shares/USER/<pipeline_name>
    job_name="rifampicin_interactions", #Unique job folder in /shares/USER/<pipeline_name>/job_name_NNN
    job_description="Protein-ligand interaction analysis using PLIP on Boltz2 structures")

pipeline.resources(
    cpu=4,
    time="4:00:00",
    memory="8GB"
)

# Load existing Boltz2 results
boltz_results = pipeline.add(LoadOutput("/shares/locbp.chem.uzh/public/BioPipelines/Boltz/rifampicin_001/ToolOutputs/1_Boltz2_output.json"))

# Analyze protein-ligand interactions with PLIP
plip_analysis = pipeline.add(PLIP(
    structures=boltz_results,  # Input structures from Boltz2
    ligand="",  # Analyze all ligands (or specify specific ligand ID)
    output_format=['xml', 'txt', 'pymol'],  # Generate XML reports, text summaries, and PyMOL sessions
    create_pymol=True,  # Generate PyMOL visualization files
    create_images=False,  # Don't generate ray-traced images (takes longer)
    analyze_peptides=False,  # Focus on small molecule ligands
    analyze_intra=False,  # Don't analyze intra-chain interactions
    analyze_dna=False,  # Don't analyze DNA interactions
    max_threads=4,  # Use parallel processing
    verbose=True  # Enable detailed output
))

# Save pipeline configuration and generate execution scripts
pipeline.save()
pipeline.slurm()