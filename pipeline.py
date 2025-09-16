#!/usr/bin/env python3
"""
Dummy Pipeline Script for Testing SLURM Submission

This script creates a simple BioPipelines workflow for testing
the pipeline → SLURM → job submission process without complex dependencies.
"""

import os, sys
sys.path.insert(0, os.getcwd())  # to see scripts in current folder

from PipelineScripts.pipeline import Pipeline
from PipelineScripts.ligand_mpnn import LigandMPNN
from PipelineScripts.slice_datasheet import SliceDatasheet

def main():
    """Create and execute a dummy pipeline."""

    # Create dummy pipeline with debug=True to avoid creating cluster directories
    pipeline = Pipeline(
        pipeline_name="DummyPipeline",
        job_name="TestJob",
        job_description="Simple test pipeline for SLURM submission testing",
        debug=True  # This prevents creation of /shares directories
    )

    # Configure resources (modest for testing)
    pipeline.resources(
        gpu="V100",
        time="2:00:00",
        memory="8GB"
    )

    print("="*60)
    print("DUMMY PIPELINE CONFIGURATION")
    print("="*60)
    print(f"Pipeline Name: {pipeline.pipeline_name}")
    print(f"Job Name: {pipeline.job_name}")
    print(f"Description: {pipeline.job_description}")
    print(f"Resources: {pipeline.global_resources}")
    print()

    # Note: This is a minimal pipeline for testing purposes
    # In a real scenario, you would add actual tools with real inputs
    print("Note: This is a dummy pipeline for testing SLURM submission.")
    print("No actual computation will be performed.")
    print()

    try:
        # Example of how to add tools (commented out since we need real inputs)
        # lmpnn = pipeline.add(LigandMPNN(
        #     structures="protein.pdb",
        #     ligand="LIG",
        #     num_sequences=10
        # ))
        #
        # first_five = pipeline.add(SliceDatasheet(
        #     input=lmpnn.output,
        #     n_rows=5
        # ))

        # For testing, we'll just generate scripts without tools
        print("Generating pipeline scripts...")
        pipeline.save()

        print("\nGenerating SLURM submission script...")
        # Use a placeholder email - user should modify this
        pipeline.slurm(email="your.email@uzh.ch")

        print("\n" + "="*60)
        print("PIPELINE GENERATION COMPLETE")
        print("="*60)
        print(f"Pipeline folder: {pipeline.folders['job']}")
        print(f"SLURM script: {pipeline.slurm_script}")
        print()
        print("To submit to SLURM, run:")
        print("    ./submit.sh")
        print()
        print("Or manually submit with:")
        print(f"    sbatch {pipeline.slurm_script}")

    except Exception as e:
        print(f"Expected error (empty pipeline): {e}")
        print("\n" + "="*60)
        print("DUMMY PIPELINE DEMONSTRATION COMPLETE")
        print("="*60)
        print("This demonstrates the BioPipelines structure:")
        print(f"- Pipeline folder structure: {pipeline.folders}")
        print("- Resource configuration: working")
        print("- Tool import system: working")
        print()
        print("To create a real pipeline:")
        print("1. Add actual tools with real inputs")
        print("2. Call pipeline.save() and pipeline.slurm()")
        print("3. Use ./submit.sh to submit to SLURM")
        print()
        print("For working examples, see ExamplePipelines/ directory.")

    return 0

if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)