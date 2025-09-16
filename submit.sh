#!/bin/bash
#
# SLURM Submission Script for BioPipelines
#
# This script submits a BioPipelines job to SLURM, mimicking the functionality
# of going to job composer and running the job with the suggested job name.
#
# Usage:
#     ./submit.sh [pipeline_name] [job_name]
#
# If no arguments provided, will look for the most recent pipeline output.
#

set -e  # Exit on any error

# Function to print usage
usage() {
    echo "Usage: $0 [pipeline_script]"
    echo ""
    echo "Run a BioPipelines script and submit the generated job to SLURM"
    echo ""
    echo "Arguments:"
    echo "  pipeline_script  Path to pipeline script (optional, defaults to pipeline.py)"
    echo ""
    echo "Examples:"
    echo "  $0                                           # Run pipeline.py"
    echo "  $0 pipeline.py                               # Run pipeline.py"
    echo "  $0 ExamplePipelines/ligandmpnn_boltz2.py    # Run specific pipeline"
    echo ""
}

# Function to run pipeline and extract SLURM script path
run_pipeline_and_get_script() {
    local pipeline_script="$1"

    echo "Running pipeline: $pipeline_script"
    echo ""

    # Run the pipeline and capture output
    local pipeline_output=$(python "$pipeline_script" 2>&1)
    local exit_code=$?

    # Show pipeline output
    echo "$pipeline_output"
    echo ""

    if [ $exit_code -ne 0 ]; then
        echo "Pipeline execution completed with exit code $exit_code"
        echo "Note: This is expected for dummy/test pipelines"
        echo ""
    fi

    # Extract SLURM script path from pipeline output
    # Look for "Slurm saved to: /path/to/slurm.sh" pattern
    local slurm_script=$(echo "$pipeline_output" | grep "Slurm saved to:" | sed 's/Slurm saved to: //')

    if [ -n "$slurm_script" ] && [ -f "$slurm_script" ]; then
        echo "Found SLURM script: $slurm_script"
        echo "$slurm_script"
    else
        echo "Could not find SLURM script"
        echo ""
    fi
}

# Function to extract job info from pipeline folder
extract_job_info() {
    local job_folder="$1"
    local pipeline_name=""
    local job_name=""
    local job_id=""

    # Extract from folder structure: /shares/USER/PipelineName/JobName_NNN/RunTime/
    if [[ "$job_folder" =~ /([^/]+)/([^/]+)_([0-9]+)/RunTime$ ]]; then
        pipeline_name="${BASH_REMATCH[1]}"
        job_name="${BASH_REMATCH[2]}"
        job_id="${BASH_REMATCH[3]}"
    elif [[ "$job_folder" =~ /([^/]+)/([^/]+)/RunTime$ ]]; then
        pipeline_name="${BASH_REMATCH[1]}"
        job_name="${BASH_REMATCH[2]}"
        job_id="001"
    else
        # Fallback: extract from parent directory names
        local parent_dir=$(dirname "$job_folder")
        local grandparent_dir=$(dirname "$parent_dir")
        pipeline_name=$(basename "$grandparent_dir")
        job_name=$(basename "$parent_dir")
        job_id="001"
    fi

    echo "$pipeline_name:$job_name:$job_id"
}

# Function to submit job
submit_job() {
    local slurm_script="$1"
    local job_name="$2"
    local pipeline_name="$3"
    local job_id="$4"

    # Create SLURM job name following BioPipelines convention: "Pipeline: Job (ID)"
    local slurm_job_name="${pipeline_name}: ${job_name} (${job_id})"

    echo "Submitting job to SLURM..."
    echo "  Pipeline: $pipeline_name"
    echo "  Job: $job_name (ID: $job_id)"
    echo "  SLURM job name: $slurm_job_name"
    echo "  Script: $slurm_script"
    echo ""

    # Check if slurm script exists
    if [ ! -f "$slurm_script" ]; then
        echo "Error: SLURM script not found: $slurm_script"
        return 1
    fi

    # Make sure script is executable
    chmod +x "$slurm_script"

    # Submit with proper job name
    echo "Running: sbatch --job-name=\"$slurm_job_name\" \"$slurm_script\""

    if command -v sbatch >/dev/null 2>&1; then
        local job_output=$(sbatch --job-name="$slurm_job_name" "$slurm_script")
        local exit_code=$?

        if [ $exit_code -eq 0 ]; then
            echo ""
            echo "✓ Job submitted successfully!"
            echo "$job_output"

            # Extract job ID if available
            if [[ "$job_output" =~ Submitted\ batch\ job\ ([0-9]+) ]]; then
                local submitted_job_id="${BASH_REMATCH[1]}"
                echo ""
                echo "Monitor job with:"
                echo "  squeue -j $submitted_job_id"
                echo "  sacct -j $submitted_job_id"
            fi
        else
            echo ""
            echo "✗ Job submission failed with exit code $exit_code"
            return $exit_code
        fi
    else
        echo "Warning: sbatch command not found. This might not be a SLURM system."
        echo "The SLURM script has been generated at: $slurm_script"
        echo "You can submit it manually when on a SLURM system."
        return 0
    fi
}

# Main execution
main() {
    local pipeline_script="$1"

    echo "BioPipelines SLURM Submission Script"
    echo "===================================="
    echo ""

    # Handle help request
    if [ "$1" = "-h" ] || [ "$1" = "--help" ]; then
        usage
        return 0
    fi

    # Default to pipeline.py if no script specified
    if [ -z "$pipeline_script" ]; then
        pipeline_script="pipeline.py"
    fi

    # Check if pipeline script exists
    if [ ! -f "$pipeline_script" ]; then
        echo "Error: Pipeline script not found: $pipeline_script"
        echo ""
        echo "Available options:"
        echo "  ./submit.sh                     # Run pipeline.py"
        echo "  ./submit.sh pipeline.py         # Run specific pipeline script"
        echo "  ./submit.sh ExamplePipelines/ligandmpnn_boltz2.py"
        return 1
    fi

    # Run the pipeline and get SLURM script path
    local slurm_script=$(run_pipeline_and_get_script "$pipeline_script")

    if [ -z "$slurm_script" ]; then
        echo "Error: Could not extract SLURM script path from pipeline output"
        echo "This happens when:"
        echo "1. The pipeline is a dummy/test pipeline (like pipeline.py)"
        echo "2. The pipeline failed to generate SLURM scripts"
        echo "3. The pipeline output format is different than expected"
        echo ""
        echo "To test with a real pipeline, try:"
        echo "  ./submit.sh ExamplePipelines/ligandmpnn_boltz2.py"
        return 1
    fi

    # Extract job information from SLURM script path
    local job_folder=$(dirname "$slurm_script")
    local job_info=$(extract_job_info "$job_folder")
    IFS=':' read -r extracted_pipeline extracted_job extracted_id <<< "$job_info"

    echo ""
    echo "Job Information:"
    echo "  Pipeline script: $pipeline_script"
    echo "  Pipeline: $extracted_pipeline"
    echo "  Job: $extracted_job"
    echo "  ID: $extracted_id"
    echo "  SLURM script: $slurm_script"
    echo ""

    # Submit the job
    submit_job "$slurm_script" "$extracted_job" "$extracted_pipeline" "$extracted_id"
}

# Run main function with all arguments
main "$@"