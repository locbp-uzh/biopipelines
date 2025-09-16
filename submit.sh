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
    echo "Usage: $0 [pipeline_name] [job_name]"
    echo ""
    echo "Submit a BioPipelines job to SLURM"
    echo ""
    echo "Arguments:"
    echo "  pipeline_name    Name of the pipeline (optional)"
    echo "  job_name         Name of the job (optional)"
    echo ""
    echo "If no arguments provided, will look for the most recent pipeline output."
    echo ""
    echo "Examples:"
    echo "  $0                                    # Submit most recent pipeline"
    echo "  $0 DummyPipeline TestJob             # Submit specific pipeline/job"
    echo ""
}

# Function to find the most recent pipeline
find_recent_pipeline() {
    local base_dir="/shares/locbp.chem.uzh/$USER"

    # If base directory doesn't exist, look locally
    if [ ! -d "$base_dir" ]; then
        base_dir="$HOME"
        echo "Warning: SLURM shares directory not found, looking locally in $base_dir"
    fi

    # Find the most recent job folder with a slurm.sh script
    local recent_job=$(find "$base_dir" -name "slurm.sh" -type f -printf '%T@ %p\n' 2>/dev/null | sort -n | tail -1 | cut -d' ' -f2-)

    if [ -n "$recent_job" ]; then
        echo "$(dirname "$recent_job")"
    else
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
    local pipeline_name="$1"
    local job_name="$2"

    echo "BioPipelines SLURM Submission Script"
    echo "===================================="
    echo ""

    # Handle help request
    if [ "$1" = "-h" ] || [ "$1" = "--help" ]; then
        usage
        return 0
    fi

    local job_folder=""
    local slurm_script=""

    if [ -n "$pipeline_name" ] && [ -n "$job_name" ]; then
        # Look for specific pipeline/job
        echo "Looking for pipeline: $pipeline_name, job: $job_name"

        # Try standard locations
        local base_dirs=("/shares/locbp.chem.uzh/$USER" "$HOME")

        for base_dir in "${base_dirs[@]}"; do
            if [ -d "$base_dir" ]; then
                # Look for job folder with pattern: PipelineName/JobName_NNN or PipelineName/JobName
                local pattern1="$base_dir/$pipeline_name/${job_name}_*/RunTime"
                local pattern2="$base_dir/$pipeline_name/$job_name/RunTime"

                for pattern in "$pattern1" "$pattern2"; do
                    local found=$(find $(dirname "$pattern") -maxdepth 2 -name "RunTime" -type d 2>/dev/null | head -1)
                    if [ -n "$found" ]; then
                        job_folder="$found"
                        break 2
                    fi
                done
            fi
        done

        if [ -z "$job_folder" ]; then
            echo "Error: Could not find job folder for pipeline '$pipeline_name', job '$job_name'"
            echo "Searched in: /shares/locbp.chem.uzh/$USER and $HOME"
            return 1
        fi
    else
        # Look for most recent pipeline
        echo "No specific pipeline specified, looking for most recent..."
        job_folder=$(find_recent_pipeline)

        if [ -z "$job_folder" ]; then
            echo "Error: No recent pipeline found with slurm.sh script"
            echo "Please run a pipeline first (e.g., python pipeline.py)"
            return 1
        fi
    fi

    # Construct slurm script path
    slurm_script="$job_folder/slurm.sh"

    echo "Found job folder: $job_folder"

    # Extract job information
    local job_info=$(extract_job_info "$job_folder")
    IFS=':' read -r extracted_pipeline extracted_job extracted_id <<< "$job_info"

    # Use extracted info if not provided
    if [ -z "$pipeline_name" ]; then
        pipeline_name="$extracted_pipeline"
    fi
    if [ -z "$job_name" ]; then
        job_name="$extracted_job"
    fi

    echo ""
    echo "Job Information:"
    echo "  Pipeline: $pipeline_name"
    echo "  Job: $job_name"
    echo "  ID: $extracted_id"
    echo "  SLURM script: $slurm_script"
    echo ""

    # Submit the job
    submit_job "$slurm_script" "$job_name" "$pipeline_name" "$extracted_id"
}

# Run main function with all arguments
main "$@"