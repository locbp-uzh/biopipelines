#!/bin/bash
# Resubmit a SLURM job with correct output path
# Usage: resubmit.sh <path_to_slurm.sh>

if [ $# -eq 0 ]; then
    echo "Error: No SLURM script provided"
    echo "Usage: resubmit.sh <path_to_slurm.sh>"
    exit 1
fi

SLURM_SCRIPT="$1"

if [ ! -f "$SLURM_SCRIPT" ]; then
    echo "Error: SLURM script not found: $SLURM_SCRIPT"
    exit 1
fi

# Get the directory containing the SLURM script
SCRIPT_DIR=$(dirname "$SLURM_SCRIPT")
SCRIPT_DIR=$(cd "$SCRIPT_DIR" && pwd)  # Get absolute path

# Submit job with output redirected to the script's directory
echo "Submitting job from: $SLURM_SCRIPT"
echo "Output will be written to: $SCRIPT_DIR/job.out"

sbatch --output="$SCRIPT_DIR/job.out" "$SLURM_SCRIPT"

if [ $? -eq 0 ]; then
    echo "Job submitted successfully"
else
    echo "Error: Job submission failed"
    exit 1
fi
