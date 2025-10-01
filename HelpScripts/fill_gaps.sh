#!/bin/bash

# Script to fill gaps in numbered tool completion files with Boltz2 placeholders
# Usage: ./fill_gaps.sh "/path/to/pipeline/folder"

if [ $# -ne 1 ]; then
    echo "Usage: $0 <pipeline_folder_path>"
    echo "Example: $0 '/shares/locbp.chem.uzh/gquarg/BioPipelines/LigandMPNN-MutationComposer-MMseqs-Cycle/HT7_Cy7_C_R_W2_004'"
    exit 1
fi

PIPELINE_FOLDER="$1"

if [ ! -d "$PIPELINE_FOLDER" ]; then
    echo "Error: Directory '$PIPELINE_FOLDER' does not exist"
    exit 1
fi

echo "Processing folder: $PIPELINE_FOLDER"

# Find all files matching pattern NNN_Tool_COMPLETED
FILES=($(find "$PIPELINE_FOLDER" -name "*_*_COMPLETED" -type f | sort))

if [ ${#FILES[@]} -eq 0 ]; then
    echo "No completion files found matching pattern NNN_Tool_COMPLETED"
    exit 0
fi

echo "Found ${#FILES[@]} completion files"

# Extract numbers and find min/max
NUMBERS=()
for file in "${FILES[@]}"; do
    basename_file=$(basename "$file")
    if [[ $basename_file =~ ^([0-9]+)_.*_COMPLETED$ ]]; then
        NUMBERS+=(${BASH_REMATCH[1]})
    fi
done

if [ ${#NUMBERS[@]} -eq 0 ]; then
    echo "No valid numbered completion files found"
    exit 0
fi

# Sort numbers and find min/max
IFS=$'\n' SORTED_NUMBERS=($(sort -n <<<"${NUMBERS[*]}"))
unset IFS

MIN_NUM=${SORTED_NUMBERS[0]}
MAX_NUM=${SORTED_NUMBERS[-1]}

echo "Number range: $MIN_NUM to $MAX_NUM"

# Convert to integers for arithmetic
MIN_NUM=$((10#$MIN_NUM))
MAX_NUM=$((10#$MAX_NUM))

# Check for gaps and fill with Boltz2
CREATED_COUNT=0
for ((i=MIN_NUM; i<=MAX_NUM; i++)); do
    # Format number with leading zeros (3 digits)
    FORMATTED_NUM=$(printf "%03d" $i)

    # Check if any file with this number exists
    PATTERN="$PIPELINE_FOLDER/${FORMATTED_NUM}_*_COMPLETED"
    if ! ls $PATTERN 1> /dev/null 2>&1; then
        # Gap found - create Boltz2 placeholder
        PLACEHOLDER_FILE="$PIPELINE_FOLDER/${FORMATTED_NUM}_Boltz2_COMPLETED"
        touch "$PLACEHOLDER_FILE"
        echo "Created: $PLACEHOLDER_FILE"
        ((CREATED_COUNT++))
    fi
done

echo "Filled $CREATED_COUNT gaps with Boltz2 placeholders"
echo "Script completed successfully"