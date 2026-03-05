#!/bin/bash
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

# Resolve a single file path from a DataStream JSON.
#
# Usage:
#   source resolve_stream_item.sh
#   resolved=$(resolve_stream_item "/path/to/datastream.json" "item_id")
#
# Requires: jq

resolve_stream_item() {
    local ds_json="$1"
    local item_id="$2"

    if [[ ! -f "$ds_json" ]]; then
        echo "ERROR: DataStream JSON not found: $ds_json" >&2
        return 1
    fi

    local ids_array files_array
    readarray -t ids_array < <(jq -r '.ids[]' "$ds_json")
    readarray -t files_array < <(jq -r '.files[]' "$ds_json")

    # Find index of item_id in ids array
    local idx=-1
    for i in "${!ids_array[@]}"; do
        if [[ "${ids_array[$i]}" == "$item_id" ]]; then
            idx=$i
            break
        fi
    done

    if [[ $idx -eq -1 ]]; then
        echo "ERROR: ID '$item_id' not found in DataStream '$ds_json'" >&2
        return 1
    fi

    # Determine the file pattern for this ID
    local pattern
    if [[ ${#files_array[@]} -eq ${#ids_array[@]} ]]; then
        pattern="${files_array[$idx]}"
    elif [[ ${#files_array[@]} -eq 1 ]]; then
        pattern="${files_array[0]}"
    else
        echo "ERROR: files/ids length mismatch in '$ds_json'" >&2
        return 1
    fi

    # Detect wildcards by checking for '*' in the path
    if [[ "$pattern" == *"*"* ]]; then
        local expanded=()
        for f in $pattern; do  # intentional: unquoted to allow glob expansion
            [[ -e "$f" ]] && expanded+=("$f")
        done

        if [[ ${#expanded[@]} -eq 0 ]]; then
            echo "ERROR: No files match pattern '$pattern' for ID '$item_id'" >&2
            return 1
        fi

        # Best match: prefer exact basename (without extension) == item_id
        local base
        for f in "${expanded[@]}"; do
            base=$(basename "$f")
            base="${base%.*}"
            if [[ "$base" == "$item_id" ]]; then
                echo "$f"
                return 0
            fi
        done

        # Fallback: first match
        echo "${expanded[0]}"
        return 0
    else
        echo "$pattern"
        return 0
    fi
}
