#!/bin/bash
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

# Resolve a single item from a DataStream JSON.
#
# For a file-based stream this resolves the item's file path. For a value-based
# stream, pass a column name as the 3rd argument to resolve that map_table
# column's value for the id instead.
#
# Usage:
#   source resolve_stream_item.sh
#   path=$(resolve_stream_item "/path/to/datastream.json" "item_id")
#   value=$(resolve_stream_item "/path/to/compounds.json" "item_id" "code")

# Capture the directory of this script at source-time so the inline
# Python can find the biopipelines package (parent of pipe_scripts/).
_RESOLVE_HELPSCRIPTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

resolve_stream_item() {
    local ds_json="$1"
    local item_id="$2"
    local column="$3"

    python -c "
import sys, os
sys.path.insert(0, os.path.join(sys.argv[3], '..'))
from biopipelines.biopipelines_io import load_datastream, resolve_file, get_value
ds = load_datastream(sys.argv[1])
column = sys.argv[4] if len(sys.argv) > 4 and sys.argv[4] else None
print(get_value(ds, sys.argv[2], column=column) if column else resolve_file(ds, sys.argv[2]))
" "$ds_json" "$item_id" "$_RESOLVE_HELPSCRIPTS_DIR" "$column"
}
