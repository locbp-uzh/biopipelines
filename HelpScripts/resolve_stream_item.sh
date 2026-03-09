#!/bin/bash
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

# Resolve a single file path from a DataStream JSON.
#
# Usage:
#   source resolve_stream_item.sh
#   resolved=$(resolve_stream_item "/path/to/datastream.json" "item_id")

# Capture the directory of this script at source-time so the inline
# Python can find the biopipelines package (parent of HelpScripts/).
_RESOLVE_HELPSCRIPTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

resolve_stream_item() {
    local ds_json="$1"
    local item_id="$2"

    python -c "
import sys, os
sys.path.insert(0, os.path.join(sys.argv[3], '..'))
from biopipelines.biopipelines_io import load_datastream, resolve_file
ds = load_datastream(sys.argv[1])
print(resolve_file(ds, sys.argv[2]))
" "$ds_json" "$item_id" "$_RESOLVE_HELPSCRIPTS_DIR"
}
