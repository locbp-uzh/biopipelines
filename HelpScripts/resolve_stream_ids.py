#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Print expanded IDs from a DataStream JSON, one per line.

Usage:
    python resolve_stream_ids.py <ds_json>

Handles lazy patterns by loading the DataStream in runtime mode,
which expands patterns against the map_table.
"""

import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from biopipelines.biopipelines_io import load_datastream

ds = load_datastream(sys.argv[1])
for item_id in ds.ids_expanded:
    print(item_id)
