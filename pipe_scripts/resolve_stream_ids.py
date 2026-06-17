#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Print expanded IDs from a DataStream JSON, one per line.

Usage:
    python resolve_stream_ids.py <ds_json> [--valid-set]

Handles lazy patterns by loading the DataStream in runtime mode, which expands
patterns against the map_table. With ``--valid-set``, restrict to ids whose file
is present on disk (a filtered Pool/Panda declares every original id but
materializes only survivors), via iterate_files which skips absent ids.
"""

import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from biopipelines.biopipelines_io import load_datastream, iterate_files

ds = load_datastream(sys.argv[1])
if "--valid-set" in sys.argv[2:]:
    # iterate_files warns to stdout for absent ids; capture so only ids reach stdout.
    import contextlib
    ids = []
    with contextlib.redirect_stdout(sys.stderr):
        for item_id, _file in iterate_files(ds):
            ids.append(item_id)
    for item_id in ids:
        print(item_id)
else:
    for item_id in ds.ids_expanded:
        print(item_id)
