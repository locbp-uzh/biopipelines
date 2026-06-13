#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Look up resolved HBDesigner constraints for a single input id.

Reads the JSON written by pipe_hbdesigner_constraints.py and prints three lines
for the requested id, so a bash loop can capture them with sed -n:

    Line 1: guide_res   (empty if unset)
    Line 2: guide_seq   (empty if unset)
    Line 3: anchor_res  (empty if unset)

Usage:
    python resolve_hbdesigner_constraints.py <constraints_json> <struct_id>
"""

import sys
import json

with open(sys.argv[1]) as f:
    data = json.load(f)

entry = data.get(sys.argv[2], {})
print(entry.get("guide_res", ""))
print(entry.get("guide_seq", ""))
print(entry.get("anchor_res", ""))
