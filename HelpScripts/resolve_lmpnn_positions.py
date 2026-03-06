#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Look up LigandMPNN position options for a single structure ID.

Usage:
    python resolve_lmpnn_positions.py <positions_json> <struct_id>

Prints two lines:
    Line 1: fixed_option (e.g. --fixed_residues "A10 A11")
    Line 2: redesigned_option (e.g. --redesigned_residues "A20 A21")
"""

import sys
import json

with open(sys.argv[1]) as f:
    data = json.load(f)

entry = data.get(sys.argv[2], {})
print(entry.get("fixed_option", ""))
print(entry.get("redesigned_option", ""))
