#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Look up resolved RFdiffusion contig options for a single PDB id.

Reads the JSON written by pipe_rfdiffusion_contigs.py and prints three lines for
the requested id, so a bash loop can capture them with sed -n:

    Line 1: contigs
    Line 2: inpaint      (empty if unset)
    Line 3: inpaint_str  (empty if unset)

Usage:
    python resolve_rfdiffusion_contigs.py <contigs_json> <pdb_id>
"""

import sys
import json

with open(sys.argv[1]) as f:
    data = json.load(f)

entry = data.get(sys.argv[2], {})
print(entry.get("contigs", ""))
print(entry.get("inpaint", ""))
print(entry.get("inpaint_str", ""))
