#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Look up LigandMPNN position options for a single structure ID.

Usage:
    python resolve_lmpnn_positions.py <positions_json> <struct_id>

Prints the argv tokens NUL-separated, so a value containing spaces stays one token and the caller reads them into a bash array with `mapfile -d ''` instead of `eval`.
"""

import shlex
import sys
import json

with open(sys.argv[1]) as f:
    data = json.load(f)

entry = data.get(sys.argv[2], {})
tokens = []
for key in ("fixed_option", "redesigned_option"):
    tokens.extend(shlex.split(entry.get(key, "") or ""))
sys.stdout.write("".join(token + chr(0) for token in tokens))
