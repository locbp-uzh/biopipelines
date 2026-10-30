#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Print the residue code RFdiffusion-AllAtom designs one input structure against.

A structure named after a ligand of the compounds stream (a prepared ligand keeps its ligand's id) gets that ligand's `code`; any other structure gets the first ligand's, which is how one ligand is broadcast over several proteins.

Usage:
    python resolve_rfdaa_ligand_code.py <ligand_json> <struct_id>
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from biopipelines.biopipelines_io import load_datastream, get_value  # noqa: E402


def main():
    ligands = load_datastream(sys.argv[1])
    ids = list(ligands.ids_expanded)
    if not ids:
        print("Error: the ligand stream holds no ligand", file=sys.stderr)
        sys.exit(1)
    struct_id = sys.argv[2]
    print(get_value(ligands, struct_id if struct_id in ids else ids[0], column="code"))


if __name__ == "__main__":
    main()
