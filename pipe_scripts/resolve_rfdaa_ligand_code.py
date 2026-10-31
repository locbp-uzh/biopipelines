#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Print the residue code RFdiffusion-AllAtom designs one input structure against.

A structure named after a ligand of the compounds stream (a prepared ligand keeps its ligand's id) gets that ligand's `code`. When the stream holds one ligand it is the ligand of every structure, which is how one ligand is broadcast over several proteins. A structure that matches none of several ligands has no ligand to be designed against: the script says so on stderr and exits 1, and the wrapper records that structure as failed.

Usage:
    python resolve_rfdaa_ligand_code.py <ligand_json> <struct_id>
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from biopipelines.biopipelines_io import load_datastream, get_value  # noqa: E402

SHOWN_LIGANDS = 5


def main():
    ligands = load_datastream(sys.argv[1])
    ids = list(ligands.ids_expanded)
    if not ids:
        print("Error: the ligand stream holds no ligand", file=sys.stderr)
        sys.exit(1)
    struct_id = sys.argv[2]
    if struct_id in ids:
        chosen = struct_id
    elif len(ids) == 1:
        chosen = ids[0]
    else:
        listed = ", ".join(ids[:SHOWN_LIGANDS]) + (", ..." if len(ids) > SHOWN_LIGANDS else "")
        print(f"Error: structure '{struct_id}' is named after none of the {len(ids)} ligands ({listed}), "
              "so which one to design it against is not defined", file=sys.stderr)
        sys.exit(1)
    print(get_value(ligands, chosen, column="code"))


if __name__ == "__main__":
    main()
