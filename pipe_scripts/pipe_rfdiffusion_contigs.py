#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Resolve per-PDB RFdiffusion contig options at execution time.

RFdiffusion (and the AllAtom / RFdiffusion3 variants) can run over multiple
input PDBs in one step. The `contigs`, `inpaint`, and `inpaint_str` arguments
are each either a plain literal (broadcast to every input PDB) or a
``TABLE_REFERENCE:<path>:<column>`` token whose per-PDB cell holds the value
(resolved by id match). This script reads the input structures DataStream,
resolves each argument for every pdb id in one pass, and writes a JSON mapping
``{pdb_id: {contigs, inpaint, inpaint_str}}`` that the generated bash loop reads
one id at a time. It writes data only (a JSON) — never bash.

Usage:
    python pipe_rfdiffusion_contigs.py <structures_json> <contigs> <inpaint> <inpaint_str> <output_json>

Each of <contigs>/<inpaint>/<inpaint_str> is a literal string, a
``TABLE_REFERENCE:...`` token, or ``-`` (treated as empty / unset).
"""

import sys
import os
import json
import pandas as pd

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from biopipelines.biopipelines_io import (
    load_datastream, iterate_files, load_table, lookup_table_value,
)


def resolve_arg(arg, pdb_id):
    """Resolve one selection argument for a single pdb id.

    `arg` is ``'-'`` (unset -> ""), a literal (broadcast — same value for every
    id), or a ``TABLE_REFERENCE:<path>:<column>`` token whose per-id cell holds
    the value. Returns the resolved string ("" when unset / empty cell).
    """
    if arg == "-" or arg == "":
        return ""
    if arg.startswith("TABLE_REFERENCE:"):
        table, column = load_table(arg)
        cell = lookup_table_value(table, pdb_id, column)
        if cell is None or (isinstance(cell, float) and pd.isna(cell)):
            return ""
        return str(cell)
    return arg  # literal, broadcast


def main():
    if len(sys.argv) != 6:
        print("Usage: python pipe_rfdiffusion_contigs.py "
              "<structures_json> <contigs> <inpaint> <inpaint_str> <output_json>",
              file=sys.stderr)
        sys.exit(1)

    structures_json = sys.argv[1]
    contigs_arg = sys.argv[2]
    inpaint_arg = sys.argv[3]
    inpaint_str_arg = sys.argv[4]
    output_json = sys.argv[5]

    ds = load_datastream(structures_json)
    pdb_ids = [pdb_id for pdb_id, _ in iterate_files(ds)]
    if not pdb_ids:
        raise ValueError(f"No structures found in DataStream: {structures_json}")

    options = {}
    for pdb_id in pdb_ids:
        contigs = resolve_arg(contigs_arg, pdb_id)
        if not contigs:
            raise ValueError(f"contigs resolved to empty for {pdb_id}")
        options[pdb_id] = {
            "contigs": contigs,
            "inpaint": resolve_arg(inpaint_arg, pdb_id),
            "inpaint_str": resolve_arg(inpaint_str_arg, pdb_id),
        }

    os.makedirs(os.path.dirname(output_json), exist_ok=True)
    with open(output_json, "w") as f:
        json.dump(options, f, indent=2)
    print(f"Resolved RFdiffusion contig options for {len(pdb_ids)} structure(s): {output_json}")


if __name__ == "__main__":
    main()
