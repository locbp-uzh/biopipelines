#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Resolve per-input HBDesigner constraints at execution time.

HBDesigner runs over multiple input structures in one step. The `guide_res`,
`guide_seq`, and `anchor_res` arguments are each either a plain literal
(broadcast to every input) or a ``TABLE_REFERENCE:<path>:<column>`` token whose
per-input cell holds the value (resolved by id match). This script reads the
input structures DataStream, resolves each argument for every id in one pass,
and writes a JSON mapping ``{id: {guide_res, guide_seq, anchor_res}}`` that the
generated bash loop reads one id at a time. It writes data only (a JSON).

Usage:
    python pipe_hbdesigner_constraints.py <config_json>

<config_json> is a JSON object with keys ``structures_json``, ``guide_res``,
``guide_seq``, ``anchor_res``, ``output_json``. Each of ``guide_res`` /
``guide_seq`` / ``anchor_res`` is a literal string, a ``TABLE_REFERENCE:...``
token, or ``-`` (treated as empty / unset).
"""

import sys
import os
import json
import pandas as pd

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from biopipelines.biopipelines_io import (
    load_datastream, iterate_files, load_table, lookup_table_value,
)


def resolve_arg(arg, struct_id):
    """Resolve one constraint argument for a single input id.

    `arg` is ``'-'`` (unset -> ""), a literal (broadcast — same value for every
    id), or a ``TABLE_REFERENCE:<path>:<column>`` token whose per-id cell holds
    the value. Returns the resolved string ("" when unset / empty cell).

    A table cell is arbitrary user data that ends up in a bash ``eval``; reject
    shell metacharacters here (literals were already validated at construction
    time). This is the runtime counterpart to ``_validate_freeform_string``.
    """
    if arg == "-" or arg == "":
        return ""
    if arg.startswith("TABLE_REFERENCE:"):
        table, column = load_table(arg)
        cell = lookup_table_value(table, struct_id, column)
        if cell is None or (isinstance(cell, float) and pd.isna(cell)):
            return ""
        value = str(cell)
        bad = [c for c in ('"', "`", "$", "\\") if c in value]
        if bad:
            raise ValueError(
                f"constraint value for {struct_id!r} (column {column!r}) contains "
                f"unsafe shell character(s) {bad}: {value!r}"
            )
        return value
    return arg  # literal, broadcast


def main():
    if len(sys.argv) != 2:
        print("Usage: python pipe_hbdesigner_constraints.py <config_json>",
              file=sys.stderr)
        sys.exit(1)

    with open(sys.argv[1]) as f:
        cfg = json.load(f)
    structures_json = cfg["structures_json"]
    guide_res_arg = cfg["guide_res"]
    guide_seq_arg = cfg["guide_seq"]
    anchor_res_arg = cfg["anchor_res"]
    output_json = cfg["output_json"]

    ds = load_datastream(structures_json)
    struct_ids = [sid for sid, _ in iterate_files(ds)]
    if not struct_ids:
        raise ValueError(f"No structures found in DataStream: {structures_json}")

    options = {}
    for struct_id in struct_ids:
        options[struct_id] = {
            "guide_res": resolve_arg(guide_res_arg, struct_id),
            "guide_seq": resolve_arg(guide_seq_arg, struct_id),
            "anchor_res": resolve_arg(anchor_res_arg, struct_id),
        }

    os.makedirs(os.path.dirname(output_json), exist_ok=True)
    with open(output_json, "w") as f:
        json.dump(options, f, indent=2)
    print(f"Resolved HBDesigner constraints for {len(struct_ids)} structure(s): {output_json}")


if __name__ == "__main__":
    main()
