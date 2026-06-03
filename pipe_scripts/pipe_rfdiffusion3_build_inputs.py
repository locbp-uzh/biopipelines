#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Build the RFdiffusion3 (foundry) inputs JSON at execution time.

RFdiffusion3 reads one JSON whose keys are design names and whose values carry
the per-design spec (contig, length, ligand, selects, and an `input` PDB path).
The input PDB paths are lazy (resolved at runtime), so this script assembles the
final JSON from:

  - a config-time TEMPLATE entry holding the shared design params,
  - the input structures DataStream (one entry per pdb id, keyed by that id,
    with its resolved `input` path) — or, with no PDB, a single de-novo entry
    keyed by the given prefix and no `input`,
  - the ligand residue `code` resolved from the compounds stream (broadcast to
    every entry), when a ligand was provided.

It writes data only (the JSON) — never bash.

Usage:
    python pipe_rfdiffusion3_build_inputs.py \\
        --template <template.json> --output <inputs.json> \\
        [--structures-json <ds.json>] [--ligand-json <compounds.json>] \\
        [--denovo-prefix <name>]
"""

import argparse
import json
import os
import sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from biopipelines.biopipelines_io import load_datastream, iterate_files, get_value


def main():
    ap = argparse.ArgumentParser(description="Build RFdiffusion3 inputs JSON")
    ap.add_argument("--template", required=True, help="Config-time template entry JSON")
    ap.add_argument("--output", required=True, help="Path to write the final inputs JSON")
    ap.add_argument("--structures-json", default=None,
                    help="Input structures DataStream JSON (one entry per pdb id)")
    ap.add_argument("--ligand-json", default=None,
                    help="Compounds-stream JSON; the residue `code` is read from it")
    ap.add_argument("--denovo-prefix", default=None,
                    help="Design key for the no-PDB (de novo) single-entry case")
    args = ap.parse_args()

    with open(args.template) as f:
        template = json.load(f)

    # Ligand code (broadcast to every entry) resolved once from the compounds
    # stream's first id — read its `code` value, not ids[0] at config time.
    ligand_code = None
    if args.ligand_json:
        lig = load_datastream(args.ligand_json)
        first_id = lig.ids_expanded[0]
        ligand_code = get_value(lig, first_id, column="code")

    config = {}

    def make_entry():
        entry = dict(template)
        if ligand_code:
            entry["ligand"] = ligand_code
        return entry

    if args.structures_json:
        ds = load_datastream(args.structures_json)
        for pdb_id, pdb_path in iterate_files(ds):
            entry = make_entry()
            entry["input"] = pdb_path
            config[pdb_id] = entry
        if not config:
            raise ValueError(f"No structures found in DataStream: {args.structures_json}")
    else:
        if not args.denovo_prefix:
            raise ValueError("--denovo-prefix is required when --structures-json is absent")
        config[args.denovo_prefix] = make_entry()

    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    with open(args.output, "w") as f:
        json.dump(config, f, indent=2)
    print(f"Built RFdiffusion3 inputs JSON with {len(config)} entr(ies): {args.output}")


if __name__ == "__main__":
    main()
