#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Resolve fixed/redesigned positions and stamp LASErMPNN B-factors.

LASErMPNN has no residue-list flag: fixed vs designed positions are driven by
the input PDB's B-factor column with --fix_beta (B=1.0 fixed, B=0.0 designed).
This script has two modes.

resolve:
    One pass over the input DataStream. For each input id, resolve the fixed or
    redesigned selection (a broadcast PyMOL string or a per-input table-column
    reference) into a designed/fixed (chain, resnum) set and write a positions
    JSON: {id: {"fixed": "A1-3+B7", "designed": "..."}}.

stamp:
    Rewrite the B-factor column of one input PDB in place (column-exact edit, so
    HETATM ligand records and all other fields are preserved) so designed protein
    residues get B=0.0 and everything else B=1.0. Only ATOM/HETATM lines are
    touched; the designed set covers protein residues only.

Usage:
    python pipe_lasermpnn_bfactor.py resolve <args_json>
    python pipe_lasermpnn_bfactor.py stamp <positions_json> <id> <in_pdb> <out_pdb>
"""

import sys
import os
import json

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from biopipelines.biopipelines_io import load_datastream, iterate_files, load_table, lookup_table_value
from biopipelines.pdb_parser import parse_pdb_file, STANDARD_RESIDUES, field_chain, field_res_seq, field_res_name
from biopipelines.sele_utils import sele_to_list, chain_aware_sele


def _resolve_selection_for_ids(reference, design_ids, default_chain):
    """Map each design id to a list of (chain, resnum) tuples.

    `reference` is "-" (empty), a broadcast PyMOL string, or a
    TABLE_REFERENCE:path:col resolved per id.
    """
    if not reference or reference == "-":
        return {sid: [] for sid in design_ids}

    if not reference.startswith("TABLE_REFERENCE:"):
        shared = [_with_chain(t, default_chain) for t in sele_to_list(reference)]
        return {sid: shared for sid in design_ids}

    table, column = load_table(reference)
    out = {}
    for sid in design_ids:
        try:
            value = lookup_table_value(table, sid, column)
            out[sid] = [_with_chain(t, default_chain) for t in sele_to_list(value)]
        except KeyError:
            print(f"Warning: no table entry for {sid} in column {column}", file=sys.stderr)
            out[sid] = []
    return out


def _with_chain(tup, default_chain):
    chain, resnum = tup
    return (chain if chain else default_chain, resnum)


def _protein_residues(pdb_path):
    """All protein (chain, resnum) tuples in a PDB."""
    residues = []
    seen = set()
    for atom in parse_pdb_file(pdb_path):
        if atom.res_name in STANDARD_RESIDUES:
            key = (atom.chain, atom.res_num)
            if key not in seen:
                seen.add(key)
                residues.append(key)
    return residues


def do_resolve(args_json):
    with open(args_json) as f:
        cfg = json.load(f)

    ds = load_datastream(cfg["structures_json"])
    entries = list(iterate_files(ds))
    if not entries:
        raise ValueError(f"No structures in DataStream: {cfg['structures_json']}")
    design_ids = [sid for sid, _ in entries]
    pdb_by_id = {sid: path for sid, path in entries}
    default_chain = cfg.get("default_chain", "A")

    fixed_ref = cfg["fixed"]
    redesigned_ref = cfg["redesigned"]

    fixed_per_id = _resolve_selection_for_ids(fixed_ref, design_ids, default_chain)
    redesigned_per_id = _resolve_selection_for_ids(redesigned_ref, design_ids, default_chain)

    result = {}
    for sid in design_ids:
        all_res = _protein_residues(pdb_by_id[sid])
        all_set = set(all_res)
        if redesigned_ref and redesigned_ref != "-":
            designed = [r for r in redesigned_per_id[sid] if r in all_set]
        elif fixed_ref and fixed_ref != "-":
            fixed_set = set(fixed_per_id[sid])
            designed = [r for r in all_res if r not in fixed_set]
        else:
            designed = list(all_res)

        designed_set = set(designed)
        fixed = [r for r in all_res if r not in designed_set]
        result[sid] = {
            "designed": chain_aware_sele(designed),
            "fixed": chain_aware_sele(fixed),
        }

    os.makedirs(os.path.dirname(args_json), exist_ok=True)
    with open(cfg["output_json"], "w") as f:
        json.dump(result, f, indent=2)
    print(f"Wrote positions JSON: {cfg['output_json']}")


def do_stamp(positions_json, design_id, in_pdb, out_pdb):
    with open(positions_json) as f:
        positions = json.load(f)

    entry = positions.get(design_id)
    if entry is None:
        raise KeyError(f"No positions entry for id {design_id}")
    designed_set = set(sele_to_list(entry["designed"]))

    out_lines = []
    with open(in_pdb) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                chain = field_chain(line)
                res_seq = field_res_seq(line)
                res_name = field_res_name(line)
                # Protein residue in the designed set -> B=0.0; otherwise B=1.0.
                # Ligand HETATM (non-standard resname) is never in the set -> B=1.0.
                is_designed = res_name in STANDARD_RESIDUES and res_seq.isdigit() and \
                    (chain, int(res_seq)) in designed_set
                beta = "  0.00" if is_designed else "  1.00"
                line = line[:60].ljust(60) + beta + line[66:]
                if not line.endswith("\n"):
                    line += "\n"
            out_lines.append(line)

    os.makedirs(os.path.dirname(out_pdb), exist_ok=True)
    with open(out_pdb, "w") as f:
        f.writelines(out_lines)


def main():
    if len(sys.argv) < 2:
        print("Usage: pipe_lasermpnn_bfactor.py <resolve|stamp> ...", file=sys.stderr)
        sys.exit(1)
    mode = sys.argv[1]
    if mode == "resolve":
        do_resolve(sys.argv[2])
    elif mode == "stamp":
        _, _, positions_json, design_id, in_pdb, out_pdb = sys.argv
        do_stamp(positions_json, design_id, in_pdb, out_pdb)
    else:
        raise ValueError(f"Unknown mode: {mode}")


if __name__ == "__main__":
    main()
