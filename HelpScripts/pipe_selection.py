#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Structure-aware selection pipeline for BioPipelines.

Applies a sequence of operations (add, subtract, expand, shrink, shift, invert)
to build and modify PyMOL-formatted selection strings per design ID.

Usage:
    python pipe_selection.py <config_json>

Config JSON:
    - operations: list of {op, refs?, value?} dicts
    - structures_json: (optional) path to DataStream JSON with PDB files
    - output_csv: path to output CSV file
"""

import sys
import os
import json
import pandas as pd
from typing import List, Tuple, Set

# Import unified I/O utilities and selection helpers
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from biopipelines.biopipelines_io import load_datastream, iterate_files
from biopipelines.sele_utils import sele_to_list, list_to_sele
from biopipelines.pdb_parser import parse_pdb_file, STANDARD_RESIDUES


# ── PDB-aware helpers ──

def get_protein_residues(pdb_path: str) -> Set[Tuple[str, int]]:
    """Get set of (chain, resnum) tuples from a PDB/CIF file."""
    atoms = parse_pdb_file(pdb_path)
    residues = set()
    for atom in atoms:
        if atom.res_name in STANDARD_RESIDUES:
            residues.add((atom.chain, atom.res_num))
    if not residues:
        raise ValueError(f"No protein residues found in {pdb_path}")
    return residues


def expand_selection(residues: Set[Tuple[str, int]], expand_by: int,
                     valid_residues: Set[Tuple[str, int]]) -> Set[Tuple[str, int]]:
    """Expand selection by adding *expand_by* residues on each side per chain."""
    if expand_by == 0:
        return residues

    # Group valid residues by chain for fast lookup
    by_chain = {}
    for chain, rnum in valid_residues:
        by_chain.setdefault(chain, sorted(set())).append(rnum)
    for chain in by_chain:
        by_chain[chain] = sorted(set(by_chain[chain]))

    expanded = set(residues)
    for chain, rnum in residues:
        chain_nums = by_chain.get(chain, [])
        try:
            idx = chain_nums.index(rnum)
        except ValueError:
            continue
        lo = max(0, idx - expand_by)
        hi = min(len(chain_nums) - 1, idx + expand_by)
        for i in range(lo, hi + 1):
            expanded.add((chain, chain_nums[i]))
    return expanded


def shrink_selection(residues: Set[Tuple[str, int]], shrink_by: int,
                     valid_residues: Set[Tuple[str, int]]) -> Set[Tuple[str, int]]:
    """Shrink selection by removing *shrink_by* residues from each boundary per chain."""
    if shrink_by == 0:
        return residues

    by_chain = {}
    for chain, rnum in valid_residues:
        by_chain.setdefault(chain, []).append(rnum)
    for chain in by_chain:
        by_chain[chain] = sorted(set(by_chain[chain]))

    # Group selected residues by chain into contiguous intervals
    sel_by_chain = {}
    for chain, rnum in residues:
        sel_by_chain.setdefault(chain, []).append(rnum)

    shrunk = set()
    for chain, nums in sel_by_chain.items():
        chain_nums = by_chain.get(chain, [])
        nums_sorted = sorted(nums)

        # Find contiguous intervals
        intervals = []
        start = nums_sorted[0]
        prev = start
        for n in nums_sorted[1:]:
            # Check if n is the next valid residue after prev
            try:
                prev_idx = chain_nums.index(prev)
            except ValueError:
                prev = n
                start = n
                continue
            if prev_idx + 1 < len(chain_nums) and chain_nums[prev_idx + 1] == n:
                prev = n
            else:
                intervals.append((start, prev))
                start = n
                prev = n
        intervals.append((start, prev))

        # Shrink each interval
        for s, e in intervals:
            try:
                s_idx = chain_nums.index(s)
                e_idx = chain_nums.index(e)
            except ValueError:
                continue
            new_s_idx = s_idx + shrink_by
            new_e_idx = e_idx - shrink_by
            if new_s_idx > new_e_idx:
                continue
            for i in range(new_s_idx, new_e_idx + 1):
                shrunk.add((chain, chain_nums[i]))
    return shrunk


def shift_selection(residues: Set[Tuple[str, int]], shift_by: int,
                    valid_residues: Set[Tuple[str, int]]) -> Set[Tuple[str, int]]:
    """Shift selection by *shift_by* positions in the valid residue sequence per chain."""
    if shift_by == 0:
        return residues

    by_chain = {}
    for chain, rnum in valid_residues:
        by_chain.setdefault(chain, []).append(rnum)
    for chain in by_chain:
        by_chain[chain] = sorted(set(by_chain[chain]))

    shifted = set()
    for chain, rnum in residues:
        chain_nums = by_chain.get(chain, [])
        try:
            idx = chain_nums.index(rnum)
        except ValueError:
            continue
        new_idx = idx + shift_by
        if 0 <= new_idx < len(chain_nums):
            shifted.add((chain, chain_nums[new_idx]))
    return shifted


def invert_selection(residues: Set[Tuple[str, int]],
                     valid_residues: Set[Tuple[str, int]]) -> Set[Tuple[str, int]]:
    """Replace selection with its complement."""
    return valid_residues - residues


# ── Table lookup ──

def lookup_column_values(table_path: str, column: str) -> dict:
    """Load a CSV table and return {id: value} mapping for the given column."""
    df = pd.read_csv(table_path)
    if column not in df.columns:
        raise ValueError(f"Column '{column}' not found in {table_path}. Available: {list(df.columns)}")
    return dict(zip(df["id"].astype(str), df[column].astype(str)))


# ── Main pipeline ──

def process_selections(config_path: str):
    """Apply sequential operations to build selections per ID."""
    with open(config_path, 'r') as f:
        config = json.load(f)

    operations = config["operations"]
    output_csv = config["output_csv"]

    # Load structures if available
    struct_map = {}
    if "structures_json" in config:
        structures_ds = load_datastream(config["structures_json"])
        for sid, spath in iterate_files(structures_ds):
            struct_map[sid] = spath

    # Collect all IDs from the first add/subtract operation's table refs
    all_ids = None
    for op in operations:
        if op["op"] in ("add", "subtract"):
            for ref in op["refs"]:
                vals = lookup_column_values(ref["table"], ref["column"])
                if all_ids is None:
                    all_ids = list(vals.keys())
                break
            if all_ids is not None:
                break

    if all_ids is None:
        raise ValueError("No add/subtract operation found — cannot determine IDs")

    print(f"Processing {len(all_ids)} IDs through {len(operations)} operations")

    # Pre-load all referenced tables
    ref_cache = {}
    for op in operations:
        if op["op"] in ("add", "subtract"):
            for ref in op["refs"]:
                key = (ref["table"], ref["column"])
                if key not in ref_cache:
                    ref_cache[key] = lookup_column_values(ref["table"], ref["column"])

    # Cache for PDB residues (avoid re-parsing)
    pdb_residue_cache = {}

    results = []
    for design_id in all_ids:
        current = set()  # set of (chain, resnum) tuples

        for op in operations:
            op_type = op["op"]

            if op_type == "add":
                for ref in op["refs"]:
                    vals = ref_cache[(ref["table"], ref["column"])]
                    value = vals.get(design_id, "")
                    current |= set(sele_to_list(value))

            elif op_type == "subtract":
                for ref in op["refs"]:
                    vals = ref_cache[(ref["table"], ref["column"])]
                    value = vals.get(design_id, "")
                    current -= set(sele_to_list(value))

            elif op_type in ("expand", "shrink", "shift", "invert"):
                pdb_path = struct_map.get(design_id)
                if not pdb_path or not os.path.exists(pdb_path):
                    print(f"Warning: No PDB found for '{design_id}', skipping PDB-aware op '{op_type}'")
                    continue

                if pdb_path not in pdb_residue_cache:
                    pdb_residue_cache[pdb_path] = get_protein_residues(pdb_path)
                valid = pdb_residue_cache[pdb_path]

                n = op.get("value", 0)
                if op_type == "expand":
                    current = expand_selection(current, n, valid)
                elif op_type == "shrink":
                    current = shrink_selection(current, n, valid)
                elif op_type == "shift":
                    current = shift_selection(current, n, valid)
                elif op_type == "invert":
                    current = invert_selection(current, valid)

        sele_str = list_to_sele(sorted(current, key=lambda x: (x[0], x[1])))
        results.append({"id": design_id, "selection": sele_str})
        print(f"  {design_id}: {sele_str}")

    if not results:
        raise ValueError("No selections were produced")

    result_df = pd.DataFrame(results)
    os.makedirs(os.path.dirname(output_csv), exist_ok=True)
    result_df.to_csv(output_csv, index=False)

    print(f"\nSelection pipeline completed — {len(results)} entries saved to {output_csv}")


def main():
    if len(sys.argv) != 2:
        print("Usage: python pipe_selection.py <config_json>")
        sys.exit(1)

    config_path = sys.argv[1]
    if not os.path.exists(config_path):
        print(f"Error: Config file not found: {config_path}")
        sys.exit(1)

    try:
        process_selections(config_path)
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
