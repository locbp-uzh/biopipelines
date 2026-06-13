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
    - operations: list of {op, refs?, value?, stream_json?, filter_expr?, filter_column?} dicts
    - structures_json: (optional) path to DataStream JSON with PDB files
    - output_csv: path to output CSV file

Stream-based operations (add/subtract with stream_json):
    The stream must have format resi-csv with a `resi` column
    and a numeric column specified in the filter expression.
    filter_expr syntax: "column op value", e.g. "propensity>0.5", "rmsf<=1.2"
"""

import sys
import os
import re
import json
import operator
import pandas as pd
from typing import List, Tuple, Set

# Import unified I/O utilities and selection helpers
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from biopipelines.biopipelines_io import load_datastream, iterate_files
from biopipelines.sele_utils import sele_to_list, list_to_sele
from biopipelines.pdb_parser import parse_pdb_file, STANDARD_RESIDUES
from biopipelines.id_map_utils import get_mapped_ids


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
                     valid_residues: Set[Tuple[str, int]],
                     direction: str = "nc") -> Set[Tuple[str, int]]:
    """Expand selection by adding *expand_by* residues per chain.

    Direction: "n" = toward N-term (lower resnums), "c" = toward C-term, "nc" = both.
    """
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
        lo = max(0, idx - expand_by) if direction in ("n", "nc") else idx
        hi = min(len(chain_nums) - 1, idx + expand_by) if direction in ("c", "nc") else idx
        for i in range(lo, hi + 1):
            expanded.add((chain, chain_nums[i]))
    return expanded


def shrink_selection(residues: Set[Tuple[str, int]], shrink_by: int,
                     valid_residues: Set[Tuple[str, int]],
                     direction: str = "nc") -> Set[Tuple[str, int]]:
    """Shrink selection by removing *shrink_by* residues from boundaries per chain.

    Direction: "n" = remove from N-term side, "c" = from C-term side, "nc" = both.
    """
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
            new_s_idx = s_idx + shrink_by if direction in ("n", "nc") else s_idx
            new_e_idx = e_idx - shrink_by if direction in ("c", "nc") else e_idx
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


# ── Predefined structural pattern helpers ──

def select_n_terminus(valid_residues: Set[Tuple[str, int]], n: int) -> Set[Tuple[str, int]]:
    """Select the first n residues of each chain."""
    by_chain = {}
    for chain, resnum in valid_residues:
        by_chain.setdefault(chain, []).append(resnum)
    result = set()
    for chain, nums in by_chain.items():
        for resnum in sorted(nums)[:n]:
            result.add((chain, resnum))
    return result


def select_c_terminus(valid_residues: Set[Tuple[str, int]], n: int) -> Set[Tuple[str, int]]:
    """Select the last n residues of each chain."""
    by_chain = {}
    for chain, resnum in valid_residues:
        by_chain.setdefault(chain, []).append(resnum)
    result = set()
    for chain, nums in by_chain.items():
        for resnum in sorted(nums)[-n:]:
            result.add((chain, resnum))
    return result


def select_gaps(valid_residues: Set[Tuple[str, int]]) -> Set[Tuple[str, int]]:
    """Select one residue on each side of each structural gap (missing residues)."""
    by_chain = {}
    for chain, resnum in valid_residues:
        by_chain.setdefault(chain, []).append(resnum)
    result = set()
    for chain, nums in by_chain.items():
        sorted_nums = sorted(nums)
        for i in range(len(sorted_nums) - 1):
            if sorted_nums[i + 1] > sorted_nums[i] + 1:
                result.add((chain, sorted_nums[i]))       # C-side of gap
                result.add((chain, sorted_nums[i + 1]))   # N-side of gap
    return result


# ── Table lookup ──

def lookup_column_values(table_path: str, column: str) -> dict:
    """Load a CSV table and return {id: value} mapping for the given column."""
    df = pd.read_csv(table_path)
    if column not in df.columns:
        raise ValueError(f"Column '{column}' not found in {table_path}. Available: {list(df.columns)}")
    return dict(zip(df["id"].astype(str), df[column].astype(str)))


# ── Stream-based residue lookup ──

_OPS = {">": operator.gt, "<": operator.lt, ">=": operator.ge, "<=": operator.le,
        "==": operator.eq, "!=": operator.ne}

def _parse_filter_expr(expr: str):
    """Parse a filter expression like 'propensity>0.5' into (column, operator_fn, threshold)."""
    m = re.match(r"^(\w+)\s*(>=|<=|!=|>|<|==)\s*(.+)$", expr.strip())
    if not m:
        raise ValueError(
            f"Invalid filter expression: '{expr}'. "
            f"Expected 'column op value', e.g. 'propensity>0.5', 'rmsf<=1.2'"
        )
    col, op_str, val_str = m.group(1), m.group(2), m.group(3)
    return col, _OPS[op_str], float(val_str)


def load_stream_residues(stream_json: str, filter_expr: str) -> dict:
    """
    Load a per-residue DataStream and return {id: set of (chain, resi)} after filtering.

    Args:
        stream_json:   Path to the DataStream JSON file.
        filter_expr:   Filter expression, e.g. ">0.5".
        filter_column: Column to apply the filter on.

    Returns:
        dict mapping structure id -> set of (chain, resi) tuples.
        Note: per-residue CSVs may not have a chain column; if absent, uses chain "A".
    """
    ds = load_datastream(stream_json)
    filter_column, op_fn, threshold = _parse_filter_expr(filter_expr)

    residues_by_id = {}
    for sid, csv_path in iterate_files(ds):
        if not csv_path or not os.path.exists(csv_path):
            print(f"Warning: per-residue CSV not found for {sid}: {csv_path}", file=sys.stderr)
            residues_by_id[sid] = set()
            continue
        df = pd.read_csv(csv_path)
        df.columns = [c.strip().lower().replace(" ", "_") for c in df.columns]
        if filter_column not in df.columns:
            raise ValueError(
                f"Filter column '{filter_column}' not found in {csv_path}. "
                f"Available: {list(df.columns)}"
            )
        if "resi" not in df.columns:
            raise ValueError(f"'resi' column not found in {csv_path}")
        mask = op_fn(df[filter_column].astype(float), threshold)
        filtered = df[mask]
        chain_col = "chain" if "chain" in df.columns else None
        res_set = set()
        for _, row in filtered.iterrows():
            chain = str(row[chain_col]) if chain_col else ""
            res_set.add((chain, int(row["resi"])))
        residues_by_id[sid] = res_set
    return residues_by_id


def remap_chainless_residues(residues, valid_residues):
    """Remap ("", resnum) tuples to the correct chain using PDB data."""
    if not any(c == "" for c, _ in residues):
        return residues
    resnum_to_chain = {}
    for chain, resnum in valid_residues:
        resnum_to_chain.setdefault(resnum, chain)
    remapped = set()
    for chain, resnum in residues:
        if chain == "":
            actual_chain = resnum_to_chain.get(resnum)
            if actual_chain is not None:
                remapped.add((actual_chain, resnum))
        else:
            remapped.add((chain, resnum))
    return remapped


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

    # Output ids are the UNION of every ADD op's input ids (not just the first
    # op's): a second add() in a different id space — e.g. a crystal pocket
    # alongside a predicted-pose consensus — must contribute its own ids too,
    # otherwise its structures silently vanish. subtract ops only REMOVE residues
    # from existing ids, so they must NOT introduce new output ids (that would
    # leak the subtracted input's id space into the result).
    all_ids = []
    seen = set()
    def _extend(ids):
        for i in ids:
            if i not in seen:
                seen.add(i); all_ids.append(i)
    for op in operations:
        if op["op"] != "add":
            continue
        for ref in op.get("refs", []):
            _extend(lookup_column_values(ref["table"], ref["column"]).keys())
        if op.get("stream_json"):
            _extend(load_datastream(op["stream_json"]).ids_expanded)

    # Predefined pattern ops (n_terminus, c_terminus, etc.) have no refs/stream,
    # so fall back to the structures map for IDs.
    if not all_ids and struct_map:
        all_ids = list(struct_map.keys())

    if not all_ids:
        raise ValueError("Cannot determine IDs: no add/subtract operation and no structures provided")

    print(f"Processing {len(all_ids)} IDs through {len(operations)} operations")

    # Resolve each selection id to its structure id via framework id matching
    # (exact, else parent/child) — a selection keyed by pose ids still finds the
    # single best33AA structure it derives from.
    id_to_struct = {}
    if struct_map:
        id_to_struct = get_mapped_ids(all_ids, list(struct_map.keys()), unique=True)

    # Each op input is keyed in its OWN id space (a predictor stream may be keyed
    # by pooled/pose ids while another input is keyed by a group id). For every
    # input, map each output id to that input's matching key via framework id
    # matching, so a per-id lookup resolves across id spaces (e.g. subtract a
    # pose-keyed predictor from a group-keyed truth).
    def _key_map(input_keys):
        ik = list(input_keys)
        m = {i: i for i in all_ids if i in ik}      # exact first
        unresolved = [i for i in all_ids if i not in m]
        if unresolved and ik:
            for did, k in get_mapped_ids(unresolved, ik, unique=True).items():
                if k:
                    m[did] = k
        return m

    # Pre-load all referenced tables
    ref_cache = {}
    ref_keymap = {}
    for op in operations:
        if op["op"] in ("add", "subtract", "intersect"):
            for ref in op.get("refs", []):
                key = (ref["table"], ref["column"])
                if key not in ref_cache:
                    ref_cache[key] = lookup_column_values(ref["table"], ref["column"])
                    ref_keymap[key] = _key_map(ref_cache[key].keys())

    # Pre-load all stream-based residue sets
    stream_cache = {}
    stream_keymap = {}
    for op in operations:
        if op["op"] in ("add", "subtract", "intersect") and op.get("stream_json"):
            filter_expr = op.get("filter_expr")
            if filter_expr is None:
                raise ValueError(
                    f"Stream-based {op['op']} operation requires include= or exclude= "
                    f"(e.g. include=\"propensity>0.5\")"
                )
            key = (op["stream_json"], filter_expr)
            if key not in stream_cache:
                stream_cache[key] = load_stream_residues(op["stream_json"], filter_expr)
                stream_keymap[key] = _key_map(stream_cache[key].keys())

    # Cache for PDB residues (avoid re-parsing)
    pdb_residue_cache = {}

    results = []
    for design_id in all_ids:
        current = set()  # set of (chain, resnum) tuples

        for op in operations:
            op_type = op["op"]

            if op_type == "add":
                for ref in op.get("refs", []):
                    k = (ref["table"], ref["column"])
                    src = ref_keymap[k].get(design_id, design_id)
                    current |= set(sele_to_list(ref_cache[k].get(src, "")))
                if op.get("stream_json"):
                    key = (op["stream_json"], op["filter_expr"])
                    src = stream_keymap[key].get(design_id, design_id)
                    current |= stream_cache[key].get(src, set())

            elif op_type == "subtract":
                for ref in op.get("refs", []):
                    k = (ref["table"], ref["column"])
                    src = ref_keymap[k].get(design_id, design_id)
                    current -= set(sele_to_list(ref_cache[k].get(src, "")))
                if op.get("stream_json"):
                    key = (op["stream_json"], op["filter_expr"])
                    src = stream_keymap[key].get(design_id, design_id)
                    current -= stream_cache[key].get(src, set())

            elif op_type == "intersect":
                for ref in op.get("refs", []):
                    k = (ref["table"], ref["column"])
                    src = ref_keymap[k].get(design_id, design_id)
                    current &= set(sele_to_list(ref_cache[k].get(src, "")))
                if op.get("stream_json"):
                    key = (op["stream_json"], op["filter_expr"])
                    src = stream_keymap[key].get(design_id, design_id)
                    current &= stream_cache[key].get(src, set())

            elif op_type in ("expand", "shrink", "shift", "invert"):
                struct_id = id_to_struct.get(design_id) or design_id
                pdb_path = struct_map.get(struct_id)
                if not pdb_path or not os.path.exists(pdb_path):
                    print(f"Warning: No PDB found for '{design_id}', skipping PDB-aware op '{op_type}'")
                    continue

                if pdb_path not in pdb_residue_cache:
                    pdb_residue_cache[pdb_path] = get_protein_residues(pdb_path)
                valid = pdb_residue_cache[pdb_path]

                # Remap chainless residues (from streams without chain column)
                current = remap_chainless_residues(current, valid)

                n = op.get("value", 0)
                direction = op.get("direction", "nc")
                if op_type == "expand":
                    current = expand_selection(current, n, valid, direction)
                elif op_type == "shrink":
                    current = shrink_selection(current, n, valid, direction)
                elif op_type == "shift":
                    current = shift_selection(current, n, valid)
                elif op_type == "invert":
                    current = invert_selection(current, valid)

            elif op_type in ("n_terminus", "c_terminus", "termini", "all_residues", "gaps"):
                pdb_path = struct_map.get(design_id)
                if not pdb_path or not os.path.exists(pdb_path):
                    print(f"Warning: No PDB found for '{design_id}', skipping op '{op_type}'")
                    continue

                if pdb_path not in pdb_residue_cache:
                    pdb_residue_cache[pdb_path] = get_protein_residues(pdb_path)
                valid = pdb_residue_cache[pdb_path]

                current = remap_chainless_residues(current, valid)

                n = op.get("value", 1)
                if op_type == "n_terminus":
                    current |= select_n_terminus(valid, n)
                elif op_type == "c_terminus":
                    current |= select_c_terminus(valid, n)
                elif op_type == "termini":
                    current |= select_n_terminus(valid, n)
                    current |= select_c_terminus(valid, n)
                elif op_type == "all_residues":
                    current |= valid
                elif op_type == "gaps":
                    current |= select_gaps(valid)

        sele_str = list_to_sele(sorted(current, key=lambda x: (x[0], x[1])))
        results.append({"id": design_id, "selection": sele_str,
                        "n_residues": len(current)})
        print(f"  {design_id}: {sele_str} ({len(current)} residues)")

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
