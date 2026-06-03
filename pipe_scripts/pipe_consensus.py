#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Per-group aggregation of a resi-csv stream.

Reads an N-file resi-csv stream (one CSV per id, rows keyed by (chain, resi)),
partitions the ids into groups by matching against a groups stream, aggregates
per (chain, resi) across each group, and writes the result back as a resi-csv
stream with one file per group, keyed by the group id.

Usage:
    python pipe_consensus.py <config_json>

Config JSON:
    - stream_json:       DataStream JSON for the input resi-csv stream
    - groups_json:       DataStream JSON whose ids define the partition
    - operations:        list of {type, params} aggregation ops
    - consensus_dir:     output directory for per-id resi-csv files
    - consensus_map_csv: output map_table path (id, file)
"""

import sys
import os
import re
import json
import operator
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from biopipelines.biopipelines_io import load_datastream, iterate_files
from biopipelines.id_map_utils import get_mapped_ids


_OPS = {">": operator.gt, "<": operator.lt, ">=": operator.ge, "<=": operator.le,
        "==": operator.eq, "!=": operator.ne}


def _parse_predicate(expr):
    """Parse 'distance<=6' into (column, operator_fn, threshold)."""
    m = re.match(r"^(\w+)\s*(>=|<=|!=|>|<|==)\s*(.+)$", str(expr).strip())
    if not m:
        raise ValueError(
            f"Invalid predicate: '{expr}'. Expected 'column op value', e.g. 'distance<=6'"
        )
    return m.group(1), _OPS[m.group(2)], float(m.group(3))


def _load_rows(stream_json):
    """Return {id: DataFrame} for every per-id resi-csv that exists."""
    ds = load_datastream(stream_json)
    rows_by_id = {}
    for sid, csv_path in iterate_files(ds):
        if not csv_path or not os.path.exists(csv_path):
            print(f"Warning: resi-csv not found for {sid}: {csv_path}", file=sys.stderr)
            continue
        df = pd.read_csv(csv_path)
        df.columns = [c.strip().lower().replace(" ", "_") for c in df.columns]
        if "resi" not in df.columns:
            raise ValueError(f"'resi' column not found in {csv_path}")
        if "chain" not in df.columns:
            df["chain"] = ""
        rows_by_id[sid] = df
    return ds, rows_by_id


def _collect_group_cells(member_dfs):
    """Return {(chain, resi): [row dicts]} pooled across a group's member frames."""
    cells = {}
    for df in member_dfs:
        for _, row in df.iterrows():
            key = (str(row["chain"]), int(row["resi"]))
            cells.setdefault(key, []).append(row)
    return cells


def _numeric_cell(row, col):
    """Coerce row[col] to float, or None if absent / non-numeric / NaN."""
    val = pd.to_numeric(row.get(col), errors="coerce")
    return None if pd.isna(val) else float(val)


def _apply_op(op, rows, n_group):
    """Compute one aggregate value for a (chain, resi) cell."""
    op_type = op["type"]
    params = op["params"]
    if op_type == "fraction":
        col, op_fn, thr = _parse_predicate(params["predicate"])
        # NaN/missing cells fail the predicate; denominator stays n_group.
        n_pass = sum(1 for r in rows if (v := _numeric_cell(r, col)) is not None and op_fn(v, thr))
        return round(n_pass / n_group, 6) if n_group else 0.0
    if op_type == "count":
        return len(rows)
    col = params["column"]
    vals = [v for r in rows if (v := _numeric_cell(r, col)) is not None]
    if not vals:
        return ""
    if op_type == "mean":
        return round(sum(vals) / len(vals), 6)
    if op_type == "min":
        return round(min(vals), 6)
    if op_type == "max":
        return round(max(vals), 6)
    raise ValueError(f"unknown operation type '{op_type}'")


def main():
    if len(sys.argv) != 2:
        print("Usage: python pipe_consensus.py <config_json>")
        sys.exit(1)

    with open(sys.argv[1]) as f:
        config = json.load(f)

    operations = config["operations"]
    consensus_dir = config["consensus_dir"]
    consensus_map_csv = config["consensus_map_csv"]

    ds, rows_by_id = _load_rows(config["stream_json"])
    pose_ids = list(ds.ids_expanded)

    groups_ds = load_datastream(config["groups_json"])
    group_ids = list(groups_ds.ids_expanded)

    # {group_id: [pose ids that belong to it]} via framework id matching.
    mapped = get_mapped_ids(group_ids, pose_ids, unique=False)

    # Aggregate each group once.
    out_names = [op["params"]["name"] for op in operations]
    group_rows = {}   # group_id -> list of output row dicts
    for gid, members in mapped.items():
        member_dfs = [rows_by_id[pid] for pid in members if pid in rows_by_id]
        n_group = len(member_dfs)
        cells = _collect_group_cells(member_dfs)
        rows = []
        for (chain, resi) in sorted(cells.keys(), key=lambda x: (x[0], x[1])):
            cell_rows = cells[(chain, resi)]
            out = {"chain": chain, "resi": resi}
            for op in operations:
                out[op["params"]["name"]] = _apply_op(op, cell_rows, n_group)
            out["n_group"] = n_group
            rows.append(out)
        group_rows[gid] = rows

    os.makedirs(consensus_dir, exist_ok=True)
    columns = ["id", "chain", "resi"] + out_names + ["n_group"]
    map_rows = []
    written = 0
    # One file per GROUP, keyed by the group id (matches the groups stream ids).
    for gid in group_ids:
        rows = group_rows.get(gid, [])
        if not rows:
            print(f"Warning: group '{gid}' has no member poses; writing empty consensus", file=sys.stderr)
        out_path = os.path.join(consensus_dir, f"{gid}.csv")
        df_rows = [dict({"id": gid}, **r) for r in rows]
        pd.DataFrame(df_rows, columns=columns).to_csv(out_path, index=False)
        map_rows.append({"id": gid, "file": out_path})
        written += 1

    os.makedirs(os.path.dirname(consensus_map_csv), exist_ok=True)
    pd.DataFrame(map_rows, columns=["id", "file"]).to_csv(consensus_map_csv, index=False)

    print(f"Consensus: {written} files written; {len(group_rows)} groups")
    if not written:
        sys.exit(1)


if __name__ == "__main__":
    main()
