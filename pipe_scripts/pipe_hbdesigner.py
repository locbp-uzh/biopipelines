#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Collect HBDesigner per-input designs into the structures stream + summary table.

HBDesigner is invoked once per input structure into its own ``<exec_root>/<id>/``
dir. Each run writes ``<prefix>_HBDes_rank_<n>.pdb`` files plus a single
``<prefix>_HBDes_stats.csv`` (columns include Rank, 1-based best-first, and
Output_PDB), where ``<prefix>`` is the input PDB's basename — NOT the stream id.
This script, for every input id:
  * reads the single ``*_HBDes_stats.csv`` in that per-id dir (located by glob,
    since the prefix is the PDB basename, not the id), orders designs by
    ``Rank`` and links each to its PDB by ``Output_PDB`` (keyed on the upstream
    schema, not filename sort order),
  * copies the top_k retained designs into the structures stream folder with
    the framework ids ``<id>_1..<id>_top_k`` (rank order),
  * extracts each design's one-letter sequence (HBDesigner emits no FASTA) into
    the content-bearing sequences stream CSV,
  * accumulates the summary rows, tagged with a ``structures.id`` provenance
    column and the framework ``id``,
  * records ranks not produced — HBDesigner saves min(top_k, #valid), so a run
    may yield fewer than top_k — in the missing table, so the completion check
    excuses those ids instead of flagging a valid partial run.

The structures stream's map_table is written separately by
pipe_update_structures_map.py from the copied PDBs; this script owns the PDB
copy, the summary, sequences, and missing tables. Per-input failures are skipped
and reported; the script errors only if every input failed.
"""

import argparse
import glob
import os
import shutil
import sys

import pandas as pd

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from biopipelines.biopipelines_io import load_datastream, iterate_files
from biopipelines.pdb_parser import parse_pdb_file, get_protein_sequence


def _pdb_sequence(pdb_path):
    """One-letter sequence of a designed PDB, chains joined in order."""
    chains = get_protein_sequence(parse_pdb_file(pdb_path))
    return "".join(chains[c] for c in sorted(chains))

def _ordered_designs(run_dir, summary_df):
    """Return [(pdb_path, summary_row)] ordered best-first (rank 1 first).

    HBDesigner writes one ``<prefix>_HBDes_stats.csv`` per run with a ``Rank``
    column (1-based, best first) and an ``Output_PDB`` column naming each row's
    PDB. Order by ``Rank`` and link by ``Output_PDB`` — keying on the upstream
    schema directly, not on filename sort order (which misranks ``rank_10``
    before ``rank_2``).
    """
    pdbs = glob.glob(os.path.join(run_dir, "*.pdb"))
    by_name = {os.path.basename(p): p for p in pdbs}

    if summary_df is None or len(summary_df) == 0:
        raise ValueError(f"{run_dir}: no stats CSV / empty summary")
    if "Rank" not in summary_df.columns or "Output_PDB" not in summary_df.columns:
        raise ValueError(
            f"{run_dir}: stats CSV missing expected Rank/Output_PDB columns "
            f"(got {list(summary_df.columns)})"
        )

    ordered = []
    for _, row in summary_df.sort_values("Rank", kind="stable").iterrows():
        path = by_name.get(os.path.basename(str(row["Output_PDB"])))
        if path is not None:
            ordered.append((path, row))
    return ordered


def main():
    parser = argparse.ArgumentParser(description="Collect HBDesigner designs")
    parser.add_argument("--structures-json", required=True)
    parser.add_argument("--exec-root", required=True)
    parser.add_argument("--structures-dir", required=True)
    parser.add_argument("--summary-csv", required=True)
    parser.add_argument("--sequences-csv", required=True)
    parser.add_argument("--missing-csv", required=True)
    parser.add_argument("--top-k", type=int, required=True)
    args = parser.parse_args()

    ds = load_datastream(args.structures_json)
    os.makedirs(args.structures_dir, exist_ok=True)

    summary_rows = []
    sequence_rows = []
    missing_rows = []   # expected <struct_id>_<rank> ids that didn't materialize
    failed = []         # inputs that produced nothing at all
    n_inputs = 0

    for struct_id, _ in iterate_files(ds):
        n_inputs += 1
        run_dir = os.path.join(args.exec_root, struct_id)
        # Glob the per-id dir: upstream names the stats CSV from the PDB basename, not the stream id.
        crashed = os.path.exists(os.path.join(run_dir, ".hbdes_failed"))
        stats_glob = glob.glob(os.path.join(run_dir, "*_HBDes_stats.csv"))
        summary_df = None
        if stats_glob:
            try:
                summary_df = pd.read_csv(stats_glob[0])
            except Exception as e:
                print(f"WARNING: {struct_id}: could not read {stats_glob[0]}: {e}", file=sys.stderr)

        try:
            ordered = _ordered_designs(run_dir, summary_df)
        except ValueError as e:
            print(f"WARNING: {struct_id}: {e}", file=sys.stderr)
            ordered = []

        keep = ordered[:args.top_k]
        if not keep:
            # Zero designs is a failure, not an intentional filter.
            print(f"WARNING: {struct_id}: no designed PDBs in {run_dir}", file=sys.stderr)
            failed.append(struct_id)

        for rank, (pdb_path, row) in enumerate(keep, start=1):
            out_id = f"{struct_id}_{rank}"
            # Extract sequence first; on failure drop the whole design so the three streams stay consistent.
            try:
                seq = _pdb_sequence(pdb_path)
            except Exception as e:
                print(f"WARNING: {out_id}: could not extract sequence, dropping: {e}",
                      file=sys.stderr)
                missing_rows.append({"id": out_id, "removed_by": "HBDesigner",
                                     "kind": "failure",
                                     "cause": f"sequence extraction failed: {str(e)[:120]}"})
                continue
            dest = os.path.join(args.structures_dir, f"{out_id}.pdb")
            shutil.copyfile(pdb_path, dest)
            sequence_rows.append({"id": out_id, "structures.id": struct_id, "sequence": seq})
            entry = {"id": out_id, "structures.id": struct_id}
            for col, val in row.items():
                # Skip our own keys and pandas' unnamed index column (the stats
                # CSV is written with a leading index, read back as "Unnamed: 0").
                if col in ("id", "structures.id") or str(col).startswith("Unnamed:"):
                    continue
                entry[col] = val
            summary_rows.append(entry)

        produced = len(keep)
        # Crash / zero designs -> failure (not excused); fewer-than-top_k -> filter (excused).
        fail = crashed or produced == 0
        for rank in range(produced + 1, args.top_k + 1):
            if fail:
                missing_rows.append({"id": f"{struct_id}_{rank}", "removed_by": "HBDesigner",
                                     "kind": "failure",
                                     "cause": "run_hbdesigner crashed" if crashed
                                     else "produced no designs"})
            else:
                missing_rows.append({"id": f"{struct_id}_{rank}", "removed_by": "HBDesigner",
                                     "kind": "filter",
                                     "cause": "fewer valid designs than top_k"})

    os.makedirs(os.path.dirname(args.summary_csv), exist_ok=True)
    summary_cols = ["id", "structures.id"]
    if summary_rows:
        df = pd.DataFrame(summary_rows)
        cols = summary_cols + [c for c in df.columns if c not in summary_cols]
        df[cols].to_csv(args.summary_csv, index=False)
        print(f"Collected {len(summary_rows)} design(s) from {n_inputs - len(failed)}/{n_inputs} input(s)")
    else:
        pd.DataFrame(columns=summary_cols).to_csv(args.summary_csv, index=False)

    os.makedirs(os.path.dirname(args.sequences_csv), exist_ok=True)
    pd.DataFrame(sequence_rows, columns=["id", "structures.id", "sequence"]).to_csv(
        args.sequences_csv, index=False)

    os.makedirs(os.path.dirname(args.missing_csv), exist_ok=True)
    pd.DataFrame(missing_rows, columns=["id", "removed_by", "kind", "cause"]).to_csv(
        args.missing_csv, index=False)

    if failed:
        print(f"Failed {len(failed)}/{n_inputs}: {failed}", file=sys.stderr)
    if not summary_rows:
        sys.exit(1)


if __name__ == "__main__":
    main()
