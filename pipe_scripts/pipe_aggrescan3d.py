#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Runtime helper for Aggrescan3D (A3D). Both phases run under biopipelines
(Python 3); the A3D binary itself runs in its own Python-2.7 env in between.

  * --emit-worklist: resolve the input structures stream to id<TAB>file lines
    for the bash loop. Runs BEFORE the A3D py2.7 segment so resolve happens
    under Python 3 (resolve_stream_ids would crash under py2.7).
  * default (post-process): parse each A3D work dir into the stream/table layout:
      - A3D.csv (cols protein, chain, residue, residue_name, score) -> per-residue
        resi-csv (id, chain, resi, score) and a global per-structure summary row.
      - output.pdb (B-factor field = A3D score) -> structures stream.
      - <chainID>.png score plots -> images stream (prefixed with the structure id).

Per-item failures are skipped; partial results are still written.
"""

import os
import sys
import glob
import shutil
import argparse

import pandas as pd

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from biopipelines.biopipelines_io import load_datastream, iterate_files


def emit_worklist(structures_json, worklist_path):
    """Resolve the input stream to id<TAB>file lines for the bash loop.

    Runs under biopipelines (Python 3) before the A3D Python-2.7 segment, so
    lazy IDs expand against the map_table at runtime (never at config time).
    """
    ds = load_datastream(structures_json)
    rows = list(iterate_files(ds))
    if not rows:
        print("Error: no structures resolved from input DataStream", file=sys.stderr)
        sys.exit(1)
    with open(worklist_path, "w") as f:
        for struct_id, struct_file in rows:
            f.write(f"{struct_id}\t{struct_file}\n")
    print(f"Wrote worklist: {worklist_path} ({len(rows)} structures)")


def parse_a3d_csv(a3d_path, struct_id):
    """Parse A3D.csv into per-residue rows: id, chain, resi, score.

    A3D writes lowercase headers: protein, chain, residue, residue_name, score.
    """
    df = pd.read_csv(a3d_path)
    rows = []
    for _, r in df.iterrows():
        rows.append({
            "id": struct_id,
            "chain": str(r["chain"]),
            "resi": int(r["residue"]),
            "score": float(r["score"]),
        })
    return rows


def postprocess(struct_id, work_root, structures_dir, images_dir, aggregation_dir):
    """Collect A3D artefacts for one structure. Returns (resi_rows, summary_row)
    or (None, None) if the A3D run produced no parsable scores."""
    work_dir = os.path.join(work_root, struct_id)
    a3d_csv = os.path.join(work_dir, "A3D.csv")

    if not os.path.exists(a3d_csv):
        print(f"Warning: no A3D.csv for {struct_id} (expected {a3d_csv})", file=sys.stderr)
        return None, None

    resi_rows = parse_a3d_csv(a3d_csv, struct_id)
    if not resi_rows:
        print(f"Warning: A3D.csv empty for {struct_id}", file=sys.stderr)
        return None, None

    # Per-residue resi-csv for this structure.
    resi_df = pd.DataFrame(resi_rows, columns=["id", "chain", "resi", "score"])
    resi_out = os.path.join(aggregation_dir, f"{struct_id}_A3D.csv")
    resi_df.to_csv(resi_out, index=False)

    # Scored output PDB (B-factor = A3D score).
    out_pdb = os.path.join(work_dir, "output.pdb")
    if os.path.exists(out_pdb):
        shutil.copy2(out_pdb, os.path.join(structures_dir, f"{struct_id}.pdb"))
    else:
        print(f"Warning: no output.pdb for {struct_id}", file=sys.stderr)

    # Per-chain score plots; A3D names them <chainID>.png. Prefix with the
    # structure id so ids stay unique across inputs. The declared images stream
    # uses one <id>.png per structure, so keep the first chain plot under that
    # id and copy any extras with a chain suffix.
    pngs = sorted(glob.glob(os.path.join(work_dir, "*.png")))
    for i, png in enumerate(pngs):
        if i == 0:
            dst = os.path.join(images_dir, f"{struct_id}.png")
        else:
            chain = os.path.splitext(os.path.basename(png))[0]
            dst = os.path.join(images_dir, f"{struct_id}_{chain}.png")
        shutil.copy2(png, dst)

    scores = resi_df["score"]
    summary_row = {
        "id": struct_id,
        "structures.id": struct_id,
        "avg_score": float(scores.mean()),
        "min_score": float(scores.min()),
        "max_score": float(scores.max()),
        "n_residues": int(len(scores)),
        "n_aggregation_prone": int((scores > 0).sum()),
    }
    print(f"{struct_id}: {summary_row['n_residues']} residues, "
          f"avg {summary_row['avg_score']:.4f}, "
          f"{summary_row['n_aggregation_prone']} aggregation-prone")
    return resi_rows, summary_row


def main():
    parser = argparse.ArgumentParser(description="Aggrescan3D runtime helper")
    parser.add_argument("--structures", required=True, help="Input DataStream JSON")
    parser.add_argument("--emit-worklist", action="store_true",
                        help="Resolve the stream to id<TAB>file lines and exit (Python 3 segment)")
    parser.add_argument("--worklist", help="Output worklist TSV path (with --emit-worklist)")
    parser.add_argument("--work_root", help="Folder with per-structure A3D work dirs")
    parser.add_argument("--structures_dir", help="Destination for output PDBs")
    parser.add_argument("--images_dir", help="Destination for score PNGs")
    parser.add_argument("--aggregation_dir", help="Destination for per-residue CSVs")
    parser.add_argument("--structures_map", help="Output structures map CSV")
    parser.add_argument("--aggregation_map", help="Output aggregation map CSV")
    parser.add_argument("--scores_csv", help="Output per-structure summary CSV")
    parser.add_argument("--aggregation_all_csv", help="Output merged per-residue CSV")
    args = parser.parse_args()

    if args.emit_worklist:
        if not args.worklist:
            parser.error("--emit-worklist requires --worklist")
        emit_worklist(args.structures, args.worklist)
        return

    required = ["work_root", "structures_dir", "images_dir", "aggregation_dir",
                "structures_map", "aggregation_map", "scores_csv", "aggregation_all_csv"]
    missing = [f"--{r}" for r in required if getattr(args, r) is None]
    if missing:
        parser.error("post-processing requires: " + ", ".join(missing))

    structures_ds = load_datastream(args.structures)
    struct_ids = structures_ds.ids_expanded
    if not struct_ids:
        print("Error: no structures found in input DataStream", file=sys.stderr)
        sys.exit(1)

    print(f"Post-processing Aggrescan3D output for {len(struct_ids)} structures")

    all_resi = []
    summary_rows = []
    structures_rows = []
    aggregation_rows = []
    failed = []

    for struct_id in struct_ids:
        try:
            resi_rows, summary_row = postprocess(
                struct_id, args.work_root,
                args.structures_dir, args.images_dir, args.aggregation_dir,
            )
        except Exception as e:
            print(f"WARNING: {struct_id} failed: {e}", file=sys.stderr)
            failed.append(struct_id)
            continue

        if resi_rows is None:
            failed.append(struct_id)
            continue

        all_resi.extend(resi_rows)
        summary_rows.append(summary_row)
        structures_rows.append({
            "id": struct_id,
            "file": os.path.join(args.structures_dir, f"{struct_id}.pdb"),
            "structures.id": struct_id,
        })
        aggregation_rows.append({
            "id": struct_id,
            "file": os.path.join(args.aggregation_dir, f"{struct_id}_A3D.csv"),
            "structures.id": struct_id,
        })

    # Always write maps/tables (possibly with whatever succeeded).
    pd.DataFrame(structures_rows, columns=["id", "file", "structures.id"]).to_csv(
        args.structures_map, index=False)
    pd.DataFrame(aggregation_rows, columns=["id", "file", "structures.id"]).to_csv(
        args.aggregation_map, index=False)
    pd.DataFrame(
        summary_rows,
        columns=["id", "structures.id", "avg_score", "min_score",
                 "max_score", "n_residues", "n_aggregation_prone"],
    ).to_csv(args.scores_csv, index=False)
    pd.DataFrame(all_resi, columns=["id", "chain", "resi", "score"]).to_csv(
        args.aggregation_all_csv, index=False)

    print(f"Wrote summary for {len(summary_rows)} structures, "
          f"{len(all_resi)} residue rows total")
    if failed:
        print(f"Failed {len(failed)}/{len(struct_ids)}: {failed}", file=sys.stderr)
    if not summary_rows:
        print("Error: no structures produced A3D scores", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
