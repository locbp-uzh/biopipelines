#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Post-processing helper for CABS-Flex.

Collects output models, parses RMSF, copies SVG plots, and builds the
final structures map and merged RMSF table. CABSflex itself runs from
the bash script under Python 2.7; this script runs under biopipelines
(Python 3) for post-processing only.
"""

import os
import re
import sys
import shutil
import argparse

import pandas as pd

# Import unified I/O utilities
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from biopipelines.biopipelines_io import load_datastream, iterate_files


def emit_worklist(structures_json, worklist_path):
    """Resolve the input stream to id<TAB>file lines for the bash loop.

    Runs under biopipelines (Python 3) before the CABSflex Python-2.7 segment,
    so lazy IDs expand against the map_table at runtime (never at config time).
    """
    ds = load_datastream(structures_json)
    rows = list(iterate_files(ds))
    if not rows:
        print("Error: no structures resolved from input DataStream", file=sys.stderr)
        sys.exit(1)
    with open(worklist_path, "w") as f:
        for struct_id, struct_file in rows:
            f.write("{0}\t{1}\n".format(struct_id, struct_file))
    print("Wrote worklist: {0} ({1} structures)".format(worklist_path, len(rows)))


def parse_rmsf(rmsf_path, struct_id):
    """
    Parse CABS-Flex RMSF.csv (tab-separated, no header, e.g. 'A2\\t3.516').

    Returns list of dicts with keys: id, chain, resi, rmsf.
    """
    rows = []
    with open(rmsf_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) != 2:
                continue
            resid, rmsf_val = parts
            # Parse chain letter and residue number from e.g. "A2", "B123"
            match = re.match(r"^([A-Za-z])(\d+)$", resid)
            if match:
                chain = match.group(1)
                resi = match.group(2)
            else:
                chain = ""
                resi = resid
            rows.append({
                "id": struct_id,
                "chain": chain,
                "resi": resi,
                "rmsf": float(rmsf_val)
            })
    return rows


def postprocess(struct_id, work_root, structures_dir, images_dir, rmsf_dir, num_models):
    """Post-process CABSflex output for a single structure."""
    work_dir = os.path.join(work_root, struct_id)

    # Copy model PDBs
    pdbs_dir = os.path.join(work_dir, "output_pdbs")
    model_idx = 0
    if os.path.isdir(pdbs_dir):
        for model_pdb in sorted(
            f for f in os.listdir(pdbs_dir) if f.startswith("model_") and f.endswith(".pdb")
        ):
            model_idx += 1
            src = os.path.join(pdbs_dir, model_pdb)
            dst = os.path.join(structures_dir, f"{struct_id}_{model_idx}.pdb")
            shutil.copy2(src, dst)
    print(f"Copied {model_idx} models for {struct_id}")

    # Copy SVG plots
    plots_dir = os.path.join(work_dir, "plots")
    if os.path.isdir(plots_dir):
        for svg_file in os.listdir(plots_dir):
            if svg_file.endswith(".svg"):
                src = os.path.join(plots_dir, svg_file)
                dst = os.path.join(images_dir, f"{struct_id}_{svg_file}")
                shutil.copy2(src, dst)

    # Parse per-structure RMSF
    rmsf_path = os.path.join(plots_dir, "RMSF.csv") if os.path.isdir(plots_dir) else ""
    if rmsf_path and os.path.exists(rmsf_path):
        rmsf_rows = parse_rmsf(rmsf_path, struct_id)
        rmsf_df = pd.DataFrame(rmsf_rows)
        rmsf_out = os.path.join(rmsf_dir, f"{struct_id}_RMSF.csv")
        rmsf_df.to_csv(rmsf_out, index=False)
        print(f"RMSF saved to {rmsf_out} ({len(rmsf_rows)} residues)")
    else:
        print(f"Warning: No RMSF.csv found for {struct_id}")


def build_outputs(structures_ds, structures_dir, rmsf_dir, num_models, rmsf_all_csv, structures_map, rmsf_map):
    """Build merged RMSF table, structures map, and rmsf map from per-ID artefacts."""
    # Merge per-ID RMSF CSVs; record the rmsf-stream map (one CSV per input id)
    rmsf_frames = []
    rmsf_rows = []
    for struct_id in structures_ds.ids_expanded:
        rmsf_path = os.path.join(rmsf_dir, f"{struct_id}_RMSF.csv")
        if os.path.exists(rmsf_path):
            df = pd.read_csv(rmsf_path)
            rmsf_frames.append(df)
            rmsf_rows.append({"id": struct_id, "file": rmsf_path, "structures.id": struct_id})
        else:
            print(f"  Warning: Expected RMSF CSV not found: {rmsf_path}")

    if rmsf_frames:
        rmsf_all = pd.concat(rmsf_frames, ignore_index=True)
        rmsf_all.to_csv(rmsf_all_csv, index=False)
        print(f"Merged RMSF written: {rmsf_all_csv} ({len(rmsf_all)} rows)")
    else:
        pd.DataFrame(columns=["id", "chain", "resi", "rmsf"]).to_csv(
            rmsf_all_csv, index=False
        )
        print("Warning: No RMSF data found, wrote empty CSV")

    rows = []
    for input_id in structures_ds.ids_expanded:
        for model_idx in range(1, num_models + 1):
            output_id = f"{input_id}_{model_idx}"
            pdb_path = os.path.join(structures_dir, f"{output_id}.pdb")
            if os.path.exists(pdb_path):
                rows.append({
                    "id": output_id,
                    "file": pdb_path,
                    "structures.id": input_id
                })
            else:
                print(f"  Warning: Expected model not found: {pdb_path}")

    df = pd.DataFrame(rows)
    df.to_csv(structures_map, index=False)
    print(f"Structures map written: {structures_map} ({len(rows)} models)")

    pd.DataFrame(rmsf_rows, columns=["id", "file", "structures.id"]).to_csv(rmsf_map, index=False)
    print(f"RMSF map written: {rmsf_map} ({len(rmsf_rows)} structures)")


def main():
    parser = argparse.ArgumentParser(description="CABS-Flex post-processor")
    parser.add_argument("--structures", required=True, help="JSON file containing input DataStream")
    parser.add_argument("--emit-worklist", action="store_true",
                        help="Resolve the stream to id<TAB>file lines and exit (Python 3 segment)")
    parser.add_argument("--worklist", help="Output worklist TSV path (with --emit-worklist)")
    parser.add_argument("--output_dir", help="Tool output folder (also used as work_root fallback)")
    parser.add_argument("--rmsf_all_csv", help="Output merged RMSF CSV path")
    parser.add_argument("--structures_map", help="Output structures map CSV path")
    parser.add_argument("--rmsf_map", help="Output rmsf map CSV path")
    parser.add_argument("--num_models", type=int, default=10, help="Number of models per structure")
    parser.add_argument("--work_root", default=None,
                        help="Folder containing per-structure CABSflex work dirs (default: output_dir)")
    parser.add_argument("--structures_dir", default=None,
                        help="Destination for model PDBs (default: output_dir)")
    parser.add_argument("--images_dir", default=None,
                        help="Destination for SVG plots (default: output_dir)")
    parser.add_argument("--rmsf_dir", default=None,
                        help="Destination for per-ID RMSF CSVs (default: output_dir)")

    args = parser.parse_args()

    if args.emit_worklist:
        if not args.worklist:
            parser.error("--emit-worklist requires --worklist")
        emit_worklist(args.structures, args.worklist)
        return

    required = ["output_dir", "rmsf_all_csv", "structures_map", "rmsf_map"]
    missing = ["--{0}".format(r) for r in required if getattr(args, r) is None]
    if missing:
        parser.error("missing required arguments for post-processing: {0}".format(", ".join(missing)))

    structures_ds = load_datastream(args.structures)

    if not structures_ds.ids_expanded:
        print("Error: No structures found in input DataStream")
        sys.exit(1)

    work_root = args.work_root or args.output_dir
    structures_dir = args.structures_dir or args.output_dir
    images_dir = args.images_dir or args.output_dir
    rmsf_dir = args.rmsf_dir or args.output_dir

    print(f"Post-processing CABS-Flex output for {len(structures_ds.ids_expanded)} structures")

    for struct_id in structures_ds.ids_expanded:
        postprocess(struct_id, work_root, structures_dir, images_dir, rmsf_dir, args.num_models)

    print("\n=== Building final outputs ===")
    build_outputs(structures_ds, structures_dir, rmsf_dir, args.num_models,
                  args.rmsf_all_csv, args.structures_map, args.rmsf_map)


if __name__ == "__main__":
    main()
