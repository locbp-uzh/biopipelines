#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Runtime helper script for CABS-Flex.

Runs CABSflex for each input structure, collects output models, parses RMSF,
copies SVG plots, and builds the final structures map and merged RMSF table.
"""

import os
import re
import sys
import shutil
import argparse
import subprocess

import pandas as pd

# Import unified I/O utilities
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from biopipelines.biopipelines_io import load_datastream, iterate_files


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


def run_cabsflex(struct_id, struct_file, output_dir, cabsflex_flags):
    """
    Run CABSflex for a single structure.

    Returns True on success, False on failure.
    """
    work_dir = os.path.join(output_dir, struct_id)
    os.makedirs(work_dir, exist_ok=True)

    cmd = ["CABSflex", "-i", struct_file, "-w", work_dir] + cabsflex_flags
    print(f"=== Processing {struct_id} ===")
    print(f"Command: {' '.join(cmd)}")

    result = subprocess.run(cmd)
    if result.returncode != 0:
        print(f"Error: CABS-Flex failed for {struct_id}")
        return False

    # Copy model PDBs to output folder with proper naming
    pdbs_dir = os.path.join(work_dir, "output_pdbs")
    model_idx = 0
    if os.path.isdir(pdbs_dir):
        for model_pdb in sorted(
            f for f in os.listdir(pdbs_dir) if f.startswith("model_") and f.endswith(".pdb")
        ):
            model_idx += 1
            src = os.path.join(pdbs_dir, model_pdb)
            dst = os.path.join(output_dir, f"{struct_id}_{model_idx}.pdb")
            shutil.copy2(src, dst)
    print(f"Copied {model_idx} models for {struct_id}")

    # Copy SVG plots to output folder
    plots_dir = os.path.join(work_dir, "plots")
    if os.path.isdir(plots_dir):
        for svg_file in os.listdir(plots_dir):
            if svg_file.endswith(".svg"):
                src = os.path.join(plots_dir, svg_file)
                dst = os.path.join(output_dir, f"{struct_id}_{svg_file}")
                shutil.copy2(src, dst)

    # Parse and save per-structure RMSF
    rmsf_path = os.path.join(plots_dir, "RMSF.csv") if os.path.isdir(plots_dir) else ""
    if rmsf_path and os.path.exists(rmsf_path):
        rmsf_rows = parse_rmsf(rmsf_path, struct_id)
        rmsf_df = pd.DataFrame(rmsf_rows)
        rmsf_out = os.path.join(output_dir, f"{struct_id}_RMSF.csv")
        rmsf_df.to_csv(rmsf_out, index=False)
        print(f"RMSF saved to {rmsf_out} ({len(rmsf_rows)} residues)")
    else:
        print(f"Warning: No RMSF.csv found for {struct_id}")

    return True


def build_outputs(structures_ds, output_dir, num_models, rmsf_all_csv, structures_map):
    """
    Build merged RMSF table and structures map from all results.
    """
    # Merge per-ID RMSF CSVs
    rmsf_frames = []
    for struct_id in structures_ds.ids:
        rmsf_path = os.path.join(output_dir, f"{struct_id}_RMSF.csv")
        if os.path.exists(rmsf_path):
            df = pd.read_csv(rmsf_path)
            rmsf_frames.append(df)

    if rmsf_frames:
        rmsf_all = pd.concat(rmsf_frames, ignore_index=True)
        rmsf_all.to_csv(rmsf_all_csv, index=False)
        print(f"Merged RMSF written: {rmsf_all_csv} ({len(rmsf_all)} rows)")
    else:
        pd.DataFrame(columns=["id", "chain", "resi", "rmsf"]).to_csv(
            rmsf_all_csv, index=False
        )
        print("Warning: No RMSF data found, wrote empty CSV")

    # Build structures map from actual output files
    rows = []
    for input_id in structures_ds.ids:
        for model_idx in range(1, num_models + 1):
            output_id = f"{input_id}_{model_idx}"
            pdb_path = os.path.join(output_dir, f"{output_id}.pdb")
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


def main():
    parser = argparse.ArgumentParser(description="CABS-Flex runner and post-processor")
    parser.add_argument("--structures", required=True, help="JSON file containing input DataStream")
    parser.add_argument("--output_dir", required=True, help="Output directory")
    parser.add_argument("--rmsf_all_csv", required=True, help="Output merged RMSF CSV path")
    parser.add_argument("--structures_map", required=True, help="Output structures map CSV path")
    parser.add_argument("--num_models", type=int, default=10, help="Number of models per structure")
    parser.add_argument("--cabsflex_flags", nargs="*", default=[], help="Flags to pass to CABSflex")

    args = parser.parse_args()

    structures_ds = load_datastream(args.structures)

    if not structures_ds.ids:
        print("Error: No structures found in input DataStream")
        sys.exit(1)

    print(f"Running CABS-Flex on {len(structures_ds.ids)} structures")

    # Run CABSflex for each structure
    for struct_id, struct_file in iterate_files(structures_ds):
        if not os.path.exists(struct_file):
            print(f"Warning: Structure file not found: {struct_file}")
            continue
        run_cabsflex(struct_id, struct_file, args.output_dir, args.cabsflex_flags)

    # Build merged outputs
    print("\n=== Post-processing ===")
    build_outputs(structures_ds, args.output_dir, args.num_models,
                  args.rmsf_all_csv, args.structures_map)


if __name__ == "__main__":
    main()
