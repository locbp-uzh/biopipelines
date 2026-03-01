#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Runtime helper script for CABS-Flex post-processing.

Merges per-structure RMSF CSVs into a single rmsf_all.csv and rebuilds
the structures map CSV with actual output model files.
"""

import os
import sys
import argparse

import pandas as pd

# Import unified I/O utilities
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from biopipelines.biopipelines_io import load_datastream


def main():
    parser = argparse.ArgumentParser(
        description="CABS-Flex post-processing: merge RMSF tables and build structures map"
    )
    parser.add_argument("--structures", required=True, help="JSON file containing input DataStream")
    parser.add_argument("--output_dir", required=True, help="Output directory with CABS-Flex results")
    parser.add_argument("--rmsf_all_csv", required=True, help="Output merged RMSF CSV path")
    parser.add_argument("--structures_map", required=True, help="Output structures map CSV path")
    parser.add_argument("--num_models", type=int, default=10, help="Number of models per structure")

    args = parser.parse_args()

    structures_ds = load_datastream(args.structures)

    if not structures_ds.ids:
        print("Error: No structures found in input DataStream")
        sys.exit(1)

    print(f"Processing {len(structures_ds.ids)} input structures")

    # --- Merge per-ID RMSF CSVs ---
    rmsf_frames = []
    for struct_id in structures_ds.ids:
        rmsf_path = os.path.join(args.output_dir, f"{struct_id}_RMSF.csv")
        if os.path.exists(rmsf_path):
            df = pd.read_csv(rmsf_path)
            rmsf_frames.append(df)
            print(f"  Loaded RMSF: {rmsf_path} ({len(df)} residues)")
        else:
            print(f"  Warning: RMSF not found: {rmsf_path}")

    if rmsf_frames:
        rmsf_all = pd.concat(rmsf_frames, ignore_index=True)
        rmsf_all.to_csv(args.rmsf_all_csv, index=False)
        print(f"Merged RMSF written: {args.rmsf_all_csv} ({len(rmsf_all)} rows)")
    else:
        pd.DataFrame(columns=["id", "chain", "resi", "rmsf"]).to_csv(
            args.rmsf_all_csv, index=False
        )
        print("Warning: No RMSF data found, wrote empty CSV")

    # --- Build structures map from actual output files ---
    rows = []
    for input_id in structures_ds.ids:
        for model_idx in range(1, args.num_models + 1):
            output_id = f"{input_id}_{model_idx}"
            pdb_path = os.path.join(args.output_dir, f"{output_id}.pdb")
            if os.path.exists(pdb_path):
                rows.append({
                    "id": output_id,
                    "file": pdb_path,
                    "structures.id": input_id
                })
            else:
                print(f"  Warning: Expected model not found: {pdb_path}")

    df = pd.DataFrame(rows)
    df.to_csv(args.structures_map, index=False)
    print(f"Structures map written: {args.structures_map} ({len(rows)} models)")


if __name__ == "__main__":
    main()
