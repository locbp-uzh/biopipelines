#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Update structures_map.csv at runtime with actual output files.

When tools use create_map_table() at config time with lazy/bracket IDs,
the map only contains prefix-expanded IDs. At runtime, after outputs are
generated, this script rewrites the map with the actual files found on disk.

Preserves provenance columns from the existing config-time map.

Usage:
    python pipe_update_structures_map.py \\
        --structures-map /path/to/structures_map.csv \\
        --output-folder /path/to/output \\
        --extension pdb \\
        [--best-poses-dir /path/to/best_poses] \\
        [--provenance structures:prot_id,compounds:lig_id]
"""

import argparse
import os
import glob
import pandas as pd


def update_map_from_folder(structures_map, output_folder, extension="pdb",
                           best_poses_dir=None, provenance_spec=None):
    """Rewrite structures_map.csv based on actual files in output_folder.

    Args:
        structures_map: Path to structures_map.csv to rewrite.
        output_folder: Folder to scan for output files.
        extension: File extension to scan for (default: "pdb").
        best_poses_dir: If set, scan this directory instead (for Gnina best_poses).
        provenance_spec: Comma-separated "stream:col" pairs to preserve from
            the existing map (e.g. "structures:structures.id,compounds:compounds.id").
    """
    scan_dir = best_poses_dir or output_folder
    pattern = os.path.join(scan_dir, f"*.{extension}")
    found_files = sorted(glob.glob(pattern))

    if not found_files:
        print(f"Warning: No .{extension} files found in {scan_dir}, keeping existing map")
        return

    # Build rows from discovered files
    rows = []
    for fpath in found_files:
        basename = os.path.basename(fpath)
        file_id = os.path.splitext(basename)[0]
        # For Gnina best_poses: strip _best suffix
        if best_poses_dir and file_id.endswith("_best"):
            file_id = file_id[:-5]
        rows.append({"id": file_id, "file": fpath})

    new_df = pd.DataFrame(rows)

    # Try to preserve provenance columns from existing map
    if os.path.exists(structures_map):
        try:
            old_df = pd.read_csv(structures_map)
            prov_cols = [c for c in old_df.columns if c not in ("id", "file", "value")]
            if prov_cols:
                # Build lookup from old map
                old_lookup = {}
                for _, row in old_df.iterrows():
                    old_lookup[str(row["id"])] = {c: row[c] for c in prov_cols}

                # Try to match new IDs to old provenance
                for col in prov_cols:
                    new_df[col] = new_df["id"].map(
                        lambda x: old_lookup.get(x, {}).get(col, "")
                    )
        except Exception as e:
            print(f"Warning: Could not preserve provenance from existing map: {e}")

    new_df.to_csv(structures_map, index=False)
    print(f"Updated {structures_map} with {len(rows)} structures")


def main():
    parser = argparse.ArgumentParser(
        description="Update structures_map.csv with actual runtime output files"
    )
    parser.add_argument("--structures-map", required=True,
                        help="Path to structures_map.csv to rewrite")
    parser.add_argument("--output-folder", required=True,
                        help="Folder containing output structure files")
    parser.add_argument("--extension", default="pdb",
                        help="File extension to scan for (default: pdb)")
    parser.add_argument("--best-poses-dir", default=None,
                        help="Scan this directory instead of output-folder (for Gnina)")

    args = parser.parse_args()

    update_map_from_folder(
        structures_map=args.structures_map,
        output_folder=args.output_folder,
        extension=args.extension,
        best_poses_dir=args.best_poses_dir,
    )


if __name__ == "__main__":
    main()
