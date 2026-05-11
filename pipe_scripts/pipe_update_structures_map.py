#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Write structures_map.csv at runtime from the actual output files.

Tools no longer write a per-design map_table at configuration time (only the
``map_table`` *path* is declared). After the tool runs, this script scans the
output folder, builds one row per produced file, and writes the map. A
constant provenance column (e.g. every design's parent PDB id) can be set via
``--set-provenance``; fan-out provenance (``design_i_1`` ← ``design_i``) is
recoverable from the id suffix and needs no column. If a config-time map
happens to exist (legacy tools), its non-id/file/value columns are carried
forward by exact id match.

Usage:
    python pipe_update_structures_map.py \\
        --structures-map /path/to/structures_map.csv \\
        --output-folder /path/to/output \\
        --extension pdb \\
        [--best-poses-dir /path/to/best_poses] \\
        [--set-provenance structures.id=4AKE]
"""

import argparse
import os
import glob
import pandas as pd


def update_map_from_folder(structures_map, output_folder, extension="pdb",
                           best_poses_dir=None, set_provenance=None):
    """Write structures_map.csv based on the actual files in output_folder.

    Args:
        structures_map: Path to the structures_map.csv to write.
        output_folder: Folder to scan for output files.
        extension: File extension to scan for (default: "pdb").
        best_poses_dir: If set, scan this directory instead (for Gnina best_poses).
        set_provenance: Optional dict {column_name: constant_value} of provenance
            columns to add to every row (e.g. {"structures.id": "4AKE"}).
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

    # Constant provenance columns (e.g. every design shares one parent PDB).
    if set_provenance:
        for col, value in set_provenance.items():
            new_df[col] = value

    # Carry forward any provenance columns from a pre-existing (legacy) map by
    # exact id match. No-op when the config-time map was not written.
    if os.path.exists(structures_map):
        try:
            old_df = pd.read_csv(structures_map)
            prov_cols = [c for c in old_df.columns
                         if c not in ("id", "file", "value") and c not in new_df.columns]
            if prov_cols:
                old_lookup = {}
                for _, row in old_df.iterrows():
                    old_lookup[str(row["id"])] = {c: row[c] for c in prov_cols}
                for col in prov_cols:
                    new_df[col] = new_df["id"].map(
                        lambda x: old_lookup.get(x, {}).get(col, "")
                    )
        except Exception as e:
            print(f"Warning: Could not carry provenance from existing map: {e}")

    new_df.to_csv(structures_map, index=False)
    print(f"Updated {structures_map} with {len(rows)} structures")


def _parse_set_provenance(spec):
    """Parse 'col=val,col2=val2' into {col: val}. Returns None if empty."""
    if not spec:
        return None
    out = {}
    for piece in spec.split(","):
        piece = piece.strip()
        if not piece:
            continue
        if "=" not in piece:
            raise ValueError(f"--set-provenance entry must be 'col=value', got: {piece!r}")
        col, val = piece.split("=", 1)
        out[col.strip()] = val.strip()
    return out or None


def main():
    parser = argparse.ArgumentParser(
        description="Write structures_map.csv from actual runtime output files"
    )
    parser.add_argument("--structures-map", required=True,
                        help="Path to structures_map.csv to write")
    parser.add_argument("--output-folder", required=True,
                        help="Folder containing output structure files")
    parser.add_argument("--extension", default="pdb",
                        help="File extension to scan for (default: pdb)")
    parser.add_argument("--best-poses-dir", default=None,
                        help="Scan this directory instead of output-folder (for Gnina)")
    parser.add_argument("--set-provenance", default=None,
                        help="Constant provenance columns, 'col=value[,col2=value2]' "
                             "(e.g. 'structures.id=4AKE')")

    args = parser.parse_args()

    update_map_from_folder(
        structures_map=args.structures_map,
        output_folder=args.output_folder,
        extension=args.extension,
        best_poses_dir=args.best_poses_dir,
        set_provenance=_parse_set_provenance(args.set_provenance),
    )


if __name__ == "__main__":
    main()
