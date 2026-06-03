#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
DiffDock input-CSV builder.

DiffDock's inference.py reads a single CSV with columns
    complex_name, protein_path, ligand_description, protein_sequence
and runs one prediction per row. This script materialises the
cartesian product of the upstream structures and compounds streams into
that CSV.

The `complex_name` we emit matches the `{protein_id}+{ligand_id}` pair-ID
the BioPipelines wrapper predicts, so post-processing can map per-complex
folders back to BioPipelines IDs without any extra bookkeeping.

`ligand_description` is the SMILES when the compound stream is value-based
(map_table has a `smiles` column) and the absolute SDF path otherwise.
"""

import argparse
import csv
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from biopipelines.biopipelines_io import load_datastream, iterate_files, iterate_values  # noqa: E402


def collect_proteins(ds_json: str):
    """Return [(protein_id, abs_pdb_path), ...]."""
    ds = load_datastream(ds_json)
    return [(pid, os.path.abspath(path)) for pid, path in iterate_files(ds)]


def collect_ligands(ds_json: str):
    """Return [(ligand_id, ligand_description), ...].

    ligand_description is the SMILES string if available in the map_table,
    falling back to the SDF/MOL file path. DiffDock accepts either form.
    """
    ds = load_datastream(ds_json)
    out = []

    files = ds.files
    has_per_id_files = isinstance(files, list) and files and not (
        len(files) == 1 and "<id>" in files[0]
    )

    if has_per_id_files or (isinstance(files, list) and files and "<id>" in files[0]):
        for lid, fpath in iterate_files(ds):
            out.append((lid, os.path.abspath(fpath)))
        return out

    # Value-based stream: read SMILES from the map_table.
    for lid, values in iterate_values(ds, columns=["smiles"]):
        smiles = values.get("smiles", "")
        if not smiles:
            raise ValueError(
                f"Compound {lid!r} has neither a file nor a 'smiles' value in the map_table; "
                "DiffDock needs at least one."
            )
        out.append((lid, smiles))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--structures-json", required=True)
    ap.add_argument("--compounds-json", required=True)
    ap.add_argument("--output-csv", required=True)
    args = ap.parse_args()

    proteins = collect_proteins(args.structures_json)
    ligands = collect_ligands(args.compounds_json)

    if not proteins:
        print("ERROR: no proteins to dock", file=sys.stderr)
        sys.exit(1)
    if not ligands:
        print("ERROR: no ligands to dock", file=sys.stderr)
        sys.exit(1)

    os.makedirs(os.path.dirname(args.output_csv), exist_ok=True)

    # At DiffDock commit a6c5275 (the version we pin in _install_script),
    # inference.py reads only `protein_path` and `ligand`. The newer
    # `ligand_description`/`complex_name`/`protein_sequence` columns belong
    # to later upstream revisions and are silently ignored here, but the
    # legacy `ligand` column is required.
    rows = []
    for prot_id, prot_path in proteins:
        for lig_id, lig_desc in ligands:
            rows.append(
                {
                    "complex_name": f"{prot_id}+{lig_id}",
                    "protein_path": prot_path,
                    "ligand": lig_desc,
                }
            )

    with open(args.output_csv, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["complex_name", "protein_path", "ligand"],
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} (protein, ligand) pair(s) to {args.output_csv}")


if __name__ == "__main__":
    main()
