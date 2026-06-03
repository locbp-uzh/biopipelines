#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
DynamicBind input staging.

DynamicBind's run_single_protein_inference.py takes one protein and one
ligand CSV per invocation. We write a CSV per protein listing every ligand
in the compounds stream — the wrapper bash then iterates proteins.

Ligand CSV format (per upstream README): a `ligand` column with SMILES. We
add a `name` column so post-processing can map rows back to BioPipelines
compound IDs.
"""

import argparse
import csv
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from biopipelines.biopipelines_io import (  # noqa: E402
    load_datastream,
    iterate_files,
    iterate_values,
)


def collect_ligands(ds_json: str):
    """Return [(lig_id, smiles), ...]. DynamicBind needs SMILES; if the
    compounds stream is file-based (SDF), we read the first MOL block.
    """
    ds = load_datastream(ds_json)

    files = ds.files
    has_per_id_files = isinstance(files, list) and files and not (
        len(files) == 1 and "<id>" in files[0]
    )
    has_template = isinstance(files, list) and files and "<id>" in files[0]

    if has_per_id_files or has_template:
        from rdkit import Chem
        out = []
        for lid, fpath in iterate_files(ds):
            mol = next(iter(Chem.SDMolSupplier(fpath, removeHs=False)), None)
            if mol is None:
                raise ValueError(f"RDKit could not read SDF for {lid!r}: {fpath}")
            out.append((lid, Chem.MolToSmiles(mol)))
        return out

    out = []
    for lid, values in iterate_values(ds, columns=["smiles"]):
        smiles = values.get("smiles", "")
        if not smiles:
            raise ValueError(
                f"Compound {lid!r}: no file and no 'smiles' value in map_table"
            )
        out.append((lid, smiles))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--structures-json", required=True)
    ap.add_argument("--compounds-json", required=True)
    ap.add_argument("--output-dir", required=True, help="One <prot>.csv per protein lands here")
    args = ap.parse_args()

    protein_ids = list(load_datastream(args.structures_json).ids_expanded)
    ligands = collect_ligands(args.compounds_json)

    if not protein_ids:
        print("ERROR: no proteins", file=sys.stderr)
        sys.exit(1)
    if not ligands:
        print("ERROR: no ligands", file=sys.stderr)
        sys.exit(1)

    os.makedirs(args.output_dir, exist_ok=True)
    for prot_id in protein_ids:
        out_csv = os.path.join(args.output_dir, f"{prot_id}.csv")
        with open(out_csv, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["ligand", "name"])
            writer.writeheader()
            for lig_id, smiles in ligands:
                writer.writerow({"ligand": smiles, "name": lig_id})

    print(
        f"Wrote {len(protein_ids)} per-protein ligand CSV(s) "
        f"({len(ligands)} ligands each) under {args.output_dir}"
    )


if __name__ == "__main__":
    main()
