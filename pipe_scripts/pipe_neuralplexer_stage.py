#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
NeuralPLexer input staging.

For each (protein, ligand) pair, materialise:
    {staging}/{pair_id}/receptor.pdb
    {staging}/{pair_id}/ligand.sdf

The pair id is composed by `biopipelines.combinatorics.predict_single_output_id`, and `--pairs-csv` records the protein and the ligand of each pair (`id | structures.id | compounds.id`) so post-processing reads them back instead of splitting the pair id.

When the ligand stream is value-based (SMILES only), embed to 3D with RDKit.
The wrapper bash then iterates the staging folder and invokes
neuralplexer-inference per pair.
"""

import argparse
import csv
import os
import shutil
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from biopipelines.biopipelines_io import (  # noqa: E402
    load_datastream,
    iterate_files,
    iterate_values,
)
from biopipelines.combinatorics import predict_single_output_id  # noqa: E402
from biopipelines.ligand_utils import ligand_chemistry  # noqa: E402


def materialise_ligand_sdf(lig_id: str, smiles: str, dst_sdf: str) -> None:
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"RDKit could not parse SMILES for {lig_id!r}: {smiles!r}")
    mol = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mol, randomSeed=42) != 0:
        raise ValueError(f"RDKit failed to embed 3D coordinates for {lig_id!r}")
    AllChem.MMFFOptimizeMolecule(mol)
    mol.SetProp("_Name", lig_id)
    writer = Chem.SDWriter(dst_sdf)
    writer.write(mol)
    writer.close()


def collect_ligands(ds_json: str):
    """Return [(lig_id, sdf_path_or_None, smiles_or_None), ...]."""
    ds = load_datastream(ds_json)
    out = []

    files = ds.files
    has_per_id_files = isinstance(files, list) and files and not (
        len(files) == 1 and "<id>" in files[0]
    )
    has_template = isinstance(files, list) and files and "<id>" in files[0]

    if has_per_id_files or has_template:
        for lid, fpath in iterate_files(ds):
            out.append((lid, os.path.abspath(fpath), None))
        return out

    for lid, values in iterate_values(ds, columns=["smiles"]):
        smiles = ligand_chemistry(values)["smiles"]
        if not smiles:
            raise ValueError(f"Compound {lid!r}: no file and no 'smiles' value")
        out.append((lid, None, smiles))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--structures-json", required=True)
    ap.add_argument("--compounds-json", required=True)
    ap.add_argument("--staging-folder", required=True)
    ap.add_argument("--pairs-csv", required=True, help="Where each pair's protein and ligand ids are recorded")
    args = ap.parse_args()

    proteins = list(iterate_files(load_datastream(args.structures_json)))
    ligands = collect_ligands(args.compounds_json)

    if not proteins or not ligands:
        print("ERROR: empty input streams", file=sys.stderr)
        sys.exit(1)

    # Wipe stale pair folders — re-running with different inputs must not
    # leave orphaned directories that the inference loop or postprocess would
    # treat as current pairs.
    if os.path.isdir(args.staging_folder):
        shutil.rmtree(args.staging_folder)
    os.makedirs(args.staging_folder, exist_ok=True)

    protein_ids = [pid for pid, _ in proteins]
    ligand_ids = [lid for lid, _, _ in ligands]
    pair_rows = []
    for i, (prot_id, prot_path) in enumerate(proteins):
        for j, (lig_id, lig_path, lig_smiles) in enumerate(ligands):
            pair_id = predict_single_output_id(
                structures=("each", protein_ids, i, [], False),
                compounds=("each", ligand_ids, j, [], False),
            )
            pair_dir = os.path.join(args.staging_folder, pair_id)
            os.makedirs(pair_dir, exist_ok=True)

            shutil.copyfile(prot_path, os.path.join(pair_dir, "receptor.pdb"))
            dst_sdf = os.path.join(pair_dir, "ligand.sdf")
            if lig_path is not None:
                shutil.copyfile(lig_path, dst_sdf)
            else:
                materialise_ligand_sdf(lig_id, lig_smiles, dst_sdf)
            pair_rows.append({"id": pair_id, "structures.id": prot_id, "compounds.id": lig_id})

    os.makedirs(os.path.dirname(args.pairs_csv), exist_ok=True)
    with open(args.pairs_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["id", "structures.id", "compounds.id"])
        w.writeheader()
        w.writerows(pair_rows)

    print(f"Staged {len(pair_rows)} (protein, ligand) pair(s) under {args.staging_folder}")


if __name__ == "__main__":
    main()
