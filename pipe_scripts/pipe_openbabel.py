#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""OpenBabel converter. Reads a DataStream JSON (compounds or structures) and:

  - convert_3d: writes one coordinate file per id (<id>.<fmt>) into the
    structures stream folder and a structures map CSV (id, file).
  - convert_1d: computes line notations (smi, inchi, ...) once per id and merges
    them as columns into a copy of the compounds map_table.

A compounds input is read from SMILES (value-based); a structures input is read
from its per-id coordinate files.
"""

import argparse
import os
import sys

import pandas as pd
from openbabel import pybel  # type: ignore

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from biopipelines.biopipelines_io import (  # noqa: E402
    load_datastream, iterate_files, iterate_values,
)


def detect_input_format(ds) -> str:
    fmt = (ds.format or "").lower()
    if fmt in ("smiles", "smi"):
        return "smi"
    return fmt or "sdf"


def _build_mol(payload, in_fmt, from_smiles):
    if from_smiles:
        if not payload:
            raise ValueError("empty SMILES")
        return pybel.readstring("smi", payload)
    return next(pybel.readfile(in_fmt, payload))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input-json", required=True)
    p.add_argument("--input-kind", required=True, choices=["compounds", "structures"])
    p.add_argument("--convert-3d", default=None, help="single 3-D coordinate format")
    p.add_argument("--structures-dir", default=None)
    p.add_argument("--structures-map", default=None)
    p.add_argument("--convert-1d", nargs="*", default=[], help="line notations (smi, inchi, ...)")
    p.add_argument("--compounds-map", default=None)
    p.add_argument("--add-hydrogens", action="store_true")
    p.add_argument("--ph", type=float, default=None)
    p.add_argument("--gen3d", action="store_true")
    p.add_argument("--gen3d-quality", default="medium",
                   choices=["fastest", "fast", "medium", "better", "best"])
    p.add_argument("--minimize", action="store_true")
    p.add_argument("--ff", default="MMFF94")
    p.add_argument("--minimize-steps", type=int, default=500)
    p.add_argument("--coords-json", default=None,
                   help="optional structures stream (id->bound coordinate file) accompanying a "
                        "compounds input: write the SDF from THESE coords + SMILES template, no make3D")
    args = p.parse_args()

    # OpenBabel/pybel's make3D() takes a (forcefield, steps) pair rather than a
    # named quality; map the obabel --gen3d speed levels onto embed-step counts.
    GEN3D_STEPS = {"fastest": 50, "fast": 100, "medium": 250, "better": 500, "best": 1000}

    ds = load_datastream(args.input_json)
    in_fmt = detect_input_format(ds)

    # A compounds stream is value-based (SMILES); a structures stream has files.
    from_smiles = args.input_kind == "compounds"
    if from_smiles:
        # (id, smiles, code) so the 3-D writer can set the HETATM residue name.
        items = [
            (cid, vals.get("smiles", ""), vals.get("code", ""))
            for cid, vals in iterate_values(ds, columns=["smiles", "code"])
        ]
        in_fmt = "smi"
    else:
        items = [(cid, path, "") for cid, path in iterate_files(ds)]

    # Posed-ligand path: bound coordinate files accompanying a compounds input.
    # When present for an id and the target is SDF, write the SDF from those
    # coordinates + the SMILES bond-order template (shared util) instead of a
    # fresh make3D embedding — preserving the bound pose.
    coords_by_id = {}
    if args.coords_json:
        coords_by_id = {cid: path for cid, path in iterate_files(load_datastream(args.coords_json))}
    if coords_by_id:
        from biopipelines.ligand_utils import write_ligand_sdf  # noqa: E402

    ff_lower = args.ff.lower()  # pybel forcefield names are lowercase

    def prep(mol):
        if args.ph is not None:
            mol.OBMol.AddHydrogens(False, True, args.ph)
        elif args.add_hydrogens:
            mol.addh()
        if args.gen3d or args.convert_3d:
            # convert_3d implies coordinates; gen3d forces a fresh embed.
            if args.gen3d or from_smiles:
                mol.make3D(forcefield=ff_lower, steps=GEN3D_STEPS[args.gen3d_quality])
            if args.minimize:
                mol.localopt(forcefield=ff_lower, steps=args.minimize_steps)

    struct_rows = []
    notation_rows = {}  # id -> {fmt: notation}
    failed = []

    # Create the structures output dir before the conversion loop writes into it.
    if args.convert_3d:
        os.makedirs(args.structures_dir, exist_ok=True)

    for cid, payload, code in items:
        try:
            if args.convert_3d:
                out_path = os.path.join(args.structures_dir, f"{cid}.{args.convert_3d}")
                coord_file = coords_by_id.get(cid)
                if coord_file and args.convert_3d in ("sdf", "mol"):
                    # Bound pose available: write SDF from those coords, using the
                    # SMILES (payload, when input is compounds) as a bond-order
                    # template. Preserves the pose; no make3D.
                    write_ligand_sdf(coord_file, out_path, payload if from_smiles else None)
                    struct_rows.append({"id": cid, "file": out_path})
                    print(f"  {cid}: posed {coord_file} -> {out_path}")
                else:
                    mol = _build_mol(payload, in_fmt, from_smiles)
                    prep(mol)
                    if code:
                        for res in mol.OBMol.GetResidues() if hasattr(mol.OBMol, "GetResidues") else []:
                            res.SetName(code)
                    mol.write(args.convert_3d, out_path, overwrite=True)
                    struct_rows.append({"id": cid, "file": out_path})
                    print(f"  {cid}: -> {out_path}")

            if args.convert_1d:
                mol = _build_mol(payload, in_fmt, from_smiles)
                if args.ph is not None:
                    mol.OBMol.AddHydrogens(False, True, args.ph)
                elif args.add_hydrogens:
                    mol.addh()
                notation_rows[cid] = {
                    fmt: mol.write(fmt).strip() for fmt in args.convert_1d
                }
                print(f"  {cid}: 1D {list(args.convert_1d)}")
        except Exception as e:
            print(f"WARNING: {cid} conversion failed: {e}", file=sys.stderr)
            failed.append(cid)

    # Write the structures map (id, file) for the coordinate files.
    if args.convert_3d:
        os.makedirs(os.path.dirname(args.structures_map), exist_ok=True)
        pd.DataFrame(struct_rows, columns=["id", "file"]).to_csv(args.structures_map, index=False)
        print(f"Structures map: {args.structures_map} ({len(struct_rows)} rows)")

    # Merge notations as columns into a copy of the source compounds map_table.
    if args.convert_1d:
        src_map = ds.map_table
        if src_map and os.path.exists(src_map):
            df = pd.read_csv(src_map)
        else:
            df = pd.DataFrame({"id": [cid for cid, _, _ in items]})
        for fmt in args.convert_1d:
            df[fmt] = df["id"].map(lambda i: notation_rows.get(i, {}).get(fmt, ""))
        os.makedirs(os.path.dirname(args.compounds_map), exist_ok=True)
        df.to_csv(args.compounds_map, index=False)
        print(f"Compounds map (+{list(args.convert_1d)}): {args.compounds_map} ({len(df)} rows)")

    if failed:
        print(f"Failed: {len(failed)}/{len(items)}: {failed}", file=sys.stderr)

    # Error only if everything that should have produced output produced nothing.
    produced = bool(struct_rows) or bool(notation_rows)
    if items and not produced:
        sys.exit(1)


if __name__ == "__main__":
    main()
