#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""RDKit descriptor and conformer-strain calculator.

With a compounds DataStream: looks up SMILES via biopipelines_io, computes the requested descriptors (and optionally a Morgan fingerprint), writes one row per input compound to the descriptors CSV.

With a structures DataStream: rebuilds each posed ligand with bond orders from a SMILES template and scores its torsional strain against a torsion-restrained reference state, writing one row per pose to the strain CSV.
"""

import argparse
import json
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from biopipelines.biopipelines_io import (  # noqa: E402
    load_datastream, iterate_values, iterate_files, load_table, lookup_table_value,
)
from biopipelines.id_map_utils import get_mapped_ids  # noqa: E402
from biopipelines.ligand_utils import conformer_strain, posed_ligand_mol  # noqa: E402

from rdkit import Chem  # noqa: E402
from rdkit.Chem import AllChem, Descriptors  # noqa: E402
from rdkit.Chem.QED import qed as compute_qed  # noqa: E402

STRAIN_COLUMNS = ["id", "smiles", "e_pose", "e_relaxed", "strain", "ff_engine"]


def get_descriptor_fn(name: str):
    """Resolve a descriptor function. Special-cases `qed`; otherwise delegates
    to Descriptors.<name>. Raises AttributeError if the name is unknown."""
    if name == "qed":
        return compute_qed
    return getattr(Descriptors, name)


def compute_descriptors(cfg):
    ds = load_datastream(cfg["compounds_json"])
    output_csv = cfg["descriptors_csv"]
    descriptor_names = cfg["descriptors"]
    morgan_fp = bool(cfg.get("morgan_fp", False))

    # Resolve descriptor functions up-front to fail fast on typos.
    desc_fns = {}
    for n in descriptor_names:
        try:
            desc_fns[n] = get_descriptor_fn(n)
        except AttributeError:
            print(f"ERROR: unknown RDKit descriptor: {n}", file=sys.stderr)
            sys.exit(1)

    rows = []
    failed = []
    for cid, values in iterate_values(ds, columns=["smiles"]):
        smiles = values.get("smiles", "")
        if not smiles:
            print(f"WARNING: {cid} has no SMILES", file=sys.stderr)
            failed.append(cid)
            continue
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"WARNING: {cid} SMILES failed to parse: {smiles}", file=sys.stderr)
            failed.append(cid)
            continue

        row = {"id": cid, "smiles": smiles}
        for name, fn in desc_fns.items():
            try:
                row[name] = float(fn(mol))
            except Exception as e:
                print(f"WARNING: {cid} descriptor {name} failed: {e}", file=sys.stderr)
                row[name] = None
        if morgan_fp:
            bv = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
            for i in range(1024):
                row[f"fp_{i}"] = int(bv.GetBit(i))
        rows.append(row)

    if not rows:
        print("ERROR: no compounds produced descriptors", file=sys.stderr)
        return False

    os.makedirs(os.path.dirname(output_csv), exist_ok=True)
    columns = ["id", "smiles"] + descriptor_names
    if morgan_fp:
        columns += [f"fp_{i}" for i in range(1024)]
    pd.DataFrame(rows, columns=columns).to_csv(output_csv, index=False)
    print(f"Descriptors: {output_csv} ({len(rows)} rows, {len(columns)} columns)")

    if failed:
        print(f"Failed: {len(failed)}/{len(failed)+len(rows)}: {failed}", file=sys.stderr)
    return True


def _smiles_resolver(cfg, structures_ds):
    """Return a callable structure_id -> template SMILES.

    The template comes from an explicit `smiles` (a literal or a TABLE_REFERENCE)
    or, failing that, the `smiles` column of the compounds stream.

    A structure id is not a compound id: a docked/co-folded pose is `prot+lig`, or a
    fan-out child, while the compound it carries is `lig`. The link is the structures
    map's `compounds.id` provenance, resolved with get_mapped_ids exactly as Gnina
    pairs a complex with the one compound it holds.
    """
    spec = cfg.get("smiles")
    if spec and str(spec).startswith("TABLE_REFERENCE:"):
        table, column = load_table(spec)
        return lambda sid: lookup_table_value(table, sid, column)
    if spec:
        return lambda sid: spec

    compounds_json = cfg.get("compounds_json")
    if not compounds_json:
        print("ERROR: strain requires a bond-order template (smiles= or a compounds stream)",
              file=sys.stderr)
        sys.exit(1)

    cds = load_datastream(compounds_json)
    by_id = {cid: values.get("smiles", "") for cid, values in iterate_values(cds, columns=["smiles"])}

    structure_ids = list(structures_ds.ids_expanded)
    compound_ids = list(by_id.keys())
    if len(compound_ids) == 1:
        structure_to_compound = {sid: compound_ids[0] for sid in structure_ids}
    else:
        structure_to_compound = get_mapped_ids(
            source_ids=structure_ids,
            target_ids=compound_ids,
            map_table_paths=[structures_ds.map_table] if structures_ds.map_table else None,
            unique=True,
        )

    def resolve(sid):
        cid = structure_to_compound.get(sid)
        if cid is None:
            raise KeyError(
                f"no compound paired with pose {sid}; the structures map carries no "
                f"compounds.id provenance and the compounds stream has {len(compound_ids)} "
                f"compounds. Pass smiles= explicitly to supply the bond-order template")
        return by_id[cid]

    return resolve


def compute_strain(cfg):
    ds = load_datastream(cfg["structures_json"])
    output_csv = cfg["strain_csv"]
    restrain_bonds = cfg.get("restrain_bonds")
    restrain_bonds = [tuple(b) for b in restrain_bonds] if restrain_bonds else None
    ff = cfg.get("ff", "auto")

    resolve_smiles = _smiles_resolver(cfg, ds)

    rows = []
    failed = []
    for sid, path in iterate_files(ds):
        try:
            smiles = resolve_smiles(sid)
            if not smiles:
                raise ValueError("no bond-order template SMILES for this id")
            mol = posed_ligand_mol(path, smiles)
            e_pose, e_relaxed, strain, engine = conformer_strain(
                mol, restrain_bonds=restrain_bonds, ff=ff)
            rows.append({"id": sid, "smiles": smiles, "e_pose": e_pose,
                         "e_relaxed": e_relaxed, "strain": strain, "ff_engine": engine})
        except Exception as e:
            print(f"WARNING: {sid} strain failed: {e}", file=sys.stderr)
            failed.append(sid)

    if not rows:
        print("ERROR: no poses produced strain values", file=sys.stderr)
        return False

    os.makedirs(os.path.dirname(output_csv), exist_ok=True)
    pd.DataFrame(rows, columns=STRAIN_COLUMNS).to_csv(output_csv, index=False)
    print(f"Strain: {output_csv} ({len(rows)} rows)")

    if failed:
        print(f"Failed: {len(failed)}/{len(failed)+len(rows)}: {failed}", file=sys.stderr)
    return True


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--config", required=True)
    args = p.parse_args()

    with open(args.config, encoding="utf-8") as f:
        cfg = json.load(f)

    ok = True
    if cfg.get("compounds_json") and cfg.get("descriptors_csv"):
        ok = compute_descriptors(cfg) and ok
    if cfg.get("structures_json"):
        ok = compute_strain(cfg) and ok

    if not ok:
        sys.exit(1)


if __name__ == "__main__":
    main()
