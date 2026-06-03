#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""RDKit descriptor calculator. Reads a compounds DataStream JSON, looks up
SMILES via biopipelines_io, computes the requested descriptors (and optionally
a Morgan fingerprint), writes one row per input compound to the output CSV.
"""

import argparse
import json
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from biopipelines.biopipelines_io import load_datastream, iterate_values  # noqa: E402

from rdkit import Chem  # noqa: E402
from rdkit.Chem import AllChem, Descriptors  # noqa: E402
from rdkit.Chem.QED import qed as compute_qed  # noqa: E402


def get_descriptor_fn(name: str):
    """Resolve a descriptor function. Special-cases `qed`; otherwise delegates
    to Descriptors.<name>. Raises AttributeError if the name is unknown."""
    if name == "qed":
        return compute_qed
    return getattr(Descriptors, name)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--config", required=True)
    args = p.parse_args()

    with open(args.config, encoding="utf-8") as f:
        cfg = json.load(f)

    compounds_json = cfg["compounds_json"]
    output_csv = cfg["output_csv"]
    descriptor_names = cfg["descriptors"]
    morgan_fp = bool(cfg.get("morgan_fp", False))

    ds = load_datastream(compounds_json)

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
        sys.exit(1)

    os.makedirs(os.path.dirname(output_csv), exist_ok=True)
    columns = ["id", "smiles"] + descriptor_names
    if morgan_fp:
        columns += [f"fp_{i}" for i in range(1024)]
    pd.DataFrame(rows, columns=columns).to_csv(output_csv, index=False)
    print(f"Descriptors: {output_csv} ({len(rows)} rows, {len(columns)} columns)")

    if failed:
        print(f"Failed: {len(failed)}/{len(failed)+len(rows)}: {failed}", file=sys.stderr)


if __name__ == "__main__":
    main()
