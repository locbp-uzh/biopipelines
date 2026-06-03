#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Emit the Boltz2 compounds stream with the residue code Boltz assigns.

Boltz collapses every ligand to the residue name ``LIG`` in PDB output; in mmcif
it keeps a CCD ligand's code (``ATP`` stays ``ATP``) but still writes ``LIG`` for
a SMILES ligand. This script reads every ligand source listed in the combinatorics
config (so bundled ligands are all represented), preserves each id/SMILES, and
sets the ``code`` column to the code Boltz actually uses, so downstream
HETATM-selector tools read the right code automatically.

The map is the sole source of the emitted compounds stream's content (map_table
contract): the wrapper only declared its path at config time.
"""

import argparse
import json
import os
import sys

import pandas as pd

_biopipelines_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'biopipelines')
sys.path.insert(0, _biopipelines_dir)
from ligand_utils import auth_ligand_field

SMILES_RESIDUE_CODE = "LIG"

OUTPUT_COLUMNS = ["id", "format", "code", "smiles", "ccd"]


def assigned_code(row, output_format: str) -> str:
    """The residue code Boltz writes for this ligand in `output_format`."""
    field, value = auth_ligand_field(row)
    # mmcif preserves a CCD ligand's code; everything else collapses to LIG.
    if field == "ccd" and output_format == "mmcif":
        return value
    return SMILES_RESIDUE_CODE


def ligand_source_paths(combinatorics_config: str):
    """Every ligand source CSV listed in the combinatorics config, in order."""
    with open(combinatorics_config) as fh:
        cfg = json.load(fh)
    axes = cfg.get("axes", {})
    paths = []
    for axis in axes.values():
        if axis.get("entity_type") != "ligand":
            continue
        for src in axis.get("sources", []):
            p = src.get("path")
            if p and p not in paths:
                paths.append(p)
    return paths


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--combinatorics-config", required=True, help="combinatorics config JSON")
    p.add_argument("--output-format", required=True, help="boltz output_format (pdb or mmcif)")
    p.add_argument("--output-compounds", required=True, help="emitted compounds map CSV")
    args = p.parse_args()

    if not os.path.exists(args.combinatorics_config):
        print(f"ERROR: combinatorics config not found: {args.combinatorics_config}", file=sys.stderr)
        sys.exit(1)

    source_paths = ligand_source_paths(args.combinatorics_config)
    if not source_paths:
        print("No ligand sources in combinatorics config; nothing to emit")
        return

    frames = []
    for path in source_paths:
        if not os.path.exists(path):
            print(f"ERROR: ligand source CSV not found: {path}", file=sys.stderr)
            sys.exit(1)
        frames.append(pd.read_csv(path))
    df = pd.concat(frames, ignore_index=True)
    if "id" not in df.columns:
        print("ERROR: ligand source CSV has no 'id' column", file=sys.stderr)
        sys.exit(1)

    out = pd.DataFrame()
    out["id"] = df["id"]
    out["format"] = df["format"] if "format" in df.columns else ""
    out["code"] = df.apply(lambda r: assigned_code(r, args.output_format), axis=1)
    out["smiles"] = df["smiles"] if "smiles" in df.columns else ""
    out["ccd"] = df["ccd"] if "ccd" in df.columns else ""

    os.makedirs(os.path.dirname(args.output_compounds), exist_ok=True)
    out[OUTPUT_COLUMNS].to_csv(args.output_compounds, index=False)
    print(f"Compounds map: {args.output_compounds} ({len(out)} rows)")
    for _, r in out.iterrows():
        print(f"  {r['id']}: code={r['code']}")


if __name__ == "__main__":
    main()
