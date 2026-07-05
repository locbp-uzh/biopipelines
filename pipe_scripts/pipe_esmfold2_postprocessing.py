#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
ESMFold2 post-processing.

Reads the per-id scores JSON written by the inference script and the predicted
mmCIFs, then writes:
  - structures_map.csv  (id, file, + {axis}.id provenance columns)
  - confidence.csv      (id, file, plddt, ptm, iptm[, max_pae])
  - compounds_map.csv   (ligand chemistry passthrough, when ligands were folded)

Runs under the 'biopipelines' env. Drops ids whose mmCIF is missing so the map
stays honest under partial failure.
"""

import argparse
import json
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))


def provenance_lookup(config):
    """{base_id: {axis: axis_id}} from the combinatorics config's provenance block."""
    predicted = config.get("predicted_ids", [])
    provenance = config.get("provenance", {})
    out = {}
    for i, pid in enumerate(predicted):
        out[pid] = {axis: vals[i] for axis, vals in provenance.items() if i < len(vals)}
    return out


def base_id(row_id, prov):
    """Strip a trailing _<digits> sample suffix to match provenance keys."""
    if row_id in prov:
        return row_id
    base, sep, suffix = row_id.rpartition("_")
    if sep and suffix.isdigit() and base in prov:
        return base
    return row_id


def main():
    p = argparse.ArgumentParser(description="ESMFold2 post-processing")
    p.add_argument("--combinatorics-config", required=True)
    p.add_argument("--predictions-dir", required=True)
    p.add_argument("--structures-map", required=True)
    p.add_argument("--confidence-csv", required=True)
    p.add_argument("--ligands-table", default=None)
    p.add_argument("--compounds-map", default=None)
    args = p.parse_args()

    with open(args.combinatorics_config) as f:
        config = json.load(f)
    prov = provenance_lookup(config)
    axes = list(config.get("provenance", {}).keys())

    scores_path = os.path.join(args.predictions_dir, "ESMFold2_scores.json")
    if not os.path.exists(scores_path):
        print(f"ERROR: scores JSON not found: {scores_path}", file=sys.stderr)
        sys.exit(1)
    with open(scores_path) as f:
        scores = json.load(f)

    struct_rows, conf_rows = [], []
    for s in scores:
        cid = s["id"]
        cif = os.path.join(args.predictions_dir, f"{cid}.cif")
        if not os.path.exists(cif):
            print(f"WARNING: missing mmCIF for {cid}, dropping", file=sys.stderr)
            continue
        srow = {"id": cid, "file": cif}
        for axis in axes:
            srow[f"{axis}.id"] = prov.get(base_id(cid, prov), {}).get(axis, "")
        struct_rows.append(srow)

        crow = {"id": cid, "file": cif, "plddt": s.get("plddt"),
                "ptm": s.get("ptm"), "iptm": s.get("iptm")}
        if "max_pae" in s:
            crow["max_pae"] = s["max_pae"]
        conf_rows.append(crow)

    if not struct_rows:
        print("ERROR: no predicted structures found", file=sys.stderr)
        sys.exit(1)

    os.makedirs(os.path.dirname(args.structures_map), exist_ok=True)
    pd.DataFrame(struct_rows).to_csv(args.structures_map, index=False)
    os.makedirs(os.path.dirname(args.confidence_csv), exist_ok=True)
    pd.DataFrame(conf_rows).to_csv(args.confidence_csv, index=False)
    print(f"Wrote {len(struct_rows)} structures to {args.structures_map}")

    if args.ligands_table and args.compounds_map:
        if not os.path.exists(args.ligands_table):
            print(f"ERROR: ligands table not found: {args.ligands_table}", file=sys.stderr)
            sys.exit(1)
        df = pd.read_csv(args.ligands_table)
        out = pd.DataFrame()
        out["id"] = df["id"]
        out["format"] = df["format"] if "format" in df.columns else ""
        # ESMFold2 writes mmCIF, which preserves CCD codes; pass the input code through.
        out["code"] = df["code"] if "code" in df.columns else ""
        out["smiles"] = df["smiles"] if "smiles" in df.columns else ""
        out["ccd"] = df["ccd"] if "ccd" in df.columns else ""
        os.makedirs(os.path.dirname(args.compounds_map), exist_ok=True)
        out.to_csv(args.compounds_map, index=False)
        print(f"Wrote {len(out)} compounds to {args.compounds_map}")


if __name__ == "__main__":
    main()
