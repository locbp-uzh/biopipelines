#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Runtime helper script for ADMET-AI predictions.

Reads SMILES from the compounds DataStream's map_table, runs the full batch
through ``ADMETModel.predict()`` once, and writes a CSV with one row per
input compound (columns: id, smiles, <ADMET endpoints>).
"""

import os
import sys
import json
import argparse

# Drop any inherited MPLBACKEND before matplotlib is imported transitively
# (admet_ai -> chemfunc -> matplotlib.pyplot). Colab leaks
# 'module://matplotlib_inline.backend_inline' into child processes, which the
# admet_ai env doesn't have the matplotlib_inline shim for, so matplotlib's
# import crashes. 'agg' is headless and always available.
os.environ["MPLBACKEND"] = "agg"

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from biopipelines.biopipelines_io import load_datastream, iterate_values


def main():
    parser = argparse.ArgumentParser(description="Run ADMET-AI on a compounds DataStream")
    parser.add_argument("--config", required=True, help="Path to admet_ai_config.json")
    args = parser.parse_args()

    with open(args.config) as f:
        config = json.load(f)

    compounds_ds = load_datastream(config["compounds_json"])
    output_csv = config["output_csv"]

    # Collect SMILES in input order. iterate_values reads from map_table.
    rows = []
    for comp_id, values in iterate_values(compounds_ds, columns=["smiles"]):
        smiles = values["smiles"]
        if smiles is None or (isinstance(smiles, float) and pd.isna(smiles)) or str(smiles).strip() == "":
            print(f"WARNING: skipping {comp_id} (empty smiles)", file=sys.stderr)
            continue
        rows.append({"id": comp_id, "smiles": str(smiles)})

    if not rows:
        print("ERROR: no compounds with valid SMILES", file=sys.stderr)
        sys.exit(1)

    # Lazy import — keeps tool wrapper importable in environments where
    # admet_ai isn't installed (config-time on cluster head node, etc.).
    from admet_ai import ADMETModel

    smiles_list = [r["smiles"] for r in rows]
    model = ADMETModel()
    preds = model.predict(smiles=smiles_list)

    # ADMETModel.predict returns a DataFrame indexed by SMILES with one
    # column per endpoint. Join back on input order so 'id' aligns with
    # the original compounds stream rather than the (possibly de-duplicated)
    # SMILES index.
    preds = preds.reset_index().rename(columns={"index": "smiles"})
    input_df = pd.DataFrame(rows)
    out = input_df.merge(preds, on="smiles", how="left")

    os.makedirs(os.path.dirname(output_csv), exist_ok=True)
    out.to_csv(output_csv, index=False)
    print(f"Wrote {len(out)} predictions to {output_csv}")


if __name__ == "__main__":
    main()
