#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Runtime helper script for ADMET-AI property prediction.

Loads compounds via DataStream, extracts SMILES from the map_table,
runs ADMET-AI prediction, and writes output CSVs:
  - predictions.csv: ADMET + physicochemical properties per compound
  - drugbank.csv: DrugBank percentile comparison (optional)
"""

import os
import sys
import argparse
import json
import pandas as pd
from typing import Dict, Any

# Import unified I/O utilities
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from biopipelines.biopipelines_io import load_datastream, iterate_values


def run_admet_ai(config_data: Dict[str, Any]) -> None:
    """
    Run ADMET-AI prediction on all compounds.

    Args:
        config_data: Configuration dictionary with paths and options
    """
    from admet_ai import ADMETModel

    compounds_ds = load_datastream(config_data['compounds_json'])
    predictions_csv = config_data['predictions_csv']
    drugbank_csv = config_data['drugbank_csv']
    drugbank = config_data.get('drugbank', True)

    print(f"Running ADMET-AI property prediction")
    print(f"Compounds: {len(compounds_ds.ids_expanded)}")

    # Extract SMILES from map_table
    ids = []
    smiles_list = []
    for comp_id, values in iterate_values(compounds_ds, columns=["smiles"]):
        ids.append(comp_id)
        smiles_list.append(values["smiles"])

    print(f"Extracted {len(smiles_list)} SMILES strings")

    output_dir = os.path.dirname(predictions_csv)
    os.makedirs(output_dir, exist_ok=True)

    if len(smiles_list) == 0:
        print("No compounds to process. Writing empty output files.")
        pd.DataFrame(columns=["id", "smiles"]).to_csv(predictions_csv, index=False)
        if drugbank:
            pd.DataFrame(columns=["id", "smiles"]).to_csv(drugbank_csv, index=False)
        return

    # Run predictions
    print("Loading ADMET-AI model...")
    model = ADMETModel()

    print("Predicting ADMET properties...")
    preds_df = model.predict(smiles=smiles_list)

    # Insert id and smiles columns at the front
    preds_df.insert(0, "smiles", smiles_list)
    preds_df.insert(0, "id", ids)

    preds_df.to_csv(predictions_csv, index=False)
    print(f"Predictions written to: {predictions_csv}")
    print(f"Properties predicted: {len(preds_df.columns) - 2}")
    print(preds_df)

    # DrugBank percentile comparison
    if drugbank:
        print("\nComputing DrugBank percentile comparison...")
        drugbank_df = model.predict(smiles=smiles_list, drugbank_comparison=True)

        drugbank_df.insert(0, "smiles", smiles_list)
        drugbank_df.insert(0, "id", ids)

        drugbank_df.to_csv(drugbank_csv, index=False)
        print(f"DrugBank comparison written to: {drugbank_csv}")
        print(drugbank_df)

    print(f"\nADMET-AI prediction completed successfully!")


def main():
    parser = argparse.ArgumentParser(description='ADMET-AI property prediction')
    parser.add_argument('--config', required=True, help='JSON config file with prediction parameters')

    args = parser.parse_args()

    if not os.path.exists(args.config):
        print(f"Error: Config file not found: {args.config}")
        sys.exit(1)

    try:
        with open(args.config, 'r') as f:
            config_data = json.load(f)
    except Exception as e:
        print(f"Error loading config: {e}")
        sys.exit(1)

    required_params = ['compounds_json', 'predictions_csv', 'drugbank_csv']
    for param in required_params:
        if param not in config_data:
            print(f"Error: Missing required parameter: {param}")
            sys.exit(1)

    try:
        run_admet_ai(config_data)
    except Exception as e:
        print(f"Error running ADMET-AI prediction: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
