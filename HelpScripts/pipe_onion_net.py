#!/usr/bin/env python3
"""
Runtime helper script for OnionNet and OnionNet-2 tools.

Processes protein-ligand complexes and predicts binding affinities using
OnionNet or OnionNet-2 models.
"""

import os
import sys
import argparse
import json
import pandas as pd
import subprocess
from pathlib import Path
from typing import Dict, List, Any, Tuple


def load_structures_from_table(table_path: str) -> List[Tuple[str, str]]:
    """
    Load structures from table.

    Returns:
        List of (id, file_path) tuples
    """
    df = pd.read_csv(table_path)
    if 'id' not in df.columns or 'file_path' not in df.columns:
        raise ValueError(f"Table must have 'id' and 'file_path' columns: {table_path}")

    return list(zip(df['id'], df['file_path']))


def run_onionnet_v1(structures: List[Tuple[str, str]], output_folder: str,
                    config: Dict[str, Any]) -> Tuple[List[Dict], bool]:
    """
    Run OnionNet v1 predictions.

    Args:
        structures: List of (id, file_path) tuples
        output_folder: Output directory
        config: Configuration dictionary

    Returns:
        Tuple of (predictions_data, success)
    """
    print(f"\nRunning OnionNet v1 on {len(structures)} structures")

    predictions = []

    # Get model paths
    model_weights = config.get("model_weights")
    scaler_model = config.get("scaler_model")

    if not model_weights or not scaler_model:
        print("Error: model_weights and scaler_model are required for OnionNet v1", file=sys.stderr)
        return [], False

    # Create temporary input file list
    input_list_file = os.path.join(output_folder, "input_structures.dat")
    with open(input_list_file, 'w') as f:
        for struct_id, struct_path in structures:
            f.write(f"{struct_path}\n")

    # Generate features
    features_file = os.path.join(output_folder, "features.csv")
    print("Generating features...")

    try:
        # Assume OnionNet is installed and generate_features.py is available
        cmd_features = [
            "python", "-m", "onionnet.generate_features",
            "-inp", input_list_file,
            "-out", features_file
        ]

        result = subprocess.run(
            cmd_features,
            check=True,
            capture_output=True,
            text=True
        )
        print(result.stdout)
        if result.stderr:
            print("STDERR:", result.stderr, file=sys.stderr)

    except subprocess.CalledProcessError as e:
        print(f"Error generating features: {e}", file=sys.stderr)
        print(f"STDOUT: {e.stdout}", file=sys.stderr)
        print(f"STDERR: {e.stderr}", file=sys.stderr)
        return [], False
    except Exception as e:
        print(f"Unexpected error generating features: {e}", file=sys.stderr)
        return [], False

    # Predict affinities
    predictions_file = os.path.join(output_folder, "predictions.csv")
    print("Predicting affinities...")

    try:
        cmd_predict = [
            "python", "-m", "onionnet.predict",
            "-fn", features_file,
            "-weights", model_weights,
            "-scaler", scaler_model,
            "-out", predictions_file
        ]

        result = subprocess.run(
            cmd_predict,
            check=True,
            capture_output=True,
            text=True
        )
        print(result.stdout)
        if result.stderr:
            print("STDERR:", result.stderr, file=sys.stderr)

    except subprocess.CalledProcessError as e:
        print(f"Error predicting affinities: {e}", file=sys.stderr)
        print(f"STDOUT: {e.stdout}", file=sys.stderr)
        print(f"STDERR: {e.stderr}", file=sys.stderr)
        return [], False
    except Exception as e:
        print(f"Unexpected error predicting affinities: {e}", file=sys.stderr)
        return [], False

    # Parse predictions
    try:
        pred_df = pd.read_csv(predictions_file)
        for i, (struct_id, struct_path) in enumerate(structures):
            if i < len(pred_df):
                affinity = pred_df.iloc[i]['predicted_affinity']
                predictions.append({
                    "id": struct_id,
                    "structure_path": struct_path,
                    "predicted_affinity_pKa": affinity
                })
    except Exception as e:
        print(f"Error parsing predictions: {e}", file=sys.stderr)
        return [], False

    return predictions, True


def run_onionnet_v2(structures: List[Tuple[str, str]], output_folder: str,
                    config: Dict[str, Any]) -> Tuple[List[Dict], bool]:
    """
    Run OnionNet v2 predictions.

    Args:
        structures: List of (id, file_path) tuples
        output_folder: Output directory
        config: Configuration dictionary

    Returns:
        Tuple of (predictions_data, success)
    """
    print(f"\nRunning OnionNet-2 on {len(structures)} structures")

    predictions = []

    # Get model paths
    model_path = config.get("model_path")
    scaler_path = config.get("scaler_path")
    shells = config.get("shells", 62)

    if not model_path or not scaler_path:
        print("Error: model_path and scaler_path are required for OnionNet-2", file=sys.stderr)
        return [], False

    # Process each structure
    for struct_id, struct_path in structures:
        print(f"Processing {struct_id}...")

        # For OnionNet-2, need separate protein and ligand files
        # Assuming structure is a complex, we need to split it
        # This is a placeholder - actual implementation would need structure parsing

        output_file = os.path.join(output_folder, f"{struct_id}_prediction.txt")

        try:
            # Assume OnionNet-2 is installed and predict.py is available
            cmd = [
                "python", "-m", "onionnet2.predict",
                "-rec_fpath", struct_path,  # Would need to extract protein
                "-lig_fpath", struct_path,  # Would need to extract ligand
                "-model", model_path,
                "-scaler", scaler_path,
                "-shells", str(shells),
                "-out_fpath", output_file
            ]

            result = subprocess.run(
                cmd,
                check=True,
                capture_output=True,
                text=True
            )
            print(result.stdout)
            if result.stderr:
                print("STDERR:", result.stderr, file=sys.stderr)

            # Read prediction
            if os.path.exists(output_file):
                with open(output_file, 'r') as f:
                    affinity = float(f.read().strip())

                predictions.append({
                    "id": struct_id,
                    "structure_path": struct_path,
                    "predicted_affinity_pKa": affinity
                })

        except subprocess.CalledProcessError as e:
            print(f"Error predicting affinity for {struct_id}: {e}", file=sys.stderr)
            print(f"STDOUT: {e.stdout}", file=sys.stderr)
            print(f"STDERR: {e.stderr}", file=sys.stderr)
            continue
        except Exception as e:
            print(f"Unexpected error processing {struct_id}: {e}", file=sys.stderr)
            continue

    return predictions, len(predictions) > 0


def process_onionnet_predictions(config_data: Dict[str, Any], version: int) -> int:
    """
    Main processing function for OnionNet predictions.

    Args:
        config_data: Configuration dictionary
        version: OnionNet version (1 or 2)

    Returns:
        0 if successful, 1 if failed
    """
    output_folder = config_data['output_folder']
    os.makedirs(output_folder, exist_ok=True)

    # Load protein-ligand structures
    structure_input = config_data['structures']
    structures = []

    if structure_input:
        if structure_input.get('type') == 'tool_output':
            # Load from tables
            tables = config_data.get('input_tables', {})
            if 'structures' in tables:
                structures = load_structures_from_table(tables['structures'])
            else:
                print("Error: No structures table found in input", file=sys.stderr)
                return 1
        elif structure_input.get('type') == 'string':
            # Single structure file
            structures = [("structure_1", structure_input['value'])]
        elif structure_input.get('type') == 'list':
            # List of structure files
            structures = [(f"structure_{i+1}", path) for i, path in enumerate(structure_input['values'])]

    if not structures:
        print("Error: No structures to process", file=sys.stderr)
        return 1

    print(f"Processing {len(structures)} structures with OnionNet v{version}")

    # Run appropriate version
    if version == 1:
        predictions_data, success = run_onionnet_v1(structures, output_folder, config_data)
    elif version == 2:
        predictions_data, success = run_onionnet_v2(structures, output_folder, config_data)
    else:
        print(f"Error: Unknown OnionNet version: {version}", file=sys.stderr)
        return 1

    if not success or not predictions_data:
        print("OnionNet prediction failed", file=sys.stderr)
        return 1

    # Save predictions
    predictions_csv = os.path.join(output_folder, "affinity_predictions.csv")

    if predictions_data:
        df_predictions = pd.DataFrame(predictions_data)
        df_predictions.to_csv(predictions_csv, index=False)
        print(f"\nSaved predictions: {predictions_csv} ({len(predictions_data)} complexes)")
    else:
        print("Warning: No predictions generated", file=sys.stderr)
        return 1

    print(f"\nOnionNet v{version} processing completed successfully")
    return 0


def main():
    parser = argparse.ArgumentParser(description='OnionNet affinity prediction runtime script')
    parser.add_argument('--config', required=True, help='JSON config file with OnionNet parameters')
    parser.add_argument('--version', type=int, required=True, choices=[1, 2],
                       help='OnionNet version (1 or 2)')

    args = parser.parse_args()

    # Load configuration
    if not os.path.exists(args.config):
        print(f"Error: Config file not found: {args.config}", file=sys.stderr)
        sys.exit(1)

    try:
        with open(args.config, 'r') as f:
            config_data = json.load(f)
    except Exception as e:
        print(f"Error loading config: {e}", file=sys.stderr)
        sys.exit(1)

    # Process OnionNet predictions
    exit_code = process_onionnet_predictions(config_data, args.version)
    sys.exit(exit_code)


if __name__ == "__main__":
    main()
