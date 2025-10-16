#!/usr/bin/env python3
"""
Runtime helper script for RF3 tool.

Processes RF3 configuration, generates JSON input files, runs RF3 predictions in batch,
and post-processes results into standardized datasheets.
"""

import os
import sys
import argparse
import json
import pandas as pd
import subprocess
from pathlib import Path
from typing import Dict, List, Any, Tuple


def load_sequences_from_datasheet(datasheet_path: str) -> List[Tuple[str, str]]:
    """
    Load sequences from datasheet.

    Returns:
        List of (id, sequence) tuples
    """
    df = pd.read_csv(datasheet_path)
    if 'id' not in df.columns or 'sequence' not in df.columns:
        raise ValueError(f"Datasheet must have 'id' and 'sequence' columns: {datasheet_path}")

    return list(zip(df['id'], df['sequence']))


def load_compounds_from_datasheet(datasheet_path: str) -> List[Tuple[str, str]]:
    """
    Load compounds from datasheet.

    Returns:
        List of (id, smiles) tuples
    """
    df = pd.read_csv(datasheet_path)
    if 'id' not in df.columns or 'smiles' not in df.columns:
        raise ValueError(f"Datasheet must have 'id' and 'smiles' columns: {datasheet_path}")

    return list(zip(df['id'], df['smiles']))


def create_rf3_input_json(sequence_id: str, sequence: str, ligand_smiles: str = None,
                          output_path: str = None) -> str:
    """
    Create RF3 input JSON file.

    Args:
        sequence_id: Identifier for the sequence
        sequence: Protein sequence
        ligand_smiles: Optional SMILES string for ligand
        output_path: Path to save JSON file

    Returns:
        Path to created JSON file
    """
    input_data = {
        "job": sequence_id,
        "sequences": [
            {
                "protein": {
                    "id": ["A"],
                    "sequence": [sequence]
                }
            }
        ]
    }

    # Add ligand if provided
    if ligand_smiles:
        input_data["sequences"][0]["sm_string"] = [ligand_smiles]

    if output_path is None:
        output_path = f"{sequence_id}.json"

    with open(output_path, 'w') as f:
        json.dump(input_data, f, indent=2)

    return output_path


def run_rf3_batch(input_jsons: List[str], output_folder: str, config: Dict[str, Any]) -> bool:
    """
    Run RF3 predictions in batch mode.

    Args:
        input_jsons: List of JSON input file paths
        output_folder: Output directory
        config: Configuration dictionary

    Returns:
        True if successful, False otherwise
    """
    print(f"\nRunning RF3 on {len(input_jsons)} inputs")

    # Build RF3 command
    cmd = ["rf3", "fold"]

    # Add input files
    cmd.extend([f"inputs={json_file}" for json_file in input_jsons])

    # Add optional parameters
    if config.get("checkpoint_path"):
        cmd.append(f"ckpt_path={config['checkpoint_path']}")

    if config.get("early_stopping_plddt"):
        cmd.append(f"early_stopping_plddt_threshold={config['early_stopping_plddt']}")

    # Set output directory via environment or parameter
    env = os.environ.copy()
    env["RF3_OUTPUT_DIR"] = output_folder

    # Ensure temp directory exists
    # RF3 uses /sctmp/<user> on clusters, but we need to ensure it exists or use a fallback
    temp_base_dirs = ["/sctmp", os.environ.get("TMPDIR"), "/tmp"]
    for temp_base in temp_base_dirs:
        if temp_base and os.path.exists(temp_base):
            user_temp_dir = os.path.join(temp_base, os.environ.get("USER", "unknown"))
            os.makedirs(user_temp_dir, exist_ok=True)
            env["TMPDIR"] = user_temp_dir
            print(f"Using temp directory: {user_temp_dir}")
            break

    print(f"Running command: {' '.join(cmd)}")
    print(f"Output directory: {output_folder}")

    try:
        result = subprocess.run(
            cmd,
            cwd=output_folder,
            env=env,
            check=True,
            capture_output=True,
            text=True
        )
        print(result.stdout)
        if result.stderr:
            print("STDERR:", result.stderr, file=sys.stderr)
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error running RF3: {e}", file=sys.stderr)
        print(f"STDOUT: {e.stdout}", file=sys.stderr)
        print(f"STDERR: {e.stderr}", file=sys.stderr)
        return False
    except Exception as e:
        print(f"Unexpected error running RF3: {e}", file=sys.stderr)
        return False


def collect_rf3_outputs(output_folder: str, sequence_ids: List[str],
                       output_format: str) -> Tuple[List[Dict], List[Dict]]:
    """
    Collect RF3 output files and metrics.

    Args:
        output_folder: Directory containing RF3 outputs
        sequence_ids: List of sequence IDs that were processed
        output_format: Expected output format ("pdb" or "cif")

    Returns:
        Tuple of (structures_data, confidence_data)
    """
    structures = []
    confidence = []

    extension = ".cif.gz" if output_format == "cif" else ".pdb"

    for seq_id in sequence_ids:
        # Look for model files
        # RF3 output naming: {job}_model_{n}.cif.gz
        model_files = list(Path(output_folder).glob(f"{seq_id}_model_*.cif.gz"))

        if not model_files and output_format == "pdb":
            # Try PDB format if CIF not found
            model_files = list(Path(output_folder).glob(f"{seq_id}_model_*.pdb"))

        # Look for metrics file
        metrics_file = os.path.join(output_folder, f"{seq_id}_metrics.csv")

        metrics_df = None
        if os.path.exists(metrics_file):
            try:
                metrics_df = pd.read_csv(metrics_file)
            except Exception as e:
                print(f"Warning: Could not read metrics for {seq_id}: {e}")

        for i, model_file in enumerate(sorted(model_files), 1):
            # Extract pLDDT score from metrics if available
            plddt_score = None
            ptm_score = None

            if metrics_df is not None and i <= len(metrics_df):
                model_metrics = metrics_df.iloc[i-1]
                plddt_score = model_metrics.get('plddt', None)
                ptm_score = model_metrics.get('ptm', None)

            structures.append({
                "id": seq_id,
                "model_id": i,
                "file_path": str(model_file),
                "plddt_score": plddt_score
            })

            confidence.append({
                "id": seq_id,
                "model_id": i,
                "plddt_score": plddt_score,
                "ptm_score": ptm_score
            })

    return structures, confidence


def process_rf3_predictions(config_data: Dict[str, Any]) -> int:
    """
    Main processing function for RF3 predictions.

    Args:
        config_data: Configuration dictionary

    Returns:
        0 if successful, 1 if failed
    """
    output_folder = config_data['output_folder']
    os.makedirs(output_folder, exist_ok=True)

    # Load protein sequences
    protein_input = config_data['proteins']
    sequences = []

    if protein_input:
        if protein_input.get('type') == 'tool_output':
            # Load from datasheets
            datasheets = config_data.get('input_datasheets', {})
            if 'sequences' in datasheets:
                sequences = load_sequences_from_datasheet(datasheets['sequences'])
            else:
                print("Error: No sequences datasheet found in input", file=sys.stderr)
                return 1
        elif protein_input.get('type') == 'string':
            # Single sequence
            sequences = [("protein_1", protein_input['value'])]
        elif protein_input.get('type') == 'list':
            # List of sequences
            sequences = [(f"protein_{i+1}", seq) for i, seq in enumerate(protein_input['values'])]

    if not sequences:
        print("Error: No sequences to process", file=sys.stderr)
        return 1

    print(f"Processing {len(sequences)} sequences")

    # Load ligands if provided
    ligands = []
    ligand_input = config_data.get('ligands')
    if ligand_input:
        if ligand_input.get('type') == 'tool_output':
            datasheets = config_data.get('input_datasheets', {})
            if 'compounds' in datasheets:
                ligands = load_compounds_from_datasheet(datasheets['compounds'])
        elif ligand_input.get('type') == 'string':
            ligands = [("ligand_1", ligand_input['value'])]

    # Create RF3 input JSON files
    json_inputs = []
    sequence_ids = []

    for seq_id, sequence in sequences:
        if ligands:
            # Create protein-ligand complexes
            for lig_id, smiles in ligands:
                complex_id = f"{seq_id}_{lig_id}"
                json_file = create_rf3_input_json(
                    complex_id, sequence, smiles,
                    os.path.join(output_folder, f"{complex_id}.json")
                )
                json_inputs.append(json_file)
                sequence_ids.append(complex_id)
        else:
            # Apo prediction
            json_file = create_rf3_input_json(
                seq_id, sequence, None,
                os.path.join(output_folder, f"{seq_id}.json")
            )
            json_inputs.append(json_file)
            sequence_ids.append(seq_id)

    print(f"Created {len(json_inputs)} input JSON files")

    # Run RF3 batch prediction
    success = run_rf3_batch(json_inputs, output_folder, config_data)

    if not success:
        print("RF3 prediction failed", file=sys.stderr)
        return 1

    # Collect and process outputs
    print("\nCollecting RF3 outputs...")
    structures_data, confidence_data = collect_rf3_outputs(
        output_folder,
        sequence_ids,
        config_data.get('output_format', 'pdb')
    )

    # Save datasheets
    structures_csv = os.path.join(output_folder, "structures.csv")
    confidence_csv = os.path.join(output_folder, "confidence_scores.csv")

    if structures_data:
        df_structures = pd.DataFrame(structures_data)
        df_structures.to_csv(structures_csv, index=False)
        print(f"Saved structures: {structures_csv} ({len(structures_data)} models)")
    else:
        print("Warning: No structures collected", file=sys.stderr)
        return 1

    if confidence_data:
        df_confidence = pd.DataFrame(confidence_data)
        df_confidence.to_csv(confidence_csv, index=False)
        print(f"Saved confidence scores: {confidence_csv}")

    print("\nRF3 processing completed successfully")
    return 0


def main():
    parser = argparse.ArgumentParser(description='RF3 structure prediction runtime script')
    parser.add_argument('--config', required=True, help='JSON config file with RF3 parameters')

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

    # Process RF3 predictions
    exit_code = process_rf3_predictions(config_data)
    sys.exit(exit_code)


if __name__ == "__main__":
    main()
