#Copyright © 2024 LOCBP @ University of Zürich
#Distributed under MIT license
"""
Extract confidence metrics from GhostFold/ColabFold output files.

GhostFold uses ColabFold internally, so the output format follows ColabFold
conventions. This script searches for JSON score files and extracts pLDDT,
max_pae, and ptm metrics from the best-ranked models.
"""
import argparse
import json
import pandas as pd
import os
import glob
import numpy as np

parser = argparse.ArgumentParser(description='Extract confidence metrics from GhostFold/ColabFold output')
parser.add_argument('output_folder', type=str, help='Path to GhostFold output folder')
parser.add_argument('confidence_csv', type=str, help='Output CSV file path for confidence metrics')

args = parser.parse_args()


def find_score_files(output_folder):
    """
    Find all ColabFold score JSON files in the output folder.

    GhostFold may organize outputs in different ways:
    - predictions/<seq_id>/<seq_id>_scores_rank_001_*.json
    - Folding/<seq_id>_scores_rank_001_*.json
    - <seq_id>_scores_rank_001_*.json

    Returns:
        List of (seq_id, json_path, pdb_path) tuples
    """
    score_files = []

    # Search patterns for ColabFold score files
    search_patterns = [
        os.path.join(output_folder, "predictions", "*", "*_scores_rank_001_*.json"),
        os.path.join(output_folder, "Folding", "*_scores_rank_001_*.json"),
        os.path.join(output_folder, "*_scores_rank_001_*.json"),
    ]

    found_files = set()
    for pattern in search_patterns:
        for json_file in glob.glob(pattern):
            if json_file not in found_files:
                found_files.add(json_file)

                # Extract sequence ID from filename
                basename = os.path.basename(json_file)
                seq_id = basename.split('_scores_rank_')[0]

                # Find corresponding PDB file (in main output folder)
                pdb_path = os.path.join(output_folder, f"{seq_id}.pdb")

                score_files.append((seq_id, json_file, pdb_path))

    return score_files


def extract_ghostfold_confidence(output_folder, confidence_csv):
    """
    Extract confidence metrics from GhostFold/ColabFold JSON files.

    Processes all rank_001 structures (best models) and extracts:
    - pLDDT (average per-residue confidence)
    - max_pae (maximum predicted aligned error)
    - ptm (predicted template modeling score)

    Args:
        output_folder: Path to GhostFold output folder
        confidence_csv: Output CSV file path
    """
    confidence_data = []
    score_files = find_score_files(output_folder)

    if not score_files:
        print(f"Warning: No score files found in {output_folder}")
        # Create empty CSV with headers
        df = pd.DataFrame(columns=['id', 'structure', 'plddt', 'max_pae', 'ptm'])
        df.to_csv(confidence_csv, index=False)
        return

    for seq_id, json_file, pdb_path in score_files:
        try:
            with open(json_file, 'r') as f:
                metrics = json.load(f)
        except Exception as e:
            print(f"Error reading JSON file {json_file}: {e}")
            continue

        # Extract metrics
        plddt_array = metrics.get('plddt', [])
        max_pae = metrics.get('max_pae', None)
        ptm = metrics.get('ptm', None)

        # Calculate average pLDDT
        avg_plddt = np.mean(plddt_array) if plddt_array else None

        confidence_row = {
            'id': seq_id,
            'structure': pdb_path,
            'plddt': round(avg_plddt, 2) if avg_plddt is not None else None,
            'max_pae': round(max_pae, 2) if max_pae is not None else None,
            'ptm': round(ptm, 4) if ptm is not None else None
        }
        confidence_data.append(confidence_row)

    # Save to CSV
    if confidence_data:
        df = pd.DataFrame(confidence_data)
        df.to_csv(confidence_csv, index=False)
        print(f"Extracted confidence metrics for {len(confidence_data)} structures")
    else:
        print("Warning: No confidence data extracted")
        df = pd.DataFrame(columns=['id', 'structure', 'plddt', 'max_pae', 'ptm'])
        df.to_csv(confidence_csv, index=False)


if __name__ == "__main__":
    extract_ghostfold_confidence(args.output_folder, args.confidence_csv)
