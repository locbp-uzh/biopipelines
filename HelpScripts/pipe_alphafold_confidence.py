# Copyright (c) 2026 Gianluca Quargnali @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.


import argparse
import json
import pandas as pd
import os
import glob
import numpy as np

parser = argparse.ArgumentParser(description='Extract confidence metrics from AlphaFold/ColabFold JSON files')
parser.add_argument('folding_folder', type=str, help='Path to AlphaFold Folding folder containing output files')
parser.add_argument('output_folder', type=str, help='Path to output folder for main structures')
parser.add_argument('confidence_csv', type=str, help='Output CSV file path for confidence metrics')

# Parse the arguments
args = parser.parse_args()

def extract_alphafold_confidence(folding_folder, output_folder, confidence_csv):
    """
    Extract confidence metrics from AlphaFold/ColabFold JSON files.

    Processes all rank_001 structures (best models) and extracts:
    - pLDDT (average per-residue confidence)
    - max_pae (maximum predicted aligned error)
    - ptm (predicted template modeling score)

    Args:
        folding_folder: Path to Folding subfolder with raw ColabFold outputs
        output_folder: Path to main output folder where simplified structures are stored
        confidence_csv: Output CSV file path
    """
    confidence_data = []

    # Look for unrelaxed or relaxed rank_001 PDB files (best models)
    for pdb_file in glob.glob(os.path.join(folding_folder, "*_rank_001_*.pdb")):
        # Extract base name and corresponding JSON file
        # ColabFold naming:
        # PDB:  <id>_unrelaxed_rank_001_alphafold2_ptm_model_2_seed_000.pdb
        # JSON: <id>_scores_rank_001_alphafold2_ptm_model_2_seed_000.json
        # Replace _unrelaxed_ or _relaxed_ with _scores_ to get JSON filename
        ##### ERROR THREE DEFINITIONS
        json_file = pdb_file.replace('.pdb', '.json')
        json_file = json_file.replace('_unrelaxed_', '_scores_')
        json_file = json_file.replace('_relaxed_', '_scores_')

        if not os.path.exists(json_file):
            print(f"Warning: JSON file not found for {pdb_file}")
            continue

        # Extract sequence ID from filename
        basename = os.path.basename(pdb_file)
        if '_relaxed_rank_' in basename:
            seq_id = basename.split('_relaxed_rank_')[0]
        elif '_unrelaxed_rank_' in basename:
            seq_id = basename.split('_unrelaxed_rank_')[0]
        else:
            print(f"Warning: Could not parse sequence ID from {basename}")
            continue

        # Read JSON file
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

        # Structure path (simplified name in main folder)
        structure_path = os.path.join(output_folder, f"{seq_id}.pdb")

        confidence_row = {
            'id': seq_id,
            'structure': structure_path,
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

# Run the extraction
if __name__ == "__main__":
    extract_alphafold_confidence(args.folding_folder, args.output_folder, args.confidence_csv)
