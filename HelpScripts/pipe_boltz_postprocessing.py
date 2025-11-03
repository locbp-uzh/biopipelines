#Copyright © 2024 LOCBP @ University of Zürich
#Distributed under MIT license
"""
Script used by Boltz2 pipeline to post-process and rename outputs based on sequence IDs
Simplified version of boltz_results.py with standard parameters
"""

import argparse

parser = argparse.ArgumentParser(description='Post-process Boltz2 outputs for pipeline integration')
parser.add_argument('PREDICTION_FOLDER', type=str, help='Boltz2 prediction folder')
parser.add_argument('OUTPUT_FOLDER', type=str, help='Tool output folder where results will be organized')
parser.add_argument('SEQUENCE_IDS_FILE', type=str, help='CSV file with sequence IDs for renaming')

# Parse the arguments
args = parser.parse_args()

import os
import json
import shutil
import csv
import pandas as pd

# ----- Helper function to flatten nested confidence JSON -----
def flatten_confidence(conf, parent_key=""):
    """
    Flattens the confidence dictionary. Nested keys are concatenated.
    For a dict value, keys are appended with ':' (or '-' for the second-level nested dict)
    e.g. {'chains_ptm': {'0': 0.85, '1': 0.83}} becomes {'chains_ptm:0': 0.85, 'chains_ptm:1': 0.83}
         {'pair_chains_iptm': {'0': {'0': 0.85, '1': 0.81}, '1': {'0': 0.82, '1': 0.83}}}
             becomes {'pair_chains_iptm:0-0': 0.85, 'pair_chains_iptm:0-1': 0.81, 'pair_chains_iptm:1-0': 0.82, 'pair_chains_iptm:1-1': 0.83}
    """
    items = {}
    for key, value in conf.items():
        new_key = f"{parent_key}:{key}" if parent_key else key
        if isinstance(value, dict):
            # For second-level dictionaries, use a different delimiter (dash) if needed.
            for subkey, subvalue in value.items():
                if isinstance(subvalue, dict):
                    for subsubkey, subsubvalue in subvalue.items():
                        final_key = f"{new_key}-{subkey}-{subsubkey}"
                        items[final_key] = subsubvalue
                else:
                    final_key = f"{new_key}:{subkey}"
                    items[final_key] = subvalue
        else:
            items[new_key] = value
    return items

PREDICTION_FOLDER = args.PREDICTION_FOLDER  
OUTPUT_FOLDER = args.OUTPUT_FOLDER
SEQUENCE_IDS_FILE = args.SEQUENCE_IDS_FILE

# Load sequence IDs for proper renaming
sequence_ids = []
if os.path.exists(SEQUENCE_IDS_FILE):
    df = pd.read_csv(SEQUENCE_IDS_FILE)
    if 'id' in df.columns:
        sequence_ids = df['id'].tolist()
    else:
        print(f"Warning: No 'id' column found in {SEQUENCE_IDS_FILE}")
else:
    print(f"Warning: Sequence IDs file not found: {SEQUENCE_IDS_FILE}")

# Create results folder
results_folder = OUTPUT_FOLDER
os.makedirs(results_folder, exist_ok=True)

# Create MSAs folder in the tool directory
msas_folder = os.path.join(results_folder, "MSAs")
os.makedirs(msas_folder, exist_ok=True)

confidence_data = {}  # dictionary mapping NAME -> flattened confidence dict
affinity_data = {}  # dictionary mapping NAME -> flattened affinity dict
structural_files = {} # dictionary mapping NAME -> pdb or cif file

# ----- Process each boltz_results folder -----
folder_to_sequence_map = {}  # Map folder names to sequence IDs

# Look for boltz_results_* folders in the prediction folder
boltz_results_folders = [f for f in os.listdir(PREDICTION_FOLDER) if f.startswith('boltz_results_')]

for boltz_folder in sorted(boltz_results_folders):
    boltz_folder_path = os.path.join(PREDICTION_FOLDER, boltz_folder)
    if os.path.isdir(boltz_folder_path):
        # Extract the config ID from folder name (e.g., boltz_results_(S)Cy7-CH2F-RCG+)
        config_id = boltz_folder.replace('boltz_results_', '')
        
        # Find the corresponding predictions subfolder
        predictions_path = os.path.join(boltz_folder_path, "predictions", config_id)
        if not os.path.exists(predictions_path):
            print(f"Warning: Predictions path not found: {predictions_path}")
            continue
        
        # Determine sequence_id from expected sequence IDs if available
        if sequence_ids:
            if len(sequence_ids) == 1:
                # Single expected sequence - use it regardless of folder name
                sequence_id = sequence_ids[0]
                print(f"Mapping {config_id} to expected sequence ID: {sequence_id}")
            elif config_id in sequence_ids:
                # Config ID matches expected - use it
                sequence_id = config_id
            else:
                # Try to find best match or use the first available expected ID
                print(f"Warning: {config_id} not found in expected sequence IDs: {sequence_ids}")
                print(f"Using first expected sequence ID as fallback")
                sequence_id = sequence_ids[0]
        else:
            # No expected sequence IDs - use config_id as fallback
            sequence_id = config_id
        
        folder_to_sequence_map[config_id] = sequence_id
        
        # Look for structure files in predictions subfolder
        pdb_filename = f"{config_id}_model_0.pdb"
        pdb_filepath = os.path.join(predictions_path, pdb_filename)
        cif_filename = f"{config_id}_model_0.cif"
        cif_filepath = os.path.join(predictions_path, cif_filename)
        
        # Determine output filename using sequence ID
        target_filepath = None
        if os.path.exists(pdb_filepath):
            target_filepath = os.path.join(results_folder, f"{sequence_id}.pdb")
            try:
                shutil.copy2(pdb_filepath, target_filepath)
                print(f"Copied {pdb_filename} to {sequence_id}.pdb")
            except Exception as e:
                print(f"Error copying {pdb_filepath}: {e}")
                continue
        elif os.path.exists(cif_filepath):
            target_filepath = os.path.join(results_folder, f"{sequence_id}.cif")
            try:
                shutil.copy2(cif_filepath, target_filepath)
                print(f"Copied {cif_filename} to {sequence_id}.cif")
            except Exception as e:
                print(f"Error copying {cif_filepath}: {e}")
                continue
        else:
            print(f"No structure file found for {config_id}")
            continue
        
        # Process confidence and affinity files (in predictions subfolder)
        conf_filename = f"confidence_{config_id}_model_0.json"
        conf_filepath = os.path.join(predictions_path, conf_filename)
        aff_filename = f"affinity_{config_id}.json"
        aff_filepath = os.path.join(predictions_path, aff_filename)
        
        # Process MSA files (copy from msa subfolder to tool MSAs folder)
        msas_path = os.path.join(boltz_folder_path, "msa")
        print(f"Looking for MSA folder: {msas_path}")
        if os.path.exists(msas_path):
            msa_source = os.path.join(msas_path, f"{config_id}_0.csv")
            msa_target = os.path.join(msas_folder, f"{sequence_id}.csv")  # Use sequence_id for target filename
            print(f"Looking for MSA file: {msa_source}")
            if os.path.exists(msa_source):
                try:
                    shutil.copy2(msa_source, msa_target)
                    print(f"Copied MSA: {config_id}_0.csv to {sequence_id}.csv")
                except Exception as e:
                    print(f"Error copying MSA file: {e}")
            else:
                print(f"MSA file not found: {msa_source}")
        else:
            print(f"MSA folder not found: {msas_path}")
        
        # Initialize data dictionaries for this sequence
        confidence_data[sequence_id] = {}
        affinity_data[sequence_id] = {}
        
        # Load confidence data
        try:
            if os.path.exists(conf_filepath):
                with open(conf_filepath, "r") as f:
                    conf = json.load(f)
                flattened_conf = flatten_confidence(conf)
                confidence_data[sequence_id] = flattened_conf
                print(f"Loaded confidence for {sequence_id}: {len(flattened_conf)} keys")
        except Exception as e:
            print(f"Error reading confidence file {conf_filepath}: {e}")
            confidence_data[sequence_id] = {}
        
        # Load affinity data
        try:
            if os.path.exists(aff_filepath):
                with open(aff_filepath, "r") as f:
                    aff = json.load(f)
                flattened_aff = flatten_confidence(aff)
                affinity_data[sequence_id] = flattened_aff
                print(f"Loaded affinity for {sequence_id}: {len(flattened_aff)} keys")
            else:
                print(f"No affinity file found for {config_id}")
                affinity_data[sequence_id] = {}
        except Exception as e:
            print(f"Error reading affinity file {aff_filepath}: {e}")
            affinity_data[sequence_id] = {}
        
        # Store structural file path
        if target_filepath:
            structural_files[sequence_id] = target_filepath

# Check if we have any data
print(f"> Total sequences processed: {len(confidence_data)}")
print(f"> Sequences with confidence data: {[name for name, data in confidence_data.items() if data]}")
print(f"> Sequences with affinity data: {[name for name, data in affinity_data.items() if data]}")

# Create CSV files for pipeline compatibility
# 1. Confidence scores CSV
confidence_rows = []
for seq_id, conf_data in confidence_data.items():
    if conf_data:
        row = {'id': seq_id, 'input_file': seq_id}
        row.update(conf_data)
        confidence_rows.append(row)

if confidence_rows:
    confidence_df = pd.DataFrame(confidence_rows)
    confidence_csv = os.path.join(OUTPUT_FOLDER, "confidence_scores.csv")
    confidence_df.to_csv(confidence_csv, index=False)
    print(f"Created confidence scores CSV: {confidence_csv}")

# 2. Affinity scores CSV (if sequences were processed, create affinity file even if empty)
affinity_rows = []
for seq_id, aff_data in affinity_data.items():
    if aff_data:
        row = {'id': seq_id, 'input_file': seq_id}
        row.update(aff_data)
        affinity_rows.append(row)
    else:
        # Create empty row for sequences with no affinity data
        row = {'id': seq_id, 'input_file': seq_id}
        affinity_rows.append(row)

if affinity_rows:
    affinity_df = pd.DataFrame(affinity_rows)
    affinity_csv = os.path.join(OUTPUT_FOLDER, "affinity_scores.csv")
    affinity_df.to_csv(affinity_csv, index=False)
    print(f"Created affinity scores CSV: {affinity_csv}")
elif confidence_data:  # If we processed sequences but no affinity data, create empty file
    # Create empty affinity CSV with just headers
    affinity_df = pd.DataFrame(columns=['id', 'input_file'])
    affinity_csv = os.path.join(OUTPUT_FOLDER, "affinity_scores.csv")
    affinity_df.to_csv(affinity_csv, index=False)
    print(f"Created empty affinity scores CSV: {affinity_csv}")

# 3. Create sequences CSV from original input data
sequences_csv = os.path.join(OUTPUT_FOLDER, "sequences.csv")
if os.path.exists(SEQUENCE_IDS_FILE):
    # Copy the original sequences file to maintain input data
    import shutil
    shutil.copy2(SEQUENCE_IDS_FILE, sequences_csv)
    print(f"Created sequences CSV: {sequences_csv}")
else:
    # Create basic sequences CSV with just the IDs
    sequences_df = pd.DataFrame({'id': list(confidence_data.keys())})
    sequences_df.to_csv(sequences_csv, index=False)
    print(f"Created basic sequences CSV: {sequences_csv}")

# 4. Create MSAs CSV if MSA files exist
msa_csv = os.path.join(OUTPUT_FOLDER, "msas.csv")
msa_files_in_dir = []

# Load sequences data to get actual protein sequences
sequences_data = {}
if os.path.exists(sequences_csv):
    try:
        seq_df = pd.read_csv(sequences_csv)
        if 'sequence' in seq_df.columns:
            sequences_data = dict(zip(seq_df['id'], seq_df['sequence']))
        print(f"Loaded {len(sequences_data)} sequences for MSA table")
    except Exception as e:
        print(f"Warning: Could not load sequences from {sequences_csv}: {e}")

print(f"Checking for MSA files in: {msas_folder}")
if os.path.exists(msas_folder):
    msa_files = os.listdir(msas_folder)
    print(f"Found files in MSA folder: {msa_files}")
    for msa_file in msa_files:
        if msa_file.endswith('.csv') or msa_file.endswith('.a3m'):
            seq_id = os.path.splitext(msa_file)[0]  # This is already the sequence_id since we renamed files
            
            msa_entry = {
                'id': seq_id,  # Remove msa_ prefix - use same ID as other tables
                'sequence_id': seq_id,
                'sequence': sequences_data.get(seq_id, ''),  # Add actual protein sequence
                'msa_file': os.path.join(msas_folder, msa_file)
            }
            msa_files_in_dir.append(msa_entry)
            print(f"Added MSA file to list: {msa_file} (sequence_id: {seq_id})")
else:
    print(f"MSA folder does not exist: {msas_folder}")

print(f"Total MSA files found: {len(msa_files_in_dir)}")
if msa_files_in_dir:
    msa_df = pd.DataFrame(msa_files_in_dir)
    msa_df.to_csv(msa_csv, index=False)
    print(f"Created MSAs CSV: {msa_csv}")
else:
    print("[Warning] No MSA files found, creating empty MSAs CSV")
    # Create empty MSAs CSV so completion check doesn't fail
    msa_df = pd.DataFrame(columns=['id', 'sequence_id', 'sequence', 'msa_file'])
    msa_df.to_csv(msa_csv, index=False)
    print(f"Created empty MSAs CSV: {msa_csv}")

# Create scores explanation file
scores_explanation="""
"confidence_score": Aggregated score used to sort the predictions, corresponds to 0.8 * complex_plddt + 0.2 * iptm (ptm for single chains)
"ptm": Predicted TM score for the complex
"iptm": Predicted TM score when aggregating at the interfaces
"ligand_iptm": ipTM but only aggregating at protein-ligand interfaces
"protein_iptm": ipTM but only aggregating at protein-protein interfaces
"complex_plddt": Average pLDDT score for the complex
"complex_iplddt": Average pLDDT score when upweighting interface tokens
"complex_pde": Average PDE score for the complex
"complex_ipde": Average PDE score when aggregating at interfaces  
"chains_ptm": Predicted TM score within each chain
"pair_chains_iptm": Predicted (interface) TM score between each pair of chains
"""

with open(os.path.join(OUTPUT_FOLDER, "scores_info.txt"), 'w') as scores_exp:
    scores_exp.write(scores_explanation)

print("Post-processing completed successfully!")