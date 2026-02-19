# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.


import argparse
import os
import pandas as pd

parser = argparse.ArgumentParser(description='Create MSAs CSV table from AlphaFold MSA files')
parser.add_argument('msas_folder', type=str, help='Path to MSAs folder containing .a3m files')
parser.add_argument('queries_csv', type=str, help='Path to queries CSV with sequences')
parser.add_argument('msa_csv', type=str, help='Output CSV file path for MSA table')

# Parse the arguments
args = parser.parse_args()

def create_alphafold_msas_table(msas_folder, queries_csv, msa_csv):
    """
    Create MSAs CSV table from AlphaFold MSA files.

    Args:
        msas_folder: Path to MSAs subfolder with .a3m files
        queries_csv: Path to queries CSV containing sequences
        msa_csv: Output CSV file path
    """
    # Load sequences data to get actual protein sequences
    sequences_data = {}
    if os.path.exists(queries_csv):
        try:
            seq_df = pd.read_csv(queries_csv)
            if 'sequence' in seq_df.columns:
                sequences_data = dict(zip(seq_df['id'], seq_df['sequence']))
            print(f"Loaded {len(sequences_data)} sequences for MSA table")
        except Exception as e:
            print(f"Warning: Could not load sequences from {queries_csv}: {e}")

    msa_entries = []

    # Check if MSAs folder exists
    if not os.path.exists(msas_folder):
        print(f"Warning: MSAs folder does not exist: {msas_folder}")
        # Create empty MSAs CSV
        msa_df = pd.DataFrame(columns=['id', 'sequence_id', 'sequence', 'msa_file'])
        msa_df.to_csv(msa_csv, index=False)
        print(f"Created empty MSAs CSV: {msa_csv}")
        return

    # Process all .a3m files in the MSAs folder
    msa_files = [f for f in os.listdir(msas_folder) if f.endswith('.a3m')]
    print(f"Found {len(msa_files)} MSA files in {msas_folder}")

    for msa_file in msa_files:
        # Extract sequence ID from filename (remove .a3m extension)
        seq_id = os.path.splitext(msa_file)[0]

        msa_entry = {
            'id': seq_id,
            'sequence_id': seq_id,
            'sequence': sequences_data.get(seq_id, ''),  # Add actual protein sequence
            'msa_file': os.path.join(msas_folder, msa_file)
        }
        msa_entries.append(msa_entry)
        print(f"Added MSA file: {msa_file} (sequence_id: {seq_id})")

    # Save to CSV
    if msa_entries:
        msa_df = pd.DataFrame(msa_entries)
        msa_df.to_csv(msa_csv, index=False)
        print(f"Created MSAs CSV with {len(msa_entries)} entries: {msa_csv}")
    else:
        print("Warning: No MSA files found, creating empty MSAs CSV")
        msa_df = pd.DataFrame(columns=['id', 'sequence_id', 'sequence', 'msa_file'])
        msa_df.to_csv(msa_csv, index=False)
        print(f"Created empty MSAs CSV: {msa_csv}")

# Run the table creation
if __name__ == "__main__":
    create_alphafold_msas_table(args.msas_folder, args.queries_csv, args.msa_csv)
