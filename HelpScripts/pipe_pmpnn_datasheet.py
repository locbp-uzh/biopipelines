#Copyright © 2024 LOCBP @ University of Zürich
#Distributed under MIT license
import argparse
import os
import pandas as pd
import glob

parser = argparse.ArgumentParser(description='Create datasheet from ProteinMPNN FASTA output files')
parser.add_argument('seqs_folder', type=str, help='Path to ProteinMPNN seqs folder containing FASTA files')
parser.add_argument('pipeline_name', type=str, help='Pipeline name for generating sequence IDs')
parser.add_argument('input_datasheet', type=str, help='Input datasheet to inherit columns from, or "-" if none')
parser.add_argument('output_datasheet', type=str, help='Output CSV datasheet path')
parser.add_argument('-d', '--duplicates', action='store_true',
                    help='Allow duplicate sequences in output')

# Parse the arguments
args = parser.parse_args()

def process_proteinmpnn_fasta(fasta_file, pipeline_name, allow_duplicates=False):
    """Process a single ProteinMPNN FASTA file and extract sequence data."""
    sequences = []
    
    try:
        with open(fasta_file, 'r') as f:
            lines = [line.strip() for line in f.readlines()]
        
        # Parse FASTA file - ProteinMPNN format
        # Lines alternate between headers (>) and sequences
        i = 0
        seq_index = 0
        source_pdb = os.path.splitext(os.path.basename(fasta_file))[0]
        
        while i < len(lines):
            if lines[i].startswith('>'):
                if i + 1 < len(lines):
                    header = lines[i][1:]  # Remove '>'
                    sequence = lines[i + 1]
                    
                    # Parse header parameters (e.g., "T=0.1, sample=0, score=...")
                    seq_data = {
                        'id': f"{pipeline_name}_{source_pdb}_seq_{seq_index}",
                        'source_pdb': source_pdb,
                        'sequence': sequence,
                        'sample_index': seq_index
                    }
                    
                    # Parse additional parameters from header
                    if ', ' in header:
                        params = header.split(', ')
                        for param in params:
                            if '=' in param:
                                key, value = param.split('=', 1)
                                key = key.strip()
                                value = value.strip()
                                # Handle special keys
                                if key == 'sample':
                                    seq_data['sample_index'] = int(value)
                                else:
                                    seq_data[key] = value
                    
                    sequences.append(seq_data)
                    seq_index += 1
                    
                i += 2
            else:
                i += 1
                
    except Exception as e:
        print(f"Error processing {fasta_file}: {e}")
        
    return sequences

def main():
    print(f"Processing ProteinMPNN sequences from: {args.seqs_folder}")
    
    # Load input datasheet if provided
    input_df = None
    if args.input_datasheet != "-" and os.path.exists(args.input_datasheet):
        try:
            input_df = pd.read_csv(args.input_datasheet)
            print(f"Loaded input datasheet: {args.input_datasheet} ({len(input_df)} rows)")
        except Exception as e:
            print(f"Error loading input datasheet: {e}")
            input_df = None
    
    # Find all FASTA files in seqs folder
    fasta_pattern = os.path.join(args.seqs_folder, "*.fa")
    fasta_files = glob.glob(fasta_pattern)
    
    if not fasta_files:
        print("No FASTA files found in seqs folder")
        # Create empty datasheet with expected columns
        columns = ['id', 'source_pdb', 'sequence']
        if input_df is not None:
            # Include columns from input datasheet (excluding 'id' and 'pdb')
            inherit_cols = [col for col in input_df.columns if col not in ['id', 'pdb']]
            columns.extend(inherit_cols)
        df = pd.DataFrame(columns=columns)
        df.to_csv(args.output_datasheet, index=False)
        return
    
    print(f"Found {len(fasta_files)} FASTA files")
    
    all_sequences = []
    seen_sequences = set() if not args.duplicates else None
    
    for fasta_file in sorted(fasta_files):
        print(f"Processing: {os.path.basename(fasta_file)}")
        sequences = process_proteinmpnn_fasta(fasta_file, args.pipeline_name, args.duplicates)
        
        for seq_data in sequences:
            if not args.duplicates:
                if seq_data['sequence'] in seen_sequences:
                    print(f"Skipped duplicate sequence: {seq_data['id']}")
                    continue
                seen_sequences.add(seq_data['sequence'])
            
            # Inherit columns from input datasheet if available
            if input_df is not None:
                source_pdb = seq_data['source_pdb']
                
                # Try to find matching row in input datasheet by PDB name
                matching_rows = input_df[input_df['pdb'].str.contains(source_pdb, na=False)]
                
                if not matching_rows.empty:
                    # Use first matching row
                    input_row = matching_rows.iloc[0]
                    
                    # Inherit all columns except 'id' and 'pdb' 
                    for col in input_df.columns:
                        if col not in ['id', 'pdb'] and col not in seq_data:
                            seq_data[col] = input_row[col]
                else:
                    print(f"Warning: No matching input datasheet row found for {source_pdb}")
            
            all_sequences.append(seq_data)
    
    # Create DataFrame
    df = pd.DataFrame(all_sequences)
    
    if not df.empty:
        print(f"Total sequences processed: {len(df)}")
        print(f"Unique structures: {df['source_pdb'].nunique()}")
        print(f"Average sequences per structure: {len(df) / df['source_pdb'].nunique():.1f}")
        
        # Show inherited columns
        if input_df is not None:
            inherited_cols = [col for col in df.columns if col in input_df.columns and col not in ['id', 'source_pdb', 'sequence']]
            if inherited_cols:
                print(f"Inherited columns: {', '.join(inherited_cols)}")
    else:
        print("No sequences found")
        # Create empty datasheet with expected columns
        columns = ['id', 'source_pdb', 'sequence']
        if input_df is not None:
            inherit_cols = [col for col in input_df.columns if col not in ['id', 'pdb']]
            columns.extend(inherit_cols)
        df = pd.DataFrame(columns=columns)
    
    # Save datasheet
    df.to_csv(args.output_datasheet, index=False)
    print(f"Datasheet saved to: {args.output_datasheet}")
    
    # Print sample of first few entries
    if not df.empty:
        print("\nFirst few entries:")
        display_cols = ['id', 'source_pdb'] + ([col for col in ['fixed', 'designed'] if col in df.columns])
        if display_cols:
            print(df[display_cols].head())
        else:
            print(df[['id', 'source_pdb']].head())

if __name__ == "__main__":
    main()