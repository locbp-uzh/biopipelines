#!/usr/bin/env python3
"""
Process sequences through MMseqs2 server and generate MSA CSV files.

This script reads a sequences CSV file, submits each sequence to the MMseqs2 server,
and generates an MSA CSV file compatible with Boltz2 format.
"""

import argparse
import pandas as pd
import os
import sys
import subprocess
import tempfile
import time

def log(message):
    """Log with timestamp."""
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] {message}")

def submit_sequence_to_server(sequence, sequence_id, client_script, output_format="csv"):
    """Submit a single sequence to MMseqs2 server and return the result file path."""
    # Create temporary output file
    with tempfile.NamedTemporaryFile(mode='w', suffix=f'.{output_format}', delete=False) as temp_file:
        temp_output = temp_file.name

    try:
        # Submit to MMseqs2 server
        cmd = [
            client_script,
            "--sequence", sequence,
            "--type", output_format,
            "--output", temp_output
        ]

        log(f"Submitting sequence {sequence_id} to MMseqs2 server")
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)

        if result.returncode != 0:
            log(f"ERROR: MMseqs2 server failed for sequence {sequence_id}")
            log(f"STDOUT: {result.stdout}")
            log(f"STDERR: {result.stderr}")
            return None

        if not os.path.exists(temp_output):
            log(f"ERROR: Output file not created for sequence {sequence_id}")
            return None

        return temp_output

    except subprocess.TimeoutExpired:
        log(f"ERROR: Timeout waiting for sequence {sequence_id}")
        return None
    except Exception as e:
        log(f"ERROR: Exception processing sequence {sequence_id}: {str(e)}")
        return None

def convert_a3m_to_csv_format(a3m_file, sequence_id, output_csv_file):
    """Convert A3M file to CSV format and return summary row."""
    try:
        with open(a3m_file, 'r') as f:
            sequences = []
            current_seq = ""

            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_seq:
                        sequences.append(current_seq)
                        current_seq = ""
                else:
                    current_seq += line

            # Add last sequence
            if current_seq:
                sequences.append(current_seq)

        # Create CSV file with all MSA sequences
        msa_data = [{'sequence': seq} for seq in sequences]
        msa_df = pd.DataFrame(msa_data)
        msa_df.to_csv(output_csv_file, index=False)
        log(f"Converted A3M to CSV: {output_csv_file}")

        # Get query sequence (first sequence in A3M)
        query_sequence = sequences[0] if sequences else ""

        # Return single summary row pointing to the converted CSV file
        return [{
            'id': f"{sequence_id}_msa",
            'sequence_id': sequence_id,
            'sequence': query_sequence,
            'msa_file': output_csv_file  # Reference to converted CSV file
        }]

    except Exception as e:
        log(f"ERROR: Failed to process A3M file {a3m_file}: {str(e)}")
        return []

def process_csv_output(csv_file, sequence_id):
    """Create summary row for CSV MSA file."""
    try:
        df = pd.read_csv(csv_file)

        # Get query sequence (first row)
        query_sequence = df.iloc[0]['sequence'] if len(df) > 0 else ""

        # Return single summary row pointing to the MSA CSV file
        return [{
            'id': f"{sequence_id}_msa",
            'sequence_id': sequence_id,
            'sequence': query_sequence,
            'msa_file': csv_file  # Reference to source CSV file
        }]

    except Exception as e:
        log(f"ERROR: Failed to process CSV file {csv_file}: {str(e)}")
        return []

def main():
    parser = argparse.ArgumentParser(description='Process sequences through MMseqs2 server')
    parser.add_argument('sequences_csv', help='Input sequences CSV file (id,sequence format)')
    parser.add_argument('output_msa_csv', help='Output MSA CSV file')
    parser.add_argument('client_script', help='Path to MMseqs2 client script')
    parser.add_argument('--output_format', default='csv', choices=['csv', 'a3m'],
                       help='Output format from server (default: csv)')

    args = parser.parse_args()

    log(f"Processing sequences from {args.sequences_csv}")
    log(f"Output format: {args.output_format}")
    log(f"Client script: {args.client_script}")

    # Read input sequences
    try:
        sequences_df = pd.read_csv(args.sequences_csv)
        log(f"Loaded {len(sequences_df)} sequences")
    except Exception as e:
        log(f"ERROR: Failed to read sequences CSV: {str(e)}")
        sys.exit(1)

    # Validate input format
    if 'id' not in sequences_df.columns or 'sequence' not in sequences_df.columns:
        log("ERROR: Input CSV must have 'id' and 'sequence' columns")
        sys.exit(1)

    # Process each sequence
    all_msa_rows = []

    for _, row in sequences_df.iterrows():
        sequence_id = row['id']
        sequence = row['sequence']

        log(f"Processing sequence: {sequence_id}")

        # Submit to server
        result_file = submit_sequence_to_server(
            sequence, sequence_id, args.client_script, args.output_format
        )

        if result_file is None:
            log(f"WARNING: Skipping sequence {sequence_id} due to server error")
            continue

        # Save individual MSA file with sequence ID name
        output_dir = os.path.dirname(args.output_msa_csv)

        # Process result based on format
        if args.output_format == 'a3m':
            # For A3M: convert to CSV format
            individual_msa_file = os.path.join(output_dir, f"{sequence_id}.csv")
            msa_rows = convert_a3m_to_csv_format(result_file, sequence_id, individual_msa_file)
        else:  # csv
            # For CSV: just copy the file
            individual_msa_file = os.path.join(output_dir, f"{sequence_id}.csv")
            try:
                import shutil
                shutil.copy2(result_file, individual_msa_file)
                log(f"Saved individual MSA file: {individual_msa_file}")
            except Exception as e:
                log(f"WARNING: Failed to save individual MSA file: {e}")

            msa_rows = process_csv_output(result_file, sequence_id)

            # Update msa_file references to point to individual files
            for row in msa_rows:
                row['msa_file'] = individual_msa_file

        all_msa_rows.extend(msa_rows)

        # Clean up temporary file
        try:
            os.remove(result_file)
        except:
            pass

        log(f"Completed sequence {sequence_id}: {len(msa_rows)} MSA entries")

    # Create output MSA dataframe
    if all_msa_rows:
        msa_df = pd.DataFrame(all_msa_rows)

        # Ensure output directory exists
        os.makedirs(os.path.dirname(args.output_msa_csv), exist_ok=True)

        # Save MSA CSV
        msa_df.to_csv(args.output_msa_csv, index=False)
        log(f"Created MSA CSV with {len(msa_df)} entries: {args.output_msa_csv}")
    else:
        log("WARNING: No MSA entries generated")
        # Create empty file with correct columns
        empty_df = pd.DataFrame(columns=['id', 'sequence_id', 'sequence', 'msa_file'])
        empty_df.to_csv(args.output_msa_csv, index=False)
        log(f"Created empty MSA CSV: {args.output_msa_csv}")

if __name__ == "__main__":
    main()