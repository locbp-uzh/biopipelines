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

def convert_a3m_to_csv_format(a3m_file, sequence_id):
    """Convert A3M file to CSV format compatible with Boltz2."""
    msa_rows = []

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

        # Create MSA rows in Boltz2 format
        for i, seq in enumerate(sequences):
            msa_rows.append({
                'id': f"{sequence_id}_msa_{i}",
                'sequence_id': sequence_id,
                'sequence': seq,
                'msa_file': a3m_file  # Reference to source A3M file
            })

    except Exception as e:
        log(f"ERROR: Failed to process A3M file {a3m_file}: {str(e)}")

    return msa_rows

def process_csv_output(csv_file, sequence_id):
    """Process CSV output from MMseqs2 server."""
    msa_rows = []

    try:
        df = pd.read_csv(csv_file)

        for i, row in df.iterrows():
            msa_rows.append({
                'id': f"{sequence_id}_msa_{i}",
                'sequence_id': sequence_id,
                'sequence': row['sequence'],
                'msa_file': csv_file  # Reference to source CSV file
            })

    except Exception as e:
        log(f"ERROR: Failed to process CSV file {csv_file}: {str(e)}")

    return msa_rows

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
        individual_msa_file = os.path.join(output_dir, f"{sequence_id}.{args.output_format}")

        try:
            # Copy the result file to the output directory with proper name
            import shutil
            shutil.copy2(result_file, individual_msa_file)
            log(f"Saved individual MSA file: {individual_msa_file}")
        except Exception as e:
            log(f"WARNING: Failed to save individual MSA file: {e}")

        # Process result based on format for combined CSV
        if args.output_format == 'a3m':
            msa_rows = convert_a3m_to_csv_format(result_file, sequence_id)
        else:  # csv
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