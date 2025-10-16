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

def map_sequence_id_to_datasheet_id(seq_id, id_map):
    """
    Map sequence ID to datasheet ID using id_map pattern.

    Args:
        seq_id: Sequence ID (e.g., "rifampicin_1_2")
        id_map: ID mapping dictionary (e.g., {"*": "*_<N>"})

    Returns:
        Mapped datasheet ID (e.g., "rifampicin_1")
    """
    import re

    # Check if id_map uses the standard pattern
    if "*" in id_map and "*_<N>" in id_map.values():
        # Strip last _<number> from sequence ID
        match = re.match(r'^(.+)_\d+$', seq_id)
        if match:
            return match.group(1)

    # No mapping or pattern doesn't match, use as-is
    return seq_id

def sele_to_list(s):
    """
    Convert selection string to list of residue numbers (1-indexed).
    Format: "10-20+30-40" or "10 15 20"
    """
    a = []
    if not s or s == "":
        return a

    # Handle both '+' separated ranges and space-separated legacy format
    if '+' in s:
        parts = s.split('+')
    else:
        # Legacy space-separated format
        parts = s.split()

    for part in parts:
        part = part.strip()
        if not part:
            continue

        if '-' in part and not part.startswith('-'):
            # Range like "10-15"
            range_parts = part.split('-')
            if len(range_parts) == 2:
                min_val, max_val = range_parts
                for ri in range(int(min_val), int(max_val) + 1):
                    a.append(ri)
            else:
                log(f"Warning: Malformed range '{part}', skipping")
        else:
            # Single residue
            try:
                a.append(int(part))
            except ValueError:
                log(f"Warning: Could not parse '{part}' as integer, skipping")

    return sorted(a)

def apply_mask_to_msa(msa_sequences, mask_positions):
    """
    Apply mask to MSA sequences by replacing positions with '-'.

    Args:
        msa_sequences: List of sequences from MSA (first is query)
        mask_positions: List of 1-indexed positions to mask

    Returns:
        List of masked sequences (query unchanged, others masked)
    """
    if not mask_positions or len(msa_sequences) == 0:
        return msa_sequences

    # First sequence is the query - keep it unchanged
    masked_sequences = [msa_sequences[0]]

    # Mask all other sequences
    for seq in msa_sequences[1:]:
        masked_seq = list(seq)
        for pos in mask_positions:
            # Convert 1-indexed to 0-indexed, check bounds
            idx = pos - 1
            if 0 <= idx < len(masked_seq):
                masked_seq[idx] = '-'
        masked_sequences.append(''.join(masked_seq))

    return masked_sequences

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

def convert_a3m_to_csv_format(a3m_file, sequence_id, output_csv_file, mask_positions=None):
    """Convert A3M file to CSV format with optional masking, and return summary row."""
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

        # Apply masking if specified
        if mask_positions:
            log(f"Applying mask to {len(sequences)} MSA sequences (positions: {mask_positions[:10]}...)")
            sequences = apply_mask_to_msa(sequences, mask_positions)

        # Create CSV file with all MSA sequences (key=-1 for all rows)
        msa_data = [{'key': -1, 'sequence': seq} for seq in sequences]
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

def process_csv_output(csv_file, sequence_id, output_csv_file, mask_positions=None):
    """Create summary row for CSV MSA file with optional masking."""
    try:
        df = pd.read_csv(csv_file)

        # Get sequences
        sequences = df['sequence'].tolist() if 'sequence' in df.columns else []

        # Apply masking if specified
        if mask_positions and sequences:
            log(f"Applying mask to {len(sequences)} MSA sequences (positions: {mask_positions[:10]}...)")
            sequences = apply_mask_to_msa(sequences, mask_positions)

            # Save masked sequences to output file (key=-1 for all rows)
            msa_data = [{'key': -1, 'sequence': seq} for seq in sequences]
            msa_df = pd.DataFrame(msa_data)
            msa_df.to_csv(output_csv_file, index=False)
            log(f"Saved masked MSA to: {output_csv_file}")
        else:
            # No masking - need to ensure key column exists
            df_original = pd.read_csv(csv_file)

            # Add key column if it doesn't exist
            if 'key' not in df_original.columns:
                df_original.insert(0, 'key', -1)

            df_original.to_csv(output_csv_file, index=False)
            log(f"Saved MSA (with key column) to: {output_csv_file}")

        # Get query sequence (first row)
        query_sequence = sequences[0] if sequences else ""

        # Return single summary row pointing to the MSA CSV file
        return [{
            'id': f"{sequence_id}_msa",
            'sequence_id': sequence_id,
            'sequence': query_sequence,
            'msa_file': output_csv_file  # Reference to output CSV file
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
    parser.add_argument('--mask_datasheet', default=None,
                       help='Path to datasheet CSV with mask information per sequence')
    parser.add_argument('--mask_column', default=None,
                       help='Column name in mask_datasheet containing selection strings')
    parser.add_argument('--mask_selection', default=None,
                       help='Direct selection string to apply to all sequences (e.g., "10-20+30-40")')
    parser.add_argument('--id_map', default='{"*": "*_<N>"}',
                       help='JSON string for ID mapping pattern (default: {"*": "*_<N>"})')

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

    # Parse id_map from JSON
    import json
    try:
        id_map = json.loads(args.id_map)
        log(f"ID mapping pattern: {id_map}")
    except Exception as e:
        log(f"ERROR: Failed to parse id_map JSON: {str(e)}")
        sys.exit(1)

    # Load mask data if provided
    mask_data = {}  # Maps sequence_id -> list of positions to mask
    datasheet_mask_data = {}  # Maps datasheet_id -> list of positions (from datasheet)

    if args.mask_datasheet and args.mask_column:
        # Per-sequence masking from datasheet
        try:
            mask_df = pd.read_csv(args.mask_datasheet)
            log(f"Loaded mask datasheet from {args.mask_datasheet}")

            if 'id' not in mask_df.columns:
                log("ERROR: Mask datasheet must have 'id' column")
                sys.exit(1)

            if args.mask_column not in mask_df.columns:
                log(f"ERROR: Mask column '{args.mask_column}' not found in datasheet")
                sys.exit(1)

            # Parse mask selections for each datasheet entry
            for _, row in mask_df.iterrows():
                datasheet_id = row['id']
                mask_selection = row[args.mask_column]
                if pd.notna(mask_selection) and str(mask_selection).strip():
                    mask_positions = sele_to_list(str(mask_selection))
                    if mask_positions:
                        datasheet_mask_data[datasheet_id] = mask_positions
                        log(f"Mask for datasheet ID '{datasheet_id}': {len(mask_positions)} positions")

            # Now map sequence IDs to datasheet IDs and populate mask_data
            for seq_id in sequences_df['id']:
                datasheet_id = map_sequence_id_to_datasheet_id(seq_id, id_map)
                if datasheet_id in datasheet_mask_data:
                    mask_data[seq_id] = datasheet_mask_data[datasheet_id]
                    if datasheet_id != seq_id:
                        log(f"Mapped sequence ID '{seq_id}' -> datasheet ID '{datasheet_id}'")

        except Exception as e:
            log(f"ERROR: Failed to load mask datasheet: {str(e)}")
            sys.exit(1)

    elif args.mask_selection:
        # Same mask for all sequences
        mask_positions = sele_to_list(args.mask_selection)
        if mask_positions:
            log(f"Using global mask: {len(mask_positions)} positions")
            # Will apply to all sequences
            for seq_id in sequences_df['id']:
                mask_data[seq_id] = mask_positions
        else:
            log("WARNING: Empty mask selection provided")

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

        # Get mask positions for this sequence
        mask_positions = mask_data.get(sequence_id, None)

        # Save individual MSA file with sequence ID name
        output_dir = os.path.dirname(args.output_msa_csv)
        individual_msa_file = os.path.join(output_dir, f"{sequence_id}.csv")

        # Process result based on format
        if args.output_format == 'a3m':
            # For A3M: convert to CSV format with optional masking
            msa_rows = convert_a3m_to_csv_format(result_file, sequence_id, individual_msa_file, mask_positions)
        else:  # csv
            # For CSV: process with optional masking
            msa_rows = process_csv_output(result_file, sequence_id, individual_msa_file, mask_positions)

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