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
import random

# Import unified ID mapping utilities
from id_map_utils import map_table_ids_to_ids

def log(message):
    """Log with timestamp."""
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] {message}", flush=True)

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

def handle_server_error(server_dir):
    """
    Handle server error by checking if server is running and resubmitting if needed.

    Logic:
    1. Check if server is running (via check_and_resubmit_server)
    2. If random < 0.2, resubmit anyway (this handles error due to server overload by many jobs)
    """

    # First, check and resubmit if server is not running
    check_and_resubmit_server(server_dir)

    # 15% probability to resubmit anyway
    if random.random() < 0.2:
        log("Server overload: Resubmitting anyway")

        MMSEQS_SERVER_DIR = server_dir
        GPU_SUBMITTING = f"{MMSEQS_SERVER_DIR}/GPU_SUBMITTING"

        os.makedirs(MMSEQS_SERVER_DIR, exist_ok=True)
        with open(GPU_SUBMITTING, 'w') as f:
            f.write(time.strftime("%H:%M:%S"))

        script_dir = os.path.dirname(os.path.abspath(__file__))
        repo_root = os.path.dirname(script_dir)

        try:
            result = subprocess.run(
                ["./submit", "ExamplePipelines/mmseqs2_server.py"],
                cwd=repo_root,
                capture_output=True,
                text=True
            )

            import re
            match = re.search(r'Submitted batch job (\d+)', result.stdout)
            if match:
                log(f"MMseqs2 server resubmitted with job ID: {match.group(1)}")
            else:
                try:
                    os.remove(GPU_SUBMITTING)
                except:
                    pass
        except Exception as e:
            log(f"Error resubmitting server: {e}")
            try:
                os.remove(GPU_SUBMITTING)
            except:
                pass

def check_and_resubmit_server(server_dir):
    """Check if MMseqs2 server is running, and resubmit if not."""
    log("Entering check_and_resubmit_server()")
    # Use server directory passed from pipeline
    MMSEQS_SERVER_DIR = server_dir
    GPU_TIMESTAMP = f"{MMSEQS_SERVER_DIR}/GPU_SERVER"
    CPU_TIMESTAMP = f"{MMSEQS_SERVER_DIR}/CPU_SERVER"
    GPU_SUBMITTING = f"{MMSEQS_SERVER_DIR}/GPU_SUBMITTING"
    CPU_SUBMITTING = f"{MMSEQS_SERVER_DIR}/CPU_SUBMITTING"
    MAX_AGE_HOURS = 3
    log(f"Server dir: {MMSEQS_SERVER_DIR}, MAX_AGE: {MAX_AGE_HOURS}h") 

    def server_is_valid(timestamp_file):
        """Check if server timestamp file is valid and not too old."""
        if not os.path.exists(timestamp_file):
            return False

        try:
            with open(timestamp_file, 'r') as f:
                server_time_str = f.read().strip()

            # Parse server timestamp (format: HH:MM:SS)
            from datetime import datetime
            current_time = datetime.now()
            server_time = datetime.strptime(server_time_str, "%H:%M:%S").replace(
                year=current_time.year,
                month=current_time.month,
                day=current_time.day
            )

            # Handle midnight wraparound
            time_diff = (current_time - server_time).total_seconds()
            if time_diff < 0:
                # Server started yesterday
                time_diff += 86400

            max_age_seconds = MAX_AGE_HOURS * 3600
            if time_diff < max_age_seconds:
                hours_old = int(time_diff / 3600)
                minutes_old = int((time_diff % 3600) / 60)
                log(f"Server started at {server_time_str} ({hours_old}h {minutes_old}m ago)")
                return True
            else:
                log(f"Server timestamp too old (started at {server_time_str})")
                return False
        except Exception as e:
            log(f"Error checking timestamp {timestamp_file}: {e}")
            return False

    def submission_in_progress(submit_file):
        """Check if server submission is in progress."""
        if not os.path.exists(submit_file):
            return False

        try:
            submit_time = os.path.getmtime(submit_file)
            current_time = time.time()
            age = int(current_time - submit_time)
            log(f"Server submission in progress (submitted {age}s ago)")
            return True
        except Exception as e:
            log(f"Error checking submission file {submit_file}: {e}")
            return False

    # Check if GPU or CPU server is running and valid
    log("Checking GPU server...")
    server_running = False
    if server_is_valid(GPU_TIMESTAMP):
        log("MMseqs2 GPU server is running and valid")
        server_running = True
    elif server_is_valid(CPU_TIMESTAMP):
        log("MMseqs2 CPU server is running and valid")
        server_running = True
    elif submission_in_progress(GPU_SUBMITTING):
        log("GPU server submission in progress, waiting...")
        server_running = True  # Consider it running if submission in progress
    elif submission_in_progress(CPU_SUBMITTING):
        log("CPU server submission in progress, waiting...")
        server_running = True  # Consider it running if submission in progress

    log(f"Server running status: {server_running}")

    if not server_running:
        log("No valid MMseqs2 server found, resubmitting server...")

        # Create submission timestamp to prevent other processes from also submitting
        os.makedirs(MMSEQS_SERVER_DIR, exist_ok=True)
        with open(GPU_SUBMITTING, 'w') as f:
            f.write(time.strftime("%H:%M:%S"))
        log(f"Created submission timestamp at {GPU_SUBMITTING}")

        # Find notebooks directory (go up from HelpScripts to repo root)
        script_dir = os.path.dirname(os.path.abspath(__file__))
        repo_root = os.path.dirname(script_dir)
        notebooks_dir = repo_root
        log(f"Repo root: {repo_root}")

        try:
            # Submit server job
            log(f"Submitting server from: {notebooks_dir}")
            result = subprocess.run(
                ["./submit", "ExamplePipelines/mmseqs2_server.py"],
                cwd=notebooks_dir,
                capture_output=True,
                text=True
            )
            log(f"Submit command completed with return code: {result.returncode}")

            # Extract job ID from output
            import re
            match = re.search(r'Submitted batch job (\d+)', result.stdout)
            if match:
                server_job_id = match.group(1)
                log(f"MMseqs2 server submitted with job ID: {server_job_id}")
                log("Proceeding to submit queries - they will queue until server is ready")
            else:
                log("Warning: Failed to submit server job, but will attempt to submit queries anyway")
                log(f"Submit output: {result.stdout}")
                log(f"Submit error: {result.stderr}")
                # Clean up submission timestamp on failure
                try:
                    os.remove(GPU_SUBMITTING)
                except:
                    pass
        except Exception as e:
            log(f"Error submitting server: {e}")
            # Clean up submission timestamp on failure
            try:
                os.remove(GPU_SUBMITTING)
            except:
                pass

    log("Exiting check_and_resubmit_server()")

def submit_sequence_to_server(sequence, sequence_id, client_script, output_format="csv", output_path=None, timeout=3600, server_dir=None):
    """Submit a single sequence to MMseqs2 server with specified output path."""
    try:
        # Build command
        cmd = [
            client_script,
            "--sequence", sequence,
            "--type", output_format
        ]

        # Add output path if specified
        if output_path:
            cmd.extend(["--output", output_path])
            log(f"Submitting sequence {sequence_id} to MMseqs2 server (output: {output_path})")
        else:
            log(f"Submitting sequence {sequence_id} to MMseqs2 server")

        result = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)

        if result.returncode != 0:
            log(f"ERROR: MMseqs2 server failed for sequence {sequence_id}")
            log(f"STDOUT: {result.stdout}")
            log(f"STDERR: {result.stderr}")
            return None

        if output_path and not os.path.exists(output_path):
            log(f"ERROR: Output file not created at {output_path} for sequence {sequence_id}")
            return None

        return output_path if output_path else None

    except subprocess.TimeoutExpired:
        log(f"ERROR: Timeout waiting for sequence {sequence_id}")

        # Handle timeout error by checking server and potentially resubmitting
        if server_dir:
            handle_server_error(server_dir)

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

            # Add key column if it doesn't exist - we already know it exists but like to waste CPU and IO
            if 'key' not in df_original.columns:
                df_original.insert(0, 'key', -1)

            df_original.to_csv(output_csv_file, index=False)
            log(f"Saved MSA to: {output_csv_file}")

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
    parser.add_argument('--server_dir', required=True,
                       help='Path to MMseqs2 server directory for timestamp files')
    parser.add_argument('--mask_table', default=None,
                       help='Path to table CSV with mask information per sequence')
    parser.add_argument('--mask_column', default=None,
                       help='Column name in mask_table containing selection strings')
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
    table_mask_data = {}  # Maps table_id -> list of positions (from table)

    if args.mask_table and args.mask_column:
        # Per-sequence masking from table
        try:
            mask_df = pd.read_csv(args.mask_table)
            log(f"Loaded mask table from {args.mask_table}")

            if 'id' not in mask_df.columns:
                log("ERROR: Mask table must have 'id' column")
                sys.exit(1)

            if args.mask_column not in mask_df.columns:
                log(f"ERROR: Mask column '{args.mask_column}' not found in table")
                sys.exit(1)

            # Parse mask selections for each table entry
            for _, row in mask_df.iterrows():
                table_id = row['id']
                mask_selection = row[args.mask_column]
                if pd.notna(mask_selection) and str(mask_selection).strip():
                    mask_positions = sele_to_list(str(mask_selection))
                    if mask_positions:
                        table_mask_data[table_id] = mask_positions
                        log(f"Mask for table ID '{table_id}': {len(mask_positions)} positions")

            # Now map sequence IDs to table IDs and populate mask_data
            for seq_id in sequences_df['id']:
                candidate_ids = map_table_ids_to_ids(seq_id, id_map)
                # Try all candidate IDs in priority order
                found = False
                for candidate_id in candidate_ids:
                    if candidate_id in table_mask_data:
                        mask_data[seq_id] = table_mask_data[candidate_id]
                        if candidate_id != seq_id:
                            log(f"Mapped sequence ID '{seq_id}' -> table ID '{candidate_id}'")
                        found = True
                        break

                if not found:
                    log(f"Warning: No mask data for sequence ID '{seq_id}'. Tried: {', '.join(candidate_ids)}")

        except Exception as e:
            log(f"ERROR: Failed to load mask table: {str(e)}")
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

    # Check server status before starting
    log("Checking server status before starting submissions...")
    check_and_resubmit_server(args.server_dir)

    # Process each sequence
    all_msa_rows = []
    submission_counter = 0

    for _, row in sequences_df.iterrows():
        sequence_id = row['id']
        sequence = row['sequence']

        # Check server status every 25 submissions
        if submission_counter > 0 and submission_counter % 25 == 0:
            log(f"Checking server status after {submission_counter} submissions...")
            check_and_resubmit_server(args.server_dir)

        # Save individual MSA file with sequence ID name
        output_dir = os.path.dirname(args.output_msa_csv)
        individual_msa_file = os.path.join(output_dir, f"{sequence_id}.csv")

        # Check if MSA already exists (for retry scenarios)
        if os.path.exists(individual_msa_file):
            log(f"MSA already exists for {sequence_id}, skipping submission")

            # Read existing MSA file to get query sequence
            try:
                existing_df = pd.read_csv(individual_msa_file)
                if 'sequence' in existing_df.columns and len(existing_df) > 0:
                    query_sequence = existing_df['sequence'].iloc[0]
                    msa_rows = [{
                        'id': f"{sequence_id}_msa",
                        'sequence_id': sequence_id,
                        'sequence': query_sequence,
                        'msa_file': individual_msa_file
                    }]
                    all_msa_rows.extend(msa_rows)
                    log(f"Reused existing MSA for {sequence_id}")
                    continue
                else:
                    log(f"WARNING: Existing MSA file is malformed, will reprocess")
                    # Fall through to reprocess
            except Exception as e:
                log(f"WARNING: Could not read existing MSA file: {e}, will reprocess")
                # Fall through to reprocess

        log(f"Processing sequence: {sequence_id}")

        # Get mask positions for this sequence
        mask_positions = mask_data.get(sequence_id, None)

        # If masking is needed, use temporary location; otherwise use final location
        if mask_positions:
            # Server writes to temp, we apply masking and save to final location
            with tempfile.NamedTemporaryFile(mode='w', suffix=f'.{args.output_format}', delete=False) as temp_file:
                temp_output = temp_file.name

            result_file = submit_sequence_to_server(
                sequence, sequence_id, args.client_script, args.output_format, temp_output,
                server_dir=args.server_dir
            )

            if result_file is None:
                log(f"WARNING: Skipping sequence {sequence_id} due to server error")
                continue

            # Process result with masking and save to final location
            if args.output_format == 'a3m':
                msa_rows = convert_a3m_to_csv_format(result_file, sequence_id, individual_msa_file, mask_positions)
            else:  # csv
                msa_rows = process_csv_output(result_file, sequence_id, individual_msa_file, mask_positions)

            # Clean up temporary file
            try:
                os.remove(result_file)
            except:
                pass
        else:
            # No masking - server writes directly to final location
            result_file = submit_sequence_to_server(
                sequence, sequence_id, args.client_script, args.output_format, individual_msa_file,
                server_dir=args.server_dir
            )

            if result_file is None:
                log(f"WARNING: Skipping sequence {sequence_id} due to server error")
                continue

            # Process result to create summary rows (no masking needed)
            if args.output_format == 'a3m':
                msa_rows = convert_a3m_to_csv_format(individual_msa_file, sequence_id, individual_msa_file, None)
            else:  # csv
                msa_rows = process_csv_output(individual_msa_file, sequence_id, individual_msa_file, None)

        all_msa_rows.extend(msa_rows)
        submission_counter += 1

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