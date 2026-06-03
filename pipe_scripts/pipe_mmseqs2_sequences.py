#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

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

# Add repo root to path so biopipelines package is importable
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from biopipelines.biopipelines_io import lookup_table_value

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
    if not s or s == "" or pd.isna(s):
        return a

    # Convert to string to handle numeric types
    s = str(s)
    if s == "nan":
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

# How long a submission lock is honoured before it's considered stale (the
# submitted server job failed to come up / never wrote its timestamp). After
# this, another client may steal the lock and resubmit. Must comfortably exceed
# SLURM queue wait + the server's DB load. The CPU server holds this lock until
# its ~746GB index is fully mlocked (several minutes of load on top of a
# big-RAM-partition queue that can run to hours), so 30 min was too short — a
# legitimately-still-loading server would get a duplicate submitted. 3h matches
# the server-validity window (MAX_AGE_HOURS) used elsewhere.
SUBMIT_LOCK_TTL_SECONDS = 3 * 60 * 60


def handle_server_error(server_dir):
    """On a query timeout, re-check the server and resubmit only if genuinely
    absent.

    There is intentionally NO random "resubmit anyway for overload" here. That
    heuristic was the main cause of server pileup: with many concurrent clients
    timing out during the server's startup window, each rolled an independent
    chance to submit yet another server. A single server already serves all
    queued queries, so the only legitimate reason to submit is that none exists
    — which check_and_resubmit_server() decides atomically.
    """
    check_and_resubmit_server(server_dir)

def _get_server_mode():
    """Which MMseqs2 server to auto-start: 'cpu' (default) or 'gpu'.

    Read from the active config at tool_overrides.mmseqs2server.mode. Sites that
    want the GPU server auto-started set it to 'gpu'; otherwise the CPU server
    (RAM-resident DBs, no GPU) is launched. Any unset/unknown value falls back to
    'cpu'.
    """
    try:
        from biopipelines.config_manager import ConfigManager
        overrides = ConfigManager()._config.get('tool_overrides', {}) or {}
        mode = (overrides.get('mmseqs2server', {}) or {}).get('mode', 'cpu')
    except Exception as e:
        log(f"Could not read tool_overrides.mmseqs2server.mode ({e}); defaulting to cpu")
        return "cpu"
    mode = str(mode).strip().lower()
    if mode not in ("cpu", "gpu"):
        log(f"Unknown mmseqs2server.mode={mode!r}; defaulting to cpu")
        return "cpu"
    return mode

def _get_num_servers():
    """How many server jobs to launch as a pool when auto-starting (default 1).

    Read from tool_overrides.mmseqs2server.num_servers. The single client that
    wins the submit lock submits this many identical server jobs; they all drain
    the one shared queue (per-query mv-claim keeps that safe). >1 trades cluster
    allocation for throughput under many concurrent queries.
    """
    try:
        from biopipelines.config_manager import ConfigManager
        overrides = ConfigManager()._config.get('tool_overrides', {}) or {}
        n = (overrides.get('mmseqs2server', {}) or {}).get('num_servers', 1)
        n = int(n)
    except Exception as e:
        log(f"Could not read tool_overrides.mmseqs2server.num_servers ({e}); defaulting to 1")
        return 1
    return n if n >= 1 else 1

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
        """True iff a fresh (non-stale) submission lock exists.

        A lock older than SUBMIT_LOCK_TTL_SECONDS means the server it was
        submitted for never came up; it's treated as absent so a client may
        steal it and resubmit.
        """
        if not os.path.exists(submit_file):
            return False

        try:
            age = int(time.time() - os.path.getmtime(submit_file))
            if age >= SUBMIT_LOCK_TTL_SECONDS:
                log(f"Submission lock is stale ({age}s old >= {SUBMIT_LOCK_TTL_SECONDS}s), ignoring")
                return False
            log(f"Server submission in progress (submitted {age}s ago)")
            return True
        except Exception as e:
            log(f"Error checking submission file {submit_file}: {e}")
            return False

    def acquire_submit_lock(lock_file):
        """Atomically claim the right to submit a server, NFS-safe.

        Returns True iff THIS process won the claim, so only one of many
        concurrent clients (typically on different nodes, sharing the lock over
        NFS) submits a server. If the lock is stale (older than the TTL — the
        submitted server never came up) it is stolen once and re-claimed.

        Uses ``os.mkdir`` as the atomic primitive rather than ``open(O_EXCL)``:
        O_EXCL file creation is NOT atomic over NFS, so it let two clients on
        different nodes both "acquire" and each submit a server. Directory
        creation IS atomic over NFS (mkdir either creates or fails with EEXIST),
        which is the classic NFS-safe lock. On success a ``GPU_SUBMITTING``
        marker file is also written for the visibility/staleness check and for
        the server to clear on startup.
        """
        lock_dir = lock_file + ".lockdir"
        os.makedirs(os.path.dirname(lock_file), exist_ok=True)
        for attempt in (1, 2):
            try:
                os.mkdir(lock_dir)  # atomic over NFS
                # Won the lock — also drop the marker file used elsewhere.
                with open(lock_file, "w") as f:
                    f.write(time.strftime("%H:%M:%S"))
                return True
            except FileExistsError:
                # Someone holds it. Steal once if stale, else yield.
                try:
                    age = int(time.time() - os.path.getmtime(lock_dir))
                except OSError:
                    age = 0  # vanished between mkdir and stat — retry
                if attempt == 1 and age >= SUBMIT_LOCK_TTL_SECONDS:
                    log(f"Stealing stale submission lock ({age}s old)")
                    try:
                        os.rmdir(lock_dir)
                    except OSError:
                        pass
                    continue
                return False
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
        # Which server to auto-start (cpu or gpu). Read from the active config's
        # tool_overrides.mmseqs2server.mode; default "cpu". A running server of
        # *either* mode was already accepted above (the checks cover both
        # GPU_SERVER and CPU_SERVER) — this only governs what we launch when none
        # is up. CPU mode serves from the RAM-resident DBs and needs no GPU.
        server_mode = _get_server_mode()
        if server_mode == "gpu":
            server_pipeline = "example_pipelines/_mmseqs2_server.py"
            submit_lock = GPU_SUBMITTING
        else:
            server_pipeline = "example_pipelines/_mmseqs2_cpu_server.py"
            submit_lock = CPU_SUBMITTING
        log(f"Server mode (tool_overrides.mmseqs2server.mode): {server_mode}")

        # Atomically claim the right to submit. If another concurrent client
        # already holds the lock, yield — its server will serve our queries too.
        if not acquire_submit_lock(submit_lock):
            log("Another client is already submitting a server; not submitting a duplicate")
            log("Exiting check_and_resubmit_server()")
            return

        log("No valid MMseqs2 server found and lock acquired, submitting server...")

        # Find biopipelines directory (go up from pipe_scripts to repo root)
        script_dir = os.path.dirname(os.path.abspath(__file__))
        repo_root = os.path.dirname(script_dir)
        biopipelines_dir = repo_root
        log(f"Repo root: {repo_root}")

        try:
            # Submit the server pool. The single client that won the lock submits
            # num_servers identical server jobs; they all drain the one shared
            # queue (per-query mv-claim keeps that safe). Default 1.
            num_servers = _get_num_servers()
            log(f"Submitting {num_servers} {server_mode} server(s) ({server_pipeline}) from: {biopipelines_dir}")
            import re
            submitted_ids = []
            result = None
            for n in range(num_servers):
                result = subprocess.run(
                    ["./submit", server_pipeline],
                    cwd=biopipelines_dir,
                    capture_output=True,
                    text=True
                )
                m = re.search(r'Submitted batch job (\d+)', result.stdout)
                if m:
                    submitted_ids.append(m.group(1))
                    log(f"  server {n+1}/{num_servers} submitted: job {m.group(1)}")
                else:
                    log(f"  WARNING: server {n+1}/{num_servers} submit produced no job id")
                    log(f"  submit stdout: {result.stdout}")
                    log(f"  submit stderr: {result.stderr}")

            # As long as at least one server was submitted, proceed (queries queue
            # until a server is ready). Only treat it as a failure if none went in.
            if submitted_ids:
                log(f"MMseqs2 server pool submitted: {submitted_ids}")
                log("Proceeding to submit queries - they will queue until a server is ready")
            else:
                log("Warning: Failed to submit server job, but will attempt to submit queries anyway")
                log(f"Submit output: {result.stdout}")
                log(f"Submit error: {result.stderr}")
                # Release the submission lock on failure so another client can retry.
                try:
                    os.remove(submit_lock)
                except OSError:
                    pass
                try:
                    os.rmdir(submit_lock + ".lockdir")
                except OSError:
                    pass
        except Exception as e:
            log(f"Error submitting server: {e}")
            # Clean up submission timestamp on failure
            try:
                os.remove(submit_lock)
            except:
                pass

    log("Exiting check_and_resubmit_server()")

def submit_sequence_to_server(sequence, sequence_id, client_script, output_format="csv", output_path=None, timeout=1800, server_dir=None):
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
            'sequences.id': sequence_id,
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
            'sequences.id': sequence_id,
            'sequence': query_sequence,
            'msa_file': output_csv_file  # Reference to output CSV file
        }]

    except Exception as e:
        log(f"ERROR: Failed to process CSV file {csv_file}: {str(e)}")
        return []

def submit_batch_job(server_dir, seqs, output_format, output_dir):
    """Submit ALL sequences as ONE job carrying a multi-FASTA; return the job_id.

    The server runs a single colabfold_search over the whole set (many-against-
    many), writing one "<alias>.<ext>" per sequence into output_dir. `seqs` is a
    list of (seq_id, sequence). Non-blocking — just drops the .fasta + .job into
    the queue; poll the result with wait_for_status().

    Queries are aliased to sequential integer headers (0,1,2,...) rather than the
    raw sequence ids: colabfold_search names each result safe_filename(header)
    (every char outside [0-9a-zA-Z_.-] -> '_'), so raw ids containing e.g. '+'
    would land under a different name, and two distinct ids could even sanitize to
    the SAME file and clobber each other. Integer aliases are collision-free and
    pass through unchanged. The caller maps "<alias>.<ext>" back to the real id.
    """
    queue_dir = os.path.join(server_dir, "job_queue")
    os.makedirs(queue_dir, exist_ok=True)
    os.makedirs(output_dir, exist_ok=True)
    job_id = f"msa_{int(time.time())}_{os.getpid()}"

    # Multi-FASTA with sequential-integer headers; the server names each result
    # "<alias>.<ext>". Lives in the queue dir so the server can read it.
    fasta_path = os.path.join(queue_dir, f"{job_id}.fasta")
    with open(fasta_path, "w") as f:
        for alias, (sid, seq) in enumerate(seqs):
            f.write(f">{alias}\n{seq}\n")

    # Write the .job last via temp+rename so the server never reads a half-written
    # job file.
    job_path = os.path.join(queue_dir, f"{job_id}.job")
    tmp_job = job_path + ".tmp"
    with open(tmp_job, "w") as f:
        f.write(f"job_id={job_id}\n")
        f.write(f"fasta={fasta_path}\n")
        f.write(f"output_format={output_format}\n")
        f.write(f"output_dir={output_dir}\n")
        f.write(f"submitted={time.strftime('%Y-%m-%d %H:%M:%S')}\n")
    os.rename(tmp_job, job_path)
    log(f"Submitted batch job {job_id}: {len(seqs)} sequences -> {output_dir}")
    return job_id

def wait_for_status(server_dir, job_id, timeout):
    """Block until the server writes <job_id>.status; return its first line
    (SUCCESS/FAILED) or None on timeout. Removes the status file once read."""
    status_file = os.path.join(server_dir, "results", f"{job_id}.status")
    start = time.time()
    while True:
        if os.path.exists(status_file):
            try:
                with open(status_file) as f:
                    status = f.readline().strip()
            except Exception:
                status = ""
            try:
                os.remove(status_file)
            except OSError:
                pass
            return status or "FAILED"
        if time.time() - start > timeout:
            log(f"ERROR: timeout ({timeout}s) waiting for job {job_id}")
            return None
        time.sleep(5)

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
    parser.add_argument('--missing_csv', default=None,
                       help='Destination missing.csv (id,removed_by,kind,cause) for sequences with no MSA')
    parser.add_argument('--upstream_missing', default=None,
                       help='Upstream missing.csv to propagate forward and exclude before searching')

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

    # Load mask data if provided
    mask_data = {}  # Maps sequence_id -> list of positions to mask

    if args.mask_table and args.mask_column:
        # Per-sequence masking from table: lookup_table_value handles ID matching
        # via pdb/id columns, id_map suffix stripping, and {stream}.id provenance.
        try:
            mask_df = pd.read_csv(args.mask_table)
            log(f"Loaded mask table from {args.mask_table}")

            if args.mask_column not in mask_df.columns:
                log(f"ERROR: Mask column '{args.mask_column}' not found in table")
                sys.exit(1)

            for seq_id in sequences_df['id']:
                try:
                    mask_selection = lookup_table_value(mask_df, seq_id, args.mask_column)
                except Exception:
                    mask_selection = None

                if mask_selection is not None and pd.notna(mask_selection) and str(mask_selection).strip():
                    mask_positions = sele_to_list(str(mask_selection))
                    if mask_positions:
                        mask_data[seq_id] = mask_positions
                        log(f"Mask for sequence ID '{seq_id}': {len(mask_positions)} positions")
                else:
                    log(f"Warning: No mask data for sequence ID '{seq_id}'")

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

    # Missing-table rows (framework convention: id,removed_by,kind,cause). Seed
    # with any sequences the upstream tool already dropped, and exclude those ids
    # from the search set so they propagate forward instead of being re-queried.
    missing_rows = []
    upstream_missing_ids = set()
    if args.upstream_missing and os.path.exists(args.upstream_missing):
        try:
            up_df = pd.read_csv(args.upstream_missing)
            if 'id' in up_df.columns:
                upstream_missing_ids = set(up_df['id'].astype(str))
                missing_rows.extend(up_df.to_dict('records'))
                log(f"Propagating {len(upstream_missing_ids)} missing id(s) from upstream")
        except Exception as e:
            log(f"WARNING: could not read upstream missing.csv {args.upstream_missing}: {e}")

    if upstream_missing_ids:
        before = len(sequences_df)
        sequences_df = sequences_df[~sequences_df['id'].astype(str).isin(upstream_missing_ids)]
        log(f"Excluded {before - len(sequences_df)} upstream-missing sequence(s) from search")

    # Check server status before starting
    log("Checking server status before starting submissions...")
    check_and_resubmit_server(args.server_dir)

    output_dir = os.path.dirname(args.output_msa_csv)
    os.makedirs(output_dir, exist_ok=True)

    # --- Phase 1: collect work. Reuse any per-sequence MSA that already exists
    # (resume after a partial/failed run), and gather the rest into ONE batch. ---
    all_msa_rows = []
    pending = []   # (sequence_id, sequence) still needing an MSA
    for _, row in sequences_df.iterrows():
        sequence_id = row['id']
        sequence = row['sequence']
        individual_msa_file = os.path.join(output_dir, f"{sequence_id}.csv")

        if os.path.exists(individual_msa_file):
            try:
                existing_df = pd.read_csv(individual_msa_file)
                if 'sequence' in existing_df.columns and len(existing_df) > 0:
                    all_msa_rows.append({
                        'id': f"{sequence_id}_msa",
                        'sequences.id': sequence_id,
                        'sequence': existing_df['sequence'].iloc[0],
                        'msa_file': individual_msa_file,
                    })
                    log(f"Reused existing MSA for {sequence_id}")
                    continue
                log(f"WARNING: existing MSA for {sequence_id} malformed, will reprocess")
            except Exception as e:
                log(f"WARNING: could not read existing MSA for {sequence_id}: {e}, will reprocess")
        pending.append((sequence_id, sequence))

    # --- Phase 2: submit ALL pending sequences as ONE multi-FASTA job and wait.
    # mmseqs2 is many-against-many, so the server runs a single colabfold_search
    # over the whole set (the expensive index sweep is amortised across every
    # sequence) instead of one search per sequence. ---
    if pending:
        ext = args.output_format
        # Unique per-invocation staging dir so concurrent client pipelines (same
        # seq ids from different runs) don't clobber each other's results.
        staging = os.path.join(args.server_dir, "results",
                               f"client_out_{int(time.time())}_{os.getpid()}")
        os.makedirs(staging, exist_ok=True)

        log(f"Submitting batch of {len(pending)} sequence(s) as one job")
        job_id = submit_batch_job(args.server_dir, pending, ext, staging)
        # Generous timeout: one amortised search, but a large set on a busy node
        # can still take a while; scale with batch size.
        timeout = max(7200, 120 * len(pending))
        status = wait_for_status(args.server_dir, job_id, timeout)

        if status is None or status.startswith("FAILED"):
            log(f"ERROR: batch job {job_id} did not succeed (status: {status})")
            cause = "batch timeout" if status is None else status.replace("\n", " ")
            for sequence_id, _seq in pending:
                missing_rows.append({'id': sequence_id, 'removed_by': 'MMseqs2',
                                     'kind': 'failure', 'cause': cause})
        else:
            # Collect each sequence's result, apply masking, write the tool output.
            # Results are named by the sequential alias (the FASTA header), so map
            # alias index -> real sequence id by enumerating `pending` in order.
            for alias, (sequence_id, _seq) in enumerate(pending):
                produced = os.path.join(staging, f"{alias}.{ext}")
                individual_msa_file = os.path.join(output_dir, f"{sequence_id}.csv")
                if not os.path.exists(produced):
                    log(f"WARNING: no MSA produced for {sequence_id} (alias {alias})")
                    missing_rows.append({'id': sequence_id, 'removed_by': 'MMseqs2',
                                         'kind': 'failure', 'cause': 'no MSA produced'})
                    continue
                mask_positions = mask_data.get(sequence_id, None)
                if ext == 'a3m':
                    msa_rows = convert_a3m_to_csv_format(produced, sequence_id, individual_msa_file, mask_positions)
                else:
                    msa_rows = process_csv_output(produced, sequence_id, individual_msa_file, mask_positions)
                all_msa_rows.extend(msa_rows)
                try:
                    os.remove(produced)
                except OSError:
                    pass
                log(f"Completed sequence {sequence_id}: {len(msa_rows)} MSA entries")

        # Remove the now-empty per-invocation staging dir.
        try:
            os.rmdir(staging)
        except OSError:
            pass

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
        empty_df = pd.DataFrame(columns=['id', 'sequences.id', 'sequence', 'msa_file'])
        empty_df.to_csv(args.output_msa_csv, index=False)
        log(f"Created empty MSA CSV: {args.output_msa_csv}")

    # Write the missing table (always, even when empty, so downstream Load can
    # read it). Columns follow the framework convention: id,removed_by,kind,cause.
    if args.missing_csv:
        os.makedirs(os.path.dirname(args.missing_csv), exist_ok=True)
        cols = ['id', 'removed_by', 'kind', 'cause']
        pd.DataFrame(missing_rows, columns=cols).to_csv(args.missing_csv, index=False)
        log(f"Wrote missing table with {len(missing_rows)} entr(y/ies): {args.missing_csv}")

if __name__ == "__main__":
    main()