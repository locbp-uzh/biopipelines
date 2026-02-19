#!/bin/bash
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# This software is freely available for use, modification, and redistribution.

set -e

# timestamped logging
timestamp() { date '+%Y-%m-%d %H:%M:%S'; }
log() { echo "[$(timestamp)] $*"; }

# Resolve a required folder: env var first, then biopipelines-config, else fail
require_folder() {
  local val="$1" name="$2" key="$3"
  if [[ -n "$val" ]]; then echo "$val"; return; fi
  val=$(biopipelines-config folder "$key" 2>/dev/null) && [[ -n "$val" ]] && { echo "$val"; return; }
  echo "ERROR: $name is not set and biopipelines-config could not resolve '$key'." >&2
  echo "Set the environment variable or configure the folder in config.yaml." >&2
  exit 1
}

# MMseqs2 client working directories - must match server script structure
USER=${USER:-$(whoami)}
MMSEQS2_SHARED_FOLDER=$(require_folder "${MMSEQS2_SHARED_FOLDER:-}" "MMSEQS2_SHARED_FOLDER" "MMseqs2Server")
JOB_QUEUE_DIR="$MMSEQS2_SHARED_FOLDER/job_queue"
RESULTS_DIR="$MMSEQS2_SHARED_FOLDER/results"

usage() {
    echo "Usage: $0 -s <sequence> [-f <fasta_file>] [-t <a3m|csv>] [-o <output_file>]"
}

# Defaults
STATUS=0
OUTPUT_FORMAT="csv"  # Changed default to csv

while [[ $# -gt 0 ]]; do
  case $1 in
    -s|--sequence)    SEQUENCE="$2"; shift 2 ;;
    -f|--fasta)       FASTA_FILE="$2"; shift 2 ;;
    -t|--type)        OUTPUT_FORMAT="$2"; shift 2 ;;
    -o|--output)      OUTPUT_PATH="$2"; shift 2 ;;
    -u|--status)      STATUS=1; shift ;;
    *) usage; exit 1 ;;
  esac
done

# -- status check --
if [[ "$STATUS" -eq 1 ]]; then
  # Check servers using timestamp files (per-user)
  MMSEQS_SERVER_DIR="$MMSEQS2_SHARED_FOLDER"
  GPU_TIMESTAMP="$MMSEQS_SERVER_DIR/GPU_SERVER"
  CPU_TIMESTAMP="$MMSEQS_SERVER_DIR/CPU_SERVER"
  MAX_AGE_HOURS=3

  check_server_timestamp() {
    local timestamp_file=$1
    local server_type=$2

    if [[ ! -f "$timestamp_file" ]]; then
      log "MMseqs2 $server_type server: NOT RUNNING (no timestamp file)"
      return 1
    fi

    # Read the timestamp from the file
    local server_time=$(cat "$timestamp_file")

    # Get current time and calculate the difference
    local current_seconds=$(date +%s)
    local server_seconds=$(date -d "today $server_time" +%s 2>/dev/null || echo "0")

    # Handle case where server started yesterday (time wrapped around midnight)
    local time_diff=$((current_seconds - server_seconds))
    if [[ $time_diff -lt 0 ]]; then
      # Server timestamp is from yesterday
      time_diff=$((time_diff + 86400))
    fi

    local max_age_seconds=$((MAX_AGE_HOURS * 3600))
    local hours_old=$((time_diff / 3600))
    local minutes_old=$(((time_diff % 3600) / 60))

    if [[ $time_diff -lt $max_age_seconds ]]; then
      log "MMseqs2 $server_type server: RUNNING (started at $server_time, ${hours_old}h ${minutes_old}m ago)"
      return 0
    else
      log "MMseqs2 $server_type server: EXPIRED (started at $server_time, ${hours_old}h ${minutes_old}m ago, max age: ${MAX_AGE_HOURS}h)"
      return 1
    fi
  }

  # Check both servers
  gpu_running=false
  cpu_running=false

  if check_server_timestamp "$GPU_TIMESTAMP" "GPU"; then
    gpu_running=true
  fi

  if check_server_timestamp "$CPU_TIMESTAMP" "CPU"; then
    cpu_running=true
  fi

  if [[ "$gpu_running" = false && "$cpu_running" = false ]]; then
    log "No valid MMseqs2 server is running"
  fi

  # Count pending jobs (.job files in queue - not yet picked up by server)
  pending_count=$(ls -1 "$JOB_QUEUE_DIR"/*.job 2>/dev/null | wc -l)

  # Count processing jobs (temporary directories created when server picks up job)
  # Jobs are being processed when they have tmp dirs but no .status file yet
  processing_count=0
  if [[ -d "$RESULTS_DIR/../tmp" ]]; then
    for tmp_dir in "$RESULTS_DIR/../tmp"/msa_*; do
      [[ -d "$tmp_dir" ]] || continue
      job_id=$(basename "$tmp_dir")
      # Check if this job doesn't have a status file yet (still processing)
      if [[ ! -f "$RESULTS_DIR/$job_id.status" ]]; then
        processing_count=$((processing_count + 1))
      fi
    done
  fi

  # Total jobs in queue
  total_jobs=$((pending_count + processing_count))

  log "Jobs in queue: $total_jobs (pending: $pending_count, processing: $processing_count)"

  exit 0
fi

if [[ -z "$SEQUENCE" && -z "$FASTA_FILE" ]]; then
  usage
  exit 1
fi

mkdir -p "$JOB_QUEUE_DIR" "$RESULTS_DIR"

# If user gave -f, use it; otherwise write SEQUENCE to a temp FASTA
if [[ -n "$FASTA_FILE" ]]; then
  fasta_file="$FASTA_FILE"
else
  job_id="msa_$(date +%s)_$$"
  fasta_file="$JOB_QUEUE_DIR/$job_id.fasta"
  echo ">query_sequence" > "$fasta_file"
  echo "$SEQUENCE"       >> "$fasta_file"
fi

if [[ -z "$job_id" ]]; then
  job_id="msa_$(date +%s)_$$"
fi

MAX_RETRIES=3
retry_count=0

while [[ $retry_count -lt $MAX_RETRIES ]]; do
  if [[ $retry_count -gt 0 ]]; then
    log "Retry attempt $retry_count/$MAX_RETRIES"
    # Wait a bit before retrying
    sleep 10
  fi

  log "Submitting MSA job: $job_id (format: $OUTPUT_FORMAT)"
  job_file="$JOB_QUEUE_DIR/$job_id.job"

  # Write out all needed params
  {
    echo "job_id=$job_id"
    echo "fasta=$fasta_file"
    echo "output_format=$OUTPUT_FORMAT"
    echo "output_path=$OUTPUT_PATH"
    echo "submitted=$(date '+%Y-%m-%d %H:%M:%S')"
  } > "$job_file"

  log "Job submitted with ID: $job_id"
  log "Waiting for results..."

  start_time=$(date +%s)
  status_file="$RESULTS_DIR/$job_id.status"

  while true; do
    if [[ -f "$status_file" ]]; then
      status=$(sed -n '1p' "$status_file")
      if [[ "$status" == "SUCCESS" ]]; then
        elapsed_time=$(($(date +%s) - start_time))
        if [[ -n "$OUTPUT_PATH" ]]; then
          # Server already copied the file to OUTPUT_PATH
          log "Results saved to: $OUTPUT_PATH (completed in ${elapsed_time}s)"

          # Cleanup: delete status file (server already cleaned up result file)
          rm -f "$status_file"
          log "Cleaned up server status file for job $job_id"
        else
          # pick up output_file=… for backward compatibility (when no OUTPUT_PATH specified)
          result_file=$(sed -n '2p' "$status_file" | cut -d'=' -f2)
          log "Results ready: $result_file (completed in ${elapsed_time}s)"
          # Cleanup status file but keep result in server location
          rm -f "$status_file"
        fi
        exit 0
      else
        # Job failed - check if it's a recoverable failure
        failure_msg=$(sed -n '2p' "$status_file")
        log "Job failed: $failure_msg"

        # Cleanup failed job status
        rm -f "$status_file"

        # Check if this is a server crash recovery (recoverable)
        if [[ "$failure_msg" =~ "Server died" ]] || [[ "$failure_msg" =~ "cannot be recovered" ]]; then
          log "Recoverable failure detected - will retry"
          retry_count=$((retry_count + 1))
          break  # Break inner loop to retry
        else
          # Permanent failure
          log "Permanent failure - not retrying"
          exit 1
        fi
      fi
    fi
    # No timeout - wait indefinitely for server to process
    # SLURM job timeout will handle any issues
    sleep 5
  done
done

# Exhausted retries
log "Exhausted all $MAX_RETRIES retry attempts"
exit 1