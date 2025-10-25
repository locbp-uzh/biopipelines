#!/bin/bash
set -e

# timestamped logging
timestamp() { date '+%Y-%m-%d %H:%M:%S'; }
log() { echo "[$(timestamp)] $*"; }

# MMseqs2 client working directories - must match server script structure
USER=${USER:-$(whoami)}
JOB_QUEUE_DIR="/shares/locbp.chem.uzh/$USER/BioPipelines/MMseqs2Server/job_queue"
RESULTS_DIR="/shares/locbp.chem.uzh/$USER/BioPipelines/MMseqs2Server/results"

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
  # Check servers using timestamp files
  MMSEQS_SERVER_DIR="/shares/locbp.chem.uzh/models/mmseqs2_server"
  GPU_TIMESTAMP="$MMSEQS_SERVER_DIR/GPU_SERVER"
  CPU_TIMESTAMP="$MMSEQS_SERVER_DIR/CPU_SERVER"
  MAX_AGE_HOURS=12

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

log "Submitting MSA job: $job_id (format: $OUTPUT_FORMAT)"
job_file="$JOB_QUEUE_DIR/$job_id.job"

# Write out all needed params
{
  echo "job_id=$job_id"
  echo "fasta=$fasta_file"
  echo "output_format=$OUTPUT_FORMAT"
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
      # pick up output_file=…
      result_file=$(sed -n '2p' "$status_file" | cut -d'=' -f2)
      elapsed_time=$(($(date +%s) - start_time))
      if [[ -n "$OUTPUT_PATH" ]]; then
        cp "$result_file" "$OUTPUT_PATH"
        log "Saved results to: $OUTPUT_PATH (completed in ${elapsed_time}s)"

        # Cleanup: delete server files after successful copy
        rm -f "$result_file" "$status_file"
        log "Cleaned up server files for job $job_id"
      else
        log "Results ready: $result_file (completed in ${elapsed_time}s)"
        # If no output path, still cleanup status file but keep result
        rm -f "$status_file"
      fi
      exit 0
    else
      log "Job failed:"
      sed -n '2p' "$status_file"
      # Cleanup failed job status
      rm -f "$status_file"
      exit 1
    fi
  fi
  # No timeout - wait indefinitely for server to process
  # SLURM job timeout will handle any issues
  sleep 5
done