#!/bin/bash
set -e

# timestamped logging
timestamp() { date '+%Y-%m-%d %H:%M:%S'; }
log() { echo "[$(timestamp)] $*"; }

# MMseqs2 client working directories - must match server script structure
USER=${USER:-$(whoami)}
JOB_QUEUE_DIR="/shares/locbp.chem.uzh/$USER/BioPipelines/MMseqs2Server/job_queue"
RESULTS_DIR="/shares/locbp.chem.uzh/$USER/BioPipelines/MMseqs2Server/results"
TIMEOUT=3600

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
  cpu_running=$(pgrep -f mmseqs2_server_cpu.sh >/dev/null && echo "yes" || echo "no")
  gpu_running=$(pgrep -f mmseqs2_server_gpu.sh >/dev/null && echo "yes" || echo "no")

  if [[ "$cpu_running" == "yes" || "$gpu_running" == "yes" ]]; then
    log "MMseqs2 server is running (CPU: $cpu_running, GPU: $gpu_running)"
  else
    log "MMseqs2 server is not running"
  fi
  job_count=$(ls -1 "$JOB_QUEUE_DIR"/*.job 2>/dev/null | wc -l)
  log "Jobs in queue: $job_count"
  if [[ -f "$RESULTS_DIR/server.log" ]]; then
    log "Last 10 lines of server.log:"
    tail -n 10 "$RESULTS_DIR/server.log"
  fi
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
      if [[ -n "$OUTPUT_PATH" ]]; then
        cp "$result_file" "$OUTPUT_PATH"
        log "Saved results to: $OUTPUT_PATH"
      else
        log "Results ready: $result_file"
      fi
      exit 0
    else
      log "Job failed:"
      sed -n '2p' "$status_file"
      exit 1
    fi
  fi
  if (( $(date +%s) - start_time > TIMEOUT )); then
    log "Timed out waiting for job."
    exit 1
  fi
  sleep 5
done