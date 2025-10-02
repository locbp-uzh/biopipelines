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
  # Check CPU server using PID file
  cpu_pid_file="$RESULTS_DIR/server_cpu.pid"
  if [[ -f "$cpu_pid_file" ]]; then
    cpu_pid=$(cat "$cpu_pid_file")
    if kill -0 "$cpu_pid" 2>/dev/null; then
      cpu_running="yes"
    else
      cpu_running="no (stale PID file)"
      rm -f "$cpu_pid_file"  # cleanup stale PID file
    fi
  else
    cpu_running="no"
  fi

  # Check GPU server using PID file
  gpu_pid_file="$RESULTS_DIR/server_gpu.pid"
  if [[ -f "$gpu_pid_file" ]]; then
    gpu_pid=$(cat "$gpu_pid_file")
    if kill -0 "$gpu_pid" 2>/dev/null; then
      gpu_running="yes"
    else
      gpu_running="no (stale PID file)"
      rm -f "$gpu_pid_file"  # cleanup stale PID file
    fi
  else
    gpu_running="no"
  fi

  if [[ "$cpu_running" == "yes" || "$gpu_running" == "yes" ]]; then
    log "MMseqs2 server is running (CPU: $cpu_running, GPU: $gpu_running)"

    # Check server activity using log timestamps
    if [[ -f "$RESULTS_DIR/server.log" ]]; then
      last_activity=$(tail -n 1 "$RESULTS_DIR/server.log" | grep -o '\[.*\]' | tr -d '[]' 2>/dev/null || echo "unknown")
      if [[ "$last_activity" != "unknown" ]]; then
        log "Last server activity: $last_activity"

        # Check if server is actively processing (recent activity within 5 minutes)
        if command -v date >/dev/null 2>&1; then
          current_time=$(date +%s)
          last_time=$(date -d "$last_activity" +%s 2>/dev/null || echo "0")
          time_diff=$((current_time - last_time))
          if [[ $time_diff -lt 300 ]]; then  # 5 minutes
            log "Server is actively processing (last activity ${time_diff}s ago)"
          else
            log "Server may be idle (last activity ${time_diff}s ago)"
          fi
        fi
      fi
    fi
  else
    log "MMseqs2 server is not running (CPU: $cpu_running, GPU: $gpu_running)"
  fi

  # Count pending jobs (.job files in queue)
  job_count=$(ls -1 "$JOB_QUEUE_DIR"/*.job 2>/dev/null | wc -l)
  log "Jobs in queue: $job_count"

  # Count processing jobs (check for .status files that are still being written)
  processing_count=$(find "$RESULTS_DIR" -name "*.status" -type f -mmin -5 2>/dev/null | wc -l)
  if [[ $processing_count -gt 0 ]]; then
    log "Jobs being processed: $processing_count"
  fi

  exit 0
fi

if [[ -z "$SEQUENCE" && -z "$FASTA_FILE" ]]; then
  usage
  exit 1
fi

mkdir -p "$JOB_QUEUE_DIR" "$RESULTS_DIR"

# Check if server is running, if not, auto-start it
check_and_start_server() {
  local cpu_pid_file="$RESULTS_DIR/server_cpu.pid"
  local gpu_pid_file="$RESULTS_DIR/server_gpu.pid"
  local cpu_running=false
  local gpu_running=false

  # Check CPU server
  if [[ -f "$cpu_pid_file" ]]; then
    local cpu_pid=$(cat "$cpu_pid_file")
    if kill -0 "$cpu_pid" 2>/dev/null; then
      cpu_running=true
    fi
  fi

  # Check GPU server
  if [[ -f "$gpu_pid_file" ]]; then
    local gpu_pid=$(cat "$gpu_pid_file")
    if kill -0 "$gpu_pid" 2>/dev/null; then
      gpu_running=true
    fi
  fi

  # If neither server is running, start one
  if [[ "$cpu_running" == false && "$gpu_running" == false ]]; then
    log "No MMseqs2 server detected, attempting to auto-start..."

    # Find the biopipelines root directory (where submit script is located)
    local script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    local biopipelines_root="$(cd "$script_dir/.." && pwd)"
    local submit_script="$biopipelines_root/submit"
    local server_script="$biopipelines_root/ExamplePipelines/mmseqs2_server.py"

    if [[ -f "$submit_script" && -f "$server_script" ]]; then
      log "Submitting MMseqs2 server: $submit_script $server_script"

      # Submit the server (run in background to not block)
      (cd "$biopipelines_root" && "$submit_script" "$server_script") &

      log "Server submitted. Waiting 10 seconds for server to start..."
      sleep 10

      # Verify server started
      local max_wait=60
      local waited=0
      while [[ $waited -lt $max_wait ]]; do
        if [[ -f "$cpu_pid_file" || -f "$gpu_pid_file" ]]; then
          log "Server started successfully"
          return 0
        fi
        sleep 5
        waited=$((waited + 5))
      done

      log "WARNING: Server may not have started properly. Proceeding anyway..."
    else
      log "WARNING: Cannot auto-start server - submit script or server script not found"
      log "  Submit script: $submit_script"
      log "  Server script: $server_script"
      log "Please start the server manually"
    fi
  fi
}

# Auto-start server if needed
check_and_start_server

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
  if (( $(date +%s) - start_time > TIMEOUT )); then
    log "Timed out waiting for job."
    exit 1
  fi
  sleep 5
done