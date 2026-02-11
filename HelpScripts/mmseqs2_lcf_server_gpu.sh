#!/bin/bash
# Copyright (c) 2026 Gianluca Quargnali @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# This software is freely available for use, modification, and redistribution.

set -e

# Check if nvidia-smi is available
if ! command -v nvidia-smi &> /dev/null
then
    echo "Could not load GPU correctly: nvidia-smi could not be found"
    exit
fi

# Display GPU information
gpu_type=$(nvidia-smi --query-gpu=gpu_name --format=csv,noheader)
echo "GPU Type: $gpu_type"
nvidia-smi --query-gpu=memory.total,memory.used,memory.free --format=csv,noheader,nounits

# MMseqs2 LCF (LocalColabFold) GPU MSA Server Script
# Uses colabfold_search for richer MSAs from UniRef30 + ColabFoldDB

# Resolve a required folder: env var first, then biopipelines-config, else fail
require_folder() {
  local val="$1" name="$2" key="$3"
  if [[ -n "$val" ]]; then echo "$val"; return; fi
  val=$(biopipelines-config folder "$key" 2>/dev/null) && [[ -n "$val" ]] && { echo "$val"; return; }
  echo "ERROR: $name is not set and biopipelines-config could not resolve '$key'." >&2
  echo "Set the environment variable or configure the folder in config.yaml." >&2
  exit 1
}

# Configuration - use environment variables from pipeline, with biopipelines-config fallback
USER=${USER:-$(whoami)}
MMSEQS2_SHARED_FOLDER=$(require_folder "${MMSEQS2_SHARED_FOLDER:-}" "MMSEQS2_SHARED_FOLDER" "MMseqs2LCFServer")
JOB_QUEUE_DIR="$MMSEQS2_SHARED_FOLDER/job_queue"
RESULTS_DIR="$MMSEQS2_SHARED_FOLDER/results"
TMP_DIR="$MMSEQS2_SHARED_FOLDER/tmp"

# Database directory (ColabFold databases setup with setup_databases.sh)
DB_DIR=$(require_folder "${COLABFOLD_DB_DIR:-}" "COLABFOLD_DB_DIR" "ColabFoldDatabases")
DATA_DIR=$(require_folder "${BIOPIPELINES_DATA_DIR:-}" "BIOPIPELINES_DATA_DIR" "data")
MMSEQS2_DIR=$(require_folder "${MMSEQS2_DIR:-}" "MMSEQS2_DIR" "MMseqs2")
ALPHAFOLD_DIR=$(require_folder "${BIOPIPELINES_ALPHAFOLD_DIR:-}" "BIOPIPELINES_ALPHAFOLD_DIR" "AlphaFold")
UNIREF_DB="uniref30_2302_db"
ENVDB="colabfold_envdb_202108_db"

THREADS=${OMP_NUM_THREADS:-16}
POLL_INTERVAL=${MMSEQS2_POLL_INTERVAL:-10}

mkdir -p "$JOB_QUEUE_DIR" "$RESULTS_DIR" "$TMP_DIR"

# Logging setup - no log file, just stdout (goes to slurm.out)
PID_FILE="$RESULTS_DIR/server_gpu.pid"

# Logging function (output to stdout only)
log() {
  echo "$(date '+%Y-%m-%d %H:%M:%S') - $*"
}

check_mmseqs_installation() {
    local mmseqs_bin="$MMSEQS2_DIR/bin/mmseqs"

    if [[ ! -f "$mmseqs_bin" ]]; then
        log "MMseqs2 not found at $mmseqs_bin, downloading..."
        local parent_dir
        parent_dir="$(dirname "$MMSEQS2_DIR")"
        cd "$parent_dir" || { log "ERROR: Cannot access $parent_dir"; exit 1; }

        # Download MMseqs2
        wget https://mmseqs.com/latest/mmseqs-linux-gpu.tar.gz

        # Extract and cleanup
        tar xvfz mmseqs-linux-gpu.tar.gz
        rm mmseqs-linux-gpu.tar.gz

        # Verify installation
        if [[ -f "$mmseqs_bin" ]]; then
            log "MMseqs2 successfully installed at $mmseqs_bin"
        else
            log "ERROR: MMseqs2 installation failed"
            exit 1
        fi
    else
        log "MMseqs2 found at $mmseqs_bin"
    fi
}

check_colabfold_search() {
    local localcolabfold_dir="$ALPHAFOLD_DIR"
    local colabfold_search="$localcolabfold_dir/colabfold-conda/bin/colabfold_search"

    if [[ ! -f "$colabfold_search" ]]; then
        log "ERROR: colabfold_search not found at $colabfold_search"
        log "Please ensure localcolabfold is installed at $ALPHAFOLD_DIR"
        exit 1
    else
        log "colabfold_search found at $colabfold_search"
    fi
}

check_databases() {
    local uniref_path="$DB_DIR/$UNIREF_DB"
    local envdb_path="$DB_DIR/$ENVDB"

    if [[ ! -f "$uniref_path" ]]; then
        log "ERROR: UniRef30 database not found at $uniref_path"
        log "Please run ColabFold's setup_databases.sh with GPU=1"
        exit 1
    fi

    if [[ ! -f "$envdb_path" ]]; then
        log "ERROR: ColabFoldDB (envdb) not found at $envdb_path"
        log "Please run ColabFold's setup_databases.sh with GPU=1"
        exit 1
    fi

    log "Databases found:"
    log "  UniRef30: $uniref_path"
    log "  EnvDB: $envdb_path"
}

echo $$ > "$PID_FILE"   # record server PID

# Create timestamp file for server detection (per-user)
MMSEQS_SERVER_DIR="$MMSEQS2_SHARED_FOLDER"
mkdir -p "$MMSEQS_SERVER_DIR"
SERVER_TIMESTAMP_FILE="$MMSEQS_SERVER_DIR/GPU_SERVER"
SUBMITTING_FILE="$MMSEQS_SERVER_DIR/GPU_SUBMITTING"

# Clean up submission timestamp if it exists (server is now starting)
if [[ -f "$SUBMITTING_FILE" ]]; then
  log "Cleaning up submission timestamp"
  rm -f "$SUBMITTING_FILE"
fi

date '+%H:%M:%S' > "$SERVER_TIMESTAMP_FILE"
log "Created server timestamp file at $SERVER_TIMESTAMP_FILE"

# Check installations and databases
check_mmseqs_installation
check_colabfold_search
check_databases

# Paths to binaries
MMSEQS_BIN="$MMSEQS2_DIR/bin/mmseqs"
COLABFOLD_SEARCH="$ALPHAFOLD_DIR/colabfold-conda/bin/colabfold_search"

# Optimized Memory Settings
export MMSEQS_MAX_MEMORY=${MMSEQS_MAX_MEMORY:-150G}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-16}
export CUDA_VISIBLE_DEVICES=0

# GPU Memory optimization
export CUDA_CACHE_MAXSIZE=2147483648  # 2GB CUDA cache
export CUDA_CACHE_DISABLE=0

log "Memory settings: MMSEQS_MAX_MEMORY=$MMSEQS_MAX_MEMORY, OMP_NUM_THREADS=$OMP_NUM_THREADS"

# Start GPU servers for both databases with optimized settings
log "Starting MMseqs2 GPU server for UniRef30: $DB_DIR/$UNIREF_DB"
CUDA_VISIBLE_DEVICES=0 $MMSEQS_BIN gpuserver "$DB_DIR/$UNIREF_DB" \
  --db-load-mode 0 \
  --prefilter-mode 1 &
GPUSERVER_UNIREF_PID=$!
log "UniRef30 GPU server PID=$GPUSERVER_UNIREF_PID"

log "Starting MMseqs2 GPU server for EnvDB: $DB_DIR/$ENVDB"
CUDA_VISIBLE_DEVICES=0 $MMSEQS_BIN gpuserver "$DB_DIR/$ENVDB" \
  --db-load-mode 0 \
  --prefilter-mode 1 &
GPUSERVER_ENVDB_PID=$!
log "EnvDB GPU server PID=$GPUSERVER_ENVDB_PID"

# Wait for databases to load into GPU memory (envdb is large, needs more time)
log "Waiting 5 minutes for databases to load into GPU memory..."
sleep 300
log "Wait complete, assuming GPU servers are ready"

# Cleanup on exit
cleanup() {
  log "MMseqs2 LCF GPU server shutting down"
  log "Stopping UniRef30 GPU server PID=$GPUSERVER_UNIREF_PID"
  kill $GPUSERVER_UNIREF_PID || true
  log "Stopping EnvDB GPU server PID=$GPUSERVER_ENVDB_PID"
  kill $GPUSERVER_ENVDB_PID || true
  # Wait for graceful shutdown
  sleep 5
  # Force kill if still running
  kill -9 $GPUSERVER_UNIREF_PID 2>/dev/null || true
  kill -9 $GPUSERVER_ENVDB_PID 2>/dev/null || true
  rm -f "$PID_FILE"
  # Remove timestamp file on shutdown
  rm -f "$SERVER_TIMESTAMP_FILE"
  log "Removed server timestamp file"
  exit 0
}
trap cleanup SIGINT SIGTERM

cleanup_old_files() {
    # Delete files older than 24 hours
    log "Running cleanup of files older than 24 hours..."

    # Clean old result files (.a3m, .status) - these are copied to client paths anyway
    local status_count=$(find "$RESULTS_DIR" -type f -name "*.status" -mtime +1 -delete -print 2>/dev/null | wc -l)
    local a3m_count=$(find "$RESULTS_DIR" -type f -name "*.a3m" -mtime +1 -delete -print 2>/dev/null | wc -l)

    if [ "$status_count" -gt 0 ]; then
        log "Removed $status_count old status files"
    fi
    if [ "$a3m_count" -gt 0 ]; then
        log "Removed $a3m_count old a3m files"
    fi

    # Clean old temporary directories (these accumulate and should be deleted)
    local tmp_dir_count=$(find "$TMP_DIR" -mindepth 1 -maxdepth 1 -type d -mtime +1 -print 2>/dev/null | wc -l)
    if [ "$tmp_dir_count" -gt 0 ]; then
        find "$TMP_DIR" -mindepth 1 -maxdepth 1 -type d -mtime +1 -exec rm -rf {} \; 2>/dev/null || true
        log "Removed $tmp_dir_count old tmp directories"
    fi

    # Clean old FASTA files in queue (orphaned submissions)
    local fasta_count=$(find "$JOB_QUEUE_DIR" -type f -name "*.fasta" -mtime +1 -delete -print 2>/dev/null | wc -l)
    if [ "$fasta_count" -gt 0 ]; then
        log "Removed $fasta_count old FASTA files"
    fi

    # Clean orphaned .job files (older than 24 hours, indicating stale/failed submissions)
    local orphaned_count=$(find "$JOB_QUEUE_DIR" -type f -name "*.job" -mtime +1 -delete -print 2>/dev/null | wc -l)
    if [ "$orphaned_count" -gt 0 ]; then
        log "Removed $orphaned_count orphaned .job files"
    fi

    log "Cleanup completed"
}

recover_orphaned_jobs() {
    # Recover jobs that were picked up by a dead server
    # If a job has been in tmp for > 30 minutes without a status file, requeue it
    local max_age_seconds=1800  # 30 minutes
    local current_time=$(date +%s)
    local recovered_count=0
    local failed_count=0

    log "Scanning for orphaned jobs (older than 30 minutes)..."

    for tmp_job_dir in "$TMP_DIR"/msa_*; do
        [[ -d "$tmp_job_dir" ]] || continue

        local job_id=$(basename "$tmp_job_dir")
        local status_file="$RESULTS_DIR/$job_id.status"

        # Skip if status file exists (job completed)
        [[ -f "$status_file" ]] && continue

        # Check age of tmp directory
        local dir_mtime=$(stat -c %Y "$tmp_job_dir" 2>/dev/null || stat -f %m "$tmp_job_dir" 2>/dev/null || echo "0")
        local age=$((current_time - dir_mtime))

        if [[ $age -gt $max_age_seconds ]]; then
            log "Recovering orphaned job: $job_id (age: $((age / 60)) minutes)"

            # Look for the .job file in the tmp/params directory
            local job_file="$tmp_job_dir/params/${job_id}.job"
            if [[ -f "$job_file" ]]; then
                # Move .job file back to queue for reprocessing
                mv "$job_file" "$JOB_QUEUE_DIR/" && log "Requeued job file: ${job_id}.job"
                recovered_count=$((recovered_count + 1))
            else
                # No .job file - cannot recover, mark as failed so client stops waiting
                log "Warning: No .job file found for orphaned job $job_id - marking as failed"
                echo -e "FAILED\nServer died during processing and job cannot be recovered" > "$status_file"
                failed_count=$((failed_count + 1))
            fi

            # Clean up the tmp directory
            rm -rf "$tmp_job_dir"
        fi
    done

    if [[ $recovered_count -gt 0 ]]; then
        log "Recovered $recovered_count orphaned jobs"
    fi
    if [[ $failed_count -gt 0 ]]; then
        log "Marked $failed_count unrecoverable orphaned jobs as failed"
    fi
}

# Track cleanup time
LAST_CLEANUP=$(date +%s)
CLEANUP_INTERVAL=3600  # Run cleanup every hour

# Run cleanup and recovery at startup (AFTER gpuserver is fully initialized)
log "Running startup cleanup to remove old files..."
cleanup_old_files
log "Running startup job recovery..."
recover_orphaned_jobs
LAST_CLEANUP=$(date +%s)  # Reset cleanup timer after startup cleanup

log "Entering job processing loop"

# Main processing loop
while true; do
  for job_meta in "$JOB_QUEUE_DIR"/*.job; do
    [[ -e "$job_meta" ]] || continue

    fname=$(basename "$job_meta")
    job_id="${fname%.job}"
    tmp="$TMP_DIR/$job_id"
    mkdir -p "$tmp"
    mkdir -p "$tmp/params"

    log "Picked up job $job_id"
    mv "$job_meta" "$tmp/params/"

    # Store job file path for cleanup after completion
    job_file_path="$tmp/params/$fname"

    # Parse params
    fasta=""; output_path=""
    while IFS='=' read -r key val; do
      case $key in
        fasta)         fasta="$val"       ;;
        output_path)   output_path="$val" ;;
      esac
    done < "$tmp/params/$(basename "$job_meta")"

    # Fallback if client didn't send fasta=
    [[ -n "$fasta" ]] || fasta="$JOB_QUEUE_DIR/${job_id}.fasta"
    if [[ ! -f "$fasta" ]]; then
      log "ERROR: FASTA not found at $fasta"
      echo -e "FAILED\nInput FASTA not found" > "$RESULTS_DIR/$job_id.status"
      continue
    fi

    # Copy fasta to tmp directory for colabfold_search
    cp "$fasta" "$tmp/query.fasta"

    # Create output directory for colabfold_search
    search_output="$tmp/search_output"
    mkdir -p "$search_output"

    # Run colabfold_search with GPU servers
    log "Running colabfold_search for job $job_id"

    # Check GPU memory before search
    gpu_mem_before=$(nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits | head -1)
    log "GPU memory before search: ${gpu_mem_before}MB"

    if ! $COLABFOLD_SEARCH "$tmp/query.fasta" "$DB_DIR" "$search_output" \
      --mmseqs "$MMSEQS_BIN" \
      --gpu 1 \
      --gpu-server 1 \
      --threads "$THREADS" \
      --use-env 1 \
      --use-templates 0 \
      --db1 "$UNIREF_DB" \
      --db3 "$ENVDB"; then
      log "colabfold_search failed for job $job_id"
      # Log GPU memory after failure
      gpu_mem_after=$(nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits | head -1)
      log "GPU memory after failed search: ${gpu_mem_after}MB"
      echo "FAILED: colabfold_search failed" > "$RESULTS_DIR/$job_id.status"
      continue
    fi

    # Log GPU memory after successful search
    gpu_mem_after=$(nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits | head -1)
    log "GPU memory after successful search: ${gpu_mem_after}MB"

    # Find the output a3m file (colabfold_search creates <query_name>.a3m)
    # The query name is typically "query_sequence" from our FASTA header
    a3m_file=$(find "$search_output" -name "*.a3m" -type f | head -1)

    if [[ -z "$a3m_file" ]] || [[ ! -f "$a3m_file" ]]; then
      log "ERROR: No A3M output found for job $job_id"
      echo "FAILED: No A3M output generated" > "$RESULTS_DIR/$job_id.status"
      continue
    fi

    log "Found A3M output: $a3m_file"

    # Copy result to client-specified path or results directory
    if [[ -n "$output_path" ]]; then
      mkdir -p "$(dirname "$output_path")"
      cp "$a3m_file" "$output_path"
      log "Copied result to client path: $output_path"
      echo "SUCCESS" > "$RESULTS_DIR/$job_id.status"
    else
      # Legacy mode: keep result in server location
      cp "$a3m_file" "$RESULTS_DIR/$job_id.a3m"
      echo -e "SUCCESS\noutput_file=$RESULTS_DIR/$job_id.a3m" > "$RESULTS_DIR/$job_id.status"
    fi

    # Count sequences in the A3M
    seq_count=$(grep -c "^>" "$a3m_file" || echo "0")
    log "Job $job_id completed successfully with $seq_count sequences in MSA"

    # Cleanup: delete job file and FASTA input
    rm -f "$job_file_path"
    [[ -f "$fasta" ]] && rm -f "$fasta"
    rm -rf "$tmp"
    log "Cleaned up job files for $job_id"
  done

  # Periodic cleanup and recovery
  CURRENT_TIME=$(date +%s)
  if (( CURRENT_TIME - LAST_CLEANUP >= CLEANUP_INTERVAL )); then
    cleanup_old_files
    recover_orphaned_jobs
    LAST_CLEANUP=$CURRENT_TIME
  fi

  sleep 5
done
