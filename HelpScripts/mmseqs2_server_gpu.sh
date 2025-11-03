#!/bin/bash
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

# Memory-Optimized MMseqs2 GPU MSA Server Script

# Configuration - must match CPU server and client directory structure
USER=${USER:-$(whoami)}
JOB_QUEUE_DIR="/shares/locbp.chem.uzh/$USER/BioPipelines/MMseqs2Server/job_queue"
RESULTS_DIR="/shares/locbp.chem.uzh/$USER/BioPipelines/MMseqs2Server/results"
DB_DIR="/shares/locbp.chem.uzh/models/mmseqs2_databases/gpu"

TMP_DIR="/shares/locbp.chem.uzh/$USER/BioPipelines/MMseqs2Server/tmp"
GPU_TMP_DIR="/shares/locbp.chem.uzh/$USER/BioPipelines/MMseqs2Server/tmp/gpu"
UNIREF_DB="uniref30_2302_db"
DB_PATH="$DB_DIR/$UNIREF_DB"  # Use databases directly from shares
THREADS=4
POLL_INTERVAL=10                     # seconds
MAX_SEQS=10000    # limit homologs per query

mkdir -p "$JOB_QUEUE_DIR" "$RESULTS_DIR" "$TMP_DIR" "$GPU_TMP_DIR"

# Logging setup - no log file, just stdout (goes to slurm.out)
PID_FILE="$RESULTS_DIR/server_gpu.pid"

# Logging function (output to stdout only)
log() {
  echo "$(date '+%Y-%m-%d %H:%M:%S') - $*"
}

check_mmseqs_installation() {
    local mmseqs_dir="/data/$USER/mmseqs"
    local mmseqs_bin="$mmseqs_dir/bin/mmseqs"

    if [[ ! -f "$mmseqs_bin" ]]; then
        echo "$(date '+%Y-%m-%d %H:%M:%S') - MMseqs2 not found at $mmseqs_bin, downloading..."
        cd "/data/$USER" || { echo "$(date '+%Y-%m-%d %H:%M:%S') - ERROR: Cannot access /data/$USER"; exit 1; }

        # Download MMseqs2
        wget https://mmseqs.com/latest/mmseqs-linux-gpu.tar.gz

        # Extract and cleanup
        tar xvfz mmseqs-linux-gpu.tar.gz
        rm mmseqs-linux-gpu.tar.gz

        # Verify installation
        if [[ -f "$mmseqs_bin" ]]; then
            echo "$(date '+%Y-%m-%d %H:%M:%S') - MMseqs2 successfully installed at $mmseqs_bin"
        else
            echo "$(date '+%Y-%m-%d %H:%M:%S') - ERROR: MMseqs2 installation failed"
            exit 1
        fi
    else
        echo "$(date '+%Y-%m-%d %H:%M:%S') - MMseqs2 found at $mmseqs_bin"
    fi
}

echo $$ > "$PID_FILE"   # record server PID

# Create timestamp file for server detection (per-user)
MMSEQS_SERVER_DIR="/shares/locbp.chem.uzh/$USER/BioPipelines/MMseqs2Server"
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

# Check and install MMseqs2 if needed
check_mmseqs_installation

# Optimized Memory Settings
export MMSEQS_MAX_MEMORY=${MMSEQS_MAX_MEMORY:-150G}  # Reduced from 200G
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}         # Reduced to 4 for GPU mode
export CUDA_VISIBLE_DEVICES=0

# GPU Memory optimization
export CUDA_CACHE_MAXSIZE=2147483648  # 2GB CUDA cache
export CUDA_CACHE_DISABLE=0

# Check available memory
log "Memory settings: MMSEQS_MAX_MEMORY=$MMSEQS_MAX_MEMORY, OMP_NUM_THREADS=$OMP_NUM_THREADS"

convert_to_a3m() {
    local result_db=$1
    local query_db=$2
    local target_db=$3
    local output_file=$4
    local original_fasta=$5
    local tmp_dir=$6

    # write raw A3M to scratch, then strip nulls
    local scratch_a3m="$tmp_dir/raw.a3m"
    local final_a3m="$output_file"

    /data/$USER/mmseqs/bin/mmseqs result2msa "$query_db" "$target_db" "$result_db" "$scratch_a3m" \
        --msa-format-mode 5 \
        --threads "$OMP_NUM_THREADS"

    # Dropping the last byte which is a null character and makes boltz crash
    log "Removing last byte from $scratch_a3m"
    head -c -1 "$scratch_a3m" > "${scratch_a3m}.clean"
    mv "${scratch_a3m}.clean" "$final_a3m"
    rm "$scratch_a3m"

    local seq_count=$(grep -c "^>" "$final_a3m" || echo "0")
    log "A3M file created with $seq_count sequences"
}

convert_a3m_to_csv() {
    local a3m_file=$1
    local output_file=$2

    # Write CSV header
    echo "key,sequence" > "$output_file"

    # Extract sequences (ignore headers) and prefix with -1
    awk '/^[^>]/ {print "-1," $0}' "$a3m_file" >> "$output_file"

    # Log the number of sequences written
    local line_count=$(wc -l < "$output_file")
    log "CSV file created with $((line_count - 1)) data rows"
}

# Start GPU server with optimized settings
log "Starting MMseqs2 GPU server for $DB_PATH"
CUDA_VISIBLE_DEVICES=0 /data/$USER/mmseqs/bin/mmseqs gpuserver "$DB_PATH" \
  --max-seqs "$MAX_SEQS" \
  --db-load-mode 0 \
  --prefilter-mode 1 &
GPUSERVER_PID=$!
log "GPU server PID=$GPUSERVER_PID"

# Test GPU server connection
log "Testing GPU server connection..."
if ! ps -p $GPUSERVER_PID > /dev/null; then
    log "ERROR: GPU server failed to start"
    exit 1
fi

# Wait until gpuserver has actually loaded the DB into GPU memory
log "Waiting for gpuserver to finish preloading DB into GPU memory"
log "Database should load to ~9-10GB of GPU memory"

# Wait for GPU memory to actually increase (database loading)
# Timeout after 10 minutes
max_wait_seconds=600
wait_start=$(date +%s)
prev_mem=-1
db_loaded=false

while true; do
  curr_mem=$(nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits | head -1)
  current_time=$(date +%s)
  elapsed=$((current_time - wait_start))

  # Check if database is loaded (GPU memory should be ~9000+ MiB)
  if [[ "$curr_mem" -gt 9000 ]]; then
    log "GPU memory at ${curr_mem} MiB - database loaded successfully"
    db_loaded=true
    break
  fi

  # Check timeout
  if [[ $elapsed -gt $max_wait_seconds ]]; then
    log "WARNING: Timeout waiting for database to load after ${max_wait_seconds}s"
    log "GPU memory only at ${curr_mem} MiB (expected >9000 MiB)"
    log "Proceeding anyway, but searches may be slow or fail"
    break
  fi

  # Log progress
  if [[ "$curr_mem" -ne "$prev_mem" ]]; then
    log "GPU memory at ${curr_mem} MiB (loading... ${elapsed}s elapsed)"
    prev_mem=$curr_mem
  fi

  sleep 5
done

if [[ "$db_loaded" == "true" ]]; then
  log "GPU server initialization complete and ready to process jobs"
else
  log "GPU server may not be fully initialized"
fi

# Cleanup on exit
cleanup() {
  log "MMseqs2 GPU server shutting down"
  log "Stopping GPU server PID=$GPUSERVER_PID"
  kill $GPUSERVER_PID || true
  # Wait for graceful shutdown
  sleep 5
  # Force kill if still running
  kill -9 $GPUSERVER_PID 2>/dev/null || true
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

    # Clean old result files (.a3m, .csv, .status) - these are copied to client paths anyway
    local status_count=$(find "$RESULTS_DIR" -type f -name "*.status" -mtime +1 -delete -print 2>/dev/null | wc -l)
    local a3m_count=$(find "$RESULTS_DIR" -type f -name "*.a3m" -mtime +1 -delete -print 2>/dev/null | wc -l)
    local csv_count=$(find "$RESULTS_DIR" -type f -name "*.csv" -mtime +1 -delete -print 2>/dev/null | wc -l)

    if [ "$status_count" -gt 0 ]; then
        log "Removed $status_count old status files"
    fi
    if [ "$a3m_count" -gt 0 ]; then
        log "Removed $a3m_count old a3m files"
    fi
    if [ "$csv_count" -gt 0 ]; then
        log "Removed $csv_count old csv files"
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
    output_format="a3m"; fasta=""; output_path=""
    while IFS='=' read -r key val; do
      case $key in
        fasta)         fasta="$val"         ;;
        output_format) output_format="$val" ;;
        output_path)   output_path="$val"   ;;
      esac
    done < "$tmp/params/$(basename "$job_meta")"

    # Fallback if client didn't send fasta=
    [[ -n "$fasta" ]] || fasta="$JOB_QUEUE_DIR/${job_id}.fasta"
    if [[ ! -f "$fasta" ]]; then
      log "ERROR: FASTA not found at $fasta"
      echo -e "FAILED\nInput FASTA not found" > "$RESULTS_DIR/$job_id.status"
      continue
    fi

    # Prepare query DB
    cp "$fasta" "$tmp/query.fasta"
    query_db="$tmp/queryDB"
    log "Creating queryDB"
    /data/$USER/mmseqs/bin/mmseqs createdb "$tmp/query.fasta" "$query_db"

    # Result database (not m8 format)
    result_db="$tmp/resultDB"

    # Set output prefix
    out_prefix="$RESULTS_DIR/$job_id"

    # Run GPU-enabled search with optimized parameters
    log "Running MMseqs2 search for job $job_id"

    # Check GPU memory before search
    gpu_mem_before=$(nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits | head -1)
    log "GPU memory before search: ${gpu_mem_before}MB"

    if ! CUDA_VISIBLE_DEVICES=0 /data/$USER/mmseqs/bin/mmseqs search "$query_db" "$DB_PATH" "$result_db" "$tmp" \
      --gpu 1 \
      --gpu-server 1 \
      --prefilter-mode 1 \
      --db-load-mode 2 \
      -a 1 \
      --alignment-mode 0 \
      --threads "$OMP_NUM_THREADS" \
      --max-seqs "$MAX_SEQS" \
      -s 7.5; then
      log "Search failed for job $job_id"
      # Log GPU memory after failure
      gpu_mem_after=$(nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits | head -1)
      log "GPU memory after failed search: ${gpu_mem_after}MB"
      echo "FAILED: Search failed" > "$RESULTS_DIR/$job_id.status"
      continue
    fi

    # Log GPU memory after successful search
    gpu_mem_after=$(nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits | head -1)
    log "GPU memory after successful search: ${gpu_mem_after}MB"

    log "Converting to A3M format"
    convert_to_a3m "$result_db" "$query_db" "$DB_PATH" "$out_prefix.a3m" "$tmp/query.fasta" "$tmp"
    # Convert output based on format
    case $output_format in
      a3m)
        if [[ -n "$output_path" ]]; then
          # Copy to client-specified path
          mkdir -p "$(dirname "$output_path")"
          cp "$out_prefix.a3m" "$output_path"
          log "Copied result to client path: $output_path"
          # Cleanup server result file
          rm -f "$out_prefix.a3m"
          echo "SUCCESS" > "$RESULTS_DIR/$job_id.status"
        else
          # Legacy mode: keep result in server location
          echo -e "SUCCESS\noutput_file=$out_prefix.a3m" > "$RESULTS_DIR/$job_id.status"
        fi
        ;;
      csv)
        log "Converting A3M to CSV format"
        convert_a3m_to_csv "$out_prefix.a3m" "$out_prefix.csv"
        if [[ -n "$output_path" ]]; then
          # Copy to client-specified path
          mkdir -p "$(dirname "$output_path")"
          cp "$out_prefix.csv" "$output_path"
          log "Copied result to client path: $output_path"
          # Cleanup server result files
          rm -f "$out_prefix.a3m" "$out_prefix.csv"
          echo "SUCCESS" > "$RESULTS_DIR/$job_id.status"
        else
          # Legacy mode: keep result in server location
          echo -e "SUCCESS\noutput_file=$out_prefix.csv" > "$RESULTS_DIR/$job_id.status"
          # Delete intermediate A3M file since client only needs CSV
          rm -f "$out_prefix.a3m"
          log "Deleted intermediate A3M file"
        fi
        ;;
      *)
        log "Unknown format $output_format"
        echo -e "FAILED\nUnknown format" > "$RESULTS_DIR/$job_id.status"
        ;;
    esac

    log "Job $job_id completed successfully"

    # Cleanup: delete job file and FASTA input
    rm -f "$job_file_path"
    [[ -f "$fasta" ]] && rm -f "$fasta"
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