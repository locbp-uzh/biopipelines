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
THREADS=32
POLL_INTERVAL=10                     # seconds
MAX_SEQS=10000    # limit homologs per query

mkdir -p "$JOB_QUEUE_DIR" "$RESULTS_DIR" "$TMP_DIR" "$GPU_TMP_DIR"

# Logging setup
LOG_FILE="$RESULTS_DIR/server.log"

# Optimized Memory Settings
export MMSEQS_MAX_MEMORY=${MMSEQS_MAX_MEMORY:-150G}  # Reduced from 200G
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-16}        # Reduced from 32
export CUDA_VISIBLE_DEVICES=0

# GPU Memory optimization
export CUDA_CACHE_MAXSIZE=2147483648  # 2GB CUDA cache
export CUDA_CACHE_DISABLE=0

# Logging function
log() {
  echo "$(date '+%Y-%m-%d %H:%M:%S') - $*" | tee -a "$LOG_FILE"
}

# Check available memory
log "Memory settings: MMSEQS_MAX_MEMORY=$MMSEQS_MAX_MEMORY, OMP_NUM_THREADS=$OMP_NUM_THREADS"
log "GPU Memory status:"
nvidia-smi --query-gpu=memory.total,memory.used,memory.free --format=csv,noheader,nounits | log

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

    /data/gquarg/mmseqs/bin/mmseqs result2msa "$query_db" "$target_db" "$result_db" "$scratch_a3m" \
        --msa-format-mode 5 \
        --threads "$OMP_NUM_THREADS" \
        2>&1 | tee -a "$LOG_FILE"

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
CUDA_VISIBLE_DEVICES=0 /data/gquarg/mmseqs/bin/mmseqs gpuserver "$DB_PATH" \
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
prev_mem=-1
while true; do
  curr_mem=$(nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits | head -1)
  if [[ "$curr_mem" -eq "$prev_mem" ]]; then
    log "GPU memory stabilized at ${curr_mem} MiB — assuming preload is done"
    break
  else
    log "GPU memory at ${curr_mem} MiB (loading…)"
    prev_mem=$curr_mem
    sleep 5
  fi
done

# Cleanup on exit
cleanup() {
  log "Stopping GPU server PID=$GPUSERVER_PID"
  kill $GPUSERVER_PID || true
  # Wait for graceful shutdown
  sleep 5
  # Force kill if still running
  kill -9 $GPUSERVER_PID 2>/dev/null || true
  exit 0
}
trap cleanup SIGINT SIGTERM

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

    # Parse params
    output_format="a3m"; fasta=""
    while IFS='=' read -r key val; do
      case $key in
        fasta)         fasta="$val"         ;;
        output_format) output_format="$val" ;;
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
    /data/gquarg/mmseqs/bin/mmseqs createdb "$tmp/query.fasta" "$query_db" \
     2>&1 | tee -a "$LOG_FILE"

    # Result database (not m8 format)
    result_db="$tmp/resultDB"

    # Set output prefix
    out_prefix="$RESULTS_DIR/$job_id"

    # Run GPU-enabled search with optimized parameters
    log "Running MMseqs2 search for job $job_id"

    # Check GPU memory before search
    gpu_mem_before=$(nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits | head -1)
    log "GPU memory before search: ${gpu_mem_before}MB"

    if ! CUDA_VISIBLE_DEVICES=0 /data/gquarg/mmseqs/bin/mmseqs search "$query_db" "$DB_PATH" "$result_db" "$tmp" \
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
        echo -e "SUCCESS\noutput_file=$out_prefix.a3m" > "$RESULTS_DIR/$job_id.status"
        ;;
      csv)
        log "Converting A3M to CSV format"
        convert_a3m_to_csv "$out_prefix.a3m" "$out_prefix.csv"
        echo -e "SUCCESS\noutput_file=$out_prefix.csv" > "$RESULTS_DIR/$job_id.status"
        ;;
      *)
        log "Unknown format $output_format"
        echo -e "FAILED\nUnknown format" > "$RESULTS_DIR/$job_id.status"
        ;;
    esac

    log "Job $job_id completed successfully"
  done
  sleep 5
done