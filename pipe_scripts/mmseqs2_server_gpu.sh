#!/bin/bash
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
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

# Memory-Optimized MMseqs2 GPU MSA Server Script

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
MMSEQS2_SHARED_FOLDER=$(require_folder "${MMSEQS2_SHARED_FOLDER:-}" "MMSEQS2_SHARED_FOLDER" "MMseqs2Server")
JOB_QUEUE_DIR="$MMSEQS2_SHARED_FOLDER/job_queue"
RESULTS_DIR="$MMSEQS2_SHARED_FOLDER/results"
DB_DIR=$(require_folder "${MMSEQS2_DB_DIR:-}" "MMSEQS2_DB_DIR" "MMseqs2Databases")
DATA_DIR=$(require_folder "${BIOPIPELINES_DATA_DIR:-}" "BIOPIPELINES_DATA_DIR" "data")
MMSEQS2_DIR=$(require_folder "${MMSEQS2_DIR:-}" "MMSEQS2_DIR" "MMseqs2")

TMP_DIR="$MMSEQS2_SHARED_FOLDER/tmp"
GPU_TMP_DIR="$MMSEQS2_SHARED_FOLDER/tmp/gpu"
# colabfold_search searches both UniRef30 (--db1) and the environmental DB
# (--db3) by default. We keep a warm gpuserver resident for each so the search
# step can use --gpu-server 1.
UNIREF_DB="uniref30_2302_db"
ENVDB="colabfold_envdb_202108_db"
UNIREF_PATH="$DB_DIR/$UNIREF_DB"  # databases used directly from shares
ENVDB_PATH="$DB_DIR/$ENVDB"
THREADS=4
POLL_INTERVAL=10                     # seconds
MAX_SEQS=10000    # limit homologs per query

# colabfold_search lives in the LocalColabFold conda env; put it on PATH.
COLABFOLD_DIR=$(require_folder "${COLABFOLD_DIR:-}" "COLABFOLD_DIR" "AlphaFold")
COLABFOLD_BIN="$COLABFOLD_DIR/colabfold-conda/bin"
if [[ ! -x "$COLABFOLD_BIN/colabfold_search" ]]; then
  echo "ERROR: colabfold_search not found at $COLABFOLD_BIN (is LocalColabFold installed?)" >&2
  exit 1
fi
export PATH="$COLABFOLD_BIN:$PATH"

mkdir -p "$JOB_QUEUE_DIR" "$RESULTS_DIR" "$TMP_DIR" "$GPU_TMP_DIR"

# Logging setup - no log file, just stdout (goes to slurm.out)
PID_FILE="$RESULTS_DIR/server_gpu.pid"

# Logging function (output to stdout only)
log() {
  echo "$(date '+%Y-%m-%d %H:%M:%S') - $*"
}

check_mmseqs_installation() {
    local mmseqs_bin="$MMSEQS2_DIR/bin/mmseqs"

    if [[ ! -f "$mmseqs_bin" ]]; then
        echo "$(date '+%Y-%m-%d %H:%M:%S') - MMseqs2 not found at $mmseqs_bin, downloading..."
        local parent_dir
        parent_dir="$(dirname "$MMSEQS2_DIR")"
        cd "$parent_dir" || { echo "$(date '+%Y-%m-%d %H:%M:%S') - ERROR: Cannot access $parent_dir"; exit 1; }

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
MMSEQS_SERVER_DIR="$MMSEQS2_SHARED_FOLDER"
mkdir -p "$MMSEQS_SERVER_DIR"
SERVER_TIMESTAMP_FILE="$MMSEQS_SERVER_DIR/GPU_SERVER"
SUBMITTING_FILE="$MMSEQS_SERVER_DIR/GPU_SUBMITTING"

# IMPORTANT: we do NOT write GPU_SERVER or clear GPU_SUBMITTING here. A client
# treats GPU_SERVER as "ready to serve" and GPU_SUBMITTING as "a server is on its
# way, don't submit another". This server is NOT ready until the DBs are loaded
# into GPU memory and the page cache is warmed (minutes after start, on top of
# SLURM queue time). Advertising readiness now would (a) let queries be picked up
# against not-yet-loaded DBs and (b) reopen the duplicate-submit window for the
# whole load period. So GPU_SUBMITTING stays held (covering queue + load) and
# GPU_SERVER is written only after warm-up (see below).

# Check and install MMseqs2 if needed
check_mmseqs_installation

# Optimized Memory Settings
export MMSEQS_MAX_MEMORY=${MMSEQS_MAX_MEMORY:-150G}  # Reduced from 200G
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}         # Reduced to 4 for GPU mode
export CUDA_VISIBLE_DEVICES=0

# Threads for colabfold_search. Its prefilter/search runs on the GPU, but the
# align / expandaln / result2msa steps run on CPU and dominate runtime, so give
# them the full SLURM CPU allocation (these GPU nodes are CPU-rich, 144 cores).
SEARCH_THREADS="${SLURM_CPUS_PER_TASK:-$(nproc)}"

# GPU Memory optimization
export CUDA_CACHE_MAXSIZE=2147483648  # 2GB CUDA cache
export CUDA_CACHE_DISABLE=0

# Check available memory
log "Memory settings: MMSEQS_MAX_MEMORY=$MMSEQS_MAX_MEMORY, OMP_NUM_THREADS=$OMP_NUM_THREADS"

# Strip a trailing NUL byte if present. mmseqs result2msa (which colabfold_search
# runs internally) can emit a final \0 that makes Boltz crash on the a3m.
strip_trailing_null() {
    local f=$1
    [[ -f "$f" ]] || return 0
    if [[ "$(tail -c 1 "$f")" == $'\0' ]]; then
        log "Removing trailing NUL byte from $f"
        head -c -1 "$f" > "${f}.clean" && mv "${f}.clean" "$f"
    fi
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

# Start a warm gpuserver for each DB colabfold_search hits (UniRef30 + envdb),
# so the search step can use --gpu-server 1 against both.
# With MMSEQS2_GPUS=2, pin the two gpuservers to separate devices (0 and 1) so
# their prefilters don't contend for one GPU; with 1, both share device 0.
NUM_GPUS="${MMSEQS2_GPUS:-1}"
if [[ "$NUM_GPUS" -ge 2 ]]; then ENVDB_DEVICE=1; else ENVDB_DEVICE=0; fi
log "Using $NUM_GPUS GPU(s): UniRef30 on device 0, envdb on device $ENVDB_DEVICE"

log "Starting MMseqs2 GPU server for $UNIREF_PATH"
CUDA_VISIBLE_DEVICES=0 $MMSEQS2_DIR/bin/mmseqs gpuserver "$UNIREF_PATH" \
  --max-seqs "$MAX_SEQS" \
  --db-load-mode 0 \
  --prefilter-mode 1 &
UNIREF_GPUSERVER_PID=$!
log "UniRef30 GPU server PID=$UNIREF_GPUSERVER_PID"

log "Starting MMseqs2 GPU server for $ENVDB_PATH"
CUDA_VISIBLE_DEVICES=$ENVDB_DEVICE $MMSEQS2_DIR/bin/mmseqs gpuserver "$ENVDB_PATH" \
  --max-seqs "$MAX_SEQS" \
  --db-load-mode 0 \
  --prefilter-mode 1 &
ENVDB_GPUSERVER_PID=$!
log "EnvDB GPU server PID=$ENVDB_GPUSERVER_PID"

# Wait for both databases to load into GPU memory
log "Waiting 3 minutes for databases to load into GPU memory..."
sleep 180
log "Wait complete, assuming GPU servers are ready"

# Warm the OS page cache with the CPU-side index files. colabfold_search's
# expandaln/align steps mmap these .idx files (--db-load-mode 2 = assume
# resident in RAM); without warming, the first read of each pages ~411GB off
# the /shares network filesystem PER QUERY, which dominated runtime (3-9 min).
# Reading them once into the page cache (held by this job's large RAM request)
# makes subsequent queries hit RAM, matching the public ColabFold server's
# ~20-30s/query. Needs a node with enough RAM to keep the indices resident
# (~512GB; the GPU nodes have ~1.5TB). The gpuserver covers the GPU side.
warm_page_cache() {
  log "Warming page cache for the DB indices (read once into RAM)..."
  local total=0
  for idx in "$UNIREF_PATH".idx "$ENVDB_PATH".idx; do
    if [[ -f "$idx" ]]; then
      local sz
      sz=$(du -h "$idx" 2>/dev/null | cut -f1)
      log "  caching $idx ($sz)"
      cat "$idx" > /dev/null 2>&1 || log "  WARNING: could not read $idx"
    fi
  done
  log "Page-cache warm-up complete; free memory now:"
  free -h | sed 's/^/    /'
}
warm_page_cache

# Now the DBs are loaded (GPU memory + warmed page cache): advertise the server
# as ready and release the submission lock. Write GPU_SERVER first, then clear
# GPU_SUBMITTING, so there is never a window where a client sees neither (which
# would let it submit a duplicate).
date '+%H:%M:%S' > "$SERVER_TIMESTAMP_FILE"
log "Created server timestamp file at $SERVER_TIMESTAMP_FILE (server is ready to serve)"
if [[ -f "$SUBMITTING_FILE" ]]; then
  log "Releasing submission lock"
  rm -f "$SUBMITTING_FILE"
fi
rmdir "${SUBMITTING_FILE}.lockdir" 2>/dev/null || true

# Cleanup on exit
cleanup() {
  log "MMseqs2 GPU server shutting down"
  for pid in "$UNIREF_GPUSERVER_PID" "$ENVDB_GPUSERVER_PID"; do
    log "Stopping GPU server PID=$pid"
    kill "$pid" 2>/dev/null || true
  done
  # Wait for graceful shutdown
  sleep 5
  # Force kill if still running
  for pid in "$UNIREF_GPUSERVER_PID" "$ENVDB_GPUSERVER_PID"; do
    kill -9 "$pid" 2>/dev/null || true
  done
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

# Idle auto-shutdown: exit once the queue has been empty for this long, so an
# idle GPU server doesn't squat the card and hurt SLURM priority. 0 disables.
IDLE_TIMEOUT="${MMSEQS2_IDLE_TIMEOUT:-600}"
LAST_JOB_TIME=$(date +%s)
log "Idle timeout: ${IDLE_TIMEOUT}s (0 = never)"

# Main processing loop
while true; do
  for job_meta in "$JOB_QUEUE_DIR"/*.job; do
    [[ -e "$job_meta" ]] || continue

    fname=$(basename "$job_meta")
    job_id="${fname%.job}"
    tmp="$TMP_DIR/$job_id"
    mkdir -p "$tmp"
    mkdir -p "$tmp/params"

    LAST_JOB_TIME=$(date +%s)   # reset idle clock on each job pickup
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

    # Prepare query FASTA
    cp "$fasta" "$tmp/query.fasta"

    # Set output prefix
    out_prefix="$RESULTS_DIR/$job_id"

    # Run colabfold_search against the warm gpuservers (UniRef30 + envdb).
    # colabfold_search owns the full search + MSA-build pipeline and writes one
    # unpacked .a3m per query into its output dir.
    log "Running colabfold_search for job $job_id"
    cf_out="$tmp/cf_out"
    mkdir -p "$cf_out"

    gpu_mem_before=$(nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits | head -1)
    log "GPU memory before search: ${gpu_mem_before}MB"

    if ! CUDA_VISIBLE_DEVICES=0 colabfold_search "$tmp/query.fasta" "$DB_DIR" "$cf_out" \
      --mmseqs "$MMSEQS2_DIR/bin/mmseqs" \
      --gpu 1 \
      --gpu-server 1 \
      --db-load-mode 2 \
      --threads "$SEARCH_THREADS"; then
      log "colabfold_search failed for job $job_id"
      gpu_mem_after=$(nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits | head -1)
      log "GPU memory after failed search: ${gpu_mem_after}MB"
      echo "FAILED: Search failed" > "$RESULTS_DIR/$job_id.status"
      continue
    fi

    gpu_mem_after=$(nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits | head -1)
    log "GPU memory after successful search: ${gpu_mem_after}MB"

    # colabfold_search writes one .a3m per query (named by query index/header).
    # Each job carries a single query, so take the sole produced a3m.
    produced_a3m=$(find "$cf_out" -maxdepth 1 -name "*.a3m" | head -1)
    if [[ -z "$produced_a3m" || ! -f "$produced_a3m" ]]; then
      log "ERROR: colabfold_search produced no .a3m for job $job_id"
      echo -e "FAILED\nNo MSA produced" > "$RESULTS_DIR/$job_id.status"
      continue
    fi
    mv "$produced_a3m" "$out_prefix.a3m"
    strip_trailing_null "$out_prefix.a3m"
    log "MSA written to $out_prefix.a3m"
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

  # Idle auto-shutdown: if the queue has been empty long enough, shut down so the
  # GPU is released. Only counts as idle when there are no pending .job files.
  if (( IDLE_TIMEOUT > 0 )); then
    pending=$(ls -1 "$JOB_QUEUE_DIR"/*.job 2>/dev/null | wc -l)
    if (( pending == 0 )) && (( CURRENT_TIME - LAST_JOB_TIME >= IDLE_TIMEOUT )); then
      log "No jobs for $((CURRENT_TIME - LAST_JOB_TIME))s (>= ${IDLE_TIMEOUT}s idle timeout); shutting down to release the GPU"
      cleanup
    fi
  fi

  sleep 5
done