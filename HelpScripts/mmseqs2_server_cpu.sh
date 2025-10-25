#!/bin/bash
set -e
# ensure world-readable output files
umask 022

# timestamped logging
timestamp() { date '+%Y-%m-%d %H:%M:%S'; }
log() { echo "[$(timestamp)] $*"; }

# -----------------------------------------------------------------------------
# CPU-Only MMseqs2 MSA Server Script (Polling, with full logging)
# -----------------------------------------------------------------------------

USER=${USER:-$(whoami)}
JOB_QUEUE_DIR="/shares/locbp.chem.uzh/$USER/BioPipelines/MMseqs2Server/job_queue"
RESULTS_DIR="/shares/locbp.chem.uzh/$USER/BioPipelines/MMseqs2Server/results"
DB_DIR="/shares/locbp.chem.uzh/models/mmseqs2_databases/cpu"

TMP_DIR="/shares/locbp.chem.uzh/$USER/BioPipelines/MMseqs2Server/tmp"
CPU_TMP_DIR="/shares/locbp.chem.uzh/$USER/BioPipelines/MMseqs2Server/tmp/cpu"
UNIREF_DB="uniref30_2302_db"
DB_PATH="$DB_DIR/$UNIREF_DB"   # Use databases directly from shares
POLL_INTERVAL=10                     # seconds
MAX_SEQS=10000    # limit homologs per query

mkdir -p "$JOB_QUEUE_DIR" "$RESULTS_DIR" "$TMP_DIR" "$CPU_TMP_DIR"

# Check and install MMseqs2 if needed
check_mmseqs_installation() {
    local mmseqs_dir="/data/$USER/mmseqs"
    local mmseqs_bin="$mmseqs_dir/bin/mmseqs"

    if [[ ! -f "$mmseqs_bin" ]]; then
        echo "[$(timestamp)] MMseqs2 not found at $mmseqs_bin, downloading..."
        cd "/data/$USER" || { echo "[$(timestamp)] ERROR: Cannot access /data/$USER"; exit 1; }

        # Download MMseqs2
        wget https://mmseqs.com/latest/mmseqs-linux-gpu.tar.gz

        # Extract and cleanup
        tar xvfz mmseqs-linux-gpu.tar.gz
        rm mmseqs-linux-gpu.tar.gz

        # Verify installation
        if [[ -f "$mmseqs_bin" ]]; then
            echo "[$(timestamp)] MMseqs2 successfully installed at $mmseqs_bin"
        else
            echo "[$(timestamp)] ERROR: MMseqs2 installation failed"
            exit 1
        fi
    else
        echo "[$(timestamp)] MMseqs2 found at $mmseqs_bin"
    fi
}

# Logging setup
LOG_FILE="$RESULTS_DIR/server.log"
PID_FILE="$RESULTS_DIR/server_cpu.pid"
: > "$LOG_FILE"   # truncate on each start
echo $$ > "$PID_FILE"   # record server PID
export OMP_NUM_THREADS=$(nproc)

# Create timestamp file for server detection
MMSEQS_SERVER_DIR="/shares/locbp.chem.uzh/models/mmseqs2_server"
mkdir -p "$MMSEQS_SERVER_DIR"
SERVER_TIMESTAMP_FILE="$MMSEQS_SERVER_DIR/CPU_SERVER"
SUBMITTING_FILE="$MMSEQS_SERVER_DIR/CPU_SUBMITTING"

# Clean up submission timestamp if it exists (server is now starting)
if [[ -f "$SUBMITTING_FILE" ]]; then
  log "Cleaning up submission timestamp"
  rm -f "$SUBMITTING_FILE"
fi

date '+%H:%M:%S' > "$SERVER_TIMESTAMP_FILE"
log "Created server timestamp file at $SERVER_TIMESTAMP_FILE"

# Check and install MMseqs2 if needed
check_mmseqs_installation

log "MMseqs2 CPU server starting (PID: $$, using $OMP_NUM_THREADS threads)"
echo # blank line

# cleanup PID file on exit
cleanup() {
    log "MMseqs2 CPU server shutting down"
    rm -f "$PID_FILE"
    # Remove timestamp file on shutdown
    rm -f "$SERVER_TIMESTAMP_FILE"
    log "Removed server timestamp file"
    exit
}
trap cleanup EXIT INT TERM

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

handle_job() {
    local job_meta="$1"
    local fname=$(basename "$job_meta")
    local job_id="${fname%.job}"
    local tmp="$TMP_DIR/$job_id"
    mkdir -p "$tmp" "$tmp/params"
    log "Picked up job $job_id"
    if [[ -f "$job_meta" ]]; then
        mv "$job_meta" "$tmp/params/"
    else
        log "Warning: job file $job_meta not found; skipping."
        return
    fi

    # Store job file path for cleanup after completion
    local job_file_path="$tmp/params/$fname"

    # parse params
    output_format="a3m"; fasta=""
    while IFS='=' read -r key val; do
        case $key in
            fasta)         fasta="$val"         ;;
            output_format) output_format="$val" ;;
        esac
    done < "$tmp/params/$fname"

    # fallback FASTA
    [[ -n "$fasta" ]] || fasta="$JOB_QUEUE_DIR/${job_id}.fasta"
    if [[ ! -f "$fasta" ]]; then
        log "ERROR: FASTA not found at $fasta"
        echo -e "FAILED\nInput FASTA not found" > "$RESULTS_DIR/$job_id.status"
        return
    fi

    cp "$fasta" "$tmp/query.fasta"
    local query_db="$tmp/queryDB"
    log "Creating queryDB"
    /data/$USER/mmseqs/bin/mmseqs createdb "$tmp/query.fasta" "$query_db" 2>&1 | tee -a "$LOG_FILE"

    local result_db="$tmp/resultDB"
    local out_prefix="$RESULTS_DIR/$job_id"

    log "Running search for $job_id"
    if ! /data/$USER/mmseqs/bin/mmseqs search \
        "$query_db" "$DB_PATH" "$result_db" "$tmp/search_tmp" \
        --db-load-mode 2 \
        -s 7.5 \
        -a 1 \
        --alignment-mode 0 \
        --threads "$OMP_NUM_THREADS" \
        --max-seqs "$MAX_SEQS" \
      2>&1 | tee -a "$LOG_FILE"; then
        log "Search failed for $job_id"
        echo -e "FAILED\nSearch error" > "$RESULTS_DIR/$job_id.status"
        return
    fi

    log "Converting to A3M format"
    convert_to_a3m "$result_db" "$query_db" "$DB_PATH" "$out_prefix.a3m" "$tmp/query.fasta" "$tmp"

    case $output_format in
        a3m)
            echo -e "SUCCESS\noutput_file=$out_prefix.a3m" > "$RESULTS_DIR/$job_id.status"
            ;;
        csv)
            log "Converting A3M to CSV format"
            convert_a3m_to_csv "$out_prefix.a3m" "$out_prefix.csv"
            echo -e "SUCCESS\noutput_file=$out_prefix.csv" > "$RESULTS_DIR/$job_id.status"
            # Delete intermediate A3M file since client only needs CSV
            rm -f "$out_prefix.a3m"
            log "Deleted intermediate A3M file"
            ;;
        *)
            log "Unknown format $output_format"
            echo -e "FAILED\nUnknown format" > "$RESULTS_DIR/$job_id.status"
            ;;
    esac

    log "Job $job_id completed successfully"

    # Cleanup: delete job file and FASTA input (keep results and tmp for debugging)
    rm -f "$job_file_path"
    [[ -f "$fasta" ]] && rm -f "$fasta"
    log "Cleaned up job files for $job_id"
}

cleanup_old_files() {
    # Delete files older than 24 hours
    log "Running cleanup of files older than 24 hours..."

    # Clean old result files (.a3m, .csv, .status)
    find "$RESULTS_DIR" -type f \( -name "*.a3m" -o -name "*.csv" -o -name "*.status" \) -mtime +1 -delete 2>/dev/null || true

    # Clean old temporary directories
    find "$TMP_DIR" -mindepth 1 -maxdepth 1 -type d -mtime +1 -exec rm -rf {} \; 2>/dev/null || true

    # Clean old FASTA files in queue (orphaned submissions)
    find "$JOB_QUEUE_DIR" -type f -name "*.fasta" -mtime +1 -delete 2>/dev/null || true

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
                # Move .job file back to queue
                mv "$job_file" "$JOB_QUEUE_DIR/" && log "Requeued job file: ${job_id}.job"
                recovered_count=$((recovered_count + 1))
            else
                log "Warning: No .job file found for orphaned job $job_id"
            fi

            # Clean up the tmp directory
            rm -rf "$tmp_job_dir"
        fi
    done

    if [[ $recovered_count -gt 0 ]]; then
        log "Recovered $recovered_count orphaned jobs"
    fi
}

# Track cleanup time
LAST_CLEANUP=$(date +%s)
CLEANUP_INTERVAL=3600  # Run cleanup every hour

# Run cleanup and recovery at startup
log "Running startup cleanup to remove old files..."
cleanup_old_files
log "Running startup job recovery..."
recover_orphaned_jobs
LAST_CLEANUP=$(date +%s)  # Reset cleanup timer after startup cleanup

while true; do
  # (1) prioritize current user's jobs
  for job_meta in "$JOB_QUEUE_DIR"/*.job; do
    [[ -e "$job_meta" ]] || continue
    if [[ $(stat -c '%U' "$job_meta") == "$USER" ]]; then
      handle_job "$job_meta"
    fi
  done

  # (2) then process other users' jobs
  for job_meta in "$JOB_QUEUE_DIR"/*.job; do
    [[ -e "$job_meta" ]] || continue
    if [[ $(stat -c '%U' "$job_meta") != "$USER" ]]; then
      handle_job "$job_meta"
    fi
  done

  # (3) periodic cleanup and recovery
  CURRENT_TIME=$(date +%s)
  if (( CURRENT_TIME - LAST_CLEANUP >= CLEANUP_INTERVAL )); then
    cleanup_old_files
    recover_orphaned_jobs
    LAST_CLEANUP=$CURRENT_TIME
  fi

  sleep $POLL_INTERVAL
done