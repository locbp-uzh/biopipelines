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

# Logging setup
LOG_FILE="$RESULTS_DIR/server.log"
PID_FILE="$RESULTS_DIR/server_cpu.pid"
: > "$LOG_FILE"   # truncate on each start
echo $$ > "$PID_FILE"   # record server PID
export OMP_NUM_THREADS=$(nproc)

log "MMseqs2 CPU server starting (PID: $$, using $OMP_NUM_THREADS threads)"
echo # blank line

# cleanup PID file on exit
cleanup() {
    log "MMseqs2 CPU server shutting down"
    rm -f "$PID_FILE"
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
            ;;
        *)
            log "Unknown format $output_format"
            echo -e "FAILED\nUnknown format" > "$RESULTS_DIR/$job_id.status"
            ;;
    esac

    log "Job $job_id completed successfully"
}

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

  sleep $POLL_INTERVAL
done