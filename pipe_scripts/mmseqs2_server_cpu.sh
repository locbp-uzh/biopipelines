#!/bin/bash
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# This software is freely available for use, modification, and redistribution.

set -e

# RAM-Optimized MMseqs2 CPU MSA Server Script
#
# CPU counterpart of mmseqs2_server_gpu.sh. It runs the SAME colabfold_search
# pipeline against the SAME ColabFold databases, but with no GPU: no gpuserver
# processes, no --gpu / --gpu-server flags. The bottleneck is RAM, not CPU —
# colabfold_search's prefilter/search/align/expandaln steps mmap the ~411GB DB
# indices, so this script warms the OS page cache once at startup (cat .idx) and
# relies on a node with enough RAM to keep the indices resident (request
# memory="512GB" via Resources()). With the indices hot in RAM, --db-load-mode 2
# (assume resident) gives query times comparable to the public ColabFold server.

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

TMP_DIR="$MMSEQS2_SHARED_FOLDER/tmp"
# colabfold_search searches both UniRef30 (--db1) and the environmental DB
# (--db3). Both are read directly from /shares; the page cache holds the
# indices resident after warm-up.
UNIREF_DB="uniref30_2302_db"
ENVDB="colabfold_envdb_202108_db"
UNIREF_PATH="$DB_DIR/$UNIREF_DB"
ENVDB_PATH="$DB_DIR/$ENVDB"
POLL_INTERVAL="${MMSEQS2_POLL_INTERVAL:-10}"   # seconds
MAX_SEQS="${MMSEQS2_MAX_SEQS:-10000}"          # limit homologs per query

# colabfold_search lives in the LocalColabFold conda env; put it on PATH.
COLABFOLD_DIR=$(require_folder "${COLABFOLD_DIR:-}" "COLABFOLD_DIR" "AlphaFold")
COLABFOLD_BIN="$COLABFOLD_DIR/colabfold-conda/bin"
if [[ ! -x "$COLABFOLD_BIN/colabfold_search" ]]; then
  echo "ERROR: colabfold_search not found at $COLABFOLD_BIN (is LocalColabFold installed?)" >&2
  exit 1
fi
export PATH="$COLABFOLD_BIN:$PATH"

# Use the mmseqs that ships with LocalColabFold (colabfold-conda/bin/mmseqs): a
# plain CPU build (no gpuserver, no NVIDIA-driver/glibc-2.29 requirement), so the
# CPU server is portable to CPU-only clusters and never depends on the GPU mmseqs
# build (which the GPU DB-build step installs under the MMseqs2 folder). It is
# v15; it reads the v16-built index fine (mmseqs is backward-compatible on the
# index format — verified).
MMSEQS_BIN="$COLABFOLD_BIN/mmseqs"
if [[ ! -x "$MMSEQS_BIN" ]]; then
  echo "ERROR: mmseqs not found at $MMSEQS_BIN (LocalColabFold should bundle it)" >&2
  exit 1
fi

mkdir -p "$JOB_QUEUE_DIR" "$RESULTS_DIR" "$TMP_DIR"

# Logging setup - no log file, just stdout (goes to slurm.out)
PID_FILE="$RESULTS_DIR/server_cpu.pid"

# Logging function (output to stdout only)
log() {
  echo "$(date '+%Y-%m-%d %H:%M:%S') - $*"
}

# vmtouch locks the DB indices into RAM so the shared-node kernel can't evict
# them (the way the public ColabFold/MMseqs2 servers keep DBs resident). It's a
# single-file C program with no deps; build it into the shared server folder if
# absent (kept off any tool-install dir so the CPU server is self-contained).
VMTOUCH_BIN="$MMSEQS2_SHARED_FOLDER/bin/vmtouch"
ensure_vmtouch() {
    if [[ -x "$VMTOUCH_BIN" ]]; then
        log "vmtouch found at $VMTOUCH_BIN"
        return 0
    fi
    log "vmtouch not found, building it..."
    local build_dir="$MMSEQS2_SHARED_FOLDER/.vmtouch_build"
    rm -rf "$build_dir"
    if git clone --depth 1 https://github.com/hoytech/vmtouch.git "$build_dir" >/dev/null 2>&1; then
        if make -C "$build_dir" >/dev/null 2>&1 && [[ -x "$build_dir/vmtouch" ]]; then
            mkdir -p "$(dirname "$VMTOUCH_BIN")"
            cp "$build_dir/vmtouch" "$VMTOUCH_BIN"
            rm -rf "$build_dir"
            log "vmtouch built at $VMTOUCH_BIN"
            return 0
        fi
    fi
    rm -rf "$build_dir"
    log "ERROR: could not build vmtouch (needed to lock DB indices into RAM)."
    log "Install git + a C compiler (cc/make) in the server env, or place a vmtouch binary at $VMTOUCH_BIN."
    exit 1
}

echo $$ > "$PID_FILE"   # record server PID

# Create timestamp file for server detection (per-user)
MMSEQS_SERVER_DIR="$MMSEQS2_SHARED_FOLDER"
mkdir -p "$MMSEQS_SERVER_DIR"
SERVER_TIMESTAMP_FILE="$MMSEQS_SERVER_DIR/CPU_SERVER"
SUBMITTING_FILE="$MMSEQS_SERVER_DIR/CPU_SUBMITTING"

# IMPORTANT: we do NOT write CPU_SERVER or clear CPU_SUBMITTING here. A client
# treats CPU_SERVER as "a server is ready to serve fast" and CPU_SUBMITTING as
# "a server is on its way, don't submit another". This server is NOT ready until
# the ~746GB index is loaded + mlocked (several minutes after start, on top of
# SLURM queue time). If we advertised readiness now:
#   - queries picked up here would run against a cold index (paging off /shares);
#   - and clearing CPU_SUBMITTING now would reopen the duplicate-submit window
#     for the whole load period.
# So CPU_SUBMITTING stays held (covering queue + install + load) and CPU_SERVER
# is written only after the index is resident (see mark_server_ready below).

log "Using mmseqs: $MMSEQS_BIN ($("$MMSEQS_BIN" version 2>/dev/null))"

# Memory / thread settings.
#
# CRITICAL for speed: mmseqs prefilter auto-splits the TARGET DB (--split 0,
# --split-mode 2) into enough chunks to fit MMSEQS_MAX_MEMORY, and processes the
# chunks SEQUENTIALLY. With a 450G cap against the ~526GB envdb index this split
# the search into ~144 serial passes — so the query ran ~70s using only ~1.6
# cores instead of one parallel sweep across all threads. But the indices are
# fully RAM-resident (vmtouch-mlocked) and read via --db-load-mode 2 (mmap, no
# copy), so there is NO reason to split for memory. Set MMSEQS_MAX_MEMORY to the
# node's real memory (well above the index size) so mmseqs picks a SINGLE split
# and uses all threads in one pass.
if [[ -n "${SLURM_MEM_PER_NODE:-}" ]]; then
  # SLURM_MEM_PER_NODE is in MB; hand mmseqs ~90% of the allocation.
  export MMSEQS_MAX_MEMORY="$(( SLURM_MEM_PER_NODE * 90 / 100 ))M"
else
  # Off-SLURM: use ~90% of total RAM (MemTotal in kB).
  _memtotal_kb=$(awk '/MemTotal/{print $2}' /proc/meminfo)
  export MMSEQS_MAX_MEMORY="$(( _memtotal_kb * 90 / 100 ))K"
fi
SEARCH_THREADS="${SLURM_CPUS_PER_TASK:-$(nproc)}"
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-$SEARCH_THREADS}
log "Memory settings: MMSEQS_MAX_MEMORY=$MMSEQS_MAX_MEMORY, OMP_NUM_THREADS=$OMP_NUM_THREADS, SEARCH_THREADS=$SEARCH_THREADS"

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

# Load the CPU-side index files into RAM and lock them there. colabfold_search's
# prefilter/search/expandaln/align steps mmap these .idx files
# (--db-load-mode 2 = assume resident in RAM). These `standard`-partition nodes
# are SHARED, so simply reading the indices into the page cache isn't enough:
# co-tenant jobs reclaim the clean file pages and the next query pages the ~746GB
# back off /shares at ~115 MB/s (mmseqs stuck in D-state), dominating runtime.
#
# Two-step load, because each tool is good at only half the job:
#  1. Parallel multi-stream `dd` populates the page cache fast. A single stream
#     (cat, or vmtouch's own serial -t) reads /shares at only ~115 MB/s → ~110min
#     for 746GB; staggered parallel reads over the 3-server NFS mount hit
#     ~1.6 GB/s → ~8 min. Override stream count with MMSEQS2_WARM_STREAMS (16).
#  2. `vmtouch -l -d` then mlocks the *already-resident* pages (no re-read, so it
#     is near-instant) and daemonises to hold the lock for the server's life —
#     the way the public ColabFold/MMseqs2 servers keep DBs pinned. mlock stops
#     co-tenant jobs from evicting them. Needs memlock=unlimited (it is on the
#     compute nodes) and a node large enough to hold the indices (~768GB+;
#     request memory="900GB" via Resources()).
#
# PIDFILE_DIR holds the per-chunk vmtouch daemon pidfiles, killed on shutdown.
WARM_STREAMS="${MMSEQS2_WARM_STREAMS:-16}"
PIDFILE_DIR="$MMSEQS2_SHARED_FOLDER/vmtouch_pids"
# Load + lock one .idx by splitting it into WARM_STREAMS byte ranges and running
# a `vmtouch -t -l -d` per range concurrently. Each chunk reads its slice off
# /shares AND mlocks it in one pass (vmtouch -p <range>), so there is no
# preload→lock eviction window and the reads parallelise across the 3-server NFS
# mount (~1.6 GB/s vs ~115 MB/s single stream). -t touch, -l mlock(2), -d
# daemonise (hold the lock), -P pidfile for clean shutdown.
lock_one_file() {
  local idx="$1" tag="$2"
  [[ -f "$idx" ]] || { log "  ERROR: index not found: $idx"; exit 1; }
  local sz_bytes total_gb chunk_gb i start_gb end_gb
  sz_bytes=$(stat -c %s "$idx")
  total_gb=$(( sz_bytes / 1073741824 + 1 ))
  chunk_gb=$(( total_gb / WARM_STREAMS + 1 ))
  log "  locking $idx (${total_gb}G) in $WARM_STREAMS ranges of ${chunk_gb}G"
  for (( i = 0; i < WARM_STREAMS; i++ )); do
    start_gb=$(( i * chunk_gb ))
    end_gb=$(( start_gb + chunk_gb ))
    [[ "$start_gb" -ge "$total_gb" ]] && break
    "$VMTOUCH_BIN" -t -l -d -P "$PIDFILE_DIR/${tag}_${i}.pid" \
      -p "${start_gb}G-${end_gb}G" "$idx" &
  done
  wait   # the parent vmtouch processes fork their -d daemons and exit immediately
}
# Wait until the two indices are (nearly) fully resident. The -d daemons load in
# the background, so without this the server would enter the job loop and serve
# queries against a cold index (I/O-bound). Poll vmtouch's resident fraction
# until >= 99% or a timeout.
wait_for_residency() {
  local target=99 timeout=1800 waited=0 pct
  while (( waited < timeout )); do
    # vmtouch -o kv prints one line of space-separated Key=Value pairs per run,
    # e.g. "... ResidentPages=N TotalPages=M ...". Sum across both files and
    # compute the combined resident percent.
    pct=$("$VMTOUCH_BIN" -o kv "$UNIREF_PATH".idx "$ENVDB_PATH".idx 2>/dev/null | awk '
      { for (i=1;i<=NF;i++) { split($i,kv,"=");
          if (kv[1]=="ResidentPages") res+=kv[2];
          if (kv[1]=="TotalPages") tot+=kv[2] } }
      END { if (tot>0) printf "%d", res*100/tot; else print 0 }')
    [[ -z "$pct" ]] && pct=0
    log "  index residency: ${pct}%"
    (( pct >= target )) && return 0
    sleep 15
    waited=$(( waited + 15 ))
  done
  log "WARNING: indices only ${pct}% resident after ${timeout}s; proceeding anyway"
}
lock_indices_in_ram() {
  ensure_vmtouch
  mkdir -p "$PIDFILE_DIR"
  rm -f "$PIDFILE_DIR"/*.pid 2>/dev/null || true
  local start end
  start=$(date +%s)
  log "Loading + locking DB indices into RAM ($WARM_STREAMS parallel vmtouch ranges/file)..."
  lock_one_file "$UNIREF_PATH".idx uniref
  lock_one_file "$ENVDB_PATH".idx envdb
  # Block until the daemons have actually loaded the pages, so the first query
  # hits a hot index rather than paging off /shares.
  wait_for_residency
  end=$(date +%s)
  log "Indices loaded + locked in $(( end - start ))s; memory now:"
  free -h | sed 's/^/    /'
  "$VMTOUCH_BIN" "$UNIREF_PATH".idx "$ENVDB_PATH".idx 2>/dev/null | sed 's/^/    /' || true
}
lock_indices_in_ram

# Now the index is resident: advertise the server as ready and release the
# submission lock. Order matters — write CPU_SERVER first, then clear
# CPU_SUBMITTING, so there is never a window where a client sees neither (which
# would let it submit a duplicate).
date '+%H:%M:%S' > "$SERVER_TIMESTAMP_FILE"
log "Created server timestamp file at $SERVER_TIMESTAMP_FILE (server is ready to serve)"
if [[ -f "$SUBMITTING_FILE" ]]; then
  log "Releasing submission lock"
  rm -f "$SUBMITTING_FILE"
fi
rmdir "${SUBMITTING_FILE}.lockdir" 2>/dev/null || true

# Cleanup on exit. Wired to EXIT (below) so a crash, `set -e` abort, or signal
# all release the mlocked indices — otherwise hundreds of GB stay pinned on the
# node. Clear the trap first so the closing `exit` can't re-enter it.
cleanup() {
  trap - EXIT SIGINT SIGTERM
  log "MMseqs2 CPU server shutting down"
  # Release the mlocked indices so the RAM is freed for other jobs: kill every
  # vmtouch daemon we started (one per byte-range, tracked via pidfiles).
  if [[ -d "$PIDFILE_DIR" ]]; then
    for pf in "$PIDFILE_DIR"/*.pid; do
      [[ -f "$pf" ]] || continue
      kill "$(cat "$pf")" 2>/dev/null || true
    done
    rm -f "$PIDFILE_DIR"/*.pid 2>/dev/null || true
  fi
  # Belt-and-braces: any stray vmtouch of ours against these indices.
  pkill -u "$USER" -f "vmtouch.*mmseqs2_databases" 2>/dev/null || true
  rm -f "$PID_FILE"
  # Remove timestamp file on shutdown
  rm -f "$SERVER_TIMESTAMP_FILE"
  log "Removed server timestamp file"
  exit 0
}
trap cleanup EXIT SIGINT SIGTERM

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
    # Recover jobs claimed by a dead server. Each claimed job lives at
    # $TMP_DIR/<job_id>/<job_id>.job. If that dir is older than 30 min and the job
    # still has no status file, the server processing it died — requeue the .job
    # so another server (or a restart) re-runs it.
    local max_age_seconds=1800  # 30 minutes
    local current_time=$(date +%s)
    local recovered_count=0
    local failed_count=0

    log "Scanning for orphaned jobs (older than 30 minutes)..."

    for job_dir in "$TMP_DIR"/*; do
        [[ -d "$job_dir" ]] || continue
        local job_id; job_id=$(basename "$job_dir")
        local job_file="$job_dir/${job_id}.job"
        [[ -f "$job_file" ]] || continue   # not a claimed-job dir

        local dir_mtime
        dir_mtime=$(stat -c %Y "$job_dir" 2>/dev/null || stat -f %m "$job_dir" 2>/dev/null || echo "0")
        local age=$((current_time - dir_mtime))
        [[ $age -gt $max_age_seconds ]] || continue

        local status_file="$RESULTS_DIR/$job_id.status"
        [[ -f "$status_file" ]] && { rm -rf "$job_dir"; continue; }   # already finished

        log "Recovering orphaned job: $job_id (age: $((age / 60)) minutes)"
        if mv "$job_file" "$JOB_QUEUE_DIR/" 2>/dev/null; then
            log "Requeued job file: ${job_id}.job"
            recovered_count=$((recovered_count + 1))
        else
            log "Warning: could not requeue ${job_id}.job - marking failed"
            echo -e "FAILED\nServer died during processing and job cannot be recovered" > "$status_file"
            failed_count=$((failed_count + 1))
        fi
        rm -rf "$job_dir"
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

# Run cleanup and recovery at startup
log "Running startup cleanup to remove old files..."
cleanup_old_files
log "Running startup job recovery..."
recover_orphaned_jobs
LAST_CLEANUP=$(date +%s)  # Reset cleanup timer after startup cleanup

# Loop-body commands can legitimately return non-zero (lost claim races, polling reads, a failed search for one job) and must not kill the daemon and leak the mlocked indices. Failures are handled per-job via explicit status files below.
set +e

log "Entering job processing loop"

# Idle auto-shutdown: exit once the queue has been empty for this long, so an
# idle server doesn't squat a large-RAM node and hurt SLURM priority. 0 disables.
IDLE_TIMEOUT="${MMSEQS2_IDLE_TIMEOUT:-600}"
LAST_JOB_TIME=$(date +%s)
log "Idle timeout: ${IDLE_TIMEOUT}s (0 = never)"

# One job = one multi-FASTA. mmseqs2 is many-against-many — it parallelises across
# QUERY sequences, not within one. So the client submits ALL its sequences as a
# single job carrying a multi-FASTA (one record per sequence, header = sequence
# id), and the server runs ONE colabfold_search over the whole set. The expensive
# ~57s envdb index sweep is then amortised across every sequence (measured
# ~11s/query for 20 vs 72s for 1). Per-sequence a3m's land in the job's output_dir.
#
# .job contract (key=value lines): job_id, fasta=<multi-FASTA>, output_format
# (a3m|csv), output_dir=<dir>. colabfold_search names each result "<header>.a3m";
# headers are the sequence ids, so result i is "<cf_out>/<seq_id>.a3m".

# Main processing loop
while true; do
  claimed=""
  for job_meta in "$JOB_QUEUE_DIR"/*.job; do
    [[ -e "$job_meta" ]] || continue
    fname=$(basename "$job_meta")
    job_id="${fname%.job}"
    job_dir="$TMP_DIR/$job_id"
    mkdir -p "$job_dir"
    # Atomic claim: rename() is atomic on NFS, so with multiple servers sharing
    # the queue exactly one wins each job; the others' mv fails and they skip it.
    jpath="$job_dir/$fname"
    if mv "$job_meta" "$jpath" 2>/dev/null; then claimed="$jpath"; break; fi
  done

  if [[ -n "$claimed" ]]; then
    LAST_JOB_TIME=$(date +%s)   # reset idle clock: we have work

    output_format="a3m"; fasta=""; output_dir=""
    while IFS='=' read -r key val; do
      case $key in
        fasta)         fasta="$val"         ;;
        output_format) output_format="$val" ;;
        output_dir)    output_dir="$val"    ;;
      esac
    done < "$claimed"

    if [[ -z "$fasta" || ! -f "$fasta" ]]; then
      log "ERROR: multi-FASTA not found for $job_id at '$fasta'"
      echo -e "FAILED\nInput FASTA not found" > "$RESULTS_DIR/$job_id.status"
      rm -rf "$job_dir"
    else
      nseq=$(grep -c '^>' "$fasta" 2>/dev/null || echo 0)
      log "Picked up job $job_id ($nseq sequences)"
      cf_out="$job_dir/cf_out"; mkdir -p "$cf_out"
      ext="a3m"; [[ "$output_format" == "csv" ]] && ext="csv"

      log "Running colabfold_search (CPU) over $nseq queries"
      if ! colabfold_search "$fasta" "$DB_DIR" "$cf_out" \
        --mmseqs "$MMSEQS_BIN" \
        --db-load-mode 2 \
        --threads "$SEARCH_THREADS"; then
        log "colabfold_search failed for job $job_id"
        echo -e "FAILED\nSearch failed" > "$RESULTS_DIR/$job_id.status"
      else
        # The client aliases queries to sequential integer headers (0,1,2,...),
        # so colabfold_search writes "<alias>.a3m" (digits are safe_filename-safe,
        # no character mangling, no collisions). Pass each result through under the
        # same alias; the client maps aliases back to real sequence ids.
        mkdir -p "$output_dir"
        produced=0; total=0
        while IFS= read -r alias; do
          total=$((total+1))
          src="$cf_out/$alias.a3m"
          if [[ ! -f "$src" ]]; then
            log "WARNING: no a3m produced for query alias $alias"
            continue
          fi
          strip_trailing_null "$src"
          if [[ "$ext" == "csv" ]]; then
            convert_a3m_to_csv "$src" "$output_dir/$alias.csv"
          else
            cp "$src" "$output_dir/$alias.a3m"
          fi
          produced=$((produced+1))
        done < <(grep '^>' "$fasta" | sed 's/^>//')

        if (( produced == total )); then
          echo "SUCCESS" > "$RESULTS_DIR/$job_id.status"
          log "Job $job_id completed: $produced/$total MSAs written to $output_dir"
        elif (( produced > 0 )); then
          echo -e "SUCCESS\npartial=$produced/$total" > "$RESULTS_DIR/$job_id.status"
          log "Job $job_id partial: $produced/$total MSAs written to $output_dir"
        else
          echo -e "FAILED\nNo MSA produced" > "$RESULTS_DIR/$job_id.status"
          log "Job $job_id FAILED: no MSAs produced"
        fi
      fi
      # Cleanup the job's working dir + queue input FASTA.
      [[ -f "$fasta" ]] && rm -f "$fasta"
      rm -rf "$job_dir"
    fi
  fi

  # Periodic cleanup and recovery
  CURRENT_TIME=$(date +%s)
  if (( CURRENT_TIME - LAST_CLEANUP >= CLEANUP_INTERVAL )); then
    cleanup_old_files
    recover_orphaned_jobs
    LAST_CLEANUP=$CURRENT_TIME
  fi

  # Idle auto-shutdown: if the queue has been empty long enough, shut down so the
  # large-RAM node is released. Only counts as idle when there are no pending .job files.
  if (( IDLE_TIMEOUT > 0 )); then
    pending=$(ls -1 "$JOB_QUEUE_DIR"/*.job 2>/dev/null | wc -l)
    if (( pending == 0 )) && (( CURRENT_TIME - LAST_JOB_TIME >= IDLE_TIMEOUT )); then
      log "No jobs for $((CURRENT_TIME - LAST_JOB_TIME))s (>= ${IDLE_TIMEOUT}s idle timeout); shutting down to release the node"
      cleanup
    fi
  fi

  sleep 5
done
