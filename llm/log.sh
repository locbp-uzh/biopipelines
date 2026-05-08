#!/usr/bin/env bash
# llm/log.sh — run a command, tee combined output to today's log.
#
# Usage: llm/log.sh <cmd> [args...]
#   e.g. llm/log.sh ssh cluster 'squeue -u $(whoami)'
#        llm/log.sh ssh cluster "cd ~/biopipelines && ./submit my_pipelines/foo.py"
#        llm/log.sh scp my_pipelines/foo.py cluster:~/biopipelines/my_pipelines/
#
# Writes a header, the (shell-quoted) command, its combined stdout/stderr, and
# the exit code to llm/logs/YYYY-MM-DD.log (relative to this script).

set -uo pipefail

if [ $# -lt 1 ]; then
  echo "Usage: $0 <cmd> [args...]" >&2
  exit 2
fi

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
LOG_DIR="$SCRIPT_DIR/logs"
mkdir -p "$LOG_DIR"
LOG_FILE="$LOG_DIR/$(date +%Y-%m-%d).log"

# Quote each arg so the log line is reproducible (copy-paste-runnable).
printf -v quoted '%q ' "$@"
{
  echo "=== $(date -Iseconds) :: ${quoted}==="
} >> "$LOG_FILE"

"$@" 2>&1 | tee -a "$LOG_FILE"
rc=${PIPESTATUS[0]}

{
  echo "=== exit=$rc ==="
  echo ""
} >> "$LOG_FILE"

exit "$rc"
