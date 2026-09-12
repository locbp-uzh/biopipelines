#!/bin/bash
# LigandMPNN's upstream requirements.txt pins x86-only wheels: a full set of
# nvidia-*-cu12==<x86 version> plus torch==2.2.1 / triton==2.2.0 that have no
# aarch64 CUDA build. On Grace this fails with "No matching distribution found
# for nvidia-cublas-cu12==12.1.3.1". Strip the explicit CUDA pins (torch pulls
# the correct aarch64 CUDA libs itself) and relax torch/triton/numpy to
# aarch64-available versions. Runs on a login node (plain file edits).
#
#   ssh daint 'bash ~/biopipelines/environments/daint/fix_ligandmpnn_aarch64.sh'
#   ssh daint 'cd ~/biopipelines && sbatch environments/daint/warm.sh LigandMPNN'
set -euo pipefail

REPO_DIR="${SCRATCH}/bp/data/LigandMPNN"
REQ="${REPO_DIR}/requirements.txt"
[ -f "$REQ" ] || { echo "no $REQ (clone LigandMPNN first)" >&2; exit 1; }

[ -f "${REQ}.orig" ] || cp "$REQ" "${REQ}.orig"
{
  grep -vE '^nvidia-[a-z-]+-cu12==' "${REQ}.orig" \
    | sed -E 's/^torch==.*/torch/; s/^triton==.*/triton/; s/^numpy==.*/numpy<2/'
  # biopipelines' runtime positions helper (pipe_lmpnn_runtime_positions.py, used
  # when redesigned= is a table reference) imports pandas, which LigandMPNN's
  # upstream requirements never list. numpy<2 keeps pandas 2.2 compatible.
  echo 'pandas<2.3'
  # prody imports pkg_resources; setuptools>=81 removed it, and conda-forge
  # ships 83 in the env. Pin <81 so LigandMPNN's run.py -> data_utils -> prody works.
  echo 'setuptools<81'
} > "$REQ"
echo "=== requirements.txt diff (orig -> aarch64) ==="
diff "${REQ}.orig" "$REQ" || true

# Drop the half-built env (micromamba envs are self-contained dirs) so the warm's
# skip-check fails and it reinstalls against the sanitized requirements.
ENV_DIR="${SCRATCH}/bp/micromamba/envs/ligandmpnn_env"
if [ -d "$ENV_DIR" ]; then
    rm -rf "$ENV_DIR"
    echo "removed half-built env $ENV_DIR"
fi
echo "done. now: sbatch environments/daint/warm.sh LigandMPNN"
