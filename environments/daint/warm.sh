#!/bin/bash
# Warm BioPipelines tool envs onto the $SCRATCH persistent root, inside the
# bp-base container on a compute node (env creation + weight downloads need
# outbound network; login nodes lack it). No GPU needed to warm.
#
#   sbatch -p debug -t 00:30:00 environments/daint/warm.sh Boltz2
#   sbatch -p debug -t 00:30:00 environments/daint/warm.sh LigandMPNN RFdiffusion3
#
# The base `biopipelines` micromamba env (framework + light tools like Panda) is
# created on first call: conda-forge deps from environments/biopipelines.yaml
# (all aarch64-native), then `pip install -e .` (NOT .[colab] -- openbabel-wheel
# is x86-only and the conda env already provides openbabel).
#SBATCH --account=uzh67
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=00:30:00
#SBATCH --job-name=bp-warm
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
set -euo pipefail

REPO="${SLURM_SUBMIT_DIR:-$HOME/biopipelines}"
[ "$#" -ge 1 ] || { echo "usage: sbatch warm.sh <Tool> [Tool...]" >&2; exit 1; }
TOOLS="$*"

srun --environment=bp-base bash -c '
set -euo pipefail
cd "'"$REPO"'"
export BIOPIPELINES_CONFIG_VARIANT=container
export BIOPIPELINES_OTF=1
export BIOPIPELINES_LOCAL_OUTPUT=0
eval "$(micromamba shell hook --shell bash)"

if ! micromamba run -n biopipelines python -c "import biopipelines" >/dev/null 2>&1; then
    echo "=== creating base biopipelines env ==="
    micromamba create -y -f environments/biopipelines.yaml
    micromamba run -n biopipelines pip install -e .
fi

echo "=== bp-warm '"$TOOLS"' ==="
micromamba run -n biopipelines bp-warm --gpu any '"$TOOLS"'
'
