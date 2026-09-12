#!/bin/bash
# Provision the `daint` config variant (scheduler: slurm, env_manager: venv).
#
# WHY. The container variant (scheduler "none") runs everything inline in one
# allocation, so Parallel(pack=N) -- which emits srun job steps -- refuses at
# config time, and every shard has to be its own OTF pipeline in its own
# top-level output folder. That is the direct cause of the 80-folder sprawl.
# The daint variant has SLURM, so pack works, the campaign is ONE pipeline in
# ONE folder, and Folder() can nest steps inside it.
#
#   sbatch -p normal -t 04:00:00 environments/daint/provision_venv.sh base
#   sbatch -p normal -t 08:00:00 environments/daint/provision_venv.sh tools
#   sbatch -p normal -t 02:00:00 environments/daint/provision_venv.sh boltz2env
#   sbatch -p debug  -t 00:20:00 environments/daint/provision_venv.sh verify
#
# THE INTERPRETER TRAP, which dictates the whole approach. A venv records the
# absolute path of the interpreter that created it. Created inside an EDF
# container, its python points at /usr/local/bin/python3.11 -- a path that does
# not exist on a bare compute node, so the venv is unusable exactly where
# base-env tools run. The host python3 is 3.6.15 with no pip, so it cannot
# create it either. The way out is the one config.daint.yaml already describes:
# build a conda env under conda_env_root (on $SCRATCH, visible everywhere) and
# create the venv with ITS python. Both roots then resolve on any node.
#SBATCH --account=uzh67
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=04:00:00
#SBATCH --job-name=daint-provision
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
set -uo pipefail

STAGE="${1:-base}"
REPO="${SLURM_SUBMIT_DIR:-$HOME/biopipelines}"
cd "$REPO"

GROUP=/capstor/store/cscs/uzh/uzh67
MAMBA="$GROUP/micromamba/micromamba"
CONDA_ROOT="$SCRATCH/mamba-envs"      # config.daint.yaml: env_manager.conda_env_root
VENV_ROOT="$SCRATCH/venvs"            # config.daint.yaml: env_manager.venv_root
BASE_CONDA="$CONDA_ROOT/biopipelines"
BASE_VENV="$VENV_ROOT/biopipelines"

export MAMBA_ROOT_PREFIX="$CONDA_ROOT"
mkdir -p "$CONDA_ROOT" "$VENV_ROOT"

echo "=== daint provisioning: stage '$STAGE' ==="
echo "  conda root : $CONDA_ROOT"
echo "  venv root  : $VENV_ROOT"
echo

# --------------------------------------------------------------- base
if [ "$STAGE" = "base" ]; then
  if [ ! -x "$BASE_CONDA/bin/python" ]; then
    echo "[base] creating the conda interpreter (python 3.12)"
    # 3.12, not 3.11: venv mode propagates THIS interpreter to every tool env
    # it creates, and rc-foundry (RFdiffusion3) pins >=3.12,<3.13.
    "$MAMBA" create -y -p "$BASE_CONDA" -c conda-forge python=3.12 pip || exit 1
  else
    echo "[base] conda interpreter already present"
  fi

  if [ ! -d "$BASE_VENV" ]; then
    echo "[base] creating the venv from that interpreter"
    # --system-site-packages so the venv inherits the conda layer, matching what
    # _env_install_block does for conda-backed envs.
    "$BASE_CONDA/bin/python" -m venv --system-site-packages "$BASE_VENV" || exit 1
  else
    echo "[base] venv already present"
  fi

  # shellcheck disable=SC1090
  source "$BASE_VENV/bin/activate"
  python -m pip install --upgrade pip >/dev/null 2>&1
  echo "[base] installing environments/biopipelines.pip.daint.txt"
  pip install -r environments/biopipelines.pip.daint.txt || exit 1
  echo "[base] installing biopipelines itself (editable, no extras)"
  # NOT .[colab]: openbabel-wheel is already pinned in the pip file and the
  # colab extra pulls x86-only wheels.
  pip install -e . || exit 1

  echo "[base] verifying on THIS node"
  python -c "import biopipelines, pandas, numpy, rdkit; print('  imports OK')" || exit 1
  echo "[base] done"
fi

# --------------------------------------------------------------- tools
# bp-warm only calls Tool.install(); it does not bootstrap the base env, which
# is why 'base' runs first. Each tool's install script is emitted by the
# framework and, for tools with an edf: entry, is automatically wrapped in
# `srun --environment=<image>` so its venv binds to that image's python
# (pipeline.py:1019). Running the driver on the bare node -- not inside another
# srun -- keeps those from nesting.
if [ "$STAGE" = "tools" ]; then
  [ -x "$BASE_VENV/bin/python" ] || { echo "run the 'base' stage first" >&2; exit 1; }
  # shellcheck disable=SC1090
  source "$BASE_VENV/bin/activate"
  export BIOPIPELINES_CONFIG_VARIANT=daint
  export BIOPIPELINES_OTF=1
  export XDG_CACHE_HOME="$SCRATCH/.cache"; export HF_HOME="$SCRATCH/.cache/huggingface"
  mkdir -p "$XDG_CACHE_HOME" "$HF_HOME"

  TOOLS="${CR_TOOLS:-RFdiffusion3 LigandMPNN Boltz2 SASA}"
  echo "[tools] warming: $TOOLS"
  bp-warm $TOOLS --time 04:00:00 --cpus 32 --memory 64GB
  echo "[tools] done"
fi


# --------------------------------------------------------------- boltz2env
# Boltz2Env is WRAPPED, not rebuilt. Three rebuild attempts all ended in
# ResolutionImpossible -- pip walks boltz back to 0.0.0 and gives up -- on both
# python 3.12 and 3.11, with and without a conda layer. The `[cuda]` extra pulls
# x86-only nvidia-* wheels, and conda-forge has no linux-aarch64 gemmi 0.6.5
# (the exact pin boltz requires), so neither resolver can close the graph here.
#
# The container variant already has a WORKING Boltz2Env (python 3.11.15, boltz
# 2.2.1, gemmi 0.6.5) built inside the bp-base image, which has the compilers.
# So create the daint venv FROM that interpreter with --system-site-packages:
# the venv inherits the working install and supplies the bin/activate the venv
# env_manager needs. This is the same "conda env wrapped in a venv shim" pattern
# used for PyMOL/mkdssp.
#
# The dependency this creates is deliberate and worth knowing: the daint
# Boltz2Env now points into the container variant's micromamba root. Both live
# on $SCRATCH and are purged together, so they age as a unit -- but deleting the
# container env silently breaks the daint one.
if [ "$STAGE" = "boltz2env" ]; then
  SRC="$SCRATCH/bp/micromamba/envs/Boltz2Env"
  B_VENV="$VENV_ROOT/Boltz2Env"

  [ -x "$SRC/bin/python" ] || {
    echo "[boltz2env] the container-variant Boltz2Env is missing at $SRC;" >&2
    echo "            build it first with: sbatch environments/daint/warm.sh Boltz2" >&2
    exit 1; }

  echo "[boltz2env] wrapping $SRC ($("$SRC/bin/python" --version 2>&1))"
  [ -d "$B_VENV" ] || "$SRC/bin/python" -m venv --system-site-packages "$B_VENV" || exit 1

  # --system-site-packages shares site-packages, NOT bin/. So `import boltz`
  # succeeds while the `boltz` CLI the tool actually calls does not exist, and
  # the step dies with "boltz: command not found" -- precisely what the first
  # native smoke hit, because the verification below checked the import rather
  # than the entry point. Symlink the console scripts across.
  echo "[boltz2env] linking console scripts from the source env"
  for exe in "$SRC"/bin/*; do
    b=$(basename "$exe")
    case "$b" in python*|pip*|activate*|conda*|*.sh) continue ;; esac
    [ -e "$B_VENV/bin/$b" ] || ln -s "$exe" "$B_VENV/bin/$b" 2>/dev/null
  done

  echo "[boltz2env] verifying on a BARE node (where the tool actually runs)"
  "$B_VENV/bin/python" -c "
import boltz, gemmi
print('  import: boltz', boltz.__version__ if hasattr(boltz,'__version__') else 'ok', '| gemmi', gemmi.__version__)
" || { echo "[boltz2env] import verification FAILED" >&2; exit 1; }
  # Verify the CLI as well as the import -- the entry point is what the
  # generated bash actually calls.
  ( . "$B_VENV/bin/activate" && command -v boltz >/dev/null && echo "  cli: $(command -v boltz)" ) \
    || { echo "[boltz2env] boltz CLI NOT on PATH after activation" >&2; exit 1; }
  echo "[boltz2env] done"
fi

# --------------------------------------------------------------- verify
if [ "$STAGE" = "verify" ]; then
  echo "=== venvs present ==="
  ls "$VENV_ROOT" 2>/dev/null | sed 's/^/  /' || echo "  none"
  echo "=== conda envs present ==="
  ls "$CONDA_ROOT" 2>/dev/null | grep -v '^pkgs$' | sed 's/^/  /' || echo "  none"
  echo
  echo "=== base env usable on a BARE node (the case that matters) ==="
  # shellcheck disable=SC1090
  source "$BASE_VENV/bin/activate" 2>/dev/null && \
    python -c "
import biopipelines, sys
from biopipelines import StructureCluster, Boltz2, RFdiffusion3, LigandMPNN, SASA
print('  python', sys.version.split()[0])
print('  biopipelines + tools import OK')
" || echo "  BASE VENV NOT USABLE"
fi
