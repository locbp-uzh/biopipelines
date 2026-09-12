#!/bin/bash
# Build the BioPipelines base image on a Daint compute node and write its EDF.
#
#   sbatch -A uzh67 -p debug -t 00:30:00 environments/daint/build_base.sh 0.1.0
#
# Must run on a COMPUTE node: CSCS requires image builds there, and the build
# needs outbound network (apt, micromamba). podman storage lives in /dev/shm
# (wiped at job end), so build + enroot import happen in one allocation.
#
# Builds into $HOME (image) and writes the EDF to ~/.edf. Nothing is written to
# the shared lab store. Tool envs are created later on $SCRATCH by bootstrap.sh.
#SBATCH --account=uzh67
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=00:30:00
#SBATCH --job-name=bp-build-base
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
set -euo pipefail

VERSION="${1:?usage: sbatch build_base.sh <version>, e.g. 0.1.0}"
NAME="biopipelines-base"

# sbatch spools the batch script, so $0 is not in the repo; locate the build
# context from the submit dir instead (submit from the repo root).
HERE="${SLURM_SUBMIT_DIR:-$(pwd)}/environments/daint"
IMAGES="${HOME}/images"
EDF_DIR="${HOME}/.edf"
# Persistent BioPipelines root on scratch: micromamba envs, caches, outputs.
BP_ROOT="${SCRATCH}/bp"

# One-time podman config: overlay store in /dev/shm (home cannot hold it; layers
# are throwaway).
mkdir -p "${HOME}/.config/containers"
if [ ! -f "${HOME}/.config/containers/storage.conf" ]; then
    cat > "${HOME}/.config/containers/storage.conf" <<EOF
[storage]
driver = "overlay"
runroot = "/dev/shm/${USER}/runroot"
graphroot = "/dev/shm/${USER}/root"
EOF
fi
mkdir -p "/dev/shm/${USER}" "${IMAGES}" "${EDF_DIR}" "${BP_ROOT}"

cd "${HERE}"
podman build -t "${NAME}:${VERSION}" -f Containerfile .

# enroot can exit non-zero on a cleanup step with the squashfs written fine;
# check for the artifact instead of trusting the exit code.
SQSH="${IMAGES}/${NAME}-${VERSION}.sqsh"
enroot import -x mount -o "${SQSH}" "podman://${NAME}:${VERSION}" || true
[ -s "${SQSH}" ] || { echo "enroot produced no image at ${SQSH}" >&2; exit 1; }

# The EDF. Every line is load-bearing (see DAINT.md traps):
#   PATH        the CE replaces the env, dropping the image's ENV PATH; without
#               this /usr/local/bin (micromamba) is not found.
#   MAMBA_*     envs live on scratch, so the first job's installs and every
#               later job's runs share one env root that survives teardown.
#   NVIDIA_*    the host driver is injected only when these are set; without
#               them there is no libcuda.so and jobs run silently on CPU.
#   CC/CXX      JIT-compiling tools (Boltz-2/Triton) need a compiler on PATH.
cat > "${EDF_DIR}/bp-base.toml" <<EOF
image = "${SQSH}"
mounts = ["\${SCRATCH}:\${SCRATCH}", "\${HOME}:\${HOME}"]
workdir = "${HOME}/biopipelines"

[env]
PATH = "/usr/local/bin:/usr/local/sbin:/usr/sbin:/usr/bin:/sbin:/bin"
MAMBA_ROOT_PREFIX = "${BP_ROOT}/micromamba"
MAMBA_EXE = "/usr/local/bin/micromamba"
NVIDIA_VISIBLE_DEVICES = "all"
NVIDIA_DRIVER_CAPABILITIES = "compute,utility"
CC = "gcc"
CXX = "g++"
MPLBACKEND = "Agg"
EOF

echo
echo "built  ${SQSH} ($(du -h "${SQSH}" | cut -f1))"
echo "edf    ${EDF_DIR}/bp-base.toml"
echo "root   ${BP_ROOT}"
echo
echo "verify on a GPU node:"
echo "  srun -A uzh67 -p debug -t 00:10:00 -N1 -n1 --environment=bp-base \\"
echo "    bash -lc 'micromamba --version; gcc --version | head -1; nvidia-smi -L'"
