#!/bin/bash
# Publish a built image to the project's shared store and write its EDF, so
# other members run it without rebuilding.
#
#   environments/scripts/share_image.sh ProteinEnv 0.1.0
#
# A file copy, not compute: runs on a login node. Building and publishing are
# separate so a rebuild does not force a republish, and an already-published
# image can have its EDF regenerated without the source copy still existing.
set -euo pipefail

ENV_NAME="${1:?usage: share_image.sh <env-name> <version>}"
VERSION="${2:?missing version, e.g. 0.1.0}"

SHARE_ROOT="${BP_SHARE_ROOT:?set BP_SHARE_ROOT to the project store, e.g. /capstor/store/cscs/uzh/<project>}"
SRC="${BP_IMAGE_DIR:-${HOME}/images}/${ENV_NAME}-${VERSION}.sqsh"
DST_DIR="${SHARE_ROOT}/images"
EDF_DIR="${SHARE_ROOT}/edf"
DST="${DST_DIR}/${ENV_NAME}-${VERSION}.sqsh"

[ -d "${SHARE_ROOT}" ] || { echo "no shared store at ${SHARE_ROOT}" >&2; exit 1; }
mkdir -p "${DST_DIR}" "${EDF_DIR}"
command -v lfs >/dev/null 2>&1 && lfs setstripe -c 4 "${DST_DIR}" 2>/dev/null || true

if [ ! -s "${SRC}" ] && [ ! -s "${DST}" ]; then
    echo "no image at ${SRC} and none at ${DST} — build it first" >&2
    exit 1
fi

if [ ! -s "${SRC}" ]; then
    echo "using the store copy; nothing to publish"
elif [ -e "${DST}" ] && [ "${BP_SHARE_FORCE:-0}" != "1" ]; then
    echo "already published: ${DST}"
    echo "images are versioned and treated as immutable — bump the version rather"
    echo "than overwriting one others may be running against."
    echo "(BP_SHARE_FORCE=1 overrides)"
else
    echo "copying $(du -h "${SRC}" | cut -f1) to ${DST}"
    cp "${SRC}" "${DST}.partial"
    mv "${DST}.partial" "${DST}"
fi
chmod g+r "${DST}"

# Every key below is load-bearing, and two of them fail silently if omitted:
#
#   PATH      the Container Engine REPLACES the environment rather than
#             extending it, so the image's own ENV PATH is dropped. Without
#             this, /opt/<env> is not on PATH and `python` is the system one.
#   NVIDIA_*  enroot injects the host driver only when these are set. Without
#             them the container has no libcuda.so and a GPU job runs on the
#             CPU — slow, and nothing in the log says so.
#   workdir   OVERRIDES the cwd of whatever called srun, so it is deliberately
#             neutral: pointing it at one person's checkout would drop everyone
#             else somewhere they do not have. Scripts cd explicitly instead.
cat > "${EDF_DIR}/${ENV_NAME}.toml" <<EOF
# Shared image for the ${ENV_NAME} environment. Written by share_image.sh.
# Point a tool at it with  edf: { <Tool>: "${EDF_DIR}/${ENV_NAME}.toml" }
image = "${DST}"
mounts = ["\${SCRATCH}:\${SCRATCH}", "\${HOME}:\${HOME}", "${SHARE_ROOT}:${SHARE_ROOT}"]
workdir = "\${HOME}"

[env]
PATH = "/opt/${ENV_NAME}/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"
MPLBACKEND = "Agg"
NVIDIA_VISIBLE_DEVICES = "all"
NVIDIA_DRIVER_CAPABILITIES = "compute,utility"
EOF
chmod g+r "${EDF_DIR}/${ENV_NAME}.toml"

echo
echo "image  ${DST}"
echo "edf    ${EDF_DIR}/${ENV_NAME}.toml"
echo
echo "verify before relying on it — a GPU env must report CUDA, not just CPU:"
echo "  srun -A <account> -p debug -t 00:10:00 --environment=${EDF_DIR}/${ENV_NAME}.toml \\"
echo "      python -c 'import sys; print(sys.prefix)'"
