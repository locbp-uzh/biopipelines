#!/bin/bash
# Build a container image for one BioPipelines environment and import it as a
# SquashFS file the CSCS Container Engine can run.
#
#   sbatch environments/scripts/build_image.sh ProteinEnv 0.1.0
#
# Images are artifacts, not repo content: this script and the Containerfile are
# what is version-controlled. Publish the result with share_image.sh.
#
#SBATCH --job-name=bp-build-image
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=02:00:00
#SBATCH --output=bp-build-image-%j.out
#
# Must run on a COMPUTE node: CSCS requires it for podman, and the build needs
# outbound network for conda-forge. Pass -A/-p on the sbatch command line, since
# the account and partition are site-specific.
set -euo pipefail

ENV_NAME="${1:?usage: sbatch build_image.sh <env-name> <version>}"
VERSION="${2:?missing version, e.g. 0.1.0}"

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENVS_DIR="$(cd "${HERE}/.." && pwd)"
VARIANT="${BIOPIPELINES_CONFIG_VARIANT:-}"
IMAGES="${BP_IMAGE_DIR:-${HOME}/images}"

# Variant-specific spec wins over the shared one, matching how the framework
# resolves <env>.<variant>.yaml before <env>.yaml. `cluster` is the last resort
# because some envs (ProteinEnv) exist only in that form, and a conda spec is
# platform-independent — the solve, not the file, is what differs per platform.
SPEC=""
for candidate in "${ENV_NAME}.${VARIANT}.yaml" "${ENV_NAME}.yaml" "${ENV_NAME}.cluster.yaml"; do
    [ -n "${VARIANT}" ] || [ "${candidate}" != "${ENV_NAME}..yaml" ] || continue
    [ -f "${ENVS_DIR}/${candidate}" ] && { SPEC="${candidate}"; break; }
done
[ -n "${SPEC}" ] || { echo "no spec for ${ENV_NAME} in ${ENVS_DIR}" >&2; exit 1; }

# Same precedence for the optional pip layer. Staged as a fixed name because
# COPY fails on a missing source, so the Containerfile always has something.
PIP_SRC=""
for candidate in "${ENV_NAME}.pip.${VARIANT}.txt" "${ENV_NAME}.pip.txt"; do
    [ -n "${VARIANT}" ] || [ "${candidate}" = "${ENV_NAME}.pip.txt" ] || continue
    [ -f "${ENVS_DIR}/${candidate}" ] && { PIP_SRC="${ENVS_DIR}/${candidate}"; break; }
done

# podman's overlay store cannot live on Lustre — a build with a home graphroot
# dies on `lsetxattr ... operation not supported`. /dev/shm is tmpfs, which
# works, and is wiped when the job ends: that is why the import below has to
# happen in this same allocation rather than a later one.
mkdir -p "${HOME}/.config/containers"
if [ ! -f "${HOME}/.config/containers/storage.conf" ]; then
    cat > "${HOME}/.config/containers/storage.conf" <<EOF
[storage]
driver = "overlay"
runroot = "/dev/shm/${USER}/runroot"
graphroot = "/dev/shm/${USER}/root"
EOF
fi
mkdir -p "/dev/shm/${USER}" "${IMAGES}"

CTX="$(mktemp -d)"
trap 'rm -rf "${CTX}"' EXIT
cp "${ENVS_DIR}/${SPEC}" "${CTX}/${ENV_NAME}.yaml"
cp "${HERE}/Containerfile" "${CTX}/"
if [ -n "${PIP_SRC}" ]; then cp "${PIP_SRC}" "${CTX}/pip-requirements.txt"; else : > "${CTX}/pip-requirements.txt"; fi

echo "building ${ENV_NAME}:${VERSION} from ${SPEC}${PIP_SRC:+ + $(basename "${PIP_SRC}")}"
podman build --build-arg "ENV_NAME=${ENV_NAME}" \
             -t "${ENV_NAME}:${VERSION}" -f "${CTX}/Containerfile" "${CTX}"

# Stripe the target directory before writing: a 1-20 GB file on a single OST is
# served by one server no matter how many jobs read it. Only affects files
# created afterwards, so it must precede the import.
command -v lfs >/dev/null 2>&1 && lfs setstripe -c 4 "${IMAGES}" 2>/dev/null || true

# enroot exits non-zero on a cleanup step even when the squashfs is written
# correctly, so trust the artifact rather than the exit code.
SQSH="${IMAGES}/${ENV_NAME}-${VERSION}.sqsh"
enroot import -x mount -o "${SQSH}" "podman://${ENV_NAME}:${VERSION}" || true
[ -s "${SQSH}" ] || { echo "enroot produced no image at ${SQSH}" >&2; exit 1; }

echo
echo "built  ${SQSH} ($(du -h "${SQSH}" | cut -f1))"
echo "publish it for the project with:"
echo "  environments/scripts/share_image.sh ${ENV_NAME} ${VERSION}"
