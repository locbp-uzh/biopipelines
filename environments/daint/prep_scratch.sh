#!/bin/bash
# One-time: stage the BioPipelines persistent root on $SCRATCH and reuse the lab
# store's read-only weights by COPYING them in (the run container mounts only
# scratch + home, never the store, so nothing is ever written to the store).
# Runs on a login node -- plain file copies, no container, no GPU.
#
#   ssh daint 'bash ~/biopipelines/environments/daint/prep_scratch.sh'
set -euo pipefail

STORE=/capstor/store/cscs/uzh/uzh67
BP_ROOT="${SCRATCH}/bp"
CFG="${HOME}/biopipelines/config.container.yaml"

mkdir -p "${BP_ROOT}/cache" "${BP_ROOT}/data"

# Boltz2 model checkpoints + CCD mols (skip the 1.8 GB redundant mols.tar).
if [ ! -e "${BP_ROOT}/cache/Boltz/boltz2_conf.ckpt" ]; then
    echo "copying Boltz2 weights from store -> scratch"
    mkdir -p "${BP_ROOT}/cache/Boltz"
    rsync -a --info=progress2 --exclude='mols.tar' \
        "${STORE}/cache/Boltz/" "${BP_ROOT}/cache/Boltz/"
else
    echo "Boltz2 weights already on scratch"
fi

# LigandMPNN repo + model_params, so its install skips the git clone + download.
if [ ! -e "${BP_ROOT}/data/LigandMPNN/model_params" ]; then
    echo "copying LigandMPNN (repo + model_params) from store -> scratch"
    rsync -a "${STORE}/data/LigandMPNN/" "${BP_ROOT}/data/LigandMPNN/"
else
    echo "LigandMPNN already on scratch"
fi

# Point the container config root at scratch (the one line meant to be edited).
# Everything else (data, cache, outputs, micromamba env root) derives from it.
sed -i 's|^\([[:space:]]*root:[[:space:]]*\).*|\1"'"${BP_ROOT}"'"|' "${CFG}"
echo -n "config root now: "; grep -E '^[[:space:]]*root:' "${CFG}" | head -1

echo
echo "scratch root ready: ${BP_ROOT}"
du -sh "${BP_ROOT}"/cache/Boltz "${BP_ROOT}"/data/LigandMPNN 2>/dev/null || true
