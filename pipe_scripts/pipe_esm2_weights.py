#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Download ESM-2 checkpoints from the Hugging Face Hub at the commit that is current now, and pin them.

A Hub id follows `refs/main`, so the weights behind `facebook/esm2_t33_650M_UR50D` can change after a run is submitted. This resolves each id to a commit hash, downloads the tokenizer files and the weights of exactly that commit into `<dest>/<model>/`, and records `<repo>@<commit>` in `REVISION` beside them. Runs load that local folder, which is offline and cannot move. A folder that already has a `REVISION` is left alone.

Weights are downloaded as safetensors when the repository has them (8M to 650M) and otherwise as the sharded `pytorch_model*.bin` files and their index (3B). The TensorFlow files of the repositories are not fetched.

Runs in the ESM2 env (needs `huggingface_hub`).

Usage:
    python pipe_esm2_weights.py --dest <weights folder> --models esm2_t33_650M_UR50D [esm2_t6_8M_UR50D ...]
"""

import argparse
import os
import re
import sys
from typing import Iterable, List

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from biopipelines.esm2 import ESM2_MODELS, HUB_ORGANIZATION  # noqa: E402

REQUIRED_FILES = ("config.json", "vocab.txt")
OPTIONAL_FILES = ("special_tokens_map.json", "tokenizer_config.json")
SHARDED_BIN = re.compile(r"pytorch_model(-\d{5}-of-\d{5})?\.bin")


def select_files(filenames: Iterable[str]) -> List[str]:
    """The repository files to download: tokenizer and config, then safetensors if any, else the pytorch .bin files."""
    names = set(filenames)
    absent = [n for n in REQUIRED_FILES if n not in names]
    if absent:
        raise ValueError(f"the repository has no {absent}")
    files = [n for n in REQUIRED_FILES + OPTIONAL_FILES if n in names]
    safetensors = sorted(n for n in names if n.endswith(".safetensors"))
    if safetensors:
        weights = safetensors + [n for n in ("model.safetensors.index.json",) if n in names]
    else:
        weights = sorted(n for n in names if SHARDED_BIN.fullmatch(n))
        weights += [n for n in ("pytorch_model.bin.index.json",) if n in names]
    if not any(n.endswith((".safetensors", ".bin")) for n in weights):
        raise ValueError("the repository has neither safetensors nor pytorch_model*.bin weights")
    return files + weights


def download(model: str, dest: str) -> None:
    folder = os.path.join(dest, model)
    if os.path.isfile(os.path.join(folder, "REVISION")):
        print(f"{model}: already pinned in {folder} ({open(os.path.join(folder, 'REVISION')).read().strip()}), skipping")
        return
    from huggingface_hub import HfApi, snapshot_download

    repo = f"{HUB_ORGANIZATION}/{model}"
    api = HfApi()
    commit = api.model_info(repo).sha
    files = select_files(s.rfilename for s in api.model_info(repo, revision=commit).siblings)
    print(f"Downloading {repo}@{commit} to {folder}: {', '.join(files)}")
    os.makedirs(folder, exist_ok=True)
    snapshot_download(repo_id=repo, revision=commit, allow_patterns=files, local_dir=folder)
    missing = [n for n in files if not os.path.isfile(os.path.join(folder, n))]
    if missing:
        raise RuntimeError(f"{repo}@{commit} did not provide {missing}")
    with open(os.path.join(folder, "REVISION"), "w") as f:
        f.write(f"{repo}@{commit}\n")
    print(f"{model}: ready in {folder} ({repo}@{commit})")


def main():
    ap = argparse.ArgumentParser(description="Download and pin ESM-2 checkpoints")
    ap.add_argument("--dest", required=True, help="Folder that holds one subfolder per checkpoint")
    ap.add_argument("--models", nargs="+", required=True, choices=ESM2_MODELS)
    args = ap.parse_args()
    for model in args.models:
        download(model, args.dest)


if __name__ == "__main__":
    main()
