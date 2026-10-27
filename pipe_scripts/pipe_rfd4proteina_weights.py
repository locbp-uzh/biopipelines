#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Download the RFD4-Proteina weights bundle at the commit the installed package pins.

A Hub id follows `refs/main`, so the weights behind `nvidia/RFD4-Proteina` can change after a run is submitted. The installed package names the repository and commit it was released against (`hf_ckpt_path`, `hf_ckpt_revision` in its `configs/inference.yaml`); this resolves that to a commit hash, downloads exactly the files upstream's own loader reads (`config.json`, `model.safetensors`) plus the LICENSE into a plain directory, and records `<repo>@<commit>` in `REVISION` beside them. Runs then load that directory (`hf_ckpt_path=<dir>`), which is offline and cannot move.

Runs in the RFD4Proteina env (needs `rfproteina`, `huggingface_hub`, `yaml`).

Usage:
    python pipe_rfd4proteina_weights.py --dest <weights dir>
"""

import argparse
import importlib.util
import os
import sys
from typing import Optional, Tuple

import yaml

BUNDLE_FILES = ["config.json", "model.safetensors", "LICENSE"]


def pinned_release(config_text: str) -> Tuple[str, Optional[str]]:
    """(repository, revision) named by upstream's inference config; the revision is None when it pins none."""
    config = yaml.safe_load(config_text) or {}
    repo = config.get("hf_ckpt_path")
    if not repo or not isinstance(repo, str):
        raise ValueError("the installed inference.yaml names no hf_ckpt_path, so there is no bundle to download")
    revision = config.get("hf_ckpt_revision")
    return repo, (str(revision) if revision else None)


def installed_inference_config() -> str:
    spec = importlib.util.find_spec("rfproteina")
    if spec is None or not spec.submodule_search_locations:
        raise RuntimeError("rfproteina is not installed in this environment")
    root = list(spec.submodule_search_locations)[0]
    with open(os.path.join(root, "configs", "inference.yaml")) as f:
        return f.read()


def main():
    ap = argparse.ArgumentParser(description="Download the pinned RFD4-Proteina weights bundle")
    ap.add_argument("--dest", required=True, help="Directory to hold config.json, model.safetensors, LICENSE, REVISION")
    args = ap.parse_args()

    from huggingface_hub import HfApi, snapshot_download

    repo, revision = pinned_release(installed_inference_config())
    commit = HfApi().model_info(repo, revision=revision).sha
    print(f"Downloading {repo}@{commit} to {args.dest}")
    os.makedirs(args.dest, exist_ok=True)
    snapshot_download(repo_id=repo, revision=commit, allow_patterns=BUNDLE_FILES, local_dir=args.dest)

    missing = [name for name in ("config.json", "model.safetensors") if not os.path.isfile(os.path.join(args.dest, name))]
    if missing:
        print(f"ERROR: {repo}@{commit} did not provide {missing}", file=sys.stderr)
        sys.exit(1)
    with open(os.path.join(args.dest, "REVISION"), "w") as f:
        f.write(f"{repo}@{commit}\n")
    print(f"Weights bundle ready: {args.dest} ({repo}@{commit})")


if __name__ == "__main__":
    main()
