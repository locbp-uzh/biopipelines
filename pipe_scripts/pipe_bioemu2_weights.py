#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Fetch the BioEmu2 checkpoint into the Hugging Face cache, at install time.

bioemu2 downloads `checkpoints/<name>/checkpoint.ckpt` and `config.yaml` from the `microsoft/bioemu2` Hub repository on the first sampling run, which on a compute node races between jobs and needs the network.
Its loader accepts a cached file from any commit and never contacts the Hub when both files are present, so the files are fetched here with the same call and cache directory.
A checkpoint name is a published, never-overwritten directory upstream, so the name pins the weights.

Exit status 0 when both files are in the cache, 1 when they could not be fetched.

Usage:
    python pipe_bioemu2_weights.py [--checkpoint bioemu2-preview] [--cache-dir <dir>]
"""

import argparse
import sys
from typing import List, Optional

HF_REPO = "microsoft/bioemu2"
DEFAULT_CHECKPOINT = "bioemu2-preview"
CHECKPOINT_FILES = ("checkpoint.ckpt", "config.yaml")


def hub_filenames(checkpoint: str) -> List[str]:
    """The repository paths upstream's loader reads for one named checkpoint."""
    return [f"checkpoints/{checkpoint}/{name}" for name in CHECKPOINT_FILES]


def fetch(checkpoint: str, cache_dir: Optional[str], local_files_only: bool = False) -> List[str]:
    """Cached paths of the checkpoint's files, downloading the missing ones."""
    from huggingface_hub import hf_hub_download

    return [hf_hub_download(repo_id=HF_REPO, filename=filename, cache_dir=cache_dir,
                            local_files_only=local_files_only)
            for filename in hub_filenames(checkpoint)]


def main():
    parser = argparse.ArgumentParser(description="Fetch the BioEmu2 checkpoint into the Hugging Face cache")
    parser.add_argument("--checkpoint", default=DEFAULT_CHECKPOINT)
    parser.add_argument("--cache-dir", default=None, help="Hub cache directory (default: the Hugging Face default)")
    parser.add_argument("--check-only", action="store_true", help="Succeed only if the files are already cached")
    args = parser.parse_args()

    try:
        paths = fetch(args.checkpoint, args.cache_dir, local_files_only=args.check_only)
    except Exception as exc:
        print(f"ERROR: could not get {HF_REPO} checkpoint {args.checkpoint}: {exc}", file=sys.stderr)
        sys.exit(1)
    for path in paths:
        print(f"Checkpoint file ready: {path}")


if __name__ == "__main__":
    main()
