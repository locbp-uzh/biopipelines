#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Name the inputs of an RFdiffusion3 batch that produced nothing, and cut one of them into a JSON of its own.

`rfd3 design` stops at the first input it cannot parse, and the inputs after it never run. The wrapper therefore runs each input that has no output alone, so a crashed one costs only itself and the error its run prints is the one recorded for it.

Usage:
    python pipe_rfdiffusion3_pending.py --json inputs.json --raw-folder DIR
        prints the key of every input with no output in DIR, one per line
    python pipe_rfdiffusion3_pending.py --json inputs.json --key KEY --single-json OUT.json
        writes the entry of KEY alone to OUT.json
"""

import argparse
import json
import os
import re


def pending_keys(config, raw_folder):
    """Keys of `config` with no `<key>_<batch>_model_<model>.*` file in `raw_folder`, in the order of the JSON."""
    names = os.listdir(raw_folder) if os.path.isdir(raw_folder) else []
    return [key for key in config
            if not any(re.match(rf"^{re.escape(key)}_\d+_model_\d+\.", name) for name in names)]


def main():
    parser = argparse.ArgumentParser(description="Inputs of an RFdiffusion3 batch that produced nothing")
    parser.add_argument("--json", required=True, help="The inputs JSON handed to `rfd3 design`")
    parser.add_argument("--raw-folder", default=None, help="Folder foundry wrote its outputs to")
    parser.add_argument("--key", default=None, help="Cut this entry into --single-json")
    parser.add_argument("--single-json", default=None, help="Where to write the single entry")
    args = parser.parse_args()

    with open(args.json) as handle:
        config = json.load(handle)

    if args.key:
        if not args.single_json:
            parser.error("--key needs --single-json")
        with open(args.single_json, "w") as handle:
            json.dump({args.key: config[args.key]}, handle, indent=2)
        return

    if not args.raw_folder:
        parser.error("--raw-folder is required to list the pending keys")
    for key in pending_keys(config, args.raw_folder):
        print(key)


if __name__ == "__main__":
    main()
