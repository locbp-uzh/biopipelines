#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Download the Protenix checkpoint and the CCD files that `protenix pred` fetches on its first run.

Runs at install time, inside the Protenix env, so the first prediction needs no network.
It uses upstream's own URL table and path configuration, so every file lands where `protenix pred` looks for it, and downloads with the standard library: upstream's `runner.inference` imports the whole model stack, which JIT-compiles a CUDA extension and fails on a node without `CUDA_HOME`.
That is under PROTENIX_ROOT_DIR, or the home directory when it is unset.

Exit status 0 means everything is in place, 1 means a CCD file could not be fetched, and 3 means only the checkpoint could not be.
The checkpoint of `protenix-v2` is a case of its own: upstream calls those weights proprietary and its download URL may refuse the request, so the install carries on and says what to do.
"""

import argparse
import os
import shutil
import sys
import urllib.request

COMMON_FILES = ("ccd_components_file", "ccd_components_rdkit_mol_file",
                "pdb_cluster_file", "obsolete_release_data_csv")
CHECKPOINT_UNAVAILABLE = 3


def targets_for(model_name, data_configs, inference_configs):
    """{upstream download key: local path} for the CCD files and the checkpoint of one model."""
    targets = {name: data_configs[name] for name in COMMON_FILES}
    targets[model_name] = os.path.join(inference_configs["load_checkpoint_dir"], f"{model_name}.pt")
    return targets


def download_from_url(url, path, check_weight):
    """Stream `url` to `path`. An HTTP refusal raises URLError, an OSError; an empty checkpoint raises RuntimeError."""
    with urllib.request.urlopen(url, timeout=120) as response, open(path, "wb") as out:
        shutil.copyfileobj(response, out, 1 << 20)
    if check_weight and os.path.getsize(path) == 0:
        raise RuntimeError("the downloaded checkpoint is empty")


def fetch_missing(model_name, targets, urls, download):
    """Download whatever is absent. Returns why the checkpoint could not be fetched, or None; a CCD file that fails raises.

    `download(url, path, check_weight)` is `download_from_url` below, passed in so tests can replace it.
    """
    unavailable = None
    for name, path in targets.items():
        if os.path.exists(path):
            print(f"Already present: {path}")
            continue
        if name not in urls:
            raise RuntimeError(f"upstream lists no download for {name}")
        os.makedirs(os.path.dirname(path), exist_ok=True)
        # A partial file at the final path would be taken for the real one on the next install.
        partial = f"{path}.part"
        print(f"Downloading {name} from {urls[name]} to {path}")
        try:
            download(urls[name], partial, name == model_name)
            os.replace(partial, path)
        except (OSError, RuntimeError) as exc:
            if os.path.exists(partial):
                os.remove(partial)
            if name != model_name:
                raise RuntimeError(f"could not download {name}: {exc}") from exc
            unavailable = f"{urls[name]}: {exc}"
    return unavailable


def main():
    parser = argparse.ArgumentParser(description="Fetch the Protenix checkpoint and CCD files")
    parser.add_argument("--model-name", required=True)
    args = parser.parse_args()

    from configs.configs_data import data_configs
    from configs.configs_inference import inference_configs
    from protenix.web_service.dependency_url import URL

    targets = targets_for(args.model_name, data_configs, inference_configs)
    try:
        unavailable = fetch_missing(args.model_name, targets, URL, download_from_url)
    except RuntimeError as exc:
        print(f"ERROR: {exc}")
        sys.exit(1)

    if unavailable:
        print(f"WARNING: the checkpoint {args.model_name} could not be downloaded from {unavailable}")
        print(f"The CCD files are in place. To run {args.model_name}, put its checkpoint at {targets[args.model_name]}.")
        print("The upstream README calls the protenix-v2 weights proprietary, so obtain them from the rights holder; any other model_name does not need them.")
        sys.exit(CHECKPOINT_UNAVAILABLE)
    print("Protenix checkpoint and CCD files are in place")


if __name__ == "__main__":
    main()
