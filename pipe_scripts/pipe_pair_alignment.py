#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Execution half of USalign and TMalign: align every pair of the pair table in every mode and write one row per (pair, mode).

An alignment that gives nothing keeps its row with NaN values and adds a kind="failure" row to the local missing file. The pair handling itself lives in biopipelines/pair_alignment.py.
"""

import argparse
import os
import sys

import pandas as pd
import yaml

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from biopipelines.biopipelines_io import (  # noqa: E402
    load_datastream,
    iterate_files,
    step_id_from_table_path,
)
from biopipelines.pair_alignment import (  # noqa: E402
    DELTA_COLUMNS,
    MISSING_COLUMNS,
    align_pairs,
    alignments_frame,
    delta_rows,
    read_pairs,
)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--config-yaml", required=True)
    p.add_argument("--structures-json", required=True)
    p.add_argument("--alignments-csv", required=True)
    p.add_argument("--delta-csv", default="")
    p.add_argument("--local-missing-csv", required=True)
    p.add_argument("--container-prefix", default="")
    args = p.parse_args()

    with open(args.config_yaml) as f:
        config = yaml.safe_load(f)

    pairs = read_pairs(config["structure_1"], config["structure_2"])
    ds = load_datastream(args.structures_json)
    files = dict(iterate_files(ds))
    step_id = step_id_from_table_path(args.local_missing_csv)

    if not pairs:
        print("ERROR: the pair table holds no pair", file=sys.stderr)
        sys.exit(1)
    print(f"{config['binary']}: {len(pairs)} pair(s), {len(config['modes'])} alignment(s) each, "
          f"{len(files)} structure(s) available")

    rows, missing = align_pairs(pairs, files, ds.ids_expanded, config["modes"], config["binary"],
                                step_id, args.container_prefix)

    outputs = [args.alignments_csv, args.local_missing_csv] + ([args.delta_csv] if args.delta_csv else [])
    for path in outputs:
        os.makedirs(os.path.dirname(path), exist_ok=True)
    alignments_frame(rows).to_csv(args.alignments_csv, index=False)
    if args.delta_csv:
        pd.DataFrame(delta_rows(rows), columns=DELTA_COLUMNS).to_csv(args.delta_csv, index=False)
    pd.DataFrame(missing, columns=MISSING_COLUMNS).to_csv(args.local_missing_csv, index=False)

    attempted, failures = len(rows), len(missing)
    print(f"Alignments: {args.alignments_csv} ({attempted} rows, {failures} without values)")

    # Nothing attempted was caught above; every alignment failing is a failed step, not a table of NaN.
    if failures == attempted:
        print(f"ERROR: all {attempted} attempted alignment(s) failed", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
