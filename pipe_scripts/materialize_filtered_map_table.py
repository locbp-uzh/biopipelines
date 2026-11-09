#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Write a map_table projected to a DataStream JSON's expanded ids at runtime.

Usage:
    python materialize_filtered_map_table.py <ds_json> <out_csv> [--columns c1,c2] [--require c1,c2]

A filtered stream restricts its ids but keeps the full upstream map_table, so a
consumer reading the raw table would reprocess filtered-out rows. This runs at
execution time (when the upstream map_table exists) and emits only the selected
ids in stream order.
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from biopipelines.biopipelines_io import (load_datastream, read_upstream_missing,
                                          write_filtered_map_table)


def main():
    parser = argparse.ArgumentParser(description="Materialize a filtered map_table from a DataStream JSON")
    parser.add_argument("ds_json")
    parser.add_argument("out_csv")
    parser.add_argument("--columns", default=None, help="comma-separated column subset (id always kept)")
    parser.add_argument("--require", default=None, help="comma-separated columns that must exist")
    # Every row excuses its id whatever its kind: these are upstream manifests, so a failure there already failed its own step.
    parser.add_argument("--excused", nargs="*", default=None,
                        help="upstream missing.csv path(s); ids listed there may be absent")
    args = parser.parse_args()

    ds = load_datastream(args.ds_json)
    ds._runtime_mode = True
    columns = args.columns.split(",") if args.columns else None
    required = args.require.split(",") if args.require else None

    excused = [r.get("id") for r in read_upstream_missing(args.excused)]

    write_filtered_map_table(ds, args.out_csv, columns=columns, required_columns=required,
                             excused_ids=excused)
    print(f"Wrote filtered map_table ({len(ds.ids_expanded)} id(s)) to {args.out_csv}")


if __name__ == "__main__":
    main()
