#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Prodigy-prot runner. For each complex PDB, instantiates
prodigy_prot.modules.prodigy.Prodigy and writes one affinity row per complex."""

import argparse
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from biopipelines.biopipelines_io import load_datastream, iterate_files, read_upstream_missing, MISSING_COLUMNS, step_id_from_table_path  # noqa: E402

from prodigy_prot.modules.parsers import parse_structure  # noqa: E402
from prodigy_prot.modules.prodigy import Prodigy as ProdigyEngine  # noqa: E402


AFF_COLS = ["id", "interface", "delta_g_kcal_mol", "kd_M",
            "n_intermol_contacts", "percent_charged_nis", "percent_apolar_nis"]


def parse_interface(s: str):
    """Convert 'A B' or 'A,B C' to the list-of-strings format prodigy expects."""
    return s.split()


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--structures-json", required=True)
    p.add_argument("--interface", required=True)
    p.add_argument("--temperature", type=float, default=25.0)
    p.add_argument("--affinity-csv", required=True)
    p.add_argument("--missing-csv", required=True)
    p.add_argument("--upstream-missing", nargs="*", default=None)
    args = p.parse_args()

    ds = load_datastream(args.structures_json)
    interface = parse_interface(args.interface)

    rows, failed = [], []
    step_id = step_id_from_table_path(args.missing_csv)
    for sid, pdb_path in iterate_files(ds):
        try:
            models, _, _ = parse_structure(pdb_path)
            engine = ProdigyEngine(models[0], name=sid, selection=interface, temp=args.temperature)
            engine.predict(temp=args.temperature, distance_cutoff=5.5, acc_threshold=0.05)
            rows.append({
                "id": sid,
                "interface": args.interface,
                "delta_g_kcal_mol": round(engine.ba_val, 3),
                "kd_M": engine.kd_val,
                "n_intermol_contacts": len(engine.ic_network),
                "percent_charged_nis": round(engine.nis_c, 2),
                "percent_apolar_nis": round(engine.nis_a, 2),
            })
            print(f"  {sid}: dG={engine.ba_val:.2f} kcal/mol  Kd={engine.kd_val:.2e} M")
        except Exception as e:
            print(f"WARNING: {sid} Prodigy failed: {e}", file=sys.stderr)
            failed.append({"id": sid, "removed_by": step_id, "kind": "failure", "cause": str(e)[:200]})

    all_missing = read_upstream_missing(args.upstream_missing) + failed

    for d in (args.affinity_csv, args.missing_csv):
        os.makedirs(os.path.dirname(d), exist_ok=True)
    pd.DataFrame(rows, columns=AFF_COLS).to_csv(args.affinity_csv, index=False)
    pd.DataFrame(all_missing, columns=MISSING_COLUMNS).to_csv(args.missing_csv, index=False)
    print(f"Affinity: {args.affinity_csv} ({len(rows)} rows)")
    print(f"Missing: {args.missing_csv} ({len(all_missing)} rows)")

    if failed:
        print(f"Failed: {len(failed)}/{len(failed)+len(rows)}", file=sys.stderr)
    if not rows:
        sys.exit(1)


if __name__ == "__main__":
    main()
