#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Execution half of PPIDdG: coordinate-constrained FastRelax followed by a bound-minus-unbound ddG, per structure.

Follows the Chaperone block `ppi_ddg` (RosettaCommons/chaperone, MIT license). The relaxed score and the unbound score use plain beta_nov16, so the constraint weight never enters the reported energies.
"""

import argparse
import os
import sys

import yaml

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from biopipelines.biopipelines_io import (  # noqa: E402
    load_datastream,
    iterate_files,
    step_id_from_table_path,
)
from biopipelines.rosetta_scoring import (  # noqa: E402
    RELAX_INIT_OPTIONS,
    init_pyrosetta,
    mean_std,
    replicate_seed,
    seed_rosetta,
    unbound_ddg,
    collect_rows,
    write_outputs,
    step_exit_status,
)

COLUMNS = ["id", "ddg", "ddg_std", "n_replicates", "relax_score"]


def relax_and_ddg(pose, binder_chain, relax_rounds, constrain, constraint_weight):
    """One FastRelax cycle on a clone of the pose, then the no-repack ddG. Returns (ddg, total_score)."""
    import pyrosetta as pyr

    relax_fxn = pyr.create_score_function("beta_nov16")
    if constrain:
        relax_fxn.set_weight(pyr.rosetta.core.scoring.coordinate_constraint, constraint_weight)
    relax = pyr.rosetta.protocols.relax.FastRelax(relax_fxn, relax_rounds)
    if constrain:
        relax.constrain_relax_to_start_coords(True)
    relaxed = pose.clone()
    relax.apply(relaxed)

    plain = pyr.create_score_function("beta_nov16")
    total = float(plain(relaxed))
    return unbound_ddg(relaxed, plain, binder_chain), total


def chains_of(pose):
    info = pose.pdb_info()
    return {info.chain(i) for i in range(1, pose.total_residue() + 1)}


def score_structure(path, config, structure_id):
    """Mean and spread of the ddG over independent relax cycles of one structure, each seeded from config["seed"]."""
    import pyrosetta as pyr

    binder, target = config["binder_chain"], config["target_chain"]
    pose = pyr.pose_from_file(path)
    present = chains_of(pose)
    for chain in (binder, target):
        if chain not in present:
            raise ValueError(f"chain {chain!r} is not in the structure (chains: {', '.join(sorted(present))})")

    ddgs, scores = [], []
    for replicate in range(config["replicates"]):
        seed_rosetta(replicate_seed(config["seed"], structure_id, replicate))
        ddg, total = relax_and_ddg(pose, binder, config["relax_rounds"],
                                   config["constrain_to_start_coords"], config["constraint_weight"])
        ddgs.append(ddg)
        scores.append(total)
    mean, std = mean_std(ddgs)
    print(f"  ddg={mean:.2f} +/- {std:.2f} over {len(ddgs)} replicate(s)")
    return {"ddg": mean, "ddg_std": std, "n_replicates": len(ddgs), "relax_score": scores[-1]}


def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument("--config-yaml", required=True)
    p.add_argument("--structures-json", required=True)
    p.add_argument("--table-csv", required=True)
    p.add_argument("--local-missing-csv", required=True)
    args = p.parse_args(argv)

    with open(args.config_yaml) as fh:
        config = yaml.safe_load(fh) or {}

    ds = load_datastream(args.structures_json)
    step_id = step_id_from_table_path(args.local_missing_csv)
    init_pyrosetta(RELAX_INIT_OPTIONS)

    rows, missing, attempted, failures = collect_rows(
        iterate_files(ds), lambda sid, path: score_structure(path, config, sid), COLUMNS, step_id)
    write_outputs(args.table_csv, rows, COLUMNS, args.local_missing_csv, missing)

    print(f"ppi_ddg: {args.table_csv} ({len(rows)} rows, {failures} failed)")
    error = step_exit_status(attempted, failures)
    if error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
