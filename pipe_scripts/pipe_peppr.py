#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Execution half of the Peppr tool: score each predicted structure against its reference with the peppr package.

Every id that enters gets a row in the table, with NaN for a metric peppr cannot compute. An id that cannot be evaluated at all (no matching reference, an unreadable file, atoms that cannot be matched, every metric failing) also gets a row of NaN, plus a kind="failure" row in the local missing file.
"""

import argparse
import math
import os
import sys
import warnings

import pandas as pd
import yaml

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from biopipelines.biopipelines_io import (  # noqa: E402
    load_datastream,
    iterate_files,
    step_id_from_table_path,
)
from biopipelines.id_map_utils import TIER_SIBLING, get_mapped_ids_with_tiers  # noqa: E402

# MonomerRMSD takes a threshold only for peppr's summary bins; it never changes the value it returns.
RMSD_SUMMARY_THRESHOLD = 2.0

# Each key builds the instance the peppr command line registers for that metric.
METRIC_FACTORIES = {
    "rmsd": lambda peppr: peppr.MonomerRMSD(RMSD_SUMMARY_THRESHOLD),
    "tm_score": lambda peppr: peppr.MonomerTMScore(),
    "lddt": lambda peppr: peppr.MonomerLDDTScore(),
    "lddt_pli": lambda peppr: peppr.LDDTPLIScore(),
    "fnat": lambda peppr: peppr.ContactFraction(),
    "irmsd": lambda peppr: peppr.InterfaceRMSD(),
    "lrmsd": lambda peppr: peppr.LigandRMSD(),
    "dockq": lambda peppr: peppr.DockQScore(),
}


def load_system(path):
    """Read a structure with the bond list peppr requires."""
    import biotite.structure.io.pdb as pdb
    import biotite.structure.io.pdbx as pdbx

    suffix = os.path.splitext(path)[1].lower()
    if suffix == ".pdb":
        return pdb.get_structure(pdb.PDBFile.read(path), model=1, include_bonds=True)
    if suffix in (".cif", ".mmcif"):
        return pdbx.get_structure(pdbx.CIFFile.read(path), model=1, include_bonds=True)
    raise ValueError(f"unsupported structure file {os.path.basename(path)!r}: expected .pdb, .cif or .mmcif")


def evaluate_pair(system_id, reference_path, pose_path, keys):
    """Return {key: float} for one pose against its reference, NaN where a metric could not be computed.

    Raises when nothing could be evaluated: a file that cannot be read, atoms that cannot be matched, or every requested metric failing.
    """
    import peppr

    metrics = [METRIC_FACTORIES[key](peppr) for key in keys]
    reference = load_system(reference_path)
    pose = load_system(pose_path)
    evaluator = peppr.Evaluator(metrics, tolerate_exceptions=True)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        evaluator.feed(system_id, reference, pose)

    for note in caught:
        print(f"  note: {note.message}", file=sys.stderr)
    unmatched = [str(w.message) for w in caught
                 if issubclass(w.category, peppr.MatchWarning)
                 or str(w.message).startswith("Failed to standardize")]
    if unmatched:
        raise ValueError(unmatched[0])
    failed = [str(w.message) for w in caught if issubclass(w.category, peppr.EvaluationWarning)]
    if len(failed) >= len(keys):
        raise ValueError(failed[0])

    results = evaluator.get_results()
    return {key: float(results[i][0][0]) for i, key in enumerate(keys)}


def pair_references(pose_ids, references, n_declared):
    """Return {pose id: reference id or None}.

    A reference stream declaring one id scores every pose, the convention of ConformationalChange and PoseChange. Several are matched to each pose by id, so a refold finds the design it derives from. A sibling id (a common ancestor, nothing more) is not a match: its structure is another molecule's, and scoring against it would look plausible and be wrong.
    """
    if n_declared == 1:
        only = next(iter(references), None)
        return {pose_id: only for pose_id in pose_ids}
    matches = get_mapped_ids_with_tiers(source_ids=list(pose_ids), target_ids=list(references), unique=True)
    return {pose_id: (None if tier == TIER_SIBLING else reference_id)
            for pose_id, (reference_id, tier) in matches.items()}


def score_all(poses, references, n_declared, keys, step_id, evaluate):
    """Score every pose; return (table rows, missing rows). poses and references map id to file path."""
    pairing = pair_references(list(poses), references, n_declared)
    rows, missing_rows = [], []
    for sid, pose_path in poses.items():
        row = {"id": sid, **{key: math.nan for key in keys}}
        try:
            reference_id = pairing.get(sid)
            if reference_id is None:
                raise ValueError("no reference with a matching id")
            print(f"{sid}: {os.path.basename(pose_path)} vs {reference_id}")
            row.update(evaluate(sid, references[reference_id], pose_path, keys))
        except Exception as e:
            print(f"WARNING: {sid} failed: {e}", file=sys.stderr)
            missing_rows.append({"id": sid, "removed_by": step_id,
                                 "kind": "failure", "cause": str(e)[:200]})
        rows.append(row)
    return rows, missing_rows


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--config-yaml", required=True)
    p.add_argument("--structures-json", required=True)
    p.add_argument("--references-json", required=True)
    p.add_argument("--peppr-csv", required=True)
    p.add_argument("--local-missing-csv", required=True)
    args = p.parse_args()

    with open(args.config_yaml) as fh:
        config = yaml.safe_load(fh) or {}
    keys = config["metrics"]
    unknown = [key for key in keys if key not in METRIC_FACTORIES]
    if unknown:
        print(f"ERROR: unknown metric(s) {unknown}; supported: {list(METRIC_FACTORIES)}", file=sys.stderr)
        sys.exit(1)

    poses_ds = load_datastream(args.structures_json)
    references_ds = load_datastream(args.references_json)
    poses = dict(iterate_files(poses_ds))
    references = dict(iterate_files(references_ds))
    step_id = step_id_from_table_path(args.local_missing_csv)

    rows, missing_rows = score_all(
        poses, references, len(references_ds.ids_expanded), keys, step_id, evaluate=evaluate_pair)

    for path in (args.peppr_csv, args.local_missing_csv):
        os.makedirs(os.path.dirname(path), exist_ok=True)
    pd.DataFrame(rows, columns=["id", *keys]).to_csv(args.peppr_csv, index=False)
    pd.DataFrame(missing_rows, columns=["id", "removed_by", "kind", "cause"]).to_csv(
        args.local_missing_csv, index=False)

    attempted, failures = len(rows), len(missing_rows)
    print(f"Peppr: {args.peppr_csv} ({attempted} rows, {failures} not evaluated)")

    if not attempted:
        print("ERROR: no input entities to process (upstream stream is empty or its files are absent)",
              file=sys.stderr)
        sys.exit(1)
    if failures == attempted:
        print(f"ERROR: all {attempted} attempted id(s) failed", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
