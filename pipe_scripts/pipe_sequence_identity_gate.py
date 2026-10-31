#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""SequenceIdentityGate runner. Aligns every design globally (Needleman-Wunsch, BLOSUM62) against every reference, and writes the maximum identity, the nearest reference and the gate verdict per design."""

import argparse
import math
import os
import sys
from typing import Dict, List, Tuple

import pandas as pd
import yaml
from Bio import Align
from Bio.Align import substitution_matrices

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from biopipelines.biopipelines_io import (  # noqa: E402
    load_datastream,
    iterate_values,
    step_id_from_table_path,
)


GATE_COLUMNS = ["id", "max_identity", "nearest_reference", "passed"]

_MATRIX = substitution_matrices.load("BLOSUM62")
_ALPHABET = set(_MATRIX.alphabet)


def make_aligner() -> Align.PairwiseAligner:
    """Global aligner with the gap scores of EMBOSS Needle (BLOSUM62, gap open -10, gap extend -0.5) and end gaps penalized, which Needle does not do by default."""
    return Align.PairwiseAligner(mode="global", substitution_matrix=_MATRIX,
                                 open_gap_score=-10.0, extend_gap_score=-0.5)


def normalize(sequence) -> str:
    """Upper-case sequence without whitespace, gap characters or stop; raises on letters BLOSUM62 lacks."""
    if sequence is None or (isinstance(sequence, float) and math.isnan(sequence)):
        raise ValueError("empty sequence")
    cleaned = "".join(str(sequence).split()).upper().replace("-", "").replace("*", "")
    if not cleaned:
        raise ValueError("empty sequence")
    unknown = sorted(set(cleaned) - _ALPHABET)
    if unknown:
        raise ValueError(f"letters outside the BLOSUM62 alphabet: {''.join(unknown)}")
    return cleaned


def global_identity(aligner: Align.PairwiseAligner, a: str, b: str) -> float:
    """Identical aligned positions over the alignment length, gaps included."""
    alignment = aligner.align(a, b)[0]
    return alignment.counts().identities / alignment.shape[1]


def load_references(references_json: str) -> List[Tuple[str, str]]:
    """Unique (id, normalized sequence) pairs; the first id of a repeated sequence is kept."""
    ds = load_datastream(references_json)
    seen: Dict[str, str] = {}
    for rid, values in iterate_values(ds, columns=["sequence"]):
        try:
            sequence = normalize(values["sequence"])
        except ValueError as e:
            raise ValueError(f"reference {rid!r}: {e}")
        seen.setdefault(sequence, rid)
    return [(rid, sequence) for sequence, rid in seen.items()]


def nearest_reference(aligner: Align.PairwiseAligner, sequence: str,
                      references: List[Tuple[str, str]]) -> Tuple[float, str]:
    """(maximum identity, id of the first reference that reaches it)."""
    best, best_id = -1.0, ""
    for rid, ref in references:
        identity = global_identity(aligner, sequence, ref)
        if identity > best:
            best, best_id = identity, rid
    return best, best_id


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--config-yaml", required=True)
    p.add_argument("--sequences-json", required=True)
    p.add_argument("--references-json", required=True)
    p.add_argument("--identity-gate-csv", required=True)
    p.add_argument("--local-missing-csv", required=True)
    args = p.parse_args()

    with open(args.config_yaml) as f:
        threshold = float(yaml.safe_load(f)["threshold"])

    try:
        references = load_references(args.references_json)
    except ValueError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)
    if not references:
        print("ERROR: the references stream holds no sequences", file=sys.stderr)
        sys.exit(1)

    ds = load_datastream(args.sequences_json)
    step_id = step_id_from_table_path(args.local_missing_csv)
    aligner = make_aligner()

    rows, local_missing = [], []
    attempted = failures = 0
    for sid, values in iterate_values(ds, columns=["sequence"]):
        attempted += 1
        try:
            identity, nearest = nearest_reference(aligner, normalize(values["sequence"]), references)
            rows.append({"id": sid, "max_identity": identity, "nearest_reference": nearest,
                         "passed": identity < threshold})
        except Exception as e:
            print(f"WARNING: {sid} failed: {e}", file=sys.stderr)
            failures += 1
            rows.append({"id": sid, "max_identity": None, "nearest_reference": None, "passed": None})
            local_missing.append({"id": sid, "removed_by": step_id, "kind": "failure",
                                  "cause": str(e)[:200]})

    for path in (args.identity_gate_csv, args.local_missing_csv):
        os.makedirs(os.path.dirname(path), exist_ok=True)
    table = pd.DataFrame(rows, columns=GATE_COLUMNS)
    table["max_identity"] = table["max_identity"].astype("float64")
    table["passed"] = table["passed"].astype("boolean")
    table.to_csv(args.identity_gate_csv, index=False)
    pd.DataFrame(local_missing, columns=["id", "removed_by", "kind", "cause"]).to_csv(
        args.local_missing_csv, index=False)
    print(f"Identity gate: {args.identity_gate_csv} ({len(rows)} rows, {failures} NaN) "
          f"against {len(references)} unique reference(s)")

    if not attempted:
        print("ERROR: no input entities to process (upstream stream is empty or its files are absent)",
              file=sys.stderr)
        sys.exit(1)
    if failures == attempted:
        print(f"ERROR: all {attempted} attempted id(s) failed", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
