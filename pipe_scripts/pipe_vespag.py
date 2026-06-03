#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""VespaG runtime helper.

Reads the input sequences CSV (columns: id, sequence), writes a FASTA, optionally
builds a VespaG mutation file from an explicit mutation request, runs

    vespag predict -i seqs.fasta -o OUTDIR --single-csv [--mutation-file muts.csv]

and parses the resulting ``OUTDIR/vespag_scores_all.csv`` (columns:
Protein, Mutation, VespaG; Mutation is 1-indexed wildtype-position-mutant, e.g.
'A42G'). Emits BioPipelines-style outputs:

  - fitness.csv : id | sequences.id | position | wildtype | mutation | fitness
                  (sorted by fitness descending across the whole table)
  - missing.csv : per-sequence failures (id | removed_by | kind | cause)

VespaG generates ESM-2 3B embeddings internally on first use (heavy download,
GPU strongly recommended). Runs under the `vespag` conda env.
"""

import argparse
import os
import re
import subprocess
import sys
import tempfile
import traceback
from typing import Dict, List, Optional, Set, Tuple

import pandas as pd

# Make biopipelines importable from a SLURM/Colab node where the package is
# not pip-installed (matches the other pipe scripts).
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from biopipelines.biopipelines_io import (  # noqa: E402
    load_datastream, load_table, lookup_table_value, step_id_from_table_path,
)

OUTPUT_COLUMNS = ["id", "sequences.id", "position", "wildtype", "mutation", "fitness"]

MUT_RE = re.compile(r"^([A-Za-z])(\d+)([A-Za-z])$")


def parse_mutation_tokens(spec: str) -> Set[Tuple[str, int, str]]:
    """Parse a '+'-joined list of wildtype-position-mutant tokens (1-indexed).

    Each token is like ``A42G``. Commas also accepted as separators. Returns a
    set of (wildtype, position, mutant) tuples. Raises ValueError on a bad token.
    """
    out: Set[Tuple[str, int, str]] = set()
    for token in spec.replace(",", "+").split("+"):
        token = token.strip()
        if not token:
            continue
        m = MUT_RE.match(token)
        if not m:
            raise ValueError(f"Bad mutation token: {token!r} (expected e.g. 'A42G')")
        wt, pos, mut = m.group(1).upper(), int(m.group(2)), m.group(3).upper()
        out.add((wt, pos, mut))
    return out


def resolve_mutations_arg(arg: str, seq_id: str) -> Optional[Set[Tuple[str, int, str]]]:
    """Resolve --mutations for a sequence. '-'/empty -> None (saturation);
    a literal token string or a TABLE_REFERENCE -> parsed token set."""
    if not arg or arg == "-":
        return None
    if arg.startswith("TABLE_REFERENCE:"):
        table, column = load_table(arg)
        cell = lookup_table_value(table, seq_id, column)
        if cell is None or (isinstance(cell, float) and pd.isna(cell)):
            return None
        return parse_mutation_tokens(str(cell))
    return parse_mutation_tokens(arg)


def read_sequences_csv(csv_path: str) -> List[Tuple[str, str]]:
    """Read (id, sequence) pairs from the content-bearing sequences CSV."""
    df = pd.read_csv(csv_path)
    if "id" not in df.columns or "sequence" not in df.columns:
        raise ValueError(
            f"sequences CSV must have 'id' and 'sequence' columns. Found: {list(df.columns)}"
        )
    return list(zip(df["id"].astype(str), df["sequence"].astype(str)))


def main():
    ap = argparse.ArgumentParser(description="VespaG fitness runtime helper")
    ap.add_argument("--sequences-json", required=True)
    ap.add_argument("--sequences-csv", required=True)
    ap.add_argument("--mutations", default="-",
                    help="'-' for saturation, a '+'-joined token list, or a TABLE_REFERENCE")
    ap.add_argument("--repo-dir", required=True,
                    help="VespaG clone dir; `vespag predict` runs with this as cwd so "
                         "load_model finds model_weights/v2/esm2.pt (resolved via Path.cwd()).")
    ap.add_argument("--fitness-csv", required=True)
    ap.add_argument("--missing-csv", required=True)
    args = ap.parse_args()

    # Loaded for validation/ID provenance; the CSV is the sequence source.
    _ = load_datastream(args.sequences_json)
    sequences = read_sequences_csv(args.sequences_csv)

    missing: List[dict] = []
    step_id = step_id_from_table_path(args.missing_csv)

    # Resolve per-sequence requested mutations up front so we can (a) build a
    # combined mutation file and (b) report any that VespaG drops.
    requested_by_seq: Dict[str, Optional[Set[Tuple[str, int, str]]]] = {}
    for seq_id, _seq in sequences:
        try:
            requested_by_seq[seq_id] = resolve_mutations_arg(args.mutations, seq_id)
        except Exception as exc:  # noqa: BLE001
            requested_by_seq[seq_id] = None
            missing.append({"id": seq_id, "removed_by": step_id, "kind": "failure",
                            "cause": f"mutation parse error: {type(exc).__name__}: {exc}"})
            traceback.print_exc()

    any_explicit = any(v is not None for v in requested_by_seq.values())

    rows: List[dict] = []

    with tempfile.TemporaryDirectory(prefix="vespag_") as work:
        fasta_path = os.path.join(work, "input.fasta")
        with open(fasta_path, "w") as fh:
            for seq_id, seq in sequences:
                fh.write(f">{seq_id}\n{seq}\n")

        cmd = [sys.executable, "-m", "vespag", "predict",
               "-i", fasta_path, "-o", work, "--single-csv"]

        # Explicit-mutation mode: write a VespaG mutation file
        # (rows of protein_id,mutation_id). VespaG defaults to 1-indexed
        # positions, matching our token convention.
        if any_explicit:
            mut_file = os.path.join(work, "mutations.csv")
            with open(mut_file, "w") as fh:
                fh.write("protein_id,mutation_id\n")
                for seq_id, muts in requested_by_seq.items():
                    if not muts:
                        continue
                    for (wt, pos, mut) in sorted(muts):
                        fh.write(f"{seq_id},{wt}{pos}{mut}\n")
            cmd += ["--mutation-file", mut_file]

        try:
            # cwd=repo_dir: vespag's load_model resolves its FNN weights as
            # Path.cwd()/model_weights/v2/esm2.pt, which live in the clone.
            # All -i/-o paths above are absolute, so cwd doesn't affect I/O.
            subprocess.run(cmd, check=True, cwd=args.repo_dir or None)
        except Exception as exc:  # noqa: BLE001 — whole-batch failure
            for seq_id, _seq in sequences:
                missing.append({"id": seq_id, "removed_by": step_id, "kind": "failure",
                                "cause": f"vespag predict failed: {type(exc).__name__}: {exc}"})
            traceback.print_exc()
            _write_outputs(rows, missing, args)
            print("VespaG: prediction failed for the whole batch.")
            return

        scores_csv = os.path.join(work, "vespag_scores_all.csv")
        if not os.path.isfile(scores_csv):
            for seq_id, _seq in sequences:
                missing.append({"id": seq_id, "removed_by": step_id, "kind": "failure",
                                "cause": f"expected output {os.path.basename(scores_csv)} not found"})
            _write_outputs(rows, missing, args)
            print("VespaG: no scores produced.")
            return

        df = pd.read_csv(scores_csv)
        # Columns: Protein, Mutation, VespaG. Mutation is 1-indexed e.g. 'A42G'.
        for _, r in df.iterrows():
            seq_id = str(r["Protein"])
            m = MUT_RE.match(str(r["Mutation"]))
            if not m:
                continue
            wt, pos, mut = m.group(1).upper(), int(m.group(2)), m.group(3).upper()
            rows.append({
                "id": f"{seq_id}_{wt}{pos}{mut}",
                "sequences.id": seq_id,
                "position": pos,
                "wildtype": wt,
                "mutation": mut,
                "fitness": float(r["VespaG"]),
            })

    # Report explicitly-requested mutations VespaG did not return.
    if any_explicit:
        produced: Dict[str, Set[Tuple[str, int, str]]] = {}
        for rr in rows:
            produced.setdefault(rr["sequences.id"], set()).add(
                (rr["wildtype"], rr["position"], rr["mutation"]))
        for seq_id, muts in requested_by_seq.items():
            if not muts:
                continue
            for (wt, pos, mut) in sorted(muts - produced.get(seq_id, set())):
                missing.append({
                    "id": f"{seq_id}_{wt}{pos}{mut}",
                    "removed_by": step_id,
                    "kind": "failure",
                    "cause": f"requested mutation {wt}{pos}{mut} not in VespaG output "
                             f"(position/wildtype mismatch against the sequence?)",
                })

    _write_outputs(rows, missing, args)
    print(f"VespaG: scored {len(rows)} mutation(s); {len(missing)} issue(s).")


def _write_outputs(rows: List[dict], missing: List[dict], args) -> None:
    fitness_df = pd.DataFrame(rows, columns=OUTPUT_COLUMNS)
    if not fitness_df.empty:
        fitness_df = fitness_df.sort_values("fitness", ascending=False).reset_index(drop=True)
    os.makedirs(os.path.dirname(args.fitness_csv), exist_ok=True)
    fitness_df.to_csv(args.fitness_csv, index=False)

    missing_df = pd.DataFrame(missing, columns=["id", "removed_by", "kind", "cause"])
    os.makedirs(os.path.dirname(args.missing_csv), exist_ok=True)
    missing_df.to_csv(args.missing_csv, index=False)


if __name__ == "__main__":
    main()
