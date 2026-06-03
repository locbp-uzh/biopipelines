#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""ThermoMPNN runtime helper.

Iterates the input structures DataStream, runs ThermoMPNN's
`analysis/custom_inference.py` per structure (site-saturation over `--chain`),
parses the per-structure ``ThermoMPNN_inference_<pdb>.csv`` it writes, optionally
filters it down to an explicit list of requested mutations, and emits
BioPipelines-style outputs:

  - ddg.csv     : id | structures.id | chain | position | wildtype | mutation | ddG_pred
                  (sorted by ddG_pred ascending across the whole table; most
                  stabilising first)
  - missing.csv : per-structure failures (id | removed_by | cause)

Upstream custom_inference.py CLI:
    python analysis/custom_inference.py --pdb X.pdb --chain A --out_dir DIR
writes DIR/ThermoMPNN_inference_<pdb_id>.csv with columns
    Model, Dataset, ddG_pred, position, wildtype, mutation, pdb, chain
The bundled checkpoint under models/ is found automatically when --model_path
is omitted.

Runs under the `thermompnn` conda env.
"""

import argparse
import os
import subprocess
import sys
import tempfile
import traceback
from typing import List, Optional, Set, Tuple

import pandas as pd

# Make biopipelines importable from a SLURM/Colab node where the package is
# not pip-installed (matches the pattern used by the other pipe scripts).
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from biopipelines.biopipelines_io import (  # noqa: E402
    load_datastream, iterate_files, load_table, lookup_table_value,
    step_id_from_table_path,
)

OUTPUT_COLUMNS = ["id", "structures.id", "chain", "position",
                  "wildtype", "mutation", "ddG_pred"]


def parse_mutation_tokens(spec: str) -> Set[Tuple[str, int, str]]:
    """Parse a '+'-joined list of wildtype-position-mutant tokens.

    Each token is like ``A42G`` (wt A at position 42 -> G). Commas are accepted
    as separators too. Returns a set of (wildtype, position, mutant) tuples.
    Raises ValueError on a malformed token.
    """
    out: Set[Tuple[str, int, str]] = set()
    for token in spec.replace(",", "+").split("+"):
        token = token.strip()
        if not token:
            continue
        # <wt><pos><mut>, e.g. A42G
        if len(token) < 3 or not token[0].isalpha() or not token[-1].isalpha():
            raise ValueError(f"Bad mutation token: {token!r} (expected e.g. 'A42G')")
        wt, mut, pos_str = token[0].upper(), token[-1].upper(), token[1:-1]
        if not pos_str.isdigit():
            raise ValueError(f"Bad mutation token: {token!r} (no integer position)")
        out.add((wt, int(pos_str), mut))
    return out


def resolve_mutations_arg(arg: str, struct_id: str) -> Optional[Set[Tuple[str, int, str]]]:
    """Resolve the --mutations argument for a given structure.

    `arg` is '-' (saturation; returns None), a literal '+'-joined token string,
    or a 'TABLE_REFERENCE:<path>:<col>' token whose per-structure cell holds
    such a string. Returns the parsed set of (wt, pos, mut) tuples, or None for
    site-saturation.
    """
    if not arg or arg == "-":
        return None
    if arg.startswith("TABLE_REFERENCE:"):
        table, column = load_table(arg)
        cell = lookup_table_value(table, struct_id, column)
        if cell is None or (isinstance(cell, float) and pd.isna(cell)):
            return None
        return parse_mutation_tokens(str(cell))
    return parse_mutation_tokens(arg)


def run_inference(repo_dir: str, pdb_file: str, chain: str, out_dir: str) -> str:
    """Run custom_inference.py and return the path to the CSV it wrote.

    Invokes the script with its own directory as cwd (its checkpoint auto-search
    is relative to analysis/ via '../models/'), then locates the single
    ThermoMPNN_inference_*.csv produced in out_dir.
    """
    analysis_dir = os.path.join(repo_dir, "analysis")
    script = os.path.join(analysis_dir, "custom_inference.py")
    if not os.path.isfile(script):
        raise FileNotFoundError(f"custom_inference.py not found at {script}")

    before = set(f for f in os.listdir(out_dir) if f.endswith(".csv"))
    cmd = [
        sys.executable, script,
        "--pdb", os.path.abspath(pdb_file),
        "--chain", chain,
        "--out_dir", os.path.abspath(out_dir),
    ]
    # cwd=analysis_dir so the '../models/' checkpoint auto-search resolves.
    subprocess.run(cmd, cwd=analysis_dir, check=True)

    after = set(f for f in os.listdir(out_dir) if f.endswith(".csv"))
    new = sorted(after - before)
    if not new:
        # Fall back to any ThermoMPNN_inference_*.csv (e.g. re-run overwriting).
        new = sorted(f for f in after if f.startswith("ThermoMPNN_inference_"))
    if not new:
        raise FileNotFoundError(
            f"ThermoMPNN produced no output CSV in {out_dir}"
        )
    return os.path.join(out_dir, new[-1])


def main():
    ap = argparse.ArgumentParser(description="ThermoMPNN ddG runtime helper")
    ap.add_argument("--structures-json", required=True)
    ap.add_argument("--repo-dir", required=True)
    ap.add_argument("--chain", default="A")
    ap.add_argument("--mutations", default="-",
                    help="'-' for site-saturation, a '+'-joined token list, or a TABLE_REFERENCE")
    ap.add_argument("--ddg-csv", required=True)
    ap.add_argument("--missing-csv", required=True)
    args = ap.parse_args()

    structures = load_datastream(args.structures_json)

    rows: List[dict] = []
    missing: List[dict] = []
    step_id = step_id_from_table_path(args.missing_csv)

    with tempfile.TemporaryDirectory(prefix="thermompnn_") as work:
        for struct_id, pdb_file in iterate_files(structures):
            try:
                requested = resolve_mutations_arg(args.mutations, struct_id)

                struct_work = os.path.join(work, struct_id)
                os.makedirs(struct_work, exist_ok=True)
                csv_path = run_inference(args.repo_dir, pdb_file, args.chain, struct_work)

                df = pd.read_csv(csv_path)
                # Upstream columns: Model, Dataset, ddG_pred, position,
                # wildtype, mutation, pdb, chain. Normalise to our schema.
                for _, r in df.iterrows():
                    pos = int(r["position"])
                    wt = str(r["wildtype"]).upper()
                    mut = str(r["mutation"]).upper()
                    if requested is not None and (wt, pos, mut) not in requested:
                        continue
                    rows.append({
                        "id": f"{struct_id}_{wt}{pos}{mut}",
                        "structures.id": struct_id,
                        "chain": args.chain,
                        "position": pos,
                        "wildtype": wt,
                        "mutation": mut,
                        "ddG_pred": float(r["ddG_pred"]),
                    })

                if requested is not None:
                    found = {(rr["wildtype"], rr["position"], rr["mutation"])
                             for rr in rows if rr["structures.id"] == struct_id}
                    for (wt, pos, mut) in sorted(requested - found):
                        missing.append({
                            "id": f"{struct_id}_{wt}{pos}{mut}",
                            "removed_by": step_id,
                            "kind": "failure",
                            "cause": f"requested mutation {wt}{pos}{mut} not in ThermoMPNN output "
                                     f"(position/wildtype mismatch on chain {args.chain}?)",
                        })
            except Exception as exc:  # noqa: BLE001 — report, don't abort the batch
                missing.append({
                    "id": struct_id,
                    "removed_by": step_id,
                    "kind": "failure",
                    "cause": f"{type(exc).__name__}: {exc}",
                })
                traceback.print_exc()

    ddg_df = pd.DataFrame(rows, columns=OUTPUT_COLUMNS)
    if not ddg_df.empty:
        ddg_df = ddg_df.sort_values("ddG_pred", ascending=True).reset_index(drop=True)
    os.makedirs(os.path.dirname(args.ddg_csv), exist_ok=True)
    ddg_df.to_csv(args.ddg_csv, index=False)

    missing_df = pd.DataFrame(missing, columns=["id", "removed_by", "kind", "cause"])
    os.makedirs(os.path.dirname(args.missing_csv), exist_ok=True)
    missing_df.to_csv(args.missing_csv, index=False)

    print(f"ThermoMPNN: scored {len(ddg_df)} mutation(s); {len(missing_df)} issue(s).")


if __name__ == "__main__":
    main()
