# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Execution half of the BFactor tool: read per-residue B-factors and summarize them over named selections.

Selections arrive in the config yaml as literal strings or as TABLE_REFERENCE:path:column. Each referenced table is loaded once here rather than once per id, so a 400-structure run costs one read per selection instead of 400.
"""

import argparse
import os
import re
import sys

import pandas as pd
import yaml

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from biopipelines.biopipelines_io import (  # noqa: E402
    load_datastream,
    iterate_files,
    load_table,
    lookup_table_value,
    step_id_from_table_path,
)

SPAN = re.compile(r"^([A-Za-z]*)(\d+)(?:-(\d+))?$")
TABLE_PREFIX = "TABLE_REFERENCE:"


def parse_selection(text):
    """'A75-77+A274' -> [(chain or None, first, last)]. A part that is not a span is an error, not an empty selection."""
    spans = []
    for part in str(text or "").replace(",", "+").split("+"):
        part = part.strip()
        if not part:
            continue
        m = SPAN.match(part)
        if not m:
            raise ValueError(f"not a residue span: {part!r} (expected e.g. A75-77 or 274)")
        chain = m.group(1) or None
        first = int(m.group(2))
        spans.append((chain, first, int(m.group(3) or first)))
    return spans


def in_selection(spans, chain, resi):
    return any((c is None or c == chain) and lo <= resi <= hi for c, lo, hi in spans)


def check_chains(name, spans, residues):
    """An unqualified span is ambiguous only when more than one chain has residues inside it.

    Counting every chain would refuse a protein with a ligand or waters, which is the ordinary predictor output.
    """
    for chain, lo, hi in spans:
        if chain is not None:
            continue
        hit = sorted({c for c, r, _i, _b in residues if lo <= r <= hi})
        if len(hit) > 1:
            raise ValueError(f"selection {name!r} names no chain but residues {lo}-{hi} exist in "
                             f"chains {', '.join(hit)}; qualify it, e.g. {hit[0]}{lo}-{hi}")


def read_residues(path):
    """-> [(chain, resi, icode, bfactor)], one entry per residue, CA-preferred; 100 and 100A stay apart."""
    ca, per_residue = {}, {}
    with open(path, errors="ignore") as fh:
        for line in fh:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            try:
                key = (line[21].strip() or "_", int(line[22:26]), line[26:27].strip())
                b = float(line[60:66])
            except ValueError:
                continue
            per_residue.setdefault(key, []).append(b)
            if line[12:16].strip() == "CA":
                ca[key] = b
    return [(chain, resi, icode, ca[key] if key in ca else sum(vs) / len(vs))
            for key, vs in sorted(per_residue.items()) for chain, resi, icode in [key]]


def stats(values):
    n = len(values)
    if not n:
        return None, None, None, 0
    mean = sum(values) / n
    sd = (sum((v - mean) ** 2 for v in values) / n) ** 0.5 if n > 1 else 0.0
    return mean, sd, min(values), n


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--config-yaml", required=True)
    p.add_argument("--structures-json", required=True)
    p.add_argument("--bfactors-dir", required=True)
    p.add_argument("--bfactors-map-csv", required=True)
    p.add_argument("--summary-csv", required=True)
    p.add_argument("--local-missing-csv", required=True)
    args = p.parse_args()

    with open(args.config_yaml) as fh:
        config = yaml.safe_load(fh) or {}
    selections = config.get("selections", {}) or {}

    # One read per referenced table, not one per structure.
    resolved = {}
    for name, spec in selections.items():
        if str(spec).startswith(TABLE_PREFIX):
            frame, column = load_table(str(spec))
            resolved[name] = ("table", frame, column)
        else:
            resolved[name] = ("literal", parse_selection(spec), None)

    os.makedirs(args.bfactors_dir, exist_ok=True)
    ds = load_datastream(args.structures_json)
    step_id = step_id_from_table_path(args.local_missing_csv)

    map_rows, summary_rows, local_missing_rows = [], [], []
    attempted = failures = 0

    for sid, structure_path in iterate_files(ds):
        attempted += 1
        row = {"id": sid}
        try:
            residues = read_residues(structure_path)
            if not residues:
                raise ValueError("no ATOM/HETATM records carrying a B-factor")

            out_csv = os.path.join(args.bfactors_dir, f"{sid}.csv")
            pd.DataFrame(
                [{"id": sid, "chain": c, "resi": r, "icode": i, "bfactor": b} for c, r, i, b in residues]
            ).to_csv(out_csv, index=False)
            map_rows.append({"id": sid, "file": out_csv})

            all_mean, all_sd, _, n_all = stats([b for _, _, _, b in residues])
            row.update({"n_residues": n_all, "all_mean": all_mean, "all_sd": all_sd})

            for name, entry in resolved.items():
                kind, first, column = entry
                if kind == "table":
                    spans = parse_selection(lookup_table_value(first, sid, column))
                else:
                    spans = first
                check_chains(name, spans, residues)
                picked = [b for c, r, _i, b in residues if in_selection(spans, c, r)]
                mean, sd, lo, n = stats(picked)
                row[f"{name}_mean"] = mean
                row[f"{name}_sd"] = sd
                row[f"{name}_min"] = lo
                row[f"{name}_n"] = n
                row[f"{name}_delta"] = (mean - all_mean) if mean is not None else None
        except Exception as e:
            print(f"WARNING: {sid} failed: {e}", file=sys.stderr)
            failures += 1
            local_missing_rows.append({"id": sid, "removed_by": step_id,
                                       "kind": "failure", "cause": str(e)[:200]})
        summary_rows.append(row)

    columns = ["id", "n_residues", "all_mean", "all_sd"]
    for name in selections:
        columns += [f"{name}_{s}" for s in ("mean", "sd", "min", "n", "delta")]

    for path in (args.summary_csv, args.bfactors_map_csv, args.local_missing_csv):
        os.makedirs(os.path.dirname(path), exist_ok=True)
    pd.DataFrame(summary_rows, columns=columns).to_csv(args.summary_csv, index=False)
    pd.DataFrame(map_rows, columns=["id", "file"]).to_csv(args.bfactors_map_csv, index=False)
    pd.DataFrame(local_missing_rows,
                 columns=["id", "removed_by", "kind", "cause"]).to_csv(args.local_missing_csv, index=False)

    print(f"Summary: {args.summary_csv} ({len(summary_rows)} rows, {failures} failed)")
    print(f"Per-residue: {args.bfactors_map_csv} ({len(map_rows)} files)")

    if not attempted:
        print("ERROR: no input entities to process (upstream stream is empty or its files are absent)",
              file=sys.stderr)
        sys.exit(1)
    if failures == attempted:
        print(f"ERROR: all {attempted} attempted id(s) failed", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
