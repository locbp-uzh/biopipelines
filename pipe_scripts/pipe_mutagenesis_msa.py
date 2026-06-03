#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Synthetic per-mutant MSA generation for the Mutagenesis tool.

Given the mutants table written by pipe_mutagenesis.py and an input MSA
map_table for the original protein(s), this derives one MSA per mutant by
copying the parent's MSA (matched on `original.id`) and substituting ONLY the
query (first) row at the mutated positions. Homolog rows pass through
unchanged, and the on-disk format (a3m / csv) is preserved.

This is the standard treatment for point substitutions: a single residue swap
does not shift alignment columns (no indels), so no realignment is needed and
the family's coevolutionary signal in the homolog rows is left intact.

Usage:
    python pipe_mutagenesis_msa.py --mutants MUTANTS_CSV --input-msas INPUT_MAP_CSV
           --output-folder OUT_DIR --output OUT_MAP_CSV

Inputs:
    --mutants: mutants table (id, original.id, sequence, mutation_positions,
        new_aa). mutation_positions is PyMOL-style ("50" or "57+58"); new_aa is
        the per-position substituted residue(s) ("A" or "AS"). Empty for
        passthrough (selection) rows.
    --input-msas: MSA map_table with columns id, sequences.id, msa_file. Keyed
        by the ORIGINAL protein id (matched against the mutants' original.id).
    --output-folder: where per-mutant MSA files are written.
    --output: output MSA map_table (id, sequences.id, original.id, sequence,
        msa_file), keyed by mutant id.

Output map columns mirror the MSA stream contract used elsewhere
(id | sequences.id | sequence | msa_file): both `id` and `sequences.id` are the
mutant id, so downstream MSA consumers (AlphaFold/Boltz copy step, which name
the copied A3M after `sequences.id`) match the mutant query correctly. An extra
`original.id` column records the parent protein the MSA was derived from.
"""

import argparse
import os
import sys

import pandas as pd

# Add repo root so biopipelines is importable on SLURM nodes / base Python.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Derive synthetic per-mutant MSAs by query-row substitution"
    )
    parser.add_argument("--mutants", required=True, help="Mutants table CSV")
    parser.add_argument("--input-msas", required=True, help="Input MSA map_table CSV (keyed by original id)")
    parser.add_argument("--output-folder", required=True, help="Output folder for per-mutant MSA files")
    parser.add_argument("--output", required=True, help="Output MSA map_table CSV")
    return parser.parse_args()


def read_a3m(path):
    """Read an A3M file → (header_lines, [(header, sequence), ...]).

    header_lines are leading `#` metadata lines (A3M count/width); records are
    (>header, sequence) pairs in file order, query first.
    """
    header_lines = []
    records = []
    cur_header = None
    cur_seq = ""
    with open(path) as f:
        for line in f:
            stripped = line.rstrip("\n")
            if cur_header is None and stripped.startswith("#"):
                header_lines.append(stripped)
                continue
            if stripped.startswith(">"):
                if cur_header is not None:
                    records.append((cur_header, cur_seq))
                cur_header = stripped
                cur_seq = ""
            elif stripped:
                cur_seq += stripped
    if cur_header is not None:
        records.append((cur_header, cur_seq))
    return header_lines, records


def write_a3m(path, header_lines, records):
    with open(path, "w") as f:
        for h in header_lines:
            f.write(h + "\n")
        for header, seq in records:
            f.write(f"{header}\n{seq}\n")


def read_csv_msa(path):
    """Read a Boltz2-style CSV MSA → list of sequences (query first)."""
    df = pd.read_csv(path)
    if "sequence" not in df.columns:
        raise ValueError(f"CSV MSA missing 'sequence' column: {path}")
    return df["sequence"].tolist()


def write_csv_msa(path, sequences):
    pd.DataFrame([{"key": -1, "sequence": s} for s in sequences]).to_csv(path, index=False)


def parse_positions(positions_str):
    """PyMOL-style positions ("50", "57+58", "60-62") → sorted unique ints.

    Mirrors pipe_mutagenesis.parse_positions_selection so position order is
    identical to how new_aa was assembled.
    """
    positions = []
    for part in str(positions_str).split("+"):
        part = part.strip()
        if not part:
            continue
        if "-" in part and not part.startswith("-"):
            lo, hi = part.split("-")
            positions.extend(range(int(lo), int(hi) + 1))
        else:
            positions.append(int(part))
    return sorted(set(positions))


def substitute_query(query_seq, positions_str, new_aa):
    """Return the query sequence with each mutated position replaced.

    positions are 1-indexed and aligned with new_aa character-by-character
    (one residue per position, same order as parse_positions). For aligned
    (A3M) sequences the leading uppercase columns correspond to query
    positions; since point substitutions don't introduce gaps in the query,
    the i-th non-gap column equals the i-th sequence position. Boltz CSV MSAs
    carry ungapped query rows, so plain index substitution applies there too.
    """
    positions = parse_positions(positions_str)
    if len(positions) != len(new_aa):
        raise ValueError(
            f"position count ({len(positions)}) != new_aa length ({len(new_aa)}): "
            f"{positions_str!r} vs {new_aa!r}"
        )
    chars = list(query_seq)
    for pos, aa in zip(positions, new_aa):
        idx = pos - 1
        if idx < 0 or idx >= len(chars):
            raise ValueError(f"position {pos} out of range for query length {len(chars)}")
        chars[idx] = aa
    return "".join(chars)


def main():
    args = parse_arguments()

    mutants = pd.read_csv(args.mutants)
    for col in ("id", "original.id"):
        if col not in mutants.columns:
            print(f"Error: mutants table missing '{col}' column", file=sys.stderr)
            sys.exit(1)

    input_msas = pd.read_csv(args.input_msas)
    if "msa_file" not in input_msas.columns:
        print("Error: input MSA table missing 'msa_file' column", file=sys.stderr)
        sys.exit(1)

    # Map original protein id -> its MSA file. The input map's own `id` column
    # is the original protein id (Mutagenesis keys MSAs by original.id).
    parent_msa = {}
    for _, row in input_msas.iterrows():
        parent_id = str(row.get("id", ""))
        msa_file = row.get("msa_file", "")
        if parent_id and isinstance(msa_file, str) and msa_file:
            parent_msa[parent_id] = msa_file

    os.makedirs(args.output_folder, exist_ok=True)

    out_rows = []
    skipped = []

    for _, row in mutants.iterrows():
        mutant_id = str(row["id"])
        original_id = str(row["original.id"])
        positions_str = "" if pd.isna(row.get("mutation_positions")) else str(row.get("mutation_positions", ""))
        new_aa = "" if pd.isna(row.get("new_aa")) else str(row.get("new_aa", ""))

        parent_file = parent_msa.get(original_id)
        if not parent_file or not os.path.exists(parent_file):
            print(f"WARNING: no input MSA for original id '{original_id}' (mutant '{mutant_id}'); skipping",
                  file=sys.stderr)
            skipped.append(mutant_id)
            continue

        is_a3m = parent_file.lower().endswith(".a3m")
        ext = ".a3m" if is_a3m else ".csv"
        out_file = os.path.join(args.output_folder, f"{mutant_id}{ext}")

        try:
            if is_a3m:
                header_lines, records = read_a3m(parent_file)
                if not records:
                    raise ValueError(f"no records in {parent_file}")
                q_header, q_seq = records[0]
                # Passthrough rows (no mutation) keep the query untouched.
                new_query = q_seq if not positions_str or not new_aa else substitute_query(q_seq, positions_str, new_aa)
                records[0] = (q_header, new_query)
                write_a3m(out_file, header_lines, records)
            else:
                sequences = read_csv_msa(parent_file)
                if not sequences:
                    raise ValueError(f"no sequences in {parent_file}")
                new_query = sequences[0] if not positions_str or not new_aa else substitute_query(sequences[0], positions_str, new_aa)
                sequences[0] = new_query
                write_csv_msa(out_file, sequences)
                new_query = sequences[0]
        except Exception as e:
            print(f"WARNING: {mutant_id} failed: {e}", file=sys.stderr)
            skipped.append(mutant_id)
            continue

        out_rows.append({
            "id": mutant_id,
            "sequences.id": mutant_id,
            "original.id": original_id,
            "sequence": new_query,
            "msa_file": out_file,
        })

    out_df = pd.DataFrame(out_rows, columns=["id", "sequences.id", "original.id", "sequence", "msa_file"])
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    out_df.to_csv(args.output, index=False)

    print(f"\nWrote {len(out_rows)} synthetic MSAs to {args.output_folder}")
    print(f"Output map: {args.output}")
    if skipped:
        print(f"Skipped {len(skipped)}/{len(out_rows) + len(skipped)}: {skipped}", file=sys.stderr)

    # Error only if every mutant failed.
    if not out_rows:
        sys.exit(1)


if __name__ == "__main__":
    main()
