#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Frame2Seq runtime helper.

Iterates the input structures DataStream, invokes Frame2seqRunner().design() per
structure, parses the combined seqs.fasta Frame2seq writes under each per-PDB
output folder, and emits BioPipelines-style outputs:

  - sequences.csv : content-bearing map_table for the `sequences` stream.
                    Columns: id, structures.id, sequence, score, recovery
  - sequences.fasta: the same records concatenated as a multi-record FASTA
  - missing.csv   : per-structure failures (id, removed_by, kind, cause)

Runs under the `frame2seq` conda env. The Frame2seq upstream API:

    from frame2seq import Frame2seqRunner
    runner = Frame2seqRunner()          # downloads weights to ~/.cache on first call
    runner.design(pdb_file, chain_id, temperature, num_samples,
                  omit_AA=None, fixed_positions=None,
                  save_indiv_seqs=False, save_indiv_neg_pll=False, verbose=True)

Frame2seq writes:
    <runner.save_dir>/seqs/seqs.fasta
with header
    >pdbid=<X> chain_id=<C> recovery=YY.YY% score=ZZ.ZZ temperature=T
"""

import argparse
import os
import re
import shutil
import sys
import traceback
from typing import List, Optional, Tuple

import pandas as pd

# Make biopipelines importable from a SLURM/Colab node where the package is
# not pip-installed (matches the pattern documented in developer_manual.md).
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from biopipelines.biopipelines_io import (  # noqa: E402
    load_datastream, iterate_files, load_table, lookup_table_value,
)
from biopipelines.pdb_parser import parse_pdb_file  # noqa: E402


HEADER_RE = re.compile(
    r">pdbid=(?P<pdbid>\S+)\s+chain_id=(?P<chain>\S+)\s+recovery=(?P<recovery>[0-9.]+)%\s+score=(?P<score>[0-9.\-]+)"
)


def parse_chain_aware_positions(spec: str, chain: str) -> Optional[List[int]]:
    """Parse a chain-aware selection string (e.g. "A10-20+A30+B5") and return
    the 1-indexed residue numbers belonging to `chain`. Returns None if `spec`
    is empty. Residues from other chains are silently dropped — Frame2seq is
    single-chain and the design call already pins `chain_id=chain`.

    Accepts either chain-aware (`A10`) or chainless (`10`) tokens; chainless
    tokens are assumed to belong to `chain`.
    """
    if not spec:
        return None
    out: List[int] = []
    for token in spec.replace(",", "+").split("+"):
        token = token.strip()
        if not token:
            continue
        # Optional leading chain letter
        m = re.match(r"^([A-Za-z])?(\d+)(?:-(\d+))?$", token)
        if not m:
            raise ValueError(f"Bad fixed_positions token: {token!r}")
        tok_chain, lo, hi = m.group(1), int(m.group(2)), m.group(3)
        if tok_chain is not None and tok_chain != chain:
            continue
        hi = int(hi) if hi else lo
        out.extend(range(lo, hi + 1))
    return sorted(set(out)) or None


def resolve_selection_arg(arg: str, struct_id: str, chain: str) -> Optional[List[int]]:
    """Resolve a fixed/redesigned selection argument to 1-indexed residue numbers.

    `arg` is either '-' (none), a chain-aware selection string, or a
    'TABLE_REFERENCE:<path>:<col>' token whose per-structure cell holds the
    selection string. Returns the sorted residue list for `chain`, or None.
    """
    if not arg or arg == "-":
        return None
    if arg.startswith("TABLE_REFERENCE:"):
        table, column = load_table(arg)
        sel = lookup_table_value(table, struct_id, column)
        return parse_chain_aware_positions(sel or "", chain)
    return parse_chain_aware_positions(arg, chain)


def chain_residue_numbers(pdb_file: str, chain: str) -> List[int]:
    """Return the sorted unique residue numbers of `chain` in the PDB."""
    atoms = parse_pdb_file(pdb_file)
    nums = {a.res_num for a in atoms if a.chain == chain}
    return sorted(nums)


def run_design_for_structure(
    runner,
    pdb_file: str,
    chain: str,
    temperature: float,
    num_samples: int,
    omit_aa: Optional[List[str]],
    fixed_positions: Optional[List[int]],
    save_dir: str,
) -> str:
    """Call Frame2seqRunner.design() with save_dir routed to `save_dir`.

    Returns the path to the combined seqs.fasta Frame2seq wrote.
    """
    # Frame2seqRunner is initialized with save_dir at construction time; we
    # rebind it per-call so each structure gets its own output folder.
    runner.save_dir = save_dir
    os.makedirs(save_dir, exist_ok=True)

    runner.design(
        pdb_file=pdb_file,
        chain_id=chain,
        temperature=temperature,
        num_samples=num_samples,
        omit_AA=omit_aa,
        fixed_positions=fixed_positions,
        save_indiv_seqs=False,
        save_indiv_neg_pll=False,
        verbose=False,
    )

    fasta_path = os.path.join(save_dir, "seqs", "seqs.fasta")
    if not os.path.isfile(fasta_path):
        raise FileNotFoundError(f"Frame2seq did not produce {fasta_path}")
    return fasta_path


def parse_seqs_fasta(fasta_path: str) -> List[Tuple[float, float, str]]:
    """Return a list of (recovery, score, sequence) tuples in the order
    Frame2seq wrote them. The order matches sample index 1..num_samples."""
    records: List[Tuple[float, float, str]] = []
    cur_header: Optional[re.Match] = None
    cur_seq: List[str] = []

    def flush():
        if cur_header is None:
            return
        seq = "".join(cur_seq).strip()
        records.append(
            (
                float(cur_header.group("recovery")),
                float(cur_header.group("score")),
                seq,
            )
        )

    with open(fasta_path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                flush()
                m = HEADER_RE.match(line)
                if not m:
                    raise ValueError(f"Unrecognised Frame2seq FASTA header: {line!r}")
                cur_header = m
                cur_seq = []
            else:
                cur_seq.append(line)
    flush()
    return records


def main():
    ap = argparse.ArgumentParser(description="Frame2Seq inverse-folding runtime helper")
    ap.add_argument("--ds-json", required=True, help="Input structures DataStream JSON")
    ap.add_argument("--output-folder", required=True, help="Per-structure raw outputs (execution_folder)")
    ap.add_argument("--sequences-csv", required=True, help="Output content-bearing sequences CSV (map_table)")
    ap.add_argument("--sequences-fasta", required=True, help="Output combined FASTA")
    ap.add_argument("--missing-csv", required=True, help="Output per-structure failure CSV")
    ap.add_argument("--num-sequences", type=int, required=True)
    ap.add_argument("--temperature", type=float, required=True)
    ap.add_argument("--chain", required=True)
    ap.add_argument("--omit-aa", default="", help="Concatenated single-letter codes to omit, e.g. 'CM'")
    ap.add_argument("--fixed", default="-",
                    help="Fixed-position selection string or TABLE_REFERENCE:path:col, or '-'")
    ap.add_argument("--redesigned", default="-",
                    help="Redesigned-position selection string or TABLE_REFERENCE:path:col, or '-'. "
                         "The complement (other residues of --chain) is held fixed.")
    args = ap.parse_args()

    # Lazy import: heavy and only needed inside the frame2seq env.
    from frame2seq import Frame2seqRunner

    runner = Frame2seqRunner()  # weights are downloaded to ~/.cache on first call

    omit_list = list(args.omit_aa) if args.omit_aa else None

    ds = load_datastream(args.ds_json)

    rows: List[dict] = []
    failed: List[Tuple[str, str]] = []

    for struct_id, pdb_file in iterate_files(ds):
        save_dir = os.path.join(args.output_folder, struct_id)
        try:
            # Resolve fixed positions per structure. `redesigned` is converted
            # to its complement against the chain's residues so Frame2seq (which
            # only takes a fixed list) designs exactly the requested positions.
            fixed_list = resolve_selection_arg(args.fixed, struct_id, args.chain)
            redesign_list = resolve_selection_arg(args.redesigned, struct_id, args.chain)
            if redesign_list is not None:
                all_resis = chain_residue_numbers(pdb_file, args.chain)
                redesign_set = set(redesign_list)
                fixed_list = sorted(r for r in all_resis if r not in redesign_set) or None

            fasta_path = run_design_for_structure(
                runner=runner,
                pdb_file=pdb_file,
                chain=args.chain,
                temperature=args.temperature,
                num_samples=args.num_sequences,
                omit_aa=omit_list,
                fixed_positions=fixed_list,
                save_dir=save_dir,
            )
            records = parse_seqs_fasta(fasta_path)
        except Exception as exc:
            cause = f"{type(exc).__name__}: {exc}".replace(",", ";").replace("\n", " ")
            print(f"  [fail] {struct_id}: {cause}", file=sys.stderr)
            traceback.print_exc()
            failed.append((struct_id, cause))
            continue

        if len(records) != args.num_sequences:
            cause = f"expected {args.num_sequences} records, got {len(records)}"
            print(f"  [warn] {struct_id}: {cause}", file=sys.stderr)
            # We still emit whatever Frame2seq produced; downstream tools work
            # with partial output. Record the discrepancy.
            failed.append((struct_id, cause))

        for n, (recovery, score, sequence) in enumerate(records, start=1):
            rows.append({
                "id": f"{struct_id}_{n}",
                "structures.id": struct_id,
                "sequence": sequence,
                "score": score,
                "recovery": recovery,
            })

    # Always write the outputs, even on partial failure — downstream tools
    # can keep working with a subset (developer_manual.md error-handling rule).
    os.makedirs(os.path.dirname(args.sequences_csv), exist_ok=True)
    pd.DataFrame(rows, columns=["id", "structures.id", "sequence", "score", "recovery"]).to_csv(
        args.sequences_csv, index=False
    )

    # FASTA twin: simple per-row records (no Frame2seq metadata, just id+seq).
    with open(args.sequences_fasta, "w") as f:
        for row in rows:
            f.write(f">{row['id']}\n{row['sequence']}\n")

    # Missing table (always write the file with the header so downstream
    # tools that read it find a well-formed CSV even when nothing failed).
    os.makedirs(os.path.dirname(args.missing_csv), exist_ok=True)
    with open(args.missing_csv, "w") as f:
        f.write("id,removed_by,kind,cause\n")
        for sid, cause in failed:
            f.write(f"{sid},Frame2Seq,failure,\"{cause}\"\n")

    if not rows:
        print("ERROR: Frame2Seq produced no sequences", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
