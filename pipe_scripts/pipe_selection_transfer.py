#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Runtime half of SelectionTransfer: map a residue selection from reference structures onto target structures.

Each selected reference chain is aligned against every target chain, with the reference sequence doubled so a circular permutation aligns as one contiguous block; the best-scoring target chain wins. Reference residues aligned to a gap are reported as unmapped.
"""

import argparse
import json
import os
import sys
from typing import Dict, List, Optional, Tuple

import pandas as pd
from Bio.Align import PairwiseAligner, substitution_matrices

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from biopipelines.biopipelines_io import (  # noqa: E402
    load_datastream, iterate_files, load_table, lookup_table_value, step_id_from_table_path,
)
from biopipelines.id_map_utils import get_mapped_ids  # noqa: E402
from biopipelines.pdb_parser import parse_pdb_file, STANDARD_RESIDUES, THREE_TO_ONE  # noqa: E402
from biopipelines.sele_utils import sele_to_list, list_to_sele  # noqa: E402

COLUMNS = ["id", "selection", "n_residues", "source_selection", "unmapped", "identity", "reference.id"]

Chains = Dict[str, List[Tuple[int, str]]]


def make_aligner() -> PairwiseAligner:
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10.0
    aligner.extend_gap_score = -0.5
    # Trimmed or extended termini, and the unused half of the doubled reference, cost nothing.
    aligner.end_gap_score = 0.0
    return aligner


def read_chains(path: str) -> Chains:
    """Ordered (resnum, one-letter) protein residues per chain."""
    chains: Chains = {}
    seen = set()
    for atom in parse_pdb_file(path):
        if atom.res_name not in STANDARD_RESIDUES:
            continue
        key = (atom.chain, atom.res_num)
        if key in seen:
            continue
        seen.add(key)
        chains.setdefault(atom.chain, []).append((atom.res_num, THREE_TO_ONE[atom.res_name]))
    if not chains:
        raise ValueError(f"no protein residues in {path}")
    return chains


def align_chain(aligner: PairwiseAligner, ref: List[Tuple[int, str]],
                tgt: List[Tuple[int, str]]) -> Tuple[float, Dict[int, int], int, int]:
    """Return (score, {ref index: target index}, aligned pairs, identical pairs)."""
    ref_seq = "".join(aa for _, aa in ref)
    tgt_seq = "".join(aa for _, aa in tgt)
    n = len(ref_seq)
    alignment = aligner.align(ref_seq + ref_seq, tgt_seq)[0]
    mapping: Dict[int, int] = {}
    identical = 0
    for (r0, r1), (t0, t1) in zip(*alignment.aligned):
        for i, j in zip(range(r0, r1), range(t0, t1)):
            if i % n in mapping:
                raise ValueError("a reference residue aligned twice across the doubled sequence")
            mapping[i % n] = j
            identical += ref_seq[i % n] == tgt_seq[j]
    return alignment.score, mapping, len(mapping), identical


def transfer(residues: List[Tuple[str, int]], ref_chains: Chains, tgt_chains: Chains,
             aligner: PairwiseAligner):
    """Map (chain, resnum) residues of the reference onto the target."""
    if any(c == "" for c, _ in residues):
        if len(ref_chains) != 1:
            raise ValueError(f"chainless selection on a reference with chains {sorted(ref_chains)}")
        only = next(iter(ref_chains))
        residues = [(only if c == "" else c, r) for c, r in residues]

    mapped, unmapped = set(), []
    aligned_total = identical_total = 0
    for ref_chain in sorted({c for c, _ in residues}):
        if ref_chain not in ref_chains:
            raise ValueError(f"chain {ref_chain} is not in the reference")
        ref = ref_chains[ref_chain]
        index = {resnum: i for i, (resnum, _) in enumerate(ref)}
        wanted = [r for c, r in residues if c == ref_chain]
        absent = [r for r in wanted if r not in index]
        if absent:
            raise ValueError(f"residues {list_to_sele([(ref_chain, r) for r in absent])} are not in the reference")

        best = None
        for tgt_chain, tgt in tgt_chains.items():
            score, mapping, n_aligned, n_identical = align_chain(aligner, ref, tgt)
            rank = (score, tgt_chain == ref_chain)
            if best is None or rank > best[0]:
                best = (rank, tgt_chain, mapping, n_aligned, n_identical)
        _, tgt_chain, mapping, n_aligned, n_identical = best
        aligned_total += n_aligned
        identical_total += n_identical

        for r in wanted:
            j = mapping.get(index[r])
            if j is None:
                unmapped.append((ref_chain, r))
            else:
                mapped.add((tgt_chain, tgt_chains[tgt_chain][j][0]))

    identity = identical_total / aligned_total if aligned_total else float("nan")
    return sorted(mapped), sorted(unmapped), identity


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--config-json", required=True)
    p.add_argument("--structures-json", required=True)
    p.add_argument("--reference-json", required=True)
    p.add_argument("--selections-csv", required=True)
    p.add_argument("--local-missing-csv", required=True)
    args = p.parse_args()

    with open(args.config_json) as f:
        selection_cfg = json.load(f)["selection"]

    structures_ds = load_datastream(args.structures_json)
    reference_ds = load_datastream(args.reference_json)
    targets = dict(iterate_files(structures_ds))
    references = dict(iterate_files(reference_ds))
    step_id = step_id_from_table_path(args.local_missing_csv)

    if not targets:
        print("ERROR: no target structures to process", file=sys.stderr)
        sys.exit(1)
    if not references:
        print("ERROR: no reference structures available", file=sys.stderr)
        sys.exit(1)

    if len(references) == 1:
        only = next(iter(references))
        pairing: Dict[str, Optional[str]] = {tid: only for tid in targets}
    else:
        maps = [structures_ds.map_table] if structures_ds.map_table else None
        pairing = get_mapped_ids(list(targets), list(references), unique=True, map_table_paths=maps)

    table = None
    if "table" in selection_cfg:
        table, _ = load_table(selection_cfg["table"])
    ref_maps = [reference_ds.map_table] if reference_ds.map_table else None

    aligner = make_aligner()
    ref_cache: Dict[str, Chains] = {}
    rows, missing_rows = [], []
    for tid, tpath in targets.items():
        rid = pairing.get(tid)
        try:
            if not rid:
                raise ValueError("no reference structure matches this id")
            if table is None:
                source = selection_cfg["literal"]
            else:
                source = lookup_table_value(table, rid, selection_cfg["column"], map_table_paths=ref_maps)
            source = "" if pd.isna(source) else str(source)
            if rid not in ref_cache:
                ref_cache[rid] = read_chains(references[rid])
            mapped, unmapped, identity = transfer(sele_to_list(source), ref_cache[rid],
                                                  read_chains(tpath), aligner)
            selection = list_to_sele(mapped)
            rows.append({"id": tid, "selection": selection, "n_residues": len(mapped),
                         "source_selection": source, "unmapped": list_to_sele(unmapped),
                         "identity": round(identity, 4), "reference.id": rid})
            note = f", unmapped {list_to_sele(unmapped)}" if unmapped else ""
            print(f"  {tid} <- {rid}: {source} -> {selection} (identity {identity:.2f}{note})")
        except Exception as e:
            print(f"WARNING: {tid} failed: {e}", file=sys.stderr)
            rows.append({"id": tid, "selection": None, "n_residues": None, "source_selection": None,
                         "unmapped": None, "identity": None, "reference.id": rid})
            missing_rows.append({"id": tid, "removed_by": step_id, "kind": "failure", "cause": str(e)[:200]})

    for path in (args.selections_csv, args.local_missing_csv):
        os.makedirs(os.path.dirname(path), exist_ok=True)
    pd.DataFrame(rows, columns=COLUMNS).to_csv(args.selections_csv, index=False)
    pd.DataFrame(missing_rows, columns=["id", "removed_by", "kind", "cause"]).to_csv(
        args.local_missing_csv, index=False)

    print(f"Selections: {args.selections_csv} ({len(rows)} rows, {len(missing_rows)} failed)")
    if len(missing_rows) == len(rows):
        print(f"ERROR: all {len(rows)} id(s) failed", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
