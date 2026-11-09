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
from typing import Dict, List, Tuple

import pandas as pd
from Bio.Align import PairwiseAligner, substitution_matrices

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from biopipelines.biopipelines_io import (  # noqa: E402
    load_datastream, iterate_files, load_table, lookup_table_value, step_id_from_table_path,
)
from biopipelines.id_map_utils import TIER_SIBLING, get_mapped_ids_with_tiers  # noqa: E402
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


def pair_references(target_ids, reference_ids, n_declared, map_table_paths=None):
    """Return {target id: (reference id or None, sibling id that was refused or None)}.

    A reference stream declaring one id is the reference of every target. Several are matched to each target by id: the same id, the id its provenance names in the target stream's map table, or a parent or child of it. A sibling id (a common ancestor, nothing more) is not a match: its selection is another design's, its loop can sit at other positions or have another length, and the alignment would still map it with plausible identity.

    `n_declared` counts the ids the stream declares, not the files present, so a reference filtered upstream is never broadcast to the other targets.
    """
    if n_declared == 1:
        only = next(iter(reference_ids), None)
        return {tid: (only, None) for tid in target_ids}
    matches = get_mapped_ids_with_tiers(list(target_ids), list(reference_ids), unique=True,
                                        map_table_paths=map_table_paths)
    return {tid: ((None, rid) if tier == TIER_SIBLING else (rid, None)) for tid, (rid, tier) in matches.items()}


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

    maps = [structures_ds.map_table] if structures_ds.map_table else None
    pairing = pair_references(list(targets), list(references), len(reference_ds.ids_expanded), maps)

    table = None
    if "table" in selection_cfg:
        table, _ = load_table(selection_cfg["table"])
    ref_maps = [reference_ds.map_table] if reference_ds.map_table else None

    aligner = make_aligner()
    ref_cache: Dict[str, Chains] = {}
    rows, missing_rows = [], []
    for tid, tpath in targets.items():
        rid, sibling = pairing.get(tid, (None, None))
        if not rid:
            reason = "no reference structure matches this id"
            if sibling:
                reason += f" ('{sibling}' is only a sibling, another design's reference, and is not used)"
            print(f"WARNING: {tid} not matched: {reason}", file=sys.stderr)
            rows.append({"id": tid, "selection": None, "n_residues": None, "source_selection": None,
                         "unmapped": None, "identity": None, "reference.id": None})
            missing_rows.append({"id": tid, "removed_by": step_id, "kind": "unmatched", "cause": reason[:200]})
            continue
        try:
            if table is None:
                source = selection_cfg["literal"]
            else:
                source = lookup_table_value(table, rid, selection_cfg["column"], map_table_paths=ref_maps)
            source = "" if pd.isna(source) else str(source)
            # sele_to_list warns and returns [] instead of raising, so an unreadable cell would pass silently.
            if not source.strip():
                raise ValueError("the selection is empty")
            residues = sele_to_list(source)
            if not residues:
                raise ValueError(f"the selection '{source}' names no residue")
            if rid not in ref_cache:
                ref_cache[rid] = read_chains(references[rid])
            mapped, unmapped, identity = transfer(residues, ref_cache[rid],
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
