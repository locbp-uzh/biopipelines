# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Runtime half of StructureCluster.

Reduces each structure to a length-normalised C-alpha trace, clusters the traces
by distance-matrix RMSD with sphere exclusion, and ranks the clusters on
per-id metrics loaded from the tables the wrapper pointed at.

Everything here is numpy + pandas so it stays in the biopipelines base env; at
10,000 structures the cost is dominated by parsing, not by the clustering.
"""

import argparse
import json
import os
import sys
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import yaml

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from biopipelines.biopipelines_io import (  # noqa: E402
    load_datastream,
    iterate_files,
    step_id_from_table_path,
)
from biopipelines.pdb_parser import parse_pdb_file  # noqa: E402


# P-SEA-style C-alpha secondary-structure criteria (Labesse et al., CABIOS 1997):
# alpha and beta are separated cleanly by the i->i+2/3/4 C-alpha distances alone.
# Bands are deliberately generous -- these run on diffusion-model output, whose
# geometry is close to but not exactly crystallographic.
_HELIX_D2 = (5.0, 6.4)
_HELIX_D3 = (4.8, 6.4)
_HELIX_D4 = (5.2, 7.2)
_STRAND_D2 = (6.2, 7.6)
_STRAND_D3 = (9.2, 10.8)
_STRAND_D4 = (11.6, 13.6)

# C-alpha contact definition for relative contact order. 8 A between C-alphas
# with |i-j| >= 3 is the standard Plaxco convention re-expressed on C-alpha.
_CONTACT_CUTOFF = 8.0
_CONTACT_MIN_SEP = 3

_LABEL_FMT = "cluster_{:03d}"

_AA3 = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "MSE", "SEC", "PYL",
}


# --------------------------------------------------------------- descriptors

def _ca_trace(path: str, want_chain: str) -> Tuple[np.ndarray, str]:
    """C-alpha coordinates of one protein chain, in residue order.

    Returns (N x 3 array, chain id). Picks the longest protein chain when
    `want_chain` is empty -- in a de novo pipeline that is always the designed
    chain, and it keeps the ligand (a one-residue HETATM chain) out of the
    trace without the caller having to name anything.
    """
    atoms = parse_pdb_file(path)
    by_chain: Dict[str, List[Tuple[int, Tuple[float, float, float]]]] = {}
    for a in atoms:
        if a.atom_name != "CA" or a.res_name not in _AA3:
            continue
        by_chain.setdefault(a.chain, []).append((a.res_num, (a.x, a.y, a.z)))
    if not by_chain:
        raise ValueError("no protein C-alpha atoms found")

    if want_chain:
        if want_chain not in by_chain:
            raise ValueError(
                f"chain {want_chain!r} has no protein C-alpha atoms "
                f"(present: {','.join(sorted(by_chain))})")
        chain = want_chain
    else:
        chain = max(by_chain, key=lambda c: len(by_chain[c]))

    # Sort by residue number, then drop duplicates (altlocs / multi-model files
    # would otherwise inject the same residue twice and distort the trace).
    seen = set()
    ordered = []
    for res_num, xyz in sorted(by_chain[chain], key=lambda t: t[0]):
        if res_num in seen:
            continue
        seen.add(res_num)
        ordered.append(xyz)
    return np.asarray(ordered, dtype=np.float64), chain


def _resample(trace: np.ndarray, n_points: int) -> np.ndarray:
    """Linearly resample a C-alpha trace to a fixed number of points.

    Length normalisation happens here, and it is what lets a 90-residue and a
    140-residue design be compared at all. Interpolation is along the residue
    INDEX, not arc length: index spacing is what preserves topology (which part
    of the chain is where), while arc length would compress a compact region
    against an extended one and call two different folds similar.
    """
    n = len(trace)
    src = np.linspace(0.0, 1.0, n)
    dst = np.linspace(0.0, 1.0, n_points)
    out = np.empty((n_points, 3), dtype=np.float64)
    for k in range(3):
        out[:, k] = np.interp(dst, src, trace[:, k])
    return out


def _descriptor(points: np.ndarray) -> np.ndarray:
    """Flattened upper-triangle distance matrix -- the invariant fingerprint."""
    d = np.linalg.norm(points[:, None, :] - points[None, :, :], axis=2)
    iu = np.triu_indices(len(points), k=1)
    return d[iu].astype(np.float32)


def _secondary_structure(trace: np.ndarray) -> Tuple[float, float, float]:
    """Helix / strand / coil fractions from C-alpha geometry (P-SEA-style).

    An estimate, not DSSP: no hydrogen bonds are involved. It is here because
    it separates an all-alpha bundle from an alpha/beta sandwich at zero extra
    dependency, which is the topology question this tool exists to answer.

    Each qualifying 5-residue window votes for its state, and every residue is
    then labelled by majority vote. Voting rather than painting is what keeps
    the three fractions a partition: a residue covered by both a helical and a
    strand window (which happens at element boundaries, since the windows
    overlap even though their distance bands do not) would otherwise be counted
    twice and the fractions would sum past 1.
    """
    n = len(trace)
    if n < 5:
        return (float("nan"),) * 3
    d = lambda k: np.linalg.norm(trace[k:] - trace[:-k], axis=1)  # noqa: E731
    d2, d3, d4 = d(2), d(3), d(4)

    h_votes = np.zeros(n, dtype=np.int32)
    s_votes = np.zeros(n, dtype=np.int32)
    for i in range(n - 4):
        if (_HELIX_D2[0] <= d2[i] <= _HELIX_D2[1]
                and _HELIX_D3[0] <= d3[i] <= _HELIX_D3[1]
                and _HELIX_D4[0] <= d4[i] <= _HELIX_D4[1]):
            h_votes[i:i + 5] += 1
        elif (_STRAND_D2[0] <= d2[i] <= _STRAND_D2[1]
                and _STRAND_D3[0] <= d3[i] <= _STRAND_D3[1]
                and _STRAND_D4[0] <= d4[i] <= _STRAND_D4[1]):
            s_votes[i:i + 5] += 1

    # Ties (equal non-zero votes) fall to coil: an ambiguous residue is exactly
    # what "coil" should absorb, and it keeps the label deterministic.
    helix = h_votes > s_votes
    strand = s_votes > h_votes
    h = float(helix.sum()) / n
    s = float(strand.sum()) / n
    return h, s, 1.0 - h - s


def _relative_contact_order(trace: np.ndarray) -> float:
    """Mean sequence separation of C-alpha contacts, divided by chain length."""
    n = len(trace)
    if n < _CONTACT_MIN_SEP + 1:
        return float("nan")
    d = np.linalg.norm(trace[:, None, :] - trace[None, :, :], axis=2)
    i, j = np.triu_indices(n, k=_CONTACT_MIN_SEP)
    hit = d[i, j] <= _CONTACT_CUTOFF
    if not hit.any():
        return float("nan")
    return float(np.abs(j[hit] - i[hit]).mean() / n)


def _radius_of_gyration(trace: np.ndarray) -> float:
    return float(np.sqrt(((trace - trace.mean(0)) ** 2).sum(1).mean()))


# ---------------------------------------------------------------- clustering

def _d0(lengths: np.ndarray) -> np.ndarray:
    """TM-score length scale, floored so short chains stay comparable."""
    return np.maximum(0.5, 1.24 * np.cbrt(np.maximum(lengths - 15.0, 1e-9)) - 1.8)


def _similarity(ref_desc: np.ndarray, ref_len: float,
                others: np.ndarray, other_lens: np.ndarray) -> np.ndarray:
    """Bounded fold similarity of one reference against many, vectorised.

    dRMSD -> 1/(1+(dRMSD/d0)^2). Computed as a batch so the whole sphere-
    exclusion pass is O(N * n_clusters) numpy work rather than a Python loop.
    """
    diff = others - ref_desc[None, :]
    drmsd = np.sqrt(np.einsum("ij,ij->i", diff, diff) / diff.shape[1])
    d0 = _d0((other_lens + ref_len) / 2.0)
    return 1.0 / (1.0 + (drmsd / d0) ** 2)


def _sphere_exclusion(desc: np.ndarray, lengths: np.ndarray, order: np.ndarray,
                      threshold: float) -> Tuple[np.ndarray, np.ndarray, List[int]]:
    """Greedy leader clustering in ranking order.

    Returns (cluster index per row, similarity to its representative,
    representative row indices in discovery order). Visiting in ranking order is
    what makes each representative the best-scoring member of its cluster.
    """
    n = len(order)
    assigned = np.full(n, -1, dtype=np.int64)
    sim_to_rep = np.zeros(n, dtype=np.float64)
    reps: List[int] = []

    unassigned = np.ones(n, dtype=bool)
    for idx in order:
        if not unassigned[idx]:
            continue
        c = len(reps)
        reps.append(int(idx))
        pool = np.flatnonzero(unassigned)
        sims = _similarity(desc[idx], lengths[idx], desc[pool], lengths[pool])
        take = pool[sims >= threshold]
        # The representative always captures itself (similarity 1.0), so `take`
        # is never empty and the loop always shrinks `unassigned`.
        assigned[take] = c
        sim_to_rep[take] = sims[sims >= threshold]
        unassigned[take] = False
    return assigned, sim_to_rep, reps


# --------------------------------------------------------------------- main

def _load_metrics(paths: List[str], columns: List[str]) -> pd.DataFrame:
    """Join the metric tables on `id`, keeping only the requested columns.

    A missing table or a missing column is a NaN column, never an exception: a
    campaign that lost one scorer should still get its clusters, with the gap
    visible in the table instead of the run dying at the last step.
    """
    merged: Optional[pd.DataFrame] = None
    for p in paths:
        if not p or not os.path.exists(p):
            print(f"WARNING: metrics table not found, skipping: {p}", file=sys.stderr)
            continue
        df = pd.read_csv(p)
        if "id" not in df.columns:
            print(f"WARNING: metrics table has no 'id' column, skipping: {p}", file=sys.stderr)
            continue
        keep = ["id"] + [c for c in columns if c in df.columns and c != "id"]
        df = df[keep].drop_duplicates(subset="id")
        merged = df if merged is None else merged.merge(df, on="id", how="outer")

    if merged is None:
        merged = pd.DataFrame({"id": pd.Series(dtype=str)})
    for c in columns:
        if c not in merged.columns:
            print(f"WARNING: metric column {c!r} not present in any metrics table; "
                  f"it will be NaN", file=sys.stderr)
            merged[c] = np.nan
    return merged


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--config-yaml", required=True)
    p.add_argument("--structures-json", required=True)
    p.add_argument("--metrics-json", required=True)
    p.add_argument("--assignments-csv", required=True)
    p.add_argument("--clusters-csv", required=True)
    p.add_argument("--local-missing-csv", required=True)
    args = p.parse_args()

    with open(args.config_yaml) as f:
        cfg = yaml.safe_load(f) or {}
    threshold = float(cfg.get("threshold", 0.5))
    n_points = int(cfg.get("n_points", 64))
    want_chain = cfg.get("chain", "") or ""
    min_residues = int(cfg.get("min_residues", 20))
    max_structures = int(cfg.get("max_structures", 0))
    rank_by: List[str] = list(cfg.get("rank_by", ["plddt", "iptm"]))
    ascending: List[bool] = list(cfg.get("ascending", [False] * len(rank_by)))

    with open(args.metrics_json) as f:
        metrics_paths = (json.load(f) or {}).get("tables", [])

    ds = load_datastream(args.structures_json)
    step_id = step_id_from_table_path(args.local_missing_csv)

    ids: List[str] = []
    rows: List[dict] = []
    desc_list: List[np.ndarray] = []
    lengths: List[float] = []
    local_missing: List[dict] = []
    attempted = 0
    failures = 0

    for sid, path in iterate_files(ds):
        attempted += 1
        if max_structures and len(ids) >= max_structures:
            local_missing.append({"id": sid, "removed_by": step_id, "kind": "filter",
                                  "cause": f"max_structures={max_structures} reached"})
            rows.append({"id": sid, "n_residues": np.nan, "chain": ""})
            continue
        try:
            trace, chain = _ca_trace(path, want_chain)
            n_res = len(trace)
            row = {
                "id": sid,
                "n_residues": n_res,
                "chain": chain,
                "radius_of_gyration": _radius_of_gyration(trace),
                "relative_contact_order": _relative_contact_order(trace),
            }
            h, s, c = _secondary_structure(trace)
            row["helix_frac"], row["strand_frac"], row["coil_frac"] = h, s, c

            if n_res < min_residues:
                # Deliberate drop from the CLUSTERING; the row survives in the
                # table, because a table records what was processed.
                rows.append(row)
                local_missing.append({"id": sid, "removed_by": step_id, "kind": "filter",
                                      "cause": f"{n_res} residues < min_residues={min_residues}"})
                continue

            desc_list.append(_descriptor(_resample(trace, n_points)))
            lengths.append(float(n_res))
            ids.append(sid)
            rows.append(row)
        except Exception as e:  # noqa: BLE001
            print(f"WARNING: {sid} failed: {e}", file=sys.stderr)
            failures += 1
            rows.append({"id": sid, "n_residues": np.nan, "chain": ""})
            local_missing.append({"id": sid, "removed_by": step_id, "kind": "failure",
                                  "cause": str(e)[:200]})

    assign = pd.DataFrame(rows)
    metrics = _load_metrics(metrics_paths, rank_by)
    assign = assign.merge(metrics, on="id", how="left")
    for c in rank_by:
        if c not in assign.columns:
            assign[c] = np.nan

    n_clustered = len(ids)
    print(f"Parsed {n_clustered} clusterable structure(s) of {attempted} attempted "
          f"({failures} failed)")

    if n_clustered:
        desc = np.vstack(desc_list)
        lens = np.asarray(lengths)

        # Visit order = ranking order, so representatives are the best models.
        # Without metrics this falls back to id order, which is arbitrary but
        # deterministic -- reproducibility matters more than the choice here.
        sub = assign.set_index("id").loc[ids]
        if metrics_paths and any(sub[c].notna().any() for c in rank_by):
            order_df = sub.sort_values(
                by=rank_by, ascending=ascending[:len(rank_by)],
                na_position="last", kind="mergesort")
            order = np.array([ids.index(i) for i in order_df.index], dtype=np.int64)
        else:
            order = np.argsort(np.array(ids, dtype=object), kind="mergesort")

        cluster_idx, sim_to_rep, reps = _sphere_exclusion(desc, lens, order, threshold)
        print(f"Formed {len(reps)} cluster(s) at similarity >= {threshold}")

        # ---- medoid: the most CENTRAL member, which is not the representative.
        # Sphere exclusion picks the best-SCORING member as representative, so it
        # can sit at the edge of its own cluster. The medoid maximises mean
        # similarity to the other members, i.e. it is the most typical fold of
        # the family -- the one to look at when asking "what does this topology
        # look like", where the representative answers "which of these should I
        # take forward". A singleton is its own medoid.
        medoid_of = {}
        for c in range(len(reps)):
            members = np.flatnonzero(cluster_idx == c)
            if len(members) == 1:
                medoid_of[c] = int(members[0])
                continue
            sub = desc[members]
            sublen = lens[members]
            # mean similarity of each member to the others (excluding itself,
            # whose similarity is 1.0 and would bias short clusters)
            tot = np.zeros(len(members))
            for k in range(len(members)):
                s = _similarity(sub[k], sublen[k], sub, sublen)
                tot[k] = (s.sum() - 1.0) / (len(members) - 1)
            medoid_of[c] = int(members[int(np.argmax(tot))])

        per_id = pd.DataFrame({
            "id": ids,
            "_cluster_idx": cluster_idx,
            "similarity_to_representative": sim_to_rep,
        })
        per_id["is_representative"] = per_id["id"].isin([ids[r] for r in reps])
        per_id["is_medoid"] = per_id["id"].isin([ids[m] for m in medoid_of.values()])
        assign = assign.merge(per_id, on="id", how="left")

        # ------------------------------------------------- per-cluster table
        joined = assign[assign["_cluster_idx"].notna()].copy()
        joined["_cluster_idx"] = joined["_cluster_idx"].astype(int)
        agg_rows = []
        for ci, grp in joined.groupby("_cluster_idx"):
            rec = {
                "_cluster_idx": ci,
                "size": len(grp),
                "representative": ids[reps[ci]],
                "medoid": ids[medoid_of[ci]],
                "mean_n_residues": grp["n_residues"].mean(),
                "mean_helix_frac": grp["helix_frac"].mean(),
                "mean_strand_frac": grp["strand_frac"].mean(),
                "mean_coil_frac": grp["coil_frac"].mean(),
                "mean_relative_contact_order": grp["relative_contact_order"].mean(),
                "mean_radius_of_gyration": grp["radius_of_gyration"].mean(),
            }
            for m in rank_by:
                col = grp[m] if m in grp.columns else pd.Series(dtype=float)
                rec[f"mean_{m}"] = col.mean()
                rec[f"median_{m}"] = col.median()
                rec[f"min_{m}"] = col.min()
                rec[f"max_{m}"] = col.max()
            agg_rows.append(rec)
        clusters = pd.DataFrame(agg_rows)

        # Rank clusters, then LABEL them -- so cluster_001 is the best cluster,
        # not the first one discovered. Size breaks ties, and is the sole key
        # when no usable metric came in.
        sort_cols = [f"mean_{m}" for m in rank_by if f"mean_{m}" in clusters.columns]
        usable = [c for c in sort_cols if clusters[c].notna().any()]
        if usable:
            asc = [ascending[sort_cols.index(c)] for c in usable]
            clusters = clusters.sort_values(by=usable + ["size"],
                                            ascending=asc + [False],
                                            na_position="last", kind="mergesort")
        else:
            clusters = clusters.sort_values(by="size", ascending=False, kind="mergesort")

        clusters = clusters.reset_index(drop=True)
        clusters["cluster_rank"] = np.arange(1, len(clusters) + 1)
        clusters["cluster"] = [_LABEL_FMT.format(r) for r in clusters["cluster_rank"]]
        clusters["fraction"] = clusters["size"] / max(1, n_clustered)

        label_of = dict(zip(clusters["_cluster_idx"], clusters["cluster"]))
        rank_of = dict(zip(clusters["_cluster_idx"], clusters["cluster_rank"]))
        assign["cluster"] = assign["_cluster_idx"].map(label_of)
        assign["cluster_rank"] = assign["_cluster_idx"].map(rank_of)
        clusters = clusters.drop(columns=["_cluster_idx"])
        assign = assign.drop(columns=["_cluster_idx"])
    else:
        for c in ("cluster", "cluster_rank", "similarity_to_representative"):
            assign[c] = np.nan
        assign["is_representative"] = False
        assign["is_medoid"] = False
        clusters = pd.DataFrame()

    # A structure that never reached the clustering is definitively not a
    # representative. Filling here keeps the column a clean bool instead of an
    # object column of True/False/NaN, so the documented downstream selection
    # -- Panda.filter("is_representative == True"), pool=folds -- matches rows
    # rather than silently returning none.
    for _c in ("is_representative", "is_medoid"):
        assign[_c] = assign[_c].fillna(False).astype(bool)

    # ------------------------------------------------------------- write out
    assign_cols = ["id", "cluster", "cluster_rank", "is_representative", "is_medoid",
                   "similarity_to_representative", "n_residues", "chain",
                   "radius_of_gyration", "helix_frac", "strand_frac", "coil_frac",
                   "relative_contact_order"] + rank_by
    cluster_cols = ["cluster", "cluster_rank", "size", "fraction", "representative", "medoid"]
    for m in rank_by:
        cluster_cols += [f"mean_{m}", f"median_{m}", f"min_{m}", f"max_{m}"]
    cluster_cols += ["mean_n_residues", "mean_helix_frac", "mean_strand_frac",
                     "mean_coil_frac", "mean_relative_contact_order",
                     "mean_radius_of_gyration"]

    for path in (args.assignments_csv, args.clusters_csv, args.local_missing_csv):
        os.makedirs(os.path.dirname(path), exist_ok=True)
    assign.reindex(columns=assign_cols).to_csv(args.assignments_csv, index=False)
    clusters.reindex(columns=cluster_cols).to_csv(args.clusters_csv, index=False)
    pd.DataFrame(local_missing, columns=["id", "removed_by", "kind", "cause"]).to_csv(
        args.local_missing_csv, index=False)

    print(f"Assignments: {args.assignments_csv} ({len(assign)} rows)")
    print(f"Clusters:    {args.clusters_csv} ({len(clusters)} rows)")

    if not attempted:
        print("ERROR: no input entities to process (upstream stream is empty or its "
              "files are absent)", file=sys.stderr)
        sys.exit(1)
    if failures == attempted:
        print(f"ERROR: all {attempted} attempted id(s) failed", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
