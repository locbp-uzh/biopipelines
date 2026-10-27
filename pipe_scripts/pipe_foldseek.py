#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Foldseek runner. Stages the structures under index names, runs `foldseek easy-cluster` (or `easy-multimercluster`) once per TM-score threshold, and writes one row per structure with its cluster number, representative and cluster count at every threshold."""

import argparse
import os
import re
import subprocess
import sys
from collections import defaultdict
from typing import Dict, List, Optional, Tuple

import pandas as pd
import yaml

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from biopipelines.biopipelines_io import (  # noqa: E402
    load_datastream,
    iterate_files,
    step_id_from_table_path,
    container_argv_prefix,
)


def threshold_tag(tm: float) -> str:
    """Column suffix of a TM-score threshold: 0.5 -> 'tm50'."""
    return f"tm{int(round(tm * 100)):02d}"


def stage_structures(files: List[Tuple[str, str]], stage_dir: str) -> Dict[int, str]:
    """Symlink each structure as <index>.pdb|.cif so Foldseek's entry names carry no id ambiguity."""
    os.makedirs(stage_dir, exist_ok=True)
    index_to_id: Dict[int, str] = {}
    for index, (sid, path) in enumerate(files):
        ext = ".cif" if path.lower().endswith((".cif", ".mmcif")) else ".pdb"
        link = os.path.join(stage_dir, f"{index}{ext}")
        if os.path.lexists(link):
            os.remove(link)
        os.symlink(os.path.abspath(path), link)
        index_to_id[index] = sid
    return index_to_id


def build_command(stage_dir: str, prefix: str, tmp_dir: str, tm: float, cfg: dict, container_prefix: str = "") -> List[str]:
    """The Foldseek argv for one threshold."""
    common = ["--alignment-type", str(cfg["alignment_type"]), "--cov-mode", "0",
              "--threads", str(cfg["threads"])]
    if cfg["multimer"]:
        cmd = ["foldseek", "easy-multimercluster", stage_dir, prefix, tmp_dir] + common + [
            "--multimer-tm-threshold", str(cfg["multimer_tm_threshold"]),
            "--chain-tm-threshold", str(tm),
            "--interface-lddt-threshold", str(cfg["interface_lddt_threshold"])]
    else:
        cmd = ["foldseek", "easy-cluster", stage_dir, prefix, tmp_dir] + common + [
            "--tmscore-threshold", str(tm)]
        if cfg["min_seq_id"] > 0:
            cmd += ["--min-seq-id", str(cfg["min_seq_id"])]
    return container_argv_prefix(container_prefix) + cmd


def run_foldseek(cmd: List[str]) -> None:
    """Run one Foldseek command; raise with the stderr tail on a non-zero exit."""
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        tail = (result.stderr.strip() or result.stdout.strip())[-300:]
        raise RuntimeError(f"foldseek exited {result.returncode}: {tail}")


def read_cluster_tsv(path: str) -> List[Tuple[str, str]]:
    """(representative, member) name pairs from a Foldseek cluster TSV."""
    pairs = []
    with open(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2 and parts[0].strip():
                pairs.append((parts[0].strip(), parts[1].strip()))
    return pairs


def entry_index(name: str) -> Optional[int]:
    """Staged index of a Foldseek entry name: '12', '12_A' and '12.pdb_A' all give 12."""
    match = re.match(r"\d+", name)
    return int(match.group()) if match else None


def assign_representatives(pairs: List[Tuple[str, str]], index_to_id: Dict[int, str],
                           multimer: bool) -> Tuple[Dict[str, str], Dict[str, str]]:
    """({id: representative id}, {id: why it has none}) for one threshold."""
    reps_of: Dict[int, List[str]] = defaultdict(list)
    for rep, member in pairs:
        index = entry_index(member)
        if index is not None:
            reps_of[index].append(rep)
    assigned: Dict[str, str] = {}
    problems: Dict[str, str] = {}
    for index, sid in index_to_id.items():
        reps = reps_of.get(index)
        if not reps:
            problems[sid] = "not reported by Foldseek (unreadable file or chain shorter than 4 residues)"
        elif len(reps) > 1 and not multimer:
            problems[sid] = (f"{len(reps)} chain or model entries were clustered separately; "
                             "use multimer=True for complexes")
        else:
            rep_id = index_to_id.get(entry_index(reps[0]))
            if rep_id is None:
                problems[sid] = f"representative {reps[0]!r} does not match a staged structure"
            else:
                assigned[sid] = rep_id
    return assigned, problems


def number_clusters(assigned: Dict[str, str], order: List[str]) -> Dict[str, int]:
    """{representative id: cluster number}, 1 for the largest cluster; ties follow input order."""
    size: Dict[str, int] = defaultdict(int)
    for rep in assigned.values():
        size[rep] += 1
    position = {sid: i for i, sid in enumerate(order)}
    ranked = sorted(size, key=lambda r: (-size[r], position[r]))
    return {rep: n for n, rep in enumerate(ranked, 1)}


def cluster_one_threshold(tm: float, stage_dir: str, index_to_id: Dict[int, str],
                          scratch_dir: str, cfg: dict, container_prefix: str):
    """Run Foldseek at one threshold; return ({id: rep}, {id: problem})."""
    work = os.path.join(scratch_dir, f"foldseek_{threshold_tag(tm)}")
    tmp_dir = os.path.join(work, "tmp")
    os.makedirs(tmp_dir, exist_ok=True)
    prefix = os.path.join(work, "res")
    run_foldseek(build_command(stage_dir, prefix, tmp_dir, tm, cfg, container_prefix))
    tsv = f"{prefix}_cluster.tsv"
    if not os.path.exists(tsv):
        raise RuntimeError(f"foldseek wrote no {os.path.basename(tsv)}")
    return assign_representatives(read_cluster_tsv(tsv), index_to_id, cfg["multimer"])


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--config-yaml", required=True)
    p.add_argument("--structures-json", required=True)
    p.add_argument("--scratch-dir", required=True)
    p.add_argument("--tm-clusters-csv", required=True)
    p.add_argument("--local-missing-csv", required=True)
    p.add_argument("--container-prefix", default="")
    args = p.parse_args()

    with open(args.config_yaml) as f:
        cfg = yaml.safe_load(f)
    thresholds = sorted(cfg["tm_thresholds"])

    ds = load_datastream(args.structures_json)
    step_id = step_id_from_table_path(args.local_missing_csv)
    files = list(iterate_files(ds))
    if not files:
        print("ERROR: no input entities to process (upstream stream is empty or its files are absent)",
              file=sys.stderr)
        sys.exit(1)

    stage_dir = os.path.join(args.scratch_dir, "foldseek_input")
    index_to_id = stage_structures(files, stage_dir)
    order = [sid for sid, _ in files]

    per_threshold: Dict[float, Tuple[Dict[str, str], Dict[str, int]]] = {}
    causes: Dict[str, List[str]] = defaultdict(list)
    failed_thresholds: List[float] = []
    for tm in thresholds:
        try:
            assigned, problems = cluster_one_threshold(
                tm, stage_dir, index_to_id, args.scratch_dir, cfg, args.container_prefix)
        except Exception as e:
            print(f"WARNING: Foldseek failed at TM {tm}: {e}", file=sys.stderr)
            failed_thresholds.append(tm)
            assigned, problems = {}, {sid: f"foldseek failed at TM {tm}: {e}" for sid in order}
        numbers = number_clusters(assigned, order)
        per_threshold[tm] = (assigned, numbers)
        for sid, cause in problems.items():
            causes[sid].append(f"TM {tm}: {cause}")
        print(f"  TM {tm}: {len(numbers)} cluster(s) over {len(assigned)}/{len(order)} structure(s)")

    columns = ["id"]
    for tm in thresholds:
        tag = threshold_tag(tm)
        columns += [f"cluster_{tag}", f"representative_{tag}", f"n_clusters_{tag}"]
    rows = []
    for sid in order:
        row = {"id": sid}
        for tm in thresholds:
            tag = threshold_tag(tm)
            assigned, numbers = per_threshold[tm]
            rep = assigned.get(sid)
            row[f"cluster_{tag}"] = numbers[rep] if rep is not None else None
            row[f"representative_{tag}"] = rep
            row[f"n_clusters_{tag}"] = len(numbers) if rep is not None else None
        rows.append(row)
    table = pd.DataFrame(rows, columns=columns)
    for col in columns:
        if col.startswith(("cluster_", "n_clusters_")):
            table[col] = table[col].astype("Int64")

    local_missing = [{"id": sid, "removed_by": step_id, "kind": "failure",
                      "cause": "; ".join(causes[sid])[:200]} for sid in order if sid in causes]

    for path in (args.tm_clusters_csv, args.local_missing_csv):
        os.makedirs(os.path.dirname(path), exist_ok=True)
    table.to_csv(args.tm_clusters_csv, index=False)
    pd.DataFrame(local_missing, columns=["id", "removed_by", "kind", "cause"]).to_csv(
        args.local_missing_csv, index=False)
    print(f"Clusters: {args.tm_clusters_csv} ({len(table)} rows)")
    print(f"Local failures: {len(local_missing)}")

    # A threshold Foldseek could not run leaves its columns NaN for every id, which is a failed step, not a per-id failure.
    if failed_thresholds:
        print(f"ERROR: Foldseek failed at TM {failed_thresholds}; their columns are NaN", file=sys.stderr)
        sys.exit(1)
    if len(local_missing) == len(order):
        print(f"ERROR: all {len(order)} attempted id(s) failed", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
