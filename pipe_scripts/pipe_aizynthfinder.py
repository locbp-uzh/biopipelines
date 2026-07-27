#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Runtime helper for AiZynthFinder.

Collects target SMILES from the compounds stream, runs one batched `aizynthcli`
search over them, then unpacks the returned route trees into a routes stream
(one JSON per route), a precursors compounds stream, and the summary tables.
"""

import argparse
import copy
import gzip
import json
import os
import subprocess
import sys

import pandas as pd
import yaml

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from biopipelines.biopipelines_io import (  # noqa: E402
    load_datastream,
    iterate_values,
    step_id_from_table_path,
    container_argv_prefix,
)

TARGET_COLUMNS = [
    "id", "smiles", "is_solved", "top_score", "number_of_routes",
    "number_of_returned_routes", "number_of_solved_routes", "number_of_nodes",
    "search_time", "first_solution_time", "first_solution_iteration",
]

ROUTE_COLUMNS = [
    "id", "compounds.id", "route", "score", "is_solved", "number_of_steps",
    "number_of_precursors", "number_of_precursors_in_stock",
]

PRECURSOR_COLUMNS = [
    "id", "format", "smiles", "ccd", "compounds.id", "route", "precursor", "in_stock",
]


def _write_search_config(cfg):
    """Copy the installed config and overlay only the limits the user set.

    The models and stock stay pointed wherever download_public_data put them;
    this run's search bounds are layered on top under the keys aizynthfinder
    documents (`search.*` for the cutoffs, `post_processing.*` for route counts).
    """
    base_config = cfg["base_config"]
    if not os.path.isfile(base_config):
        raise FileNotFoundError(
            f"AiZynthFinder config not found: {base_config}. "
            f"Run AiZynthFinder.install() to download the public models, or pass config=."
        )
    with open(base_config) as f:
        config = yaml.safe_load(f) or {}
    config = copy.deepcopy(config)

    search = config.setdefault("search", {})
    if cfg.get("time_limit") is not None:
        search["time_limit"] = cfg["time_limit"]
    if cfg.get("iteration_limit") is not None:
        search["iteration_limit"] = cfg["iteration_limit"]

    post = config.setdefault("post_processing", {})
    post["min_routes"] = cfg["min_routes"]
    post["max_routes"] = cfg["max_routes"]

    out_path = cfg["search_config_yaml"]
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as f:
        yaml.safe_dump(config, f, sort_keys=False)
    return out_path


def _run_cli(cfg, search_config, n_targets, container_prefix):
    cmd = container_argv_prefix(container_prefix) + [
        "aizynthcli",
        "--config", search_config,
        "--smiles", cfg["targets_smi"],
        "--output", cfg["cli_output"],
    ]
    for flag, key in (("--stocks", "stocks"),
                      ("--policy", "expansion_policy"),
                      ("--filter", "filter_policy")):
        if cfg.get(key):
            cmd += [flag] + list(cfg[key])
    # aizynthcli splits the target list across processes; one process per target
    # at most, or it spawns workers with nothing to do.
    if cfg.get("nproc"):
        cmd += ["--nproc", str(min(int(cfg["nproc"]), n_targets))]

    os.makedirs(os.path.dirname(cfg["cli_output"]), exist_ok=True)
    print(f"Running: {' '.join(cmd)}", flush=True)
    subprocess.run(cmd, check=True)


def _read_cli_output(path):
    """aizynthcli writes JSON (.json.gz) or HDF5 depending on the extension.

    The JSON is pandas' orient="table" envelope ({"schema":..., "data":[row,...]}),
    so the rows come from "data" rather than from the mapping itself.
    """
    if path.endswith((".json.gz", ".json")):
        opener = gzip.open if path.endswith(".gz") else open
        with opener(path, "rt") as f:
            payload = json.load(f)
        return payload["data"] if isinstance(payload, dict) and "data" in payload else payload
    return pd.read_hdf(path, "table").to_dict("records")


def _leaf_molecules(tree):
    """Walk a reaction tree and yield its leaf molecule nodes (the precursors)."""
    if not isinstance(tree, dict):
        return
    children = tree.get("children") or []
    if tree.get("type") == "mol" and not children:
        yield tree
        return
    for child in children:
        yield from _leaf_molecules(child)


def _count_reactions(tree):
    if not isinstance(tree, dict):
        return 0
    n = 1 if tree.get("type") == "reaction" else 0
    for child in tree.get("children") or []:
        n += _count_reactions(child)
    return n


def _as_bool(value):
    """Normalize the mixed bool/str/NaN the tree dicts carry for in_stock.

    NaN is falsy here by fiat: bool(float("nan")) is True, so an absent value on
    the HDF5 path would otherwise read as "in stock" and mark a route solved.
    """
    if value is None:
        return False
    if isinstance(value, float) and value != value:
        return False
    if isinstance(value, str):
        return value.strip().lower() in ("true", "1", "yes")
    return bool(value)


def main():
    p = argparse.ArgumentParser(description="Run AiZynthFinder over a compounds stream")
    p.add_argument("--config", required=True)
    p.add_argument("--container-prefix", default="")
    args = p.parse_args()

    with open(args.config) as f:
        cfg = json.load(f)

    step_id = step_id_from_table_path(cfg["local_missing_csv"])
    compounds_ds = load_datastream(cfg["compounds_json"])

    # Collect targets in stream order. A compound with no SMILES cannot be a
    # retrosynthesis target at all, so it is a deliberate drop, not a failure.
    targets, local_missing_rows = [], []
    for comp_id, values in iterate_values(compounds_ds, columns=["smiles"]):
        smiles = values.get("smiles")
        if smiles is None or (isinstance(smiles, float) and pd.isna(smiles)) or not str(smiles).strip():
            print(f"WARNING: skipping {comp_id} (no smiles)", file=sys.stderr)
            local_missing_rows.append({
                "id": comp_id, "removed_by": step_id, "kind": "filter",
                "cause": "no smiles in compounds stream",
            })
            continue
        targets.append({"id": comp_id, "smiles": str(smiles).strip()})

    if not targets:
        print("ERROR: no compounds with valid SMILES to search "
              "(upstream stream is empty or carries no smiles)", file=sys.stderr)
        sys.exit(1)

    os.makedirs(os.path.dirname(cfg["targets_smi"]), exist_ok=True)
    with open(cfg["targets_smi"], "w") as f:
        for target in targets:
            f.write(f"{target['smiles']}\n")

    search_config = _write_search_config(cfg)
    _run_cli(cfg, search_config, len(targets), args.container_prefix)

    results = _read_cli_output(cfg["cli_output"])
    if len(results) != len(targets):
        raise RuntimeError(
            f"aizynthcli returned {len(results)} rows for {len(targets)} targets; "
            f"cannot align results to input ids"
        )

    os.makedirs(cfg["routes_dir"], exist_ok=True)
    target_rows, route_rows, precursor_rows, routes_map_rows = [], [], [], []
    solved_targets = 0

    # Rows come back in the order the SMILES were written, which is stream order.
    for target, result in zip(targets, results):
        tid = target["id"]
        trees = result.get("trees") or []

        # number_of_routes is what the search FOUND; len(trees) is what
        # post_processing.max_routes let through and thus what is emitted here.
        target_rows.append({
            "id": tid,
            "smiles": target["smiles"],
            "is_solved": _as_bool(result.get("is_solved")),
            "top_score": result.get("top_score"),
            "number_of_routes": result.get("number_of_routes"),
            "number_of_returned_routes": len(trees),
            "number_of_solved_routes": result.get("number_of_solved_routes"),
            "number_of_nodes": result.get("number_of_nodes"),
            "search_time": result.get("search_time"),
            "first_solution_time": result.get("first_solution_time"),
            "first_solution_iteration": result.get("first_solution_iteration"),
        })

        if not trees:
            local_missing_rows.append({
                "id": tid, "removed_by": step_id, "kind": "filter",
                "cause": "no retrosynthetic route found",
            })
            continue
        solved_targets += 1

        for route_idx, tree in enumerate(trees, start=1):
            route_id = f"{tid}_{route_idx}"
            leaves = list(_leaf_molecules(tree))
            in_stock_flags = [_as_bool(leaf.get("in_stock")) for leaf in leaves]

            route_path = os.path.join(cfg["routes_dir"], f"{route_id}.json")
            with open(route_path, "w") as f:
                json.dump(tree, f, indent=2)
            routes_map_rows.append({
                "id": route_id, "file": route_path, "compounds.id": tid,
            })

            # Each tree carries its own scores dict, which is also the
            # authoritative source for the reaction/precursor counts.
            tree_scores = tree.get("scores") or {}
            route_rows.append({
                "id": route_id,
                "compounds.id": tid,
                "route": route_idx,
                "score": tree_scores.get("state score"),
                "is_solved": all(in_stock_flags) if leaves else False,
                "number_of_steps": tree_scores.get("number of reactions", _count_reactions(tree)),
                "number_of_precursors": tree_scores.get("number of pre-cursors", len(leaves)),
                "number_of_precursors_in_stock": tree_scores.get(
                    "number of pre-cursors in stock", sum(in_stock_flags)),
            })

            for leaf_idx, (leaf, in_stock) in enumerate(zip(leaves, in_stock_flags), start=1):
                precursor_rows.append({
                    "id": f"{route_id}_{leaf_idx}",
                    "format": "smiles",
                    "smiles": leaf.get("smiles"),
                    "ccd": "",
                    "compounds.id": tid,
                    "route": route_idx,
                    "precursor": leaf_idx,
                    "in_stock": in_stock,
                })

    for path in (cfg["targets_csv"], cfg["routes_csv"], cfg["precursors_csv"],
                 cfg["routes_map_csv"], cfg["local_missing_csv"]):
        os.makedirs(os.path.dirname(path), exist_ok=True)

    pd.DataFrame(target_rows, columns=TARGET_COLUMNS).to_csv(cfg["targets_csv"], index=False)
    pd.DataFrame(route_rows, columns=ROUTE_COLUMNS).to_csv(cfg["routes_csv"], index=False)
    # ccd is empty for every precursor (SMILES have no CCD identity); keep it an
    # empty string so downstream compounds consumers don't read a literal "nan".
    precursors_df = pd.DataFrame(precursor_rows, columns=PRECURSOR_COLUMNS)
    precursors_df["ccd"] = precursors_df["ccd"].fillna("")
    precursors_df.to_csv(cfg["precursors_csv"], index=False)
    pd.DataFrame(routes_map_rows, columns=["id", "file", "compounds.id"]).to_csv(
        cfg["routes_map_csv"], index=False)
    pd.DataFrame(local_missing_rows, columns=["id", "removed_by", "kind", "cause"]).to_csv(
        cfg["local_missing_csv"], index=False)

    print(f"Targets: {cfg['targets_csv']} ({len(target_rows)} rows, {solved_targets} with routes)")
    print(f"Routes: {cfg['routes_csv']} ({len(route_rows)} rows)")
    print(f"Precursors: {cfg['precursors_csv']} ({len(precursor_rows)} rows)")

    # Every target searched but none produced a route: the tool worked, it just
    # found nothing, so this is a filter, not a failure.
    if not route_rows:
        print(f"WARNING: no routes found for any of {len(targets)} target(s)", file=sys.stderr)


if __name__ == "__main__":
    main()
