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

FROZEN_FILTER_KEY = "frozen_substructure"

# Filter `type` values aizynthfinder resolves to its QuickKerasFilter.
QUICK_FILTER_TYPES = ("quick-filter", "feasibility", "quick_keras_filter", "QuickKerasFilter")


def _apply_filter_cutoff(filters, cutoff):
    """Set the feasibility cutoff on every quick-filter entry of the filter section."""
    found = False
    for key, spec in list(filters.items()):
        if isinstance(spec, dict):
            if spec.get("type") is not None and spec["type"] not in QUICK_FILTER_TYPES:
                continue
            spec["filter_cutoff"] = cutoff
        else:
            filters[key] = {"model": spec, "filter_cutoff": cutoff}
        found = True
    if not found:
        raise ValueError(
            "filter_cutoff is set but the config has no quick-filter in its `filter:` section")


def build_search_config(config, cfg):
    """Overlay the options the user set onto a copy of the installed config.

    The models and stock stay pointed wherever download_public_data put them; each option lands under the key aizynthfinder documents for it (`search.*` for the algorithm and cutoffs, `post_processing.*` for route selection, `scorer`, `filter` and `stock` for those collections).
    """
    config = copy.deepcopy(config or {})

    search = config.setdefault("search", {})
    for key in ("time_limit", "iteration_limit", "max_transforms"):
        if cfg.get(key) is not None:
            search[key] = cfg[key]
    if cfg.get("search_algorithm"):
        search["algorithm"] = cfg["search_algorithm"]

    post = config.setdefault("post_processing", {})
    post["min_routes"] = cfg["min_routes"]
    post["max_routes"] = cfg["max_routes"]
    if cfg.get("route_scorer"):
        post["route_scorers"] = [cfg["route_scorer"]]

    # A scorer class named with no arguments is loaded as is; one the user configured keeps its arguments.
    scorer = config.get("scorer") or {}
    for item in cfg.get("extra_scorers") or []:
        if item.get("class"):
            scorer.setdefault(item["class"], None)
    if scorer:
        config["scorer"] = scorer

    filters = config.get("filter") or {}
    if cfg.get("filter_cutoff") is not None:
        _apply_filter_cutoff(filters, cfg["filter_cutoff"])
    if cfg.get("frozen_substructures"):
        if FROZEN_FILTER_KEY in filters:
            raise ValueError(
                f"the config already defines a `{FROZEN_FILTER_KEY}` filter; rename it, "
                f"or drop frozen_substructures")
        filters[FROZEN_FILTER_KEY] = {
            "type": "frozen_substructure",
            "smarts_list": list(cfg["frozen_substructures"]),
        }
    if filters:
        config["filter"] = filters

    if cfg.get("stock_stop_criteria"):
        stock = config.get("stock") or {}
        stock["stop_criteria"] = {**(stock.get("stop_criteria") or {}), **cfg["stock_stop_criteria"]}
        config["stock"] = stock

    return config


def _write_search_config(cfg):
    base_config = cfg["base_config"]
    if not os.path.isfile(base_config):
        raise FileNotFoundError(
            f"AiZynthFinder config not found: {base_config}. "
            f"Run AiZynthFinder.install() to download the public models, or pass config=."
        )
    with open(base_config) as f:
        config = yaml.safe_load(f) or {}
    config = build_search_config(config, cfg)

    out_path = cfg["search_config_yaml"]
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as f:
        yaml.safe_dump(config, f, sort_keys=False)
    return out_path


def _filter_selection(cfg):
    """Filters passed to --filter; empty leaves aizynthcli's default of all loaded."""
    names = list(cfg.get("filter_policy") or [])
    if names and cfg.get("frozen_substructures") and FROZEN_FILTER_KEY not in names:
        names.append(FROZEN_FILTER_KEY)
    return names


def _run_cli(cfg, search_config, n_targets, container_prefix):
    cmd = container_argv_prefix(container_prefix) + [
        "aizynthcli",
        "--config", search_config,
        "--smiles", cfg["targets_smi"],
        "--output", cfg["cli_output"],
    ]
    for flag, names in (("--stocks", cfg.get("stocks")),
                        ("--policy", cfg.get("expansion_policy")),
                        ("--filter", _filter_selection(cfg))):
        if names:
            cmd += [flag] + list(names)
    # aizynthcli splits the target list across processes; one process per target
    # at most, or it spawns workers with nothing to do.
    if cfg.get("nproc"):
        cmd += ["--nproc", str(min(int(cfg["nproc"]), n_targets))]

    os.makedirs(os.path.dirname(cfg["cli_output"]), exist_ok=True)
    print(f"Running: {' '.join(cmd)}", flush=True)
    # aizynthcli writes aizynthcli<N>.log into its working directory when it runs several processes.
    subprocess.run(cmd, check=True, cwd=os.path.dirname(cfg["cli_output"]))


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


def _canonical_smiles(smiles):
    """The string aizynthcli reports as `target`, or None if it cannot set the search up.

    Mirrors aizynthfinder's Molecule.sanitize: parse unsanitized, sanitize, then write the canonical SMILES.
    """
    from rdkit import Chem

    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    if mol is None:
        return None
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        return None
    return Chem.MolToSmiles(mol)


def _validate_smarts(patterns):
    from rdkit import Chem

    for pattern in patterns or []:
        if Chem.MolFromSmarts(pattern) is None:
            raise ValueError(f"frozen_substructures: {pattern!r} is not a valid SMARTS pattern")


def _align_results(targets, results):
    """Pair each output row with its target.

    aizynthcli keeps the order of the SMILES file but skips a target it cannot set up, so the row count can fall short. Rows carry the canonical SMILES as `target`; walking both lists in order matches each row and leaves the skipped targets unpaired.
    """
    if len(results) == len(targets):
        return list(zip(targets, results)), []
    if len(results) > len(targets):
        raise RuntimeError(
            f"aizynthcli returned {len(results)} rows for {len(targets)} targets; "
            f"cannot align results to input ids")

    pairs, skipped, position = [], [], 0
    for result in results:
        key = result.get("target")
        start = position
        while position < len(targets) and targets[position]["canonical"] != key:
            position += 1
        if position == len(targets):
            raise RuntimeError(
                f"aizynthcli returned {len(results)} rows for {len(targets)} targets and the "
                f"row for {key!r} matches none of the remaining targets; cannot align results to input ids")
        skipped.extend(targets[start:position])
        pairs.append((targets[position], result))
        position += 1
    skipped.extend(targets[position:])
    return pairs, skipped


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


def _unavailable_target_row(target):
    """A targets row for an id that entered but was never searched: empty, not dropped."""
    row = {column: None for column in TARGET_COLUMNS}
    row["id"] = target["id"]
    row["smiles"] = target["smiles"]
    return row


def _convert_result(target, result, extra_scorers, routes_dir, step_id):
    """Turn one aizynthcli output row into its table rows and route files."""
    tid = target["id"]
    trees = result.get("trees") or []
    out = {"target": None, "routes": [], "precursors": [], "routes_map": [], "missing": []}

    # number_of_routes is what the search FOUND; len(trees) is what post_processing.max_routes let through.
    out["target"] = {
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
    }

    if not trees:
        out["missing"].append({
            "id": tid, "removed_by": step_id, "kind": "filter",
            "cause": "no retrosynthetic route found",
        })
        return out

    for route_idx, tree in enumerate(trees, start=1):
        route_id = f"{tid}_{route_idx}"
        leaves = list(_leaf_molecules(tree))
        in_stock_flags = [_as_bool(leaf.get("in_stock")) for leaf in leaves]

        route_path = os.path.join(routes_dir, f"{route_id}.json")
        with open(route_path, "w") as f:
            json.dump(tree, f, indent=2)
        out["routes_map"].append({
            "id": route_id, "file": route_path, "compounds.id": tid,
        })

        # Each tree carries its own scores dict, keyed by scorer name, and it is the authoritative source for the counts.
        tree_scores = tree.get("scores") or {}
        route_row = {
            "id": route_id,
            "compounds.id": tid,
            "route": route_idx,
            "score": tree_scores.get("state score"),
            "is_solved": all(in_stock_flags) if leaves else False,
            "number_of_steps": tree_scores.get("number of reactions", _count_reactions(tree)),
            "number_of_precursors": tree_scores.get("number of pre-cursors", len(leaves)),
            "number_of_precursors_in_stock": tree_scores.get(
                "number of pre-cursors in stock", sum(in_stock_flags)),
        }
        for item in extra_scorers:
            route_row[item["column"]] = tree_scores.get(item["name"])
        out["routes"].append(route_row)

        for leaf_idx, (leaf, in_stock) in enumerate(zip(leaves, in_stock_flags), start=1):
            out["precursors"].append({
                "id": f"{route_id}_{leaf_idx}",
                "format": "smiles",
                "smiles": leaf.get("smiles"),
                "ccd": "",
                "compounds.id": tid,
                "route": route_idx,
                "precursor": leaf_idx,
                "in_stock": in_stock,
            })
    return out


def _write_tables(cfg, route_columns, stream_order, target_rows, route_rows, precursor_rows,
                  routes_map_rows, local_missing_rows):
    """Write every table the step declares; also called before a step that searched nothing exits."""
    target_rows.sort(key=lambda row: stream_order[row["id"]])

    for path in (cfg["targets_csv"], cfg["routes_csv"], cfg["precursors_csv"],
                 cfg["routes_map_csv"], cfg["local_missing_csv"]):
        os.makedirs(os.path.dirname(path), exist_ok=True)

    pd.DataFrame(target_rows, columns=TARGET_COLUMNS).to_csv(cfg["targets_csv"], index=False)
    pd.DataFrame(route_rows, columns=route_columns).to_csv(cfg["routes_csv"], index=False)
    # ccd is empty for every precursor (SMILES have no CCD identity); keep it an
    # empty string so downstream compounds consumers don't read a literal "nan".
    precursors_df = pd.DataFrame(precursor_rows, columns=PRECURSOR_COLUMNS)
    precursors_df["ccd"] = precursors_df["ccd"].fillna("")
    precursors_df.to_csv(cfg["precursors_csv"], index=False)
    pd.DataFrame(routes_map_rows, columns=["id", "file", "compounds.id"]).to_csv(
        cfg["routes_map_csv"], index=False)
    pd.DataFrame(local_missing_rows, columns=["id", "removed_by", "kind", "cause"]).to_csv(
        cfg["local_missing_csv"], index=False)


def main():
    p = argparse.ArgumentParser(description="Run AiZynthFinder over a compounds stream")
    p.add_argument("--config", required=True)
    p.add_argument("--container-prefix", default="")
    args = p.parse_args()

    with open(args.config) as f:
        cfg = json.load(f)

    step_id = step_id_from_table_path(cfg["local_missing_csv"])
    compounds_ds = load_datastream(cfg["compounds_json"])
    extra_scorers = cfg.get("extra_scorers") or []
    route_columns = ROUTE_COLUMNS + [item["column"] for item in extra_scorers]

    _validate_smarts(cfg.get("frozen_substructures"))

    # No SMILES is a deliberate drop; one RDKit cannot sanitize is a failure that keeps an empty targets row.
    targets, local_missing_rows, target_rows, stream_order = [], [], [], {}
    for comp_id, values in iterate_values(compounds_ds, columns=["smiles"]):
        stream_order[comp_id] = len(stream_order)
        smiles = values.get("smiles")
        if smiles is None or (isinstance(smiles, float) and pd.isna(smiles)) or not str(smiles).strip():
            print(f"WARNING: skipping {comp_id} (no smiles)", file=sys.stderr)
            local_missing_rows.append({
                "id": comp_id, "removed_by": step_id, "kind": "filter",
                "cause": "no smiles in compounds stream",
            })
            continue
        target = {"id": comp_id, "smiles": str(smiles).strip()}
        target["canonical"] = _canonical_smiles(target["smiles"])
        if target["canonical"] is None:
            print(f"WARNING: skipping {comp_id} (SMILES cannot be sanitized)", file=sys.stderr)
            local_missing_rows.append({
                "id": comp_id, "removed_by": step_id, "kind": "failure",
                "cause": "SMILES cannot be parsed or sanitized by RDKit",
            })
            target_rows.append(_unavailable_target_row(target))
            continue
        targets.append(target)

    if not targets:
        _write_tables(cfg, route_columns, stream_order, target_rows, [], [], [], local_missing_rows)
        print("ERROR: no compounds with valid SMILES to search "
              "(upstream stream is empty or carries no usable smiles)", file=sys.stderr)
        sys.exit(1)

    os.makedirs(os.path.dirname(cfg["targets_smi"]), exist_ok=True)
    with open(cfg["targets_smi"], "w") as f:
        for target in targets:
            f.write(f"{target['smiles']}\n")

    search_config = _write_search_config(cfg)
    _run_cli(cfg, search_config, len(targets), args.container_prefix)

    results = _read_cli_output(cfg["cli_output"])
    pairs, skipped = _align_results(targets, results)
    for target in skipped:
        print(f"WARNING: aizynthcli returned no row for {target['id']}", file=sys.stderr)
        local_missing_rows.append({
            "id": target["id"], "removed_by": step_id, "kind": "failure",
            "cause": "aizynthcli could not set up the search for this target",
        })
        target_rows.append(_unavailable_target_row(target))

    os.makedirs(cfg["routes_dir"], exist_ok=True)
    route_rows, precursor_rows, routes_map_rows = [], [], []
    solved_targets = 0

    # Rows come back in the order the SMILES were written, which is stream order.
    for target, result in pairs:
        converted = _convert_result(target, result, extra_scorers, cfg["routes_dir"], step_id)
        target_rows.append(converted["target"])
        route_rows.extend(converted["routes"])
        precursor_rows.extend(converted["precursors"])
        routes_map_rows.extend(converted["routes_map"])
        local_missing_rows.extend(converted["missing"])
        solved_targets += bool(converted["routes"])

    _write_tables(cfg, route_columns, stream_order, target_rows, route_rows, precursor_rows,
                  routes_map_rows, local_missing_rows)

    print(f"Targets: {cfg['targets_csv']} ({len(target_rows)} rows, {solved_targets} with routes)")
    print(f"Routes: {cfg['routes_csv']} ({len(route_rows)} rows)")
    print(f"Precursors: {cfg['precursors_csv']} ({len(precursor_rows)} rows)")

    # Every target searched but none produced a route: the tool worked, it just
    # found nothing, so this is a filter, not a failure.
    if not route_rows:
        print(f"WARNING: no routes found for any of {len(targets)} target(s)", file=sys.stderr)


if __name__ == "__main__":
    main()
