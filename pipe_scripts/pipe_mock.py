"""Runtime helper for the Mock tool.

Reads mock_config.json and emits:
  - empty files per ID for file-based streams
  - a map_table CSV per stream (with provenance columns when available)
  - tables with either auto-filled rows or explicit rows

Handles lazy `children` patterns by applying each `produce` suffix to every
parent ID at runtime. Drops IDs listed in `missing` (matches either expanded
IDs or parents)."""

import csv
import json
import os
import sys


def _is_lazy(value):
    return isinstance(value, str) and "[" in value and "]" in value


def _strip_brackets(value):
    if not isinstance(value, str):
        return value
    out, depth = [], 0
    for ch in value:
        if ch == "[":
            depth += 1
            continue
        if ch == "]":
            depth = max(depth - 1, 0)
            continue
        if depth == 0:
            out.append(ch)
    return "".join(out)


def _expand_range_or_set(pattern):
    """Expand <a..b> or <A B C>. Returns list of strings."""
    import re

    m = re.search(r"<([^<>\[\]]+)>", pattern)
    if not m:
        return [pattern]
    token = m.group(1)
    pre, post = pattern[: m.start()], pattern[m.end():]

    if ".." in token:
        lo, hi = token.split("..")
        parts = [str(i) for i in range(int(lo), int(hi) + 1)]
    else:
        parts = token.split()

    out = []
    for p in parts:
        out.extend(_expand_range_or_set(f"{pre}{p}{post}"))
    return out


def _load_source_ids(source_json_path):
    with open(source_json_path) as f:
        data = json.load(f)
    ids = data.get("ids", [])
    map_table = data.get("map_table")
    if map_table and os.path.exists(map_table):
        with open(map_table) as mf:
            rows = list(csv.DictReader(mf))
        if rows:
            return [r["id"] for r in rows]
    return list(ids)


def _cell(table_cfg, col):
    fill = table_cfg.get("fill")
    if isinstance(fill, dict):
        return fill.get(col, "")
    if fill is not None:
        return fill
    return ""


def _resolve_parent_ids(cfg):
    """Re-resolve parent IDs at runtime from source stream JSONs, if any.

    For lazy upstream streams, this reads the runtime-populated map_table,
    picking up whatever the upstream tool actually wrote.
    """
    sources = cfg.get("source_streams") or []
    if not sources:
        return list(cfg["parent_ids"])

    # Single-axis: flat list of IDs.
    if len(sources) == 1:
        return _load_source_ids(sources[0]["path"])

    # Multi-axis: cartesian product of each source's runtime IDs, joined by '+'
    per_axis = [_load_source_ids(src["path"]) for src in sources]
    out = [""]
    for axis_ids in per_axis:
        out = [f"{prefix}+{a}" if prefix else a for prefix in out for a in axis_ids]
    return out


def _apply_children(parent_ids, cfg):
    """Expand parent IDs using cfg['children'] + cfg['produce']. Returns
    (expanded_ids, parent_per_expanded_id)."""
    children = cfg.get("children")
    if not children:
        return list(parent_ids), list(parent_ids)

    if _is_lazy(children):
        produce = cfg.get("produce") or []
        expanded, parents_mapped = [], []
        for pid in parent_ids:
            for suffix in produce:
                # lazy brackets were around the suffix (e.g. "[_<N><A V>]")
                # so we just concatenate the produce entry.
                expanded.append(f"{pid}{suffix}")
                parents_mapped.append(pid)
        return expanded, parents_mapped

    # Deterministic pattern: expand and concatenate per parent.
    suffixes = _expand_range_or_set(children)
    expanded, parents_mapped = [], []
    for pid in parent_ids:
        for s in suffixes:
            expanded.append(f"{pid}_{s}")
            parents_mapped.append(pid)
    return expanded, parents_mapped


def _apply_missing(ids, parents_mapped, missing):
    """Drop IDs whose value is in `missing` OR whose parent is in `missing`."""
    if not missing:
        return ids, parents_mapped
    miss = set(missing)
    keep_ids, keep_parents = [], []
    for i, p in zip(ids, parents_mapped):
        if i in miss or p in miss:
            continue
        keep_ids.append(i)
        keep_parents.append(p)
    return keep_ids, keep_parents


def _rebuild_runtime_provenance(cfg, parent_ids):
    """Re-project config-time provenance onto the runtime parent IDs.

    Config-time `parent_ids` may be compact patterns (e.g. "s1_<1..2>"), and
    config-time `provenance[axis]` is indexed alongside them. At runtime we
    resolve upstream map_tables to their real IDs (["s1_1", "s1_2"]). For
    provenance lookup to stay correct, we rebuild the per-axis lists so that
    their length and order match the runtime parent_ids.

    Single-axis: each axis list becomes `parent_ids` (since parent == axis id).
    Multi-axis: the runtime parent_ids are '+'-joined cartesian products, so
    we split on '+' to recover per-axis values in declared axis order.
    """
    sources = cfg.get("source_streams") or []
    axis_names = list((cfg.get("provenance") or {}).keys()) or list(cfg.get("axis_names") or [])
    if not sources or not axis_names:
        return cfg.get("provenance") or {}

    if len(sources) == 1:
        axis = axis_names[0]
        return {axis: list(parent_ids)}

    prov = {axis: [] for axis in axis_names}
    for pid in parent_ids:
        parts = pid.split("+")
        for i, axis in enumerate(axis_names):
            prov[axis].append(parts[i] if i < len(parts) else "")
    return prov


def main(config_path):
    with open(config_path) as f:
        cfg = json.load(f)

    os.makedirs(cfg["output_folder"], exist_ok=True)

    parent_ids = _resolve_parent_ids(cfg)
    # Re-project provenance onto runtime parent IDs (compact → real IDs).
    cfg["provenance"] = _rebuild_runtime_provenance(cfg, parent_ids)
    cfg["parent_ids"] = list(parent_ids)
    expanded_ids, parents_mapped = _apply_children(parent_ids, cfg)
    expanded_ids, parents_mapped = _apply_missing(
        expanded_ids, parents_mapped, cfg.get("missing") or []
    )

    # Streams: files + map_tables
    for name, spec in (cfg.get("streams") or {}).items():
        template = spec.get("file")
        values = spec.get("values")
        map_path = spec.get("map_table")
        # Per-stream folder set by the Mock config-time code. Older configs
        # without this field fall back to joining at output_folder, which
        # preserves the pre-refactor flat layout.
        stream_folder = spec.get("stream_folder") or cfg["output_folder"]

        if template:
            os.makedirs(stream_folder, exist_ok=True)
            for sid in expanded_ids:
                path = os.path.join(stream_folder, template.replace("<id>", str(sid)))
                os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
                open(path, "w").close()

        if map_path:
            os.makedirs(os.path.dirname(map_path) or ".", exist_ok=True)
            with open(map_path, "w", newline="") as mf:
                w = csv.writer(mf)
                headers = ["id", "file", "value"]
                # Axis-provenance columns from config time, if any. Provenance
                # is indexed per parent (len == len(parent_ids)); when children
                # fans out we look it up via the parent_ids index, not the
                # expanded row index.
                prov = cfg.get("provenance") or {}
                parent_index = {pid: i for i, pid in enumerate(cfg["parent_ids"])}
                axis_cols = [f"{axis}.id" for axis in prov.keys()]
                # Always emit the parent-link column when children is set, so
                # downstream tools can reconstruct the fan-out lineage. For
                # passthrough (no children), parent == id so the column is
                # redundant and omitted.
                parent_axis = f"{name}.parent" if cfg.get("children") else None
                if parent_axis is not None:
                    axis_cols.append(parent_axis)
                headers.extend(axis_cols)
                w.writerow(headers)

                for idx, sid in enumerate(expanded_ids):
                    file_val = (
                        os.path.join(stream_folder, template.replace("<id>", str(sid)))
                        if template else ""
                    )
                    parent = parents_mapped[idx]
                    pidx = parent_index.get(parent)
                    # Values may be per-parent; if so, propagate the parent's value.
                    if values is not None and pidx is not None:
                        val = values[pidx]
                    else:
                        val = ""
                    row = [sid, file_val, val]
                    for axis in prov.keys():
                        if pidx is not None and pidx < len(prov[axis]):
                            row.append(prov[axis][pidx])
                        else:
                            row.append("")
                    if parent_axis is not None:
                        row.append(parent)
                    w.writerow(row)

    # Tables
    for name, spec in (cfg.get("tables") or {}).items():
        cols = spec["columns"]
        rows = spec.get("rows")
        os.makedirs(os.path.dirname(spec["path"]) or ".", exist_ok=True)
        with open(spec["path"], "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(cols)
            if rows is not None:
                for row in rows:
                    w.writerow(row)
            else:
                for sid in expanded_ids:
                    row = [sid if c == "id" else _cell(spec, c) for c in cols]
                    w.writerow(row)

    print(
        f"Mock: wrote {len(cfg.get('streams') or {})} stream(s) and "
        f"{len(cfg.get('tables') or {})} table(s) for {len(expanded_ids)} "
        f"ID(s) (from {len(parent_ids)} parent ID(s))"
    )


if __name__ == "__main__":
    main(sys.argv[1])
