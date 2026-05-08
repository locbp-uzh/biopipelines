#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Execution-time script for the Pool tool.

Reads a pool_config.json describing N upstream runs that share the same
stream + table schema and produces:

* For each shared stream: copies (or symlinks where possible) the per-run
  files into the gather folder under composite ids ``<orig_id>_<pool_idx>``
  (1-based pool index). Writes a single map_table CSV per stream carrying
  the composite ids, file paths, and a ``pool.path`` column.
* For each shared table: vertically concatenates the per-run CSVs after
  renumbering the ``id`` column the same way and appending the
  ``pool.path`` column.

Usage:
    python pipe_pool.py <pool_config.json>
"""

import json
import os
import platform
import shutil
import sys
from typing import Optional

import pandas as pd


def _link_or_copy(source_path: str, dest_path: str) -> None:
    """Create a relative symlink on POSIX, copy on Windows or if linking fails."""
    os.makedirs(os.path.dirname(dest_path), exist_ok=True)
    if os.path.exists(dest_path) or os.path.islink(dest_path):
        os.remove(dest_path)
    if not os.path.exists(source_path):
        # Source missing — let the next stage discover it; don't crash here
        # so a partially-failed parallel run still produces partial outputs.
        print(f"  Warning: source missing, skipping: {source_path}")
        return
    if platform.system() == "Windows":
        try:
            os.symlink(os.path.abspath(source_path), dest_path)
        except OSError:
            shutil.copy2(source_path, dest_path)
    else:
        rel = os.path.relpath(source_path, os.path.dirname(dest_path))
        try:
            os.symlink(rel, dest_path)
        except OSError:
            shutil.copy2(source_path, dest_path)


def _renumber_ids(df: pd.DataFrame, pool_idx: int,
                  recount_prefix: Optional[str] = None,
                  recount_start: int = 1) -> pd.DataFrame:
    """Renumber the 'id' column.

    Two modes:

    - Default (``recount_prefix=None``): append ``_<pool_idx>`` to every
      entry in the 'id' column. ``<axis>.id`` provenance columns get the
      same suffix so cross-stream lineage stays consistent within the
      gather output. Other columns are left untouched.

    - Recount (``recount_prefix`` set): replace 'id' with a flat 1-based
      sequence ``<recount_prefix>_<recount_start + i>`` (i = 0..N-1) and
      preserve the original ids in a new ``original.id`` column. ``<axis>.id``
      provenance columns are NOT touched in this mode — they refer to
      upstream entities whose ids stay meaningful only in their own
      namespace, so renaming them to a pool-level counter would break
      downstream joins against the upstream tables.
    """
    if "id" not in df.columns:
        return df
    df = df.copy()
    if recount_prefix is None:
        df["id"] = df["id"].astype(str) + f"_{pool_idx}"
        for col in df.columns:
            if col == "id" or not col.endswith(".id"):
                continue
            df[col] = df[col].astype(str) + f"_{pool_idx}"
    else:
        df["original.id"] = df["id"].astype(str)
        df["id"] = [f"{recount_prefix}_{recount_start + i}" for i in range(len(df))]
    return df


def _gather_stream(stream_key: str, runs_streams: list, out_spec: dict,
                   recount_prefix: Optional[str] = None) -> None:
    """Copy / symlink per-run files into the gather folder and write the
    consolidated map_table for one shared stream.

    If ``recount_prefix`` is set, every emitted row gets a flat 1-based id
    of the form ``<recount_prefix>_<n>`` (n counts across all runs), the
    original id is preserved in an ``original.id`` column, and ``<axis>.id``
    provenance columns are NOT renumbered (their values still refer to
    upstream entities in their own namespace).
    """
    stream_dir = out_spec["stream_dir"]
    out_map = out_spec["map_table"]
    os.makedirs(stream_dir, exist_ok=True)

    rows = []
    counter = 0  # only used when recount_prefix is set
    for pool_idx, run_streams in enumerate(runs_streams, start=1):
        # Find this stream in this run's payload (by stream_key).
        match = next((s for s in run_streams if s["stream_key"] == stream_key), None)
        if match is None:
            continue
        ids = list(match["ids"])
        files = list(match["files"])

        # Prefer the run's existing map_table for the up-to-date id list
        # (covers tools whose map is rewritten at runtime, e.g. PDB).
        run_map = match.get("map_table") or ""
        if run_map and os.path.exists(run_map):
            try:
                run_df = pd.read_csv(run_map)
            except Exception:
                run_df = None
        else:
            run_df = None

        # File-column convention varies by upstream tool: most stream
        # map_tables use 'file', but the PDB tool's structures.csv uses
        # 'file_path'. Probe both.
        file_col = None
        if run_df is not None:
            for cand in ("file", "file_path"):
                if cand in run_df.columns:
                    file_col = cand
                    break

        if run_df is not None and "id" in run_df.columns:
            ids = [str(v) for v in run_df["id"].tolist()]
            if file_col is not None:
                files = [str(v) if pd.notna(v) else "" for v in run_df[file_col].tolist()]
            elif files and len(files) == 1 and "<id>" in files[0]:
                tmpl = files[0]
                files = [tmpl.replace("<id>", oid) for oid in ids]

        for j, oid in enumerate(ids):
            if recount_prefix is not None:
                counter += 1
                composite = f"{recount_prefix}_{counter}"
            else:
                composite = f"{oid}_{pool_idx}"
            src_file = files[j] if j < len(files) else ""
            # Skip wildcards (e.g. 'p.*' from a config-time file template
            # that was never resolved into a concrete extension).
            if src_file and (src_file.endswith(".*") or "*" in os.path.basename(src_file)):
                src_file = ""
            ext = os.path.splitext(src_file)[1] if src_file else ""
            dest_file = (
                os.path.join(stream_dir, f"{composite}{ext}") if src_file else ""
            )
            if src_file and ext:
                _link_or_copy(src_file, dest_file)

            row = {"id": composite, "file": dest_file, "value": "", "pool.path": pool_idx}
            if recount_prefix is not None:
                row["original.id"] = oid

            # Carry over any extra provenance / value columns from the input
            # map_table that aren't already covered. Renumber `<axis>.id`
            # columns for consistency with the new composite id — but only
            # in the default (non-recount) path. Under recount the upstream
            # axis ids stay verbatim because they refer to entities in
            # their own namespace.
            if run_df is not None and j < len(run_df):
                run_row = run_df.iloc[j].to_dict()
                for col, val in run_row.items():
                    if col in ("id", "file", "file_path"):
                        continue
                    if (recount_prefix is None and col.endswith(".id")
                            and pd.notna(val) and str(val) != ""):
                        row[col] = f"{val}_{pool_idx}"
                    elif col not in row:
                        row[col] = "" if pd.isna(val) else val
            rows.append(row)

    out_df = pd.DataFrame(rows)
    os.makedirs(os.path.dirname(out_map), exist_ok=True)
    out_df.to_csv(out_map, index=False)
    print(f"  Wrote stream '{stream_key}' map: {out_map} ({len(out_df)} rows)")


def _gather_table(table_name: str, runs_tables: list, out_path: str,
                  recount_prefix: Optional[str] = None) -> None:
    """Vertically concatenate per-run CSVs after renumbering the id column
    and appending pool.path. If ``recount_prefix`` is set, ids become a
    flat 1-based counter ``<recount_prefix>_<n>`` across the combined
    rows; the original id is kept in ``original.id``."""
    frames = []
    counter = 0
    for pool_idx, run_tables in enumerate(runs_tables, start=1):
        match = next((t for t in run_tables if t["name"] == table_name), None)
        if match is None:
            continue
        path = match["path"]
        if not path or not os.path.exists(path):
            print(f"  Skip table '{table_name}' run #{pool_idx}: missing {path}")
            continue
        try:
            df = pd.read_csv(path)
        except Exception as exc:
            print(f"  Skip table '{table_name}' run #{pool_idx}: read failed ({exc})")
            continue
        df = _renumber_ids(df, pool_idx,
                           recount_prefix=recount_prefix,
                           recount_start=counter + 1)
        df["pool.path"] = pool_idx
        if recount_prefix is not None:
            counter += len(df)
        frames.append(df)

    if not frames:
        # Emit an empty DataFrame with at least the pool.path column so
        # downstream consumers see the file exists.
        out_df = pd.DataFrame(columns=["id", "pool.path"])
    else:
        out_df = pd.concat(frames, ignore_index=True, sort=False)

    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    out_df.to_csv(out_path, index=False)
    print(f"  Wrote table '{table_name}': {out_path} ({len(out_df)} rows)")


def main():
    if len(sys.argv) != 2:
        print("Usage: pipe_pool.py <pool_config.json>")
        sys.exit(2)
    config_path = sys.argv[1]
    with open(config_path, "r") as f:
        config = json.load(f)

    runs_cfg = config["runs"]
    shared_streams = config["shared_streams"]
    shared_tables = config["shared_tables"]
    out_streams = config["out_streams"]
    out_tables = config["out_tables"]
    recount_prefix = config.get("recount_prefix")

    print(f"Pool: {len(runs_cfg)} input runs"
          + (f" (recount_prefix={recount_prefix!r})" if recount_prefix else ""))

    runs_streams = [r["streams"] for r in runs_cfg]
    runs_tables = [r["tables"] for r in runs_cfg]

    for stream_key in shared_streams:
        _gather_stream(stream_key, runs_streams, out_streams[stream_key],
                       recount_prefix=recount_prefix)

    for table_name in shared_tables:
        _gather_table(table_name, runs_tables, out_tables[table_name],
                      recount_prefix=recount_prefix)

    print("Pool: done")


if __name__ == "__main__":
    main()
