#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Execution-time script for the Pool tool.

Reads a pool_config.json describing N upstream runs that share the same stream + table schema and produces:

* For each shared stream: copies (or symlinks where possible) the per-run files into the gather folder under composite ids ``<orig_id>_<pool_idx>`` (1-based pool index).
  It writes one map_table CSV per stream with the composite ids, file paths, a ``pool.path`` column and an ``original.id`` column holding the id before renaming.
* For each shared table: concatenates the per-run CSVs after renumbering the ``id`` column the same way and appending ``pool.path`` and ``original.id``.
* With ``keep_ids`` no id is suffixed, and an id carried by more than one run raises a ValueError instead of overwriting.

Usage:
    python pipe_pool.py <pool_config.json>
"""

import json
import os
import platform
import shutil
import sys
import tempfile
from typing import Dict, List, Optional

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))


def _raise_on_collision(label: str, ids_by_run: Dict[int, List[str]]) -> None:
    """keep_ids: raise a ValueError naming ids that more than one run carries."""
    from biopipelines.pool import raise_on_id_collision
    raise_on_id_collision(label, ids_by_run)


def _recounted_id(recount_prefix: str, numbering: dict, pool_idx: int, original_id: str) -> str:
    """recount_prefix: the pooled id of a row, one number per (run, original id) across every stream and table."""
    from biopipelines.combinatorics import recounted_id
    return recounted_id(recount_prefix, numbering, pool_idx, original_id)


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


def _is_axis_id(column: str) -> bool:
    """An ``<axis>.id`` provenance column, which the default mode suffixes; ``original.id`` is kept."""
    return column.endswith(".id") and column not in ("id", "original.id")


def _renumber_ids(df: pd.DataFrame, pool_idx: int,
                  recount_prefix: Optional[str] = None,
                  numbering: Optional[dict] = None,
                  keep_ids: bool = False) -> pd.DataFrame:
    """Renumber the 'id' column.

    Three modes, all of which record the id before renaming in ``original.id``:

    - Default (``recount_prefix=None``): append ``_<pool_idx>`` to every entry in the 'id' column.
      ``<axis>.id`` provenance columns get the same suffix so cross-stream lineage stays consistent within the gather output.
      Other columns are left untouched, and an ``original.id`` the input already carried is kept as it was.

    - Keep (``keep_ids``): no id changes; ``original.id`` equals 'id'.

    - Recount (``recount_prefix`` set): replace 'id' with ``<recount_prefix>_<n>``, n being the number ``numbering`` gives the (run, original id) pair, so a row has the same id here as in every stream and table that names it.
      The original ids are stored in ``original.id``, replacing any such column.
      ``<axis>.id`` provenance columns are NOT touched in this mode: they name upstream entities, and a pool-level counter would break joins against the upstream tables.
    """
    if "id" not in df.columns:
        return df
    df = df.copy()
    original = df["id"].astype(str)
    if recount_prefix is not None:
        numbering = {} if numbering is None else numbering
        df["original.id"] = original
        df["id"] = [_recounted_id(recount_prefix, numbering, pool_idx, oid) for oid in original]
        return df
    carried = "original.id" in df.columns
    if not keep_ids:
        df["id"] = original + f"_{pool_idx}"
        for col in df.columns:
            if not _is_axis_id(col):
                continue
            # a blank cell (the row of a failed id) stays blank, as in _gather_stream
            named = df[col].notna() & (df[col].astype(str) != "")
            df[col] = df[col].where(~named, df[col].astype(str) + f"_{pool_idx}")
    if not carried:
        df["original.id"] = original
    return df


def _gather_shared_stream(stream_key: str, runs_streams: list, out_spec: dict,
                          recount_prefix: Optional[str] = None,
                          keep_ids: bool = False,
                          numbering: Optional[dict] = None) -> None:
    """Slice each run's shared artifact (renaming to composite ids) and
    merge the parts into one combined artifact. Writes the consolidated
    map_table with one row per composite id, all pointing at the merged
    file path.
    """
    import sys as _sys
    import os as _os
    _here = _os.path.dirname(_os.path.abspath(__file__))
    _sys.path.insert(0, _os.path.dirname(_here))
    from biopipelines.stream_slicers import get_slicer, get_merger

    stream_dir = out_spec["stream_dir"]
    out_map = out_spec["map_table"]
    fmt = out_spec.get("format", "")

    if keep_ids:
        _raise_on_collision(f"stream '{stream_key}'", {
            pool_idx: list(s["ids"])
            for pool_idx, run_streams in enumerate(runs_streams, start=1)
            for s in run_streams if s["stream_key"] == stream_key
        })

    numbering = {} if numbering is None else numbering
    parts = []
    rows = []
    tmp_root = tempfile.mkdtemp(prefix="pool_shared_")
    for pool_idx, run_streams in enumerate(runs_streams, start=1):
        match = next(
            (s for s in run_streams if s["stream_key"] == stream_key), None
        )
        if match is None:
            continue
        src = match["files"]
        # Defensive: tolerate runs whose payload predates the is_shared_file
        # flag — if files came as a 1-element list with a real path, treat
        # it as shared too.
        if isinstance(src, list):
            if len(src) == 1 and isinstance(src[0], str):
                src = src[0]
            else:
                continue
        if not src or not os.path.exists(src):
            print(f"  Skip shared stream '{stream_key}' run #{pool_idx}: missing {src}")
            continue

        ids = list(match["ids"])
        # Build composite ids + rename map for this run.
        rename: Dict[str, str] = {}
        run_composites: List[str] = []
        for oid in ids:
            if recount_prefix is not None:
                composite = _recounted_id(recount_prefix, numbering, pool_idx, oid)
            elif keep_ids:
                composite = oid
            else:
                composite = f"{oid}_{pool_idx}"
            rename[oid] = composite
            run_composites.append(composite)

        slicer = get_slicer(fmt)
        part = os.path.join(tmp_root, f"pool{pool_idx}_{os.path.basename(src)}")
        slicer(src, part, ids, rename_map=rename)
        parts.append(part)

        # Map rows: one per composite id, file = merged dest (set after merge)
        for oid, composite in zip(ids, run_composites):
            rows.append({"id": composite, "value": "", "pool.path": pool_idx,
                         "original.id": oid})

    # Merge all parts into the canonical dest under stream folder.
    if parts:
        canonical_basename = os.path.basename(
            next(
                (
                    s["files"] if isinstance(s["files"], str) else s["files"][0]
                    for run_streams in runs_streams for s in run_streams
                    if s.get("stream_key") == stream_key
                    and s.get("files")
                ),
                "shared",
            )
        )
        dest = os.path.join(stream_dir, canonical_basename)
        merger = get_merger(fmt)
        merger(parts, dest)
        for row in rows:
            row["file"] = dest
    else:
        dest = ""
        for row in rows:
            row["file"] = ""

    out_df = pd.DataFrame(rows)
    # Reorder so 'file' sits next to 'id' (matches the per-id _gather_stream
    # convention).
    if "file" in out_df.columns:
        cols = ["id", "file", "value", "pool.path"]
        if "original.id" in out_df.columns:
            cols.append("original.id")
        existing = [c for c in cols if c in out_df.columns]
        out_df = out_df[existing + [c for c in out_df.columns if c not in existing]]

    os.makedirs(os.path.dirname(out_map), exist_ok=True)
    out_df.to_csv(out_map, index=False)
    print(f"  Wrote shared stream '{stream_key}' map: {out_map} ({len(out_df)} rows -> {dest})")


def _id_columns(path: str) -> List[str]:
    """Identifier columns of a CSV: `id`, `original.id` and every `<stream>.id` provenance column."""
    try:
        header = pd.read_csv(path, nrows=0).columns
    except Exception:
        return []
    return [c for c in header if c == "id" or c == "original.id" or c.endswith(".id")]


def _read_ids_as_text(path: str, **kwargs):
    """``pd.read_csv`` with the identifier columns forced to text.

    Left to its own inference pandas turns the id `007` into 7 and the id `NA` into a
    missing value, and a later ``str()`` cannot undo either. Only the identifier columns
    are coerced, so a result column stays the number it is.
    """
    converters = {c: str for c in _id_columns(path)}
    return pd.read_csv(path, converters=converters, **kwargs)


def _load_run_stream(match: dict):
    """(ids, files, run_df) of one run's stream.

    The payload ids are the VIEW the pipeline asked for; the run's map_table is the whole
    warehouse behind it. The map is still read, because a tool that rewrites it at runtime
    (PDB) is the only one that knows its real ids — but the view is intersected with it,
    never replaced by it, or an explicit selection would silently take back every row it
    was built to drop. ``run_df`` is None when no readable map exists.
    """
    from biopipelines.id_patterns import is_pattern

    ids = [str(v) for v in match["ids"]]
    files = list(match["files"])

    run_map = match.get("map_table") or ""
    if run_map and os.path.exists(run_map):
        try:
            run_df = _read_ids_as_text(run_map)
        except Exception:
            run_df = None
    else:
        run_df = None

    if run_df is None or "id" not in run_df.columns:
        return ids, files, run_df

    map_ids = [str(v) for v in run_df["id"].tolist()]
    requested = set(ids)
    # A view is lazy when it has no ids yet, or when its ids are still patterns: only the map can resolve it.
    lazy = (not ids) or any(is_pattern(i) for i in ids)
    keep = [True] * len(map_ids) if lazy else [i in requested for i in map_ids]
    # No overlap at all means the tool renamed its ids at runtime, and the map is then the only truthful list.
    if not lazy and not any(keep):
        keep = [True] * len(map_ids)
    run_df = run_df[pd.Series(keep, index=run_df.index)].reset_index(drop=True)
    ids = [i for i, kept in zip(map_ids, keep) if kept]

    if "file" in run_df.columns:
        files = [str(v) if pd.notna(v) else "" for v in run_df["file"].tolist()]
    elif files and len(files) == 1 and "<id>" in files[0]:
        tmpl = files[0]
        files = [tmpl.replace("<id>", oid) for oid in ids]
    return ids, files, run_df


def _gather_stream(stream_key: str, runs_streams: list, out_spec: dict,
                   recount_prefix: Optional[str] = None,
                   content_file_col: Optional[str] = None,
                   keep_ids: bool = False,
                   numbering: Optional[dict] = None) -> None:
    """Copy / symlink per-run files into the gather folder and write the
    consolidated map_table for one shared stream.

    If ``recount_prefix`` is set, every emitted row gets an id of the form ``<recount_prefix>_<n>``, n being the number ``numbering`` gives its (run, original id) pair, the same in every stream and table of the pool.
    The original id is preserved in an ``original.id`` column, and ``<axis>.id`` provenance columns are NOT renumbered (their values still refer to upstream entities in their own namespace).

    If ``keep_ids`` is set, ids and ``<axis>.id`` columns are left as the runs
    wrote them, and an id carried by more than one run raises before any file
    is copied. Without either option every id gets ``_<pool_idx>``, and every
    row records the id before that in ``original.id``.

    Shared-file streams (each run owns one artifact, e.g. a multi-record
    FASTA) take a different path: each run is sliced with its per-run
    rename map, then the parts are merged into a single combined artifact
    under the stream folder. The map_table is written with one row per
    composite id, all pointing at the merged file.

    Content streams (the map_table IS the content table) own a domain
    schema; ``content_file_col`` names the column that holds the per-id file
    path (``file_path`` for PDB's structures.csv, ``file`` otherwise) so the
    pooled content table keeps the upstream schema instead of inventing a
    generic ``file``/``value`` pair. When None, this is a non-content stream
    and the generic ``file``/``value`` columns are used.
    """
    stream_dir = out_spec["stream_dir"]
    out_map = out_spec["map_table"]
    fmt = out_spec.get("format", "")
    os.makedirs(stream_dir, exist_ok=True)

    # Shared-file fast path: dispatch to slicer + merger and write the map.
    any_shared = any(
        isinstance(s, dict) and s.get("stream_key") == stream_key and s.get("is_shared_file")
        for run_streams in runs_streams for s in run_streams
    )
    numbering = {} if numbering is None else numbering
    if any_shared:
        _gather_shared_stream(stream_key, runs_streams, out_spec, recount_prefix, keep_ids, numbering)
        return

    # Find this stream in each run's payload (by stream_key).
    loaded = []
    for pool_idx, run_streams in enumerate(runs_streams, start=1):
        match = next((s for s in run_streams if s["stream_key"] == stream_key), None)
        if match is not None:
            loaded.append((pool_idx, *_load_run_stream(match)))
    if keep_ids:
        _raise_on_collision(f"stream '{stream_key}'",
                            {pool_idx: ids for pool_idx, ids, _files, _df in loaded})

    rows = []
    for pool_idx, ids, files, run_df in loaded:
        # An upstream original.id means something else (e.g. a parent id): keep it.
        upstream_original = run_df is not None and "original.id" in run_df.columns
        for j, oid in enumerate(ids):
            if recount_prefix is not None:
                composite = _recounted_id(recount_prefix, numbering, pool_idx, oid)
            elif keep_ids:
                composite = oid
            else:
                composite = f"{oid}_{pool_idx}"
            src_file = files[j] if j < len(files) else ""
            # Skip wildcards (e.g. 'p.*' from a config-time file template
            # that was never resolved into a concrete extension).
            if src_file and (src_file.endswith(".*") or "*" in os.path.basename(src_file)):
                src_file = ""
            # Skip sources that don't exist on disk. The row still gets
            # appended below with an empty file column so the id is
            # preserved in the map_table, but we never claim a dest_file
            # that was never created. This matches the wildcard-skip
            # convention immediately above.
            if src_file and not os.path.exists(src_file):
                print(f"  Warning: source missing for id {oid!r}, dropping file from map: {src_file}")
                src_file = ""
            ext = os.path.splitext(src_file)[1] if src_file else ""
            dest_file = (
                os.path.join(stream_dir, f"{composite}{ext}") if src_file else ""
            )
            if src_file and ext:
                _link_or_copy(src_file, dest_file)

            if content_file_col is not None:
                # Content table: write the file path under the upstream column
                # name and don't inject the generic file/value pair.
                row = {"id": composite, content_file_col: dest_file,
                       "pool.path": pool_idx}
            else:
                row = {"id": composite, "file": dest_file, "value": "",
                       "pool.path": pool_idx}
            if recount_prefix is not None or not upstream_original:
                row["original.id"] = oid

            # Carry the input map's other columns; only the default mode suffixes `<axis>.id`
            if run_df is not None and j < len(run_df):
                run_row = run_df.iloc[j].to_dict()
                for col, val in run_row.items():
                    if col in ("id", "file") or col in row:
                        continue
                    if (recount_prefix is None and not keep_ids and _is_axis_id(col)
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
                  recount_prefix: Optional[str] = None,
                  keep_ids: bool = False,
                  numbering: Optional[dict] = None) -> None:
    """Vertically concatenate per-run CSVs after renumbering the id column
    and appending pool.path.
    If ``recount_prefix`` is set, ids become ``<recount_prefix>_<n>`` with the number ``numbering`` gives the (run, original id) pair, as in the streams, and the original id is kept in ``original.id``.
    If ``keep_ids`` is set, ids stay as written and an id carried by more than one run raises."""
    frames = []
    ids_by_run = {}
    numbering = {} if numbering is None else numbering
    for pool_idx, run_tables in enumerate(runs_tables, start=1):
        match = next((t for t in run_tables if t["name"] == table_name), None)
        if match is None:
            continue
        path = match["path"]
        if not path or not os.path.exists(path):
            print(f"  Skip table '{table_name}' run #{pool_idx}: missing {path}")
            continue
        try:
            df = _read_ids_as_text(path)
        except Exception as exc:
            print(f"  Skip table '{table_name}' run #{pool_idx}: read failed ({exc})")
            continue
        if "id" in df.columns:
            ids_by_run[pool_idx] = df["id"].astype(str).tolist()
        df = _renumber_ids(df, pool_idx,
                           recount_prefix=recount_prefix,
                           numbering=numbering,
                           keep_ids=keep_ids)
        df["pool.path"] = pool_idx
        frames.append(df)

    if keep_ids:
        _raise_on_collision(f"table '{table_name}'", ids_by_run)

    if not frames:
        # Emit an empty DataFrame with the provenance columns so downstream
        # consumers see the file exists.
        out_df = pd.DataFrame(columns=["id", "pool.path", "original.id"])
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
    # The configuration step numbered every declared (run, id) pair; a pair it did not declare is numbered after them.
    numbering = {(int(pool_idx), str(oid)): composite
                 for pool_idx, oid, composite in config.get("recount_numbering", [])}
    keep_ids = bool(config.get("keep_ids", False))
    # Content streams: _gather_stream writes the combined file (preserving
    # the upstream domain schema + copying any per-id files); the table of the
    # same name points at that file, so skip the redundant table concat.
    content_bearing_streams = set(config.get("content_bearing_streams", []))
    content_file_cols = config.get("content_file_cols", {})

    print(f"Pool: {len(runs_cfg)} input runs"
          + (f" (recount_prefix={recount_prefix!r})" if recount_prefix else "")
          + (" (keep_ids)" if keep_ids else ""))

    runs_streams = [r["streams"] for r in runs_cfg]
    runs_tables = [r["tables"] for r in runs_cfg]

    for stream_key in shared_streams:
        _gather_stream(stream_key, runs_streams, out_streams[stream_key],
                       recount_prefix=recount_prefix,
                       content_file_col=content_file_cols.get(stream_key)
                       if stream_key in content_bearing_streams else None,
                       keep_ids=keep_ids, numbering=numbering)

    # A content table shares its file with the stream of the same name (already
    # written by _gather_stream above); don't overwrite it with the table concat.
    for table_name in shared_tables:
        if table_name in content_bearing_streams:
            continue
        _gather_table(table_name, runs_tables, out_tables[table_name],
                      recount_prefix=recount_prefix, keep_ids=keep_ids, numbering=numbering)

    print("Pool: done")


if __name__ == "__main__":
    main()
