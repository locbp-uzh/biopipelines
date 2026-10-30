#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Write structures_map.csv at runtime from the actual output files.

Tools no longer write a per-design map_table at configuration time (only the
``map_table`` *path* is declared). After the tool runs, this script scans the
output folder, builds one row per produced file, and writes the map. A
constant provenance column (e.g. every design's parent PDB id) can be set via
``--set-provenance``. A per-design parent column is filled by exact id from
the inputs that exist at runtime (``--parents-json`` + ``--suffix-pattern``,
through ``combinatorics.multiplied_id_parents``), never by stripping a suffix
off the file id. If a config-time map happens to exist (legacy tools), its
non-id/file/value columns are carried forward by exact id match.

With ``--local-missing-csv`` the script also answers for the inputs whose run
produced nothing: every id a run that did start would have made but did not
becomes a ``failure`` row there, its cause the real error line of that run's
captured stderr (``--stderr-dir``), and the script exits non-zero only when no
file at all was produced by the runs that started. A parent the input stream
carries without a file (dropped by an upstream filter) never ran and gets no row.

Usage:
    python pipe_update_structures_map.py \\
        --structures-map /path/to/structures_map.csv \\
        --output-folder /path/to/output \\
        --extension pdb \\
        [--best-poses-dir /path/to/best_poses] \\
        [--set-provenance structures.id=4AKE] \\
        [--parents-json /path/to/input_structures.json \\
         --suffix-pattern "<1..5>" --parent-column structures.id \\
         [--key-prefix op7]] \\
        [--local-missing-csv /path/to/_execution/local_missing.csv \\
         [--stderr-dir /path/to/_execution/stderr]]
"""

import argparse
import os
import glob
import re
import sys
import pandas as pd

MISSING_COLUMNS = ["id", "removed_by", "kind", "cause"]
# A foundry validator line ("Error: <reason>") states the reason better than the exception that wraps it.
_FOUNDRY_ERROR = re.compile(r"^Error: .+")
_EXCEPTION_LINE = re.compile(r"^[A-Za-z_][\w.]*(Error|Exception)\b.*:.+")
_HYDRA_HINT = "Set the environment variable HYDRA_FULL_ERROR"


def error_cause(text, limit=200):
    """The reason a run failed: its last error line, not the hint or box drawing printed after it."""
    lines = [line.strip() for line in str(text).splitlines() if line.strip()]
    for pattern in (_FOUNDRY_ERROR, _EXCEPTION_LINE):
        hits = [line for line in lines if pattern.match(line)]
        if hits:
            return (hits[-1][len("Error: "):] if pattern is _FOUNDRY_ERROR else hits[-1])[:limit]
    plain = [line for line in lines if not line.startswith(_HYDRA_HINT)]
    return plain[-1][:limit] if plain else ""


def update_map_from_folder(structures_map, output_folder, extension="pdb",
                           best_poses_dir=None, set_provenance=None,
                           parent_column=None, parent_of=None):
    """Write structures_map.csv based on the actual files in output_folder.

    Args:
        structures_map: Path to the structures_map.csv to write.
        output_folder: Folder to scan for output files.
        extension: File extension to scan for (default: "pdb").
        best_poses_dir: If set, scan this directory instead (for Gnina best_poses).
        set_provenance: Optional dict {column_name: constant_value} of provenance
            columns to add to every row (e.g. {"structures.id": "4AKE"}).
        parent_column: Optional column name for the per-row parent id.
        parent_of: {output id: parent id} for every id the step can produce. The
            parent column is looked up by exact id; a file whose id is not a key
            gets a blank parent and is named on stderr.
    """
    scan_dir = best_poses_dir or output_folder
    pattern = os.path.join(scan_dir, f"*.{extension}")
    found_files = sorted(glob.glob(pattern))

    if not found_files:
        print(f"Warning: No .{extension} files found in {scan_dir}, keeping existing map")
        return

    # Build rows from discovered files
    rows = []
    for fpath in found_files:
        basename = os.path.basename(fpath)
        file_id = os.path.splitext(basename)[0]
        # For Gnina best_poses: strip _best suffix
        if best_poses_dir and file_id.endswith("_best"):
            file_id = file_id[:-5]
        rows.append({"id": file_id, "file": fpath})

    new_df = pd.DataFrame(rows)

    # Constant provenance columns (e.g. every design shares one parent PDB).
    if set_provenance:
        for col, value in set_provenance.items():
            new_df[col] = value

    # Per-row parent id by exact lookup; an id outside the map is no input's output.
    if parent_column and parent_of is not None:
        new_df[parent_column] = new_df["id"].map(parent_of).fillna("")
        strays = new_df.loc[new_df[parent_column] == "", "id"].tolist()
        if strays:
            print(f"WARNING: {len(strays)} of {len(new_df)} file(s) in {scan_dir} are no input's "
                  f"output and carry no {parent_column}: {strays[:10]}", file=sys.stderr)
        if len(strays) == len(new_df):
            raise ValueError(f"No file in {scan_dir} matches an id the inputs and the suffix "
                             f"pattern produce (found: {strays[:10]}); {structures_map} not written")

    # Carry forward any provenance columns from a pre-existing (legacy) map by
    # exact id match. No-op when the config-time map was not written.
    if os.path.exists(structures_map):
        try:
            old_df = pd.read_csv(structures_map)
            prov_cols = [c for c in old_df.columns
                         if c not in ("id", "file", "value") and c not in new_df.columns]
            if prov_cols:
                old_lookup = {}
                for _, row in old_df.iterrows():
                    old_lookup[str(row["id"])] = {c: row[c] for c in prov_cols}
                for col in prov_cols:
                    new_df[col] = new_df["id"].map(
                        lambda x: old_lookup.get(x, {}).get(col, "")
                    )
        except Exception as e:
            print(f"Warning: Could not carry provenance from existing map: {e}")

    new_df.to_csv(structures_map, index=False)
    print(f"Updated {structures_map} with {len(rows)} structures")


def _produced_ids(output_folder, extension="pdb"):
    """Ids of the files the runs left in the output folder."""
    pattern = os.path.join(output_folder, f"*.{extension}")
    return {os.path.splitext(os.path.basename(path))[0] for path in glob.glob(pattern)}


def _parents_that_ran(parents_json):
    """Ids of the input stream whose file exists: a parent an upstream filter dropped never ran."""
    sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    import contextlib
    from biopipelines.biopipelines_io import load_datastream, iterate_files
    # iterate_files names the ids it skips on stdout; keep that out of the step log's data.
    with contextlib.redirect_stdout(sys.stderr):
        return {item_id for item_id, _ in iterate_files(load_datastream(parents_json))}


def record_failures(local_missing_csv, parent_of, ran, produced, stderr_dir=None, key_prefix=None):
    """Write a `failure` row for every id a run that started did not make, with the error that run printed.

    Args:
        local_missing_csv: This step's local missing file; its folder names the step in `removed_by`.
        parent_of: {output id: parent id} for every id the step can make (None without input structures).
        ran: Parents whose run started.
        produced: Ids of the files the runs left.
        stderr_dir: Folder of `<key>.stderr`, one capture per run; the key is the parent id, prefixed when the tool names its outputs so.
        key_prefix: The prefix of that key.
    """
    sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from biopipelines.biopipelines_io import step_id_from_table_path
    step_id = step_id_from_table_path(local_missing_csv)
    rows = []
    for out_id, parent in (parent_of or {}).items():
        if parent not in ran or out_id in produced:
            continue
        key = f"{key_prefix}_{parent}" if key_prefix else parent
        captured = os.path.join(stderr_dir, f"{key}.stderr") if stderr_dir else ""
        text = open(captured, errors="replace").read() if captured and os.path.exists(captured) else ""
        rows.append({"id": out_id, "removed_by": step_id, "kind": "failure",
                     "cause": error_cause(text) or "no design was produced and no error was captured for this input"})
    os.makedirs(os.path.dirname(os.path.abspath(local_missing_csv)), exist_ok=True)
    pd.DataFrame(rows, columns=MISSING_COLUMNS).to_csv(local_missing_csv, index=False)
    if rows:
        print(f"WARNING: {len(rows)} design(s) of {len({parent_of[r['id']] for r in rows})} input(s) were not "
              f"produced; failure rows in {local_missing_csv}", file=sys.stderr)
    return rows


def _parse_set_provenance(spec):
    """Parse 'col=val,col2=val2' into {col: val}. Returns None if empty."""
    if not spec:
        return None
    out = {}
    for piece in spec.split(","):
        piece = piece.strip()
        if not piece:
            continue
        if "=" not in piece:
            raise ValueError(f"--set-provenance entry must be 'col=value', got: {piece!r}")
        col, val = piece.split("=", 1)
        out[col.strip()] = val.strip()
    return out or None


def _parents_from_stream(parents_json, suffix_pattern, key_prefix):
    """{output id: parent id} from the input stream's real ids at runtime."""
    sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from biopipelines.biopipelines_io import load_datastream
    from biopipelines.combinatorics import multiplied_id_parents
    parent_ids = load_datastream(parents_json).ids_expanded
    return multiplied_id_parents(parent_ids, suffix_pattern, key_prefix or "")


def main():
    parser = argparse.ArgumentParser(
        description="Write structures_map.csv from actual runtime output files"
    )
    parser.add_argument("--structures-map", required=True,
                        help="Path to structures_map.csv to write")
    parser.add_argument("--output-folder", required=True,
                        help="Folder containing output structure files")
    parser.add_argument("--extension", default="pdb",
                        help="File extension to scan for (default: pdb)")
    parser.add_argument("--best-poses-dir", default=None,
                        help="Scan this directory instead of output-folder (for Gnina)")
    parser.add_argument("--set-provenance", default=None,
                        help="Constant provenance columns, 'col=value[,col2=value2]' "
                             "(e.g. 'structures.id=4AKE')")
    parser.add_argument("--parents-json", default=None,
                        help="Input DataStream JSON whose ids are the parents; fills "
                             "--parent-column by exact id (needs --suffix-pattern)")
    parser.add_argument("--suffix-pattern", default=None,
                        help="Suffix pattern appended to each parent id, e.g. '<1..5>'")
    parser.add_argument("--parent-column", default=None,
                        help="Column to fill with the parent id, e.g. 'structures.id'")
    parser.add_argument("--key-prefix", default=None,
                        help="Prefix a tool puts in front of the parent in its output "
                             "ids ('<prefix>_<parent>_<n>'); the parent stays bare")
    parser.add_argument("--local-missing-csv", default=None,
                        help="Write a failure row here for every id a run that started did not make, "
                             "and exit non-zero when those runs made no file at all")
    parser.add_argument("--stderr-dir", default=None,
                        help="Folder of <key>.stderr, the captured stderr of each run, for the cause of a failure row")

    args = parser.parse_args()
    if args.parents_json and not (args.suffix_pattern and args.parent_column):
        parser.error("--parents-json needs --suffix-pattern and --parent-column")
    parent_of = (_parents_from_stream(args.parents_json, args.suffix_pattern, args.key_prefix)
                 if args.parents_json else None)

    nothing_produced = False
    if args.local_missing_csv:
        ran = _parents_that_ran(args.parents_json) if args.parents_json else None
        produced = _produced_ids(args.output_folder, args.extension)
        record_failures(args.local_missing_csv, parent_of, ran or set(), produced,
                        stderr_dir=args.stderr_dir, key_prefix=args.key_prefix)
        # Runs that started and made nothing; a stream every parent of which an upstream filter dropped started none.
        nothing_produced = (ran is None or bool(ran)) and not produced

    update_map_from_folder(
        structures_map=args.structures_map,
        output_folder=args.output_folder,
        extension=args.extension,
        best_poses_dir=args.best_poses_dir,
        set_provenance=_parse_set_provenance(args.set_provenance),
        parent_column=args.parent_column,
        parent_of=parent_of,
    )
    if nothing_produced:
        print(f"ERROR: no design was produced in {args.output_folder}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
