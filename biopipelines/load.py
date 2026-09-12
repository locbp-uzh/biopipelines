# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Load tool for loading previously saved pipeline tool outputs.

Allows reusing results from previous pipeline runs without re-execution,
enabling incremental pipeline development and efficient result reuse.
"""

import os
import re
import json
import contextlib
from typing import Dict, List, Any, Optional, Union
from datetime import datetime

try:
    from .base_config import BaseConfig, StandardizedOutput, TableInfo
    from .file_paths import Path
    from .datastream import DataStream
    from . import id_patterns
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo
    from file_paths import Path
    from datastream import DataStream
    import id_patterns


def _under_prefix(value, prefix: str) -> bool:
    """True when a path sits inside prefix, matching whole components only.

    A bare startswith() also matches a sibling run: with prefix /jobs/Job_001,
    /jobs/Job_001_retry/pose.pdb is not under it but shares its string.
    """
    if not isinstance(value, str) or not value:
        return False
    p = os.path.normpath(value)
    root = os.path.normpath(prefix)
    return p == root or p.startswith(root + os.sep) or p.startswith(root + "/")


def _rebase_path(value, old_prefix: str, new_prefix: str):
    """Swap old_prefix for new_prefix on a path genuinely under it."""
    if not _under_prefix(value, old_prefix):
        return value
    if value.startswith(old_prefix):
        # Slice the original rather than the normalized form, so a POSIX path
        # recorded on a cluster does not come back with Windows separators.
        return new_prefix + value[len(old_prefix):]
    tail = os.path.normpath(value)[len(os.path.normpath(old_prefix)):]
    return new_prefix + tail



# Retired per-tool spellings of the canonical map_table file column. A run folder
# written before 1.4.0 still carries one; reading it as data would silently yield
# a stream with no files, so it is named and refused instead.
_RETIRED_FILE_COLUMNS = ("file_path", "msa_file")


def _check_retired_file_column(rows, map_table_path: str, ds_dict) -> None:
    """Refuse a pre-1.4.0 map_table rather than reading zero files from it.

    Only for a stream that declares files. A value-based stream carries no
    per-id file, so a path column there names a source rather than the stream's
    own artifact and keeps whatever the tool called it.
    """
    if rows is None or "file" in rows.columns:
        return
    if not ds_dict.get("files"):
        return
    found = [c for c in _RETIRED_FILE_COLUMNS if c in rows.columns]
    if not found:
        return
    stream_name = str(ds_dict.get("name", "?"))
    raise ValueError(
        f"Load: stream '{stream_name}' has a map_table whose file column is "
        f"'{found[0]}', not 'file': {map_table_path}. "
        f"That spelling was retired in 1.4.0 -- every map_table now records the "
        f"per-id file under 'file'. This run folder was written by an earlier "
        f"version. Rename the column in that CSV to 'file', or re-run the step."
    )


class Load(BaseConfig):
    """
    Load results from a previously executed pipeline tool, with optional filtering.

    This tool allows reusing outputs from previous pipeline runs without
    re-execution, enabling incremental pipeline development and result reuse.

    The Load tool:
    - Loads .expected_outputs.json from a tool's output folder
    - Rebases file paths when the folder has been moved
    - Optionally applies filters to loaded data (expression string or Filter tool output)
    - Provides filtered results through the standard pipeline interface
    - Enables subsequent tools to work with the filtered subset

    Partial outputs are a supported input, not an error. Every generated tool
    script exits 0 by design, so a step can leave a half-written folder behind
    and the run continues; Load is therefore the component that has to recover
    whatever such a step did produce. It never raises on an incomplete folder —
    absent files, a truncated or absent map_table, a ``_FAILED`` marker and a
    populated ``missing.csv`` are all recovered from — and it reports what it
    recovered and what it could not through ``recovery``, ``recovery_summary()``
    and ``recovered_nothing``.
    """

    # Tool identification
    TOOL_NAME = "Load"
    TOOL_VERSION = "1.4"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        return """echo "=== Load ==="
echo "Uses biopipelines environment (no additional installation needed)."
touch "$INSTALL_SUCCESS"
echo "=== Load ready ==="
"""

    def __init__(self, path: str, filter = None, validate_files: bool = True, **kwargs):
        """
        Initialize Load tool.

        Args:
            path: Path to the tool's output folder (containing .expected_outputs.json)
            filter: Optional filter result (ToolOutput from Filter tool) to only load specific IDs
            validate_files: Whether to validate that all referenced files exist
            **kwargs: Additional parameters for BaseConfig

        Output:
            Streams: inherits from the loaded tool (structures, sequences, compounds, etc.)
            Tables: inherits all tables from the loaded tool
        """
        self.tool_folder = os.path.abspath(path)
        self.result_file = os.path.join(self.tool_folder, ".expected_outputs.json")
        self.filter_input = filter
        self.validate_files = validate_files
        self.loaded_result = None
        self.missing_files = []
        self.unresolved_streams = []
        self.original_tool_name = None
        self.filtered_ids = None
        self.completion_status = None
        self.recovery = {}
        self.absent_tables = []
        self.absent_files_by_stream = {}
        self.excused_ids = []

        # Load and validate the result file
        self._load_and_validate_result()
        self._detect_completion_status()

        # Process filter if provided
        if self.filter_input:
            self._process_filter()

        # BaseConfig reads 'name', not 'job_name' — the old key was parked in .params
        # and discarded, so every Load ran with an unnamed job.
        if not kwargs.get('name'):
            original_job_name = self.loaded_result.get('job_name')

            if not original_job_name or original_job_name == 'unknown':
                if ('configuration' in self.loaded_result and
                    'pipeline_context' in self.loaded_result['configuration']):
                    pipeline_context = self.loaded_result['configuration']['pipeline_context']
                    original_job_name = pipeline_context.get('pipeline_job_name', 'unknown')

            if not original_job_name:
                original_job_name = 'unknown'

            kwargs['name'] = f"load_{original_job_name}"

        # Initialize base class
        super().__init__(**kwargs)

        # Set up dependency on filter input if provided
        if self.filter_input and hasattr(self.filter_input, 'config'):
            self.dependencies.append(self.filter_input.config)

        # Initialize folders to prevent AttributeError during script generation
        self.folders = {"pipe_scripts": "pipe_scripts"}  # Will be updated in configure_inputs

    def _load_and_validate_result(self):
        """Load and validate the result file, rebasing paths if the folder was moved."""
        if not os.path.exists(self.result_file):
            raise ValueError(f".expected_outputs.json not found in: {self.tool_folder}")

        try:
            with open(self.result_file, 'r') as f:
                self.loaded_result = json.load(f)
        except json.JSONDecodeError as e:
            raise ValueError(f"Invalid JSON in {self.result_file}: {e}")

        # Validate required fields
        required_fields = ['tool_name', 'tool_class', 'output_structure']
        for field in required_fields:
            if field not in self.loaded_result:
                raise ValueError(f"Missing required field '{field}' in {self.result_file}")

        self.original_tool_name = self.loaded_result['tool_name']

        # Rebase paths if the folder has been moved
        output_structure = self.loaded_result['output_structure']
        old_output_folder = output_structure.get('output_folder', '')
        if old_output_folder and os.path.normpath(old_output_folder) != os.path.normpath(self.tool_folder):
            # The folder was moved — rebase all paths
            old_prefix = os.path.dirname(old_output_folder)  # parent of old tool folder
            new_prefix = os.path.dirname(self.tool_folder)    # parent of new tool folder
            print(f"Load: Rebasing paths from {old_prefix} -> {new_prefix}")
            self._rebase_paths(output_structure, old_prefix, new_prefix)
            # A map_table is a CSV on disk: rebasing the JSON moves the pointer to it
            # but not the paths recorded inside, and those routinely point at a sibling
            # step of the same run. Downstream tools read that CSV directly at runtime,
            # so rewrite it rather than only correcting what this class returns.
            self._path_rebase = (old_prefix, new_prefix)
            self._rebased_map_tables = set()
            self._rebase_map_tables(output_structure, old_prefix, new_prefix)

        # Validate file existence if requested
        if self.validate_files:
            self._validate_file_existence()

    def _detect_completion_status(self):
        """Read the step's COMPLETED/FAILED marker, if the producer left one.

        Tool scripts exit 0 whatever happens, so the marker is the only record
        that a step ended badly. Load surfaces it rather than acting on it: a
        FAILED step's partial output is still worth loading, but the user has to
        be told which folder it came from.
        """
        parent_dir = os.path.dirname(self.tool_folder)
        folder_name = os.path.basename(self.tool_folder.rstrip(os.sep))
        tool = self.original_tool_name or 'Unknown'
        prefix = f"{folder_name.split('_')[0]}_{tool}" if (
            '_' in folder_name and folder_name.split('_')[0].isdigit()) else tool

        for status in ("COMPLETED", "FAILED"):
            if os.path.exists(os.path.join(parent_dir, f"{prefix}_{status}")):
                self.completion_status = status
                break

        if self.completion_status == "FAILED":
            print(f"Load: {prefix} carries a _FAILED marker — the step did not "
                  f"finish. Loading whatever it produced.")

    def _record_recovery(self, stream_name: str, *, declared: int, recovered: int,
                         source: str, unrecovered_ids=()):
        """Note what one stream yielded, for ``recovery_summary()``."""
        self.recovery[stream_name] = {
            'declared': declared,
            'recovered': recovered,
            'source': source,
            'unrecovered_ids': list(unrecovered_ids),
            'absent_files': list(self.absent_files_by_stream.get(stream_name, [])),
        }

    @property
    def recovered_nothing(self) -> bool:
        """True when every stream carrying declared ids yielded none of them.

        A step that reported COMPLETED having produced nothing is where exit-0
        costs the user, and an empty stream otherwise reads as success.
        """
        if not self.recovery:
            return False
        if not any(r['declared'] for r in self.recovery.values()):
            return False
        return all(r['recovered'] == 0 for r in self.recovery.values())

    def recovery_summary(self) -> List[str]:
        """Lines stating precisely what was recovered from the folder and what
        was not. Populated by ``get_output_files()``."""
        tool = self.original_tool_name or 'Unknown'
        lines = [f"Load: recovery summary — {tool} from {self.tool_folder}"]

        if self.completion_status:
            lines.append(f"  completion marker: {self.completion_status}")
        else:
            lines.append("  completion marker: none found")

        if not self.recovery:
            lines.append("  streams: none declared")
        for name, rec in self.recovery.items():
            detail = f"  {name}: recovered {rec['recovered']} of {rec['declared']} " \
                     f"declared id(s) via {rec['source']}"
            shortfall = rec['declared'] - rec['recovered']
            if shortfall > 0:
                detail += f"; {shortfall} not produced"
            lines.append(detail)
            if rec['unrecovered_ids']:
                shown = ', '.join(rec['unrecovered_ids'][:5])
                more = (f" ... and {len(rec['unrecovered_ids']) - 5} more"
                        if len(rec['unrecovered_ids']) > 5 else "")
                lines.append(f"    unmatched id(s): {shown}{more}")
            if rec['absent_files']:
                lines.append(f"    {len(rec['absent_files'])} recovered path(s) "
                             f"are absent from disk")

        if self.excused_ids:
            lines.append(f"  {len(self.excused_ids)} id(s) excused by "
                         f"tables/missing.csv")
        if self.absent_tables:
            lines.append(f"  declared table(s) absent from disk: "
                         f"{', '.join(sorted(self.absent_tables))}")
        if self.unresolved_streams:
            names = sorted({n for n, _p in self.unresolved_streams})
            lines.append(f"  stream(s) with an unresolvable file template: "
                         f"{', '.join(names)}")

        if self.recovered_nothing:
            declared = sum(r['declared'] for r in self.recovery.values())
            lines.append(f"  RECOVERED NOTHING: 0 of {declared} declared output(s) "
                         f"found. This folder carries no usable results.")

        return lines

    def _rebase_paths(self, obj, old_prefix: str, new_prefix: str):
        """
        Recursively rebase all file paths in the output structure.

        Replaces old_prefix with new_prefix in all string values that
        start with old_prefix. Modifies the structure in place.
        """
        old_norm = os.path.normpath(old_prefix)

        if isinstance(obj, dict):
            for key in obj:
                if isinstance(obj[key], str) and os.path.normpath(obj[key]).startswith(old_norm):
                    obj[key] = obj[key].replace(old_prefix, new_prefix, 1)
                elif isinstance(obj[key], (dict, list)):
                    self._rebase_paths(obj[key], old_prefix, new_prefix)
        elif isinstance(obj, list):
            for i, item in enumerate(obj):
                if isinstance(item, str) and os.path.normpath(item).startswith(old_norm):
                    obj[i] = item.replace(old_prefix, new_prefix, 1)
                elif isinstance(item, (dict, list)):
                    self._rebase_paths(item, old_prefix, new_prefix)

    def _validate_file_existence(self):
        """Validate the output structure against what is actually on disk.

        A `files` entry containing `<id>` is a template, not a path — statting it
        always fails. Resolution goes through the stream's map_table, which lists
        one row per id the producer actually wrote, so it is the authoritative
        answer to "what exists".
        """
        self.missing_files = []
        self.unresolved_streams = []
        self.absent_tables = []
        self.absent_files_by_stream = {}
        output_structure = self.loaded_result['output_structure']

        for key, value in output_structure.items():
            if key in ('tables', 'output_folder'):
                continue
            if not isinstance(value, dict) or 'files' not in value:
                continue

            # A readable map with zero rows means nothing was produced, which is an answer, not an unresolved template.
            rows = self._read_map_table(value)
            if rows is not None:
                absent = []
                if 'file' in rows.columns:
                    absent = [str(p) for p in rows['file'].tolist()
                              if isinstance(p, str) and p and not os.path.exists(p)]
                self.absent_files_by_stream[key] = absent
                self.missing_files.extend(absent)
                continue

            files = value.get('files', [])
            candidates = [files] if isinstance(files, str) else files
            absent = []
            for file_path in candidates:
                if not isinstance(file_path, str):
                    continue
                if '<id>' in file_path or any(c in file_path for c in '*?'):
                    self.unresolved_streams.append((key, file_path))
                    continue
                if not os.path.exists(file_path):
                    absent.append(file_path)
            self.absent_files_by_stream[key] = absent
            self.missing_files.extend(absent)

        # Check table files
        if 'tables' in output_structure:
            tables = output_structure['tables']
            if isinstance(tables, dict):
                for table_name, table_info in tables.items():
                    if isinstance(table_info, dict) and 'path' in table_info:
                        path = table_info['path']
                        if not os.path.exists(path):
                            self.absent_tables.append(table_name)
                            self.missing_files.append(path)

        # Check output folder
        if 'output_folder' in output_structure:
            output_folder = output_structure['output_folder']
            if not os.path.exists(output_folder):
                self.missing_files.append(output_folder)

        # Report but don't fail — files might be on a different machine.
        if self.unresolved_streams:
            print(f"Warning: {len(self.unresolved_streams)} stream(s) in {self.result_file} "
                  f"carry an unresolved file template and no map_table to resolve it against:")
            for stream_name, pattern in self.unresolved_streams[:5]:
                print(f"  - {stream_name}: {pattern}")
            if len(self.unresolved_streams) > 5:
                print(f"  ... and {len(self.unresolved_streams) - 5} more")

        if self.missing_files:
            print(f"Warning: {len(self.missing_files)} file(s) listed in {self.result_file} "
                  f"are absent from disk:")
            for missing_file in self.missing_files[:5]:
                print(f"  - {missing_file}")
            if len(self.missing_files) > 5:
                print(f"  ... and {len(self.missing_files) - 5} more")

    def _rebase_map_tables(self, output_structure: Dict[str, Any],
                           old_prefix: str, new_prefix: str):
        """Write rebased copies of every map_table and point the streams at them.

        Tools downstream of a Load read `upstream_map_table` straight off disk, so
        correcting the paths only inside this class would leave them resolving ids
        against the machine the run came from.
        """
        import pandas as pd

        def _fix(v):
            return _rebase_path(v, old_prefix, new_prefix)

        for key, value in output_structure.items():
            if not isinstance(value, dict):
                continue
            mt = value.get('map_table')
            if not mt or not os.path.exists(mt):
                continue
            try:
                # dtype=str / keep_default_na=False: this frame is written back to
                # disk, so type inference would persist its damage -- a zero-padded
                # id reduced to an int, a "NA" cell blanked -- into the CSV every
                # downstream tool then reads.
                rows = pd.read_csv(mt, dtype=str, keep_default_na=False)
            except Exception as e:
                print(f"Warning: could not rebase map_table {mt}: {e}")
                continue
            # Any column holding paths under the old root, not just the two
            # canonical names: tools name their own path column (msa_file,
            # sdf_file, best_pose_file, pocket_file, session_file...), and one
            # left un-rebased sends a downstream tool at paths on the machine the
            # run came from.
            cols = [c for c in rows.columns
                    if rows[c].map(lambda v: _under_prefix(v, old_prefix)).any()]
            if not cols:
                continue
            for c in cols:
                rows[c] = rows[c].map(_fix)
            rebased = os.path.join(os.path.dirname(mt),
                                   f".rebased_{os.path.basename(mt)}")
            try:
                rows.to_csv(rebased, index=False)
            except Exception as e:
                # Leaving the old map_table in place would hand every tool that
                # reads upstream_map_table directly a set of paths on the machine
                # the run came from, and they resolve to nothing here.
                raise RuntimeError(
                    f"Load detected a moved run and could not write the rebased "
                    f"map_table {rebased}: {e}. Copy the run somewhere writable, or "
                    f"load it from its original path.") from e
            value['map_table'] = rebased
            # The on-disk copy is already correct, so _read_map_table must not
            # substitute again -- a new prefix containing the old one (Job_001 ->
            # Job_001_v2) would otherwise rebase twice.
            self._rebased_map_tables.add(os.path.normpath(rebased))
            print(f"Load: rebased map_table for '{key}' -> {rebased}")

    def _map_table_file_paths(self, ds_dict: Dict[str, Any]) -> Optional[List[str]]:
        """Per-id file paths from a stream's map_table, or None if unusable here.

        Returns None (not an empty list) when there is no readable map_table or it
        carries no per-id files, so the caller can fall back to the declared paths.
        """
        rows = self._read_map_table(ds_dict)
        if rows is None or 'file' not in rows.columns:
            return None
        paths = [str(p) for p in rows['file'].tolist() if isinstance(p, str) and p]
        return paths or None

    def _stream_item_count(self, ds_dict: Dict[str, Any]) -> int:
        """Number of items a stream carries — from its map_table when readable.

        Falls back to expanding the declared ids, since a compact pattern like
        ``design_<1..300>`` is one list entry standing for 300 items.
        """
        rows = self._read_map_table(ds_dict)
        if rows is not None:
            return len(rows)
        try:
            return len(DataStream.from_dict(ds_dict).ids_expanded)
        except Exception:
            return len(ds_dict.get('ids', []))

    def _read_map_table(self, ds_dict: Dict[str, Any]):
        """Read a stream's map_table CSV, or None if there isn't a readable one."""
        import pandas as pd

        map_table = ds_dict.get('map_table')
        if not map_table or not os.path.exists(map_table):
            return None
        try:
            # dtype=str / keep_default_na=False: the map is authoritative for ids, and
            # inference rewrites them -- "0001" loads as 1 while its file stays
            # 0001.pdb, breaking every downstream join on id. It also eats real
            # strings, turning a "NA" or "nan" cell into a blank.
            rows = pd.read_csv(map_table, dtype=str, keep_default_na=False)
        except Exception as e:
            self._warn_once(f"read:{map_table}",
                            f"Warning: could not read map_table {map_table}: {e}")
            return None

        rebase = getattr(self, '_path_rebase', None)
        already = os.path.normpath(map_table) in getattr(self, '_rebased_map_tables', set())
        if rebase and not already:
            old_prefix, new_prefix = rebase
            # Every column, not just file/file_path — a tool's own path column
            # (msa_file, sdf_file, pocket_file...) needs rebasing too.
            for col in rows.columns:
                rows[col] = rows[col].map(
                    lambda v: _rebase_path(v, old_prefix, new_prefix))
        if 'id' not in rows.columns:
            self._warn_once(f"no_id:{map_table}",
                            f"Warning: map_table {map_table} has no 'id' column; "
                            f"ignoring it")
            return None
        _check_retired_file_column(rows, str(map_table), ds_dict)
        return rows

    def _warn_once(self, key: str, message: str):
        """Print a warning the first time only — the same map_table is read by
        validation, reconciliation and counting, and three copies of one warning
        buries the rest of the recovery report."""
        seen = self.__dict__.setdefault('_warned', set())
        if key not in seen:
            seen.add(key)
            print(message)

    def _reconcile_from_map_table(self, stream_name: str,
                                  ds_dict: Dict[str, Any]) -> Optional[Dict[str, Any]]:
        """Narrow a stream's ids and files to the rows of its map_table.

        The declared ids are the producer's prediction (and may be a compact
        pattern); the map_table rows are what it actually wrote. Returns None when
        there is no readable map_table, leaving the caller to resolve by globbing.
        """
        rows = self._read_map_table(ds_dict)
        if rows is None:
            return None

        map_ids = [str(i) for i in rows['id'].tolist()]
        reconciled = dict(ds_dict)
        reconciled['ids'] = map_ids

        try:
            declared_count = len(DataStream.from_dict(ds_dict).ids_expanded)
        except Exception:
            declared_count = len(ds_dict.get('ids', []))

        files = ds_dict.get('files', [])
        if isinstance(files, str) and files:
            # Shared-file stream: one artifact covers every id, nothing to narrow.
            self._record_recovery(stream_name, declared=declared_count,
                                  recovered=len(map_ids), source='map_table')
            print(f"Load: {stream_name} — {len(map_ids)} ids from map_table (shared file)")
            return reconciled

        if 'file' in rows.columns:
            map_files = [str(p) if isinstance(p, str) else '' for p in rows['file'].tolist()]
            if any(map_files):
                reconciled['files'] = map_files
            else:
                reconciled['files'] = []
        elif files:
            reconciled['files'] = []

        self._record_recovery(stream_name, declared=declared_count,
                              recovered=len(map_ids), source='map_table')

        note = ""
        if declared_count != len(map_ids):
            note = f" (declared {declared_count} — {declared_count - len(map_ids)} not produced)"
        print(f"Load: {stream_name} — {len(map_ids)} ids from map_table{note}")
        return reconciled

    def _resolve_file_paths(self, output_structure: Dict[str, Any]) -> Dict[str, Any]:
        """
        Resolve glob patterns and match IDs to actual file paths.

        Called at configuration time when validate_files=True.

        Args:
            output_structure: The output structure to resolve

        Returns:
            Updated output structure with resolved file paths
        """
        import glob as glob_module

        resolved = output_structure.copy()

        for file_type, ds_dict in list(resolved.items()):
            if file_type in ('tables', 'output_folder'):
                continue
            if not isinstance(ds_dict, dict) or 'ids' not in ds_dict:
                continue
            files = ds_dict.get('files', [])
            ids = ds_dict.get('ids', [])

            if not files or not ids:
                continue

            # The map_table is the record of what the producer actually wrote —
            # reconcile ids against it rather than trusting the declared ids.
            reconciled = self._reconcile_from_map_table(file_type, ds_dict)
            if reconciled is not None:
                resolved[file_type] = reconciled
                continue

            declared_count = self._stream_item_count(ds_dict)

            if isinstance(files, str):
                self._record_recovery(file_type, declared=declared_count,
                                      recovered=declared_count, source='declared')
                continue

            # Already matched (same length, no globs, no templates)
            if len(files) == len(ids):
                # A one-id stream declaring '<id>.pdb' has as many entries as ids, but a template is still not a path: taking the shortcut reported it recovered whether or not the producer ever wrote the file.
                unresolved = any(
                    '*' in f or '?' in f or '<id>' in f
                    for f in files if isinstance(f, str))
                if not unresolved:
                    self._record_recovery(file_type, declared=declared_count,
                                          recovered=len(ids), source='declared')
                    continue

            # Need to resolve
            print(f"Load: Resolving {file_type} paths ({len(files)} patterns -> {len(ids)} IDs)")

            # Collect all potential files
            all_files = []
            templates = [f for f in files if isinstance(f, str) and '<id>' in f]
            for file_pattern in files:
                if not isinstance(file_pattern, str):
                    continue

                if '<id>' in file_pattern:
                    # A template is not a path: statting it always fails, but substituting its <id> slot with * globs exactly this stream's files.
                    expanded = glob_module.glob(
                        id_patterns.glob_from_file_pattern(file_pattern))
                    all_files.extend(expanded)
                    print(f"  Expanded template {os.path.basename(file_pattern)}: "
                          f"{len(expanded)} files")
                elif '*' in file_pattern or '?' in file_pattern:
                    expanded = glob_module.glob(file_pattern)
                    all_files.extend(expanded)
                    print(f"  Expanded glob: {len(expanded)} files")
                elif os.path.isfile(file_pattern):
                    all_files.append(file_pattern)
                elif os.path.isdir(file_pattern):
                    extensions = ['.pdb', '.cif', '.mmcif', '.fasta', '.fa', '.csv', '.sdf', '.mol2']
                    for ext in extensions:
                        all_files.extend(glob_module.glob(os.path.join(file_pattern, f"*{ext}")))

            # Also check output_folder if we don't have enough files
            if len(all_files) < len(ids) and 'output_folder' in resolved:
                output_folder = resolved['output_folder']
                if os.path.isdir(output_folder):
                    extensions = ['.pdb', '.cif', '.mmcif', '.fasta', '.fa', '.csv', '.sdf', '.mol2']
                    for ext in extensions:
                        for f in glob_module.glob(os.path.join(output_folder, f"**/*{ext}"), recursive=True):
                            if f not in all_files:
                                all_files.append(f)

            # The declared template states how a name is built from an id, so inverting it recovers the id rather than guessing at it — and a name that does not fit the template is not this stream's file.
            by_template_id = {}
            template_owned = set()
            for template in templates:
                for file_path in all_files:
                    recovered = id_patterns.id_from_file_pattern(template, file_path)
                    if recovered:
                        template_owned.add(file_path)
                        by_template_id.setdefault(recovered, file_path)

            # Then exact basenames, and one file serves one id: substring matching in declared order paired design_10 with design_1.pdb and reported design_1 twice.
            slots = [None] * len(ids)
            used_files = set()
            by_basename = {}
            for file_path in all_files:
                by_basename.setdefault(
                    os.path.splitext(os.path.basename(file_path))[0], file_path)

            pending = []
            for pos, item_id in enumerate(ids):
                for exact in (by_template_id.get(item_id), by_basename.get(item_id)):
                    if exact is not None and exact not in used_files:
                        used_files.add(exact)
                        slots[pos] = (item_id, exact)
                        break
                else:
                    pending.append(pos)

            unrecovered_ids = []
            claimed_ids = {s[0] for s in slots if s is not None}
            for pos in pending:
                item_id = ids[pos]
                matched_file = None
                for file_path in all_files:
                    # A name the template already resolved belongs to the id it named, declared or not: design_20_best.pdb must stay out of design_2's reach.
                    if file_path in used_files or file_path in template_owned:
                        continue
                    basename = os.path.splitext(os.path.basename(file_path))[0]

                    if item_id in basename or basename in item_id:
                        matched_file = file_path
                        break
                    if item_id.startswith("rank") and basename.startswith("rank"):
                        id_num = ''.join(filter(str.isdigit, item_id))
                        file_rank_part = basename.split('_')[0]
                        file_num = ''.join(filter(str.isdigit, file_rank_part))
                        if id_num and file_num and id_num == file_num:
                            matched_file = file_path
                            break

                if matched_file is None:
                    unrecovered_ids.append(item_id)
                    print(f"  Warning: No file found for ID '{item_id}'")
                    continue

                actual_id = os.path.splitext(os.path.basename(matched_file))[0]
                if actual_id in claimed_ids:
                    unrecovered_ids.append(item_id)
                    print(f"  Warning: No file found for ID '{item_id}' "
                          f"(nearest match already belongs to '{actual_id}')")
                    continue

                used_files.add(matched_file)
                claimed_ids.add(actual_id)
                slots[pos] = (actual_id, matched_file)

            resolved_ids = [s[0] for s in slots if s is not None]
            resolved_files = [s[1] for s in slots if s is not None]

            resolved[file_type]['files'] = resolved_files
            resolved[file_type]['ids'] = resolved_ids
            # Wildcards resolved — files are now concrete paths

            self._record_recovery(file_type, declared=len(ids),
                                  recovered=len(resolved_ids), source='glob',
                                  unrecovered_ids=unrecovered_ids)
            print(f"  Resolved: {len(resolved_files)}/{len(ids)} {file_type}")

        return resolved

    def _process_filter(self):
        """Process filter input to determine which IDs to keep."""
        import pandas as pd

        # Find the main table to filter
        output_structure = self.loaded_result['output_structure']
        main_table_path = None

        if 'tables' in output_structure:
            tables = output_structure['tables']
            if isinstance(tables, dict):
                priority_names = ['structures', 'analysis', 'combined', 'results']
                for name in priority_names:
                    if name in tables:
                        ds_info = tables[name]
                        if isinstance(ds_info, dict) and 'path' in ds_info:
                            main_table_path = ds_info['path']
                            break

                if not main_table_path:
                    first_ds = next(iter(tables.values()))
                    if isinstance(first_ds, dict) and 'path' in first_ds:
                        main_table_path = first_ds['path']

        if not main_table_path or not os.path.exists(main_table_path):
            raise ValueError(f"Cannot find main table to filter in loaded output: {main_table_path}")

        print(f"Load: Applying filter to table: {main_table_path}")

        try:
            df = pd.read_csv(main_table_path)
            print(f"  - Loaded {len(df)} rows")
        except Exception as e:
            raise ValueError(f"Error loading table for filtering: {e}")

        # Apply filter
        if isinstance(self.filter_input, str):
            try:
                filtered_df = df.query(self.filter_input)
                print(f"  - Applied expression: {self.filter_input}")
            except Exception as e:
                raise ValueError(f"Error applying filter expression '{self.filter_input}': {e}")

        elif hasattr(self.filter_input, 'tables'):
            filter_tables = self.filter_input.tables
            if isinstance(filter_tables, dict):
                if 'filtered' in filter_tables:
                    filter_ds_info = filter_tables['filtered']
                else:
                    filter_ds_info = next(iter(filter_tables.values()))

                if isinstance(filter_ds_info, dict) and 'path' in filter_ds_info:
                    filter_csv_path = filter_ds_info['path']
                elif isinstance(filter_ds_info, str):
                    filter_csv_path = filter_ds_info
                else:
                    raise ValueError("Cannot determine filter CSV path")

                if not os.path.exists(filter_csv_path):
                    raise ValueError(f"Filter CSV file not found: {filter_csv_path}")

                try:
                    filter_df = pd.read_csv(filter_csv_path)
                    if 'id' in filter_df.columns and 'id' in df.columns:
                        filtered_df = df[df['id'].isin(filter_df['id'])]
                        print(f"  - Filtered using precomputed results: {len(filter_df)} IDs")
                    else:
                        raise ValueError("Cannot match filtered IDs - missing 'id' columns")
                except Exception as e:
                    raise ValueError(f"Error loading filter results from {filter_csv_path}: {e}")
            else:
                raise ValueError("Invalid filter input format")
        else:
            raise ValueError(f"Invalid filter input type: {type(self.filter_input)}")

        print(f"  - Result: {len(filtered_df)}/{len(df)} items kept")

        if 'id' in filtered_df.columns:
            self.filtered_ids = set(filtered_df['id'].tolist())
        else:
            self.filtered_ids = set(filtered_df.index.tolist())

        print(f"  - Filtered IDs: {sorted(list(self.filtered_ids)) if len(self.filtered_ids) <= 10 else f'{len(self.filtered_ids)} items'}")

    def _apply_filter_to_output_structure(self, output_structure: Dict[str, Any]) -> Dict[str, Any]:
        """Apply filtering to the output structure based on filtered IDs."""

        print(f"Load: Filtering output structure to {len(self.filtered_ids)} items")

        filtered_structure = output_structure.copy()

        # Filter DataStream dicts
        for stream_name, ds_dict in list(output_structure.items()):
            if stream_name in ('tables', 'output_folder'):
                continue
            if not isinstance(ds_dict, dict) or 'ids' not in ds_dict:
                continue
            ids = ds_dict.get('ids', [])
            files = ds_dict.get('files', [])

            filtered_ids = []
            filtered_files = []

            if len(files) == len(ids):
                for item_id, file_path in zip(ids, files):
                    if item_id in self.filtered_ids:
                        filtered_ids.append(item_id)
                        filtered_files.append(file_path)
            elif len(files) <= 1:
                # Single file or no files — just filter IDs
                filtered_ids = [i for i in ids if i in self.filtered_ids]
                # For single CSV files, update path to Load's output folder
                filtered_files = []
                for file_path in files:
                    if isinstance(file_path, str) and file_path.endswith('.csv'):
                        filtered_files.append(os.path.join(self.output_folder, os.path.basename(file_path)))
                    else:
                        filtered_files.append(file_path)

            filtered_structure[stream_name] = dict(ds_dict)
            filtered_structure[stream_name]['ids'] = filtered_ids
            filtered_structure[stream_name]['files'] = filtered_files
            print(f"  - Filtered {stream_name}: {len(filtered_ids)}/{len(ids)}")

        # Update tables paths to point to Load's own tables/ folder and
        # refresh counts. The source CSV stays in the upstream tool's
        # folder; Load just retargets the metadata so downstream consumers
        # read from a clean, locally-owned file after filtering.
        if 'tables' in filtered_structure:
            for ds_name, ds_info in filtered_structure['tables'].items():
                if isinstance(ds_info, dict):
                    if 'path' in ds_info:
                        original_path = ds_info['path']
                        if isinstance(original_path, str) and original_path.endswith('.csv'):
                            ds_info['path'] = self.table_path(ds_name)
                    if 'count' in ds_info and isinstance(ds_info['count'], int):
                        ds_info['count'] = len(self.filtered_ids)

        return filtered_structure

    def _exclude_missing_ids(self, output_structure: Dict[str, Any]) -> set:
        """
        Check for 'missing' table and return set of IDs to exclude.
        """
        import pandas as pd

        if 'tables' not in output_structure:
            return set()

        tables = output_structure['tables']
        if not isinstance(tables, dict) or 'missing' not in tables:
            return set()

        missing_info = tables['missing']
        if isinstance(missing_info, dict) and 'path' in missing_info:
            missing_path = missing_info['path']
        else:
            return set()

        # A declared-but-absent missing.csv is what a step killed mid-write leaves behind, and raising here would deny the user every id it did produce.
        if not os.path.exists(missing_path):
            print(f"Warning: Load: tables/missing.csv is declared but absent "
                  f"({missing_path}); no ids excused")
            return set()

        try:
            missing_df = pd.read_csv(missing_path)
        except Exception as e:
            print(f"Warning: Load: could not read {missing_path}: {e}; no ids excused")
            return set()

        if 'id' not in missing_df.columns:
            print(f"Warning: Load: {missing_path} has no 'id' column; no ids excused")
            return set()

        missing_ids = set(str(i) for i in missing_df['id'].tolist())
        self.excused_ids = sorted(missing_ids)
        print(f"Load: Excluding {len(missing_ids)} IDs from missing.csv")
        return missing_ids

    def validate_params(self):
        """Validate Load parameters."""
        if not self.loaded_result:
            raise ValueError("No result loaded - initialization failed")

        if self.validate_files and len(self.missing_files) > 0:
            print(f"Warning: {len(self.missing_files)} referenced files are absent from disk")
            print("Consider setting validate_files=False if files have been moved")

        if self.validate_files and len(self.unresolved_streams) > 0:
            print(f"Warning: {len(self.unresolved_streams)} stream(s) could not be resolved — "
                  "their file template has no map_table behind it")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure Load - no inputs needed since we're loading existing results."""
        self.folders = pipeline_folders

    def get_output_files(self) -> Dict[str, Any]:
        """
        Return the loaded output structure converted to DataStream format.

        Returns:
            Dict with DataStream objects for all streams,
            plus tables dict and output_folder string.
        """
        if not self.loaded_result:
            raise RuntimeError("No result loaded")

        output_structure = self.loaded_result['output_structure'].copy()

        # Resolve glob patterns and match IDs to files when validate_files=True
        if self.validate_files:
            output_structure = self._resolve_file_paths(output_structure)
        else:
            print("Load: validate_files=False — ids are propagated as declared, "
                  "not reconciled against the map_table. Downstream tools must "
                  "consume tables/missing.csv to excuse ids that were never produced.")

        # Apply filtering if filter was provided
        if self.filtered_ids is not None:
            output_structure = self._apply_filter_to_output_structure(output_structure)
        else:
            # Check for missing table and exclude those IDs
            missing_ids = self._exclude_missing_ids(output_structure)
            if missing_ids:
                self._filter_missing_ids(output_structure, missing_ids)

        # The framework calls get_output_files() several times per configuration, and six copies of the report is noise the real warnings hide behind.
        summary = "\n".join(self.recovery_summary())
        if summary != getattr(self, '_reported_summary', None):
            self._reported_summary = summary
            print(summary)

        return self._convert_to_datastream_format(output_structure)

    def _filter_missing_ids(self, output_structure: Dict[str, Any], missing_ids: set):
        """Filter out missing IDs from output structure in place."""
        for stream_name, ds_dict in list(output_structure.items()):
            if stream_name in ('tables', 'output_folder'):
                continue
            if not isinstance(ds_dict, dict) or 'ids' not in ds_dict:
                continue
            ids = ds_dict.get('ids', [])
            files = ds_dict.get('files', [])
            if len(files) == len(ids):
                filtered_files = []
                filtered_ids = []
                for f, i in zip(files, ids):
                    if i not in missing_ids:
                        filtered_files.append(f)
                        filtered_ids.append(i)
                ds_dict['files'] = filtered_files
                ds_dict['ids'] = filtered_ids
                if stream_name in self.recovery:
                    self.recovery[stream_name]['recovered'] = len(filtered_ids)
                if len(filtered_files) < len(files):
                    print(f"  - Filtered {stream_name}: {len(filtered_files)}/{len(files)} kept")

    def _is_datastream_dict(self, obj: Any) -> bool:
        """Check if a dict looks like a serialized DataStream."""
        return isinstance(obj, dict) and 'ids' in obj and 'files' in obj

    def _convert_to_datastream_format(self, output_structure: Dict[str, Any]) -> Dict[str, Any]:
        """
        Convert output structure to DataStream format.

        Args:
            output_structure: Dict with DataStream dicts

        Returns:
            Dict with DataStream objects for all streams, plus tables and output_folder
        """
        from .datastream import DataStream

        output_folder = output_structure.get('output_folder', '')
        reserved_keys = {'tables', 'output_folder'}

        result = {"output_folder": output_folder}

        # Convert all DataStream dicts (structures, sequences, compounds, plots, etc.)
        for key, value in output_structure.items():
            if key in reserved_keys:
                continue
            if self._is_datastream_dict(value):
                result[key] = DataStream.from_dict(value)

        # Convert tables to TableInfo objects
        tables_dict = {}
        raw_tables = output_structure.get('tables', {})
        if isinstance(raw_tables, dict):
            for table_name, table_info in raw_tables.items():
                if isinstance(table_info, dict) and 'path' in table_info:
                    tables_dict[table_name] = TableInfo(
                        name=table_name,
                        path=table_info['path'],
                        columns=table_info.get('columns', []),
                        description=table_info.get('description', ''),
                    )
                elif isinstance(table_info, TableInfo):
                    tables_dict[table_name] = table_info
        result["tables"] = tables_dict

        # Pass through tool metadata (e.g. rendering_parameters)
        metadata_keys = {k for k in output_structure
                         if k not in reserved_keys
                         and not self._is_datastream_dict(output_structure[k])}
        for key in metadata_keys:
            result[key] = output_structure[key]

        return result

    def generate_script(self, script_path: str) -> str:
        """
        Generate no-op script since files already exist.
        """
        original_tool = self.loaded_result.get('tool_name', 'Unknown')
        original_job = self.loaded_result.get('job_name')

        if not original_job or original_job == 'unknown':
            if ('configuration' in self.loaded_result and
                'pipeline_context' in self.loaded_result['configuration']):
                pipeline_context = self.loaded_result['configuration']['pipeline_context']
                original_job = pipeline_context.get('pipeline_job_name', 'unknown')

        if not original_job:
            original_job = 'unknown'

        script_content = "#!/bin/bash\n"
        script_content += f"# Load script - loading results from {original_tool}\n"
        script_content += f"# Original job: {original_job}\n"
        script_content += f"# Tool folder: {self.tool_folder}\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()
        script_content += f"""
echo "Loading output from previous {original_tool} execution"
echo "Original job: {original_job}"
echo "Tool folder: {self.tool_folder}"

# Validate that key files exist
echo "Validating loaded files..."

# Check output folder
if [ ! -d "{self.loaded_result['output_structure'].get('output_folder', '')}" ]; then
    echo "Warning: Output folder missing: {self.loaded_result['output_structure'].get('output_folder', '')}"
fi

# Check some key files (first few from each category)
"""

        output_structure = self.loaded_result['output_structure']

        for stream_name, stream_data in output_structure.items():
            if stream_name in ('tables', 'output_folder'):
                continue
            if not isinstance(stream_data, dict) or 'files' not in stream_data:
                continue

            # Templates are patterns, not paths — resolve through the map_table.
            files_list = self._map_table_file_paths(stream_data)
            if files_list is None:
                files = stream_data.get('files', [])
                candidates = [files] if isinstance(files, str) else files
                files_list = [f for f in candidates
                              if isinstance(f, str) and '<id>' not in f
                              and not any(c in f for c in '*?')]

            for file_path in files_list[:3]:
                script_content += f"""
if [ ! -f "{file_path}" ]; then
    echo "Warning: {stream_name} file missing: {file_path}"
fi"""

        # Add filtering section if filter was applied
        if self.filter_input and self.filtered_ids:
            filter_config = {
                "filtered_ids": list(self.filtered_ids),
                "output_structure": self.loaded_result['output_structure'],
                "output_folder": self.output_folder
            }

            filter_config_path = self.configuration_path("filter_config.json")
            with open(filter_config_path, 'w') as f:
                json.dump(filter_config, f, indent=2)

            script_content += f"""

# Apply filter and create filtered tables
echo "Creating filtered copies of tables..."
echo "Filter: {self.filter_input if isinstance(self.filter_input, str) else 'ToolOutput filter'}"
echo "Filtered IDs: {len(self.filtered_ids)} items"
echo "Output folder: {self.output_folder}"

python "{os.path.join(self.folders.get('pipe_scripts', 'pipe_scripts'), 'pipe_load_output_filter.py')}" \\
  --config "{filter_config_path}"

if [ $? -ne 0 ]; then
    echo "Error: Failed to create filtered tables"
    exit 1
fi

echo "Filtered tables created successfully"
"""

        script_content += f"""

echo "Load complete"
echo "Files loaded from {original_tool} are ready for use"

"""
        script_content += self.generate_completion_check_footer()

        return script_content

    def get_config_display(self) -> List[str]:
        """Get Load configuration display."""
        config_lines = super().get_config_display()

        if self.loaded_result:
            original_tool = self.loaded_result.get('tool_name', 'Unknown')
            original_job = self.loaded_result.get('job_name')

            if not original_job or original_job == 'unknown':
                if ('configuration' in self.loaded_result and
                    'pipeline_context' in self.loaded_result['configuration']):
                    pipeline_context = self.loaded_result['configuration']['pipeline_context']
                    original_job = pipeline_context.get('pipeline_job_name', 'unknown')

            if not original_job:
                original_job = 'unknown'

            execution_order = self.loaded_result.get('execution_order', 'unknown')

            config_lines.extend([
                f"Loading from: {original_tool}",
                f"Original job: {original_job}",
                f"Original order: {execution_order}",
                f"Tool folder: {self.tool_folder}",
                f"File validation: {'enabled' if self.validate_files else 'disabled'}",
                f"Completion marker: {self.completion_status or 'none found'}",
            ])

            output_structure = self.loaded_result['output_structure']
            for key, value in output_structure.items():
                if key in ('tables', 'output_folder'):
                    continue
                if isinstance(value, dict) and 'ids' in value:
                    count = self._stream_item_count(value)
                    if count > 0:
                        config_lines.append(f"Loaded {key}: {count}")

            if 'tables' in output_structure:
                table_count = len(output_structure['tables'])
                config_lines.append(f"Loaded tables: {table_count}")

            if self.missing_files:
                config_lines.append(f"Absent files: {len(self.missing_files)}")

            if self.unresolved_streams:
                config_lines.append(f"Unresolved templates: {len(self.unresolved_streams)}")

        return config_lines

    def get_loaded_metadata(self) -> Dict[str, Any]:
        """
        Get metadata about the loaded result.
        """
        if not self.loaded_result:
            return {}

        metadata = {
            'original_tool_name': self.loaded_result.get('tool_name'),
            'original_tool_class': self.loaded_result.get('tool_class'),
            'original_job_name': self.loaded_result.get('job_name'),
            'execution_order': self.loaded_result.get('execution_order'),
            'tool_folder': self.tool_folder,
            'missing_files_count': len(self.missing_files),
            'unresolved_streams_count': len(self.unresolved_streams),
            'validation_enabled': self.validate_files,
            'completion_status': self.completion_status,
            'absent_tables': list(self.absent_tables),
            'recovery': {k: dict(v) for k, v in self.recovery.items()},
        }

        if 'execution_metadata' in self.loaded_result:
            metadata['execution_metadata'] = self.loaded_result['execution_metadata']

        if 'configuration' in self.loaded_result:
            metadata['original_configuration'] = self.loaded_result['configuration']

        return metadata

    def to_dict(self) -> Dict[str, Any]:
        """Serialize Load configuration."""
        base_dict = super().to_dict()
        base_dict.update({
            "load_params": {
                "tool_folder": self.tool_folder,
                "validate_files": self.validate_files,
                "original_tool_name": self.original_tool_name,
                "missing_files_count": len(self.missing_files) if self.missing_files else 0,
                "unresolved_streams_count": len(self.unresolved_streams) if self.unresolved_streams else 0
            }
        })
        return base_dict

    def __str__(self) -> str:
        """String representation."""
        original_tool = self.original_tool_name or 'Unknown'
        return f"Load(from {original_tool}: {os.path.basename(self.tool_folder)})"


def LoadMultiple(path: str,
                tool: Optional[str] = None,
                suffix: Optional[str] = None,
                in_suffix: Optional[Union[str, List[str]]] = None,
                not_in_suffix: Optional[Union[str, List[str]]] = None,
                ascending: bool = True,
                folder: Optional[str] = "LoadMultiple",
                **load_output_kwargs) -> Dict[str, Load]:
    """
    Load multiple tool outputs from a job folder by scanning tool subfolders.

    Scans for subfolders matching the NNN_ToolName[_Suffix] pattern and loads
    .expected_outputs.json from each.

    Args:
        path: Path to the job output folder (containing NNN_ToolName/ subfolders)
        tool: Filter by tool name (e.g., "Boltz2", "RFdiffusion")
        suffix: Filter by exact folder suffix
        in_suffix: Filter requiring suffix to contain string(s)
        not_in_suffix: Filter excluding suffix containing string(s)
        ascending: Sort order (True = ascending by name, False = descending)
        folder: Name of a Folder() to nest the created Load steps under so they
                don't clutter the job root (default "LoadMultiple"). Pass None to
                keep the Load steps at the job root. Note this nests only the
                Load steps themselves; tools you chain off the returned outputs
                are constructed after LoadMultiple returns, so they fall outside
                this folder unless you open your own Folder() around them.
        **load_output_kwargs: Additional parameters passed to Load constructor

    Returns:
        Dictionary mapping folder names to Load objects

    Examples:
        # Load all tool outputs from a job
        >>> data = LoadMultiple("/path/to/Job_001")

        # Load only Boltz2 outputs
        >>> boltz = LoadMultiple("/path/to/Job_001", tool="Boltz2")

        # Load outputs with a specific suffix
        >>> cycle10 = LoadMultiple("/path/to/Job_001", suffix="Cycle10")

        # Keep the Load steps at the job root (no nesting)
        >>> data = LoadMultiple("/path/to/Job_001", folder=None)
    """
    from .pipeline import Folder, Pipeline
    job_folder = os.path.abspath(path)

    if not os.path.exists(job_folder):
        raise ValueError(f"Job folder not found: {job_folder}")

    if not os.path.isdir(job_folder):
        raise ValueError(f"Path is not a directory: {job_folder}")

    # Scan for tool subfolders matching NNN_ToolName[_Suffix] pattern
    tool_folder_pattern = re.compile(r'^(\d+)_(.+)$')

    candidates = []
    for entry in os.listdir(job_folder):
        entry_path = os.path.join(job_folder, entry)
        if not os.path.isdir(entry_path):
            continue
        match = tool_folder_pattern.match(entry)
        if not match:
            continue
        # Check that .expected_outputs.json exists at the tool-folder root.
        expected_json = os.path.join(entry_path, ".expected_outputs.json")
        if not os.path.exists(expected_json):
            continue
        candidates.append((entry, match.group(1), match.group(2)))  # (folder_name, index, rest)

    if not candidates:
        raise ValueError(f"No tool folders with .expected_outputs.json found in: {job_folder}")

    # Filter and load outputs. Nest the created Load steps under Folder(folder)
    # when a pipeline is active and a name was given, so they don't scatter at
    # the job root. nullcontext keeps the loop body identical when not nesting
    # (folder=None, or called outside a Pipeline context).
    loaded_outputs = {}
    nest = (folder is not None and Pipeline.get_active_pipeline() is not None)
    folder_ctx = Folder(folder) if nest else contextlib.nullcontext()

    with folder_ctx:
        for folder_name, index, rest in candidates:
            # Parse tool name and suffix from the rest (e.g., "Boltz2" or "Boltz2_Cycle10")
            parts = rest.split('_', 1)
            folder_tool_name = parts[0]
            folder_suffix = parts[1] if len(parts) > 1 else None

            # Apply tool name filter
            if tool is not None and folder_tool_name != tool:
                continue

            # Apply suffix filters
            if suffix is not None:
                if folder_suffix != suffix:
                    continue

            if in_suffix is not None:
                if folder_suffix is None:
                    continue
                checks = [in_suffix] if isinstance(in_suffix, str) else in_suffix
                if not all(c in folder_suffix for c in checks):
                    continue

            if not_in_suffix is not None:
                if folder_suffix is not None:
                    checks = [not_in_suffix] if isinstance(not_in_suffix, str) else not_in_suffix
                    if any(c in folder_suffix for c in checks):
                        continue

            # Passed all filters — create Load
            folder_path = os.path.join(job_folder, folder_name)
            try:
                load_output = Load(folder_path, **load_output_kwargs)
                loaded_outputs[folder_name] = load_output
            except Exception as e:
                print(f"Warning: Could not create Load for {folder_name}: {e}")
                continue

    if not loaded_outputs:
        raise ValueError(
            f"No outputs matched the specified filters:\n"
            f"  Path: {job_folder}\n"
            f"  Tool: {tool if tool else 'any'}\n"
            f"  Suffix: {suffix if suffix else 'any'}\n"
            f"  Found {len(candidates)} tool folders, but none matched filters"
        )

    # Sort outputs by name
    sorted_outputs = dict(sorted(loaded_outputs.items(), reverse=not ascending))

    print(f"LoadMultiple: Loaded {len(sorted_outputs)} outputs from {job_folder}")
    if tool:
        print(f"  - Tool filter: {tool}")
    if suffix:
        print(f"  - Suffix filter: {suffix}")
    print(f"  - Sort order: {'ascending' if ascending else 'descending'}")
    print(f"  - Keys: {list(sorted_outputs.keys())}")

    return sorted_outputs
