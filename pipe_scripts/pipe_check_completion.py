#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Generic completion checker for pipeline tools.

Checks if expected output files exist and creates COMPLETED/FAILED status files accordingly.
This allows bash scripts to skip already completed steps and provides clear status indication.
"""

import os
import sys
import glob
import argparse
import json
from pathlib import Path
from typing import Dict, List, Any, Optional, Tuple

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from biopipelines import id_patterns

def check_file_exists(file_path: str) -> bool:
    """
    Check if a file or directory exists.

    Supports glob wildcard patterns (e.g. rank0001_*.cif).

    Args:
        file_path: Path to check (may contain * or ? wildcards)

    Returns:
        True if exists, False otherwise
    """
    if '*' in file_path or '?' in file_path:
        return len(glob.glob(file_path)) > 0
    return os.path.exists(file_path) and (os.path.isfile(file_path) or os.path.isdir(file_path))

def check_files_exist(file_list: List[str]) -> Tuple[bool, List[str]]:
    """
    Check if all files in a list exist.
    
    Args:
        file_list: List of file paths to check
        
    Returns:
        Tuple of (all_exist: bool, missing_files: List[str])
    """
    missing_files = []
    for file_path in file_list:
        if not check_file_exists(file_path):
            missing_files.append(file_path)
    
    return len(missing_files) == 0, missing_files

def extract_id_file_pairs(category_data) -> List[Tuple[Optional[str], str]]:
    """
    Extract (owner_id, file_path) pairs from a category entry.

    The owner_id is the stream id that produced each path — carried through
    from the ``<id>`` substitution rather than re-derived from the filename.
    This is what lets excusal match on the id directly, independent of where
    ``<id>`` sits in the template: ``<id>.pdb``, ``rank1_<id>.pdb``, and
    ``<id>/model.pdb`` all yield the same owner_id. owner_id is ``None`` when
    the entry carries no per-id template (legacy bare lists, shared-file
    streams, raw strings) — such a path can never be id-excused.

    The category data may be:
    - A list of file paths (legacy format) -> owner_id None
    - A bare string path (shared-file stream serialized at top level) -> None
    - A dict (serialized DataStream) with a 'files' key whose value is
      either a list (per-id) or a string (shared single artifact)

    When files contains a single '<id>' template, expands it using the 'ids'
    field (which may contain compact patterns like 'name_<0..9>'); each
    expanded path is paired with the id it was built from. For lazy ids the
    bracketed slot becomes a '*' glob and the owner_id is the deterministic
    prefix the glob was built from.
    """
    if isinstance(category_data, dict):
        files = category_data.get('files', [])
        ids = category_data.get('ids', [])
        # Shared-file form: a single str path covers all ids; check it once.
        if isinstance(files, str):
            return [(None, files)] if files else []
        if len(files) == 1 and '<id>' in files[0] and ids:
            template = files[0]
            expanded_ids, is_complete = id_patterns.try_expand_ids(ids)
            if is_complete:
                return [(eid, template.replace('<id>', eid)) for eid in expanded_ids]
            else:
                # Lazy patterns: replace [...] with '*' then expand <..> slots.
                # The glob_id is both the substitution and the owner id (its
                # deterministic prefix is what an excused upstream id matches).
                glob_ids = id_patterns.glob_from_lazy_ids(ids)
                return [(gid, template.replace('<id>', gid)) for gid in glob_ids]
        return [(None, f) for f in files]
    elif isinstance(category_data, list):
        return [(None, f) for f in category_data]
    elif isinstance(category_data, str):
        return [(None, category_data)] if category_data else []
    return []


def check_pairs_exist(pairs: List[Tuple[Optional[str], str]]) -> Tuple[bool, List[Tuple[Optional[str], str]]]:
    """Like :func:`check_files_exist` but over (owner_id, path) pairs.

    Returns (all_exist, missing_pairs) preserving each missing path's owner id.
    """
    missing = [(oid, p) for (oid, p) in pairs if not check_file_exists(p)]
    return len(missing) == 0, missing


# Top-level keys in expected_outputs that are not stream-shaped and must
# not be iterated by the generic stream-coverage loop.
_RESERVED_TOP_LEVEL_KEYS = {'tables', 'output_folder'}


def _load_expected_missing_ids(expected_outputs: Dict[str, Any],
                               step_id: str, tool_name: str) -> List[str]:
    """
    Load IDs whose missing output files are *excused*.

    Reads ``tables.missing.path`` (schema: id | removed_by | kind | cause).
    An ID is excused when:
      * ``removed_by`` is neither this step's id nor its bare tool name —
        the row was propagated from an upstream tool; the file is expected
        to be absent.
      * ``kind == "filter"`` — the current tool dropped this ID on purpose
        (dedup, rename, intentional filter); the file is expected to be
        absent.

    A row whose ``removed_by`` IS this step (step id or bare tool name) and
    whose ``kind == "failure"`` is NOT excused — its missing files indicate
    a real failure and should trigger a FAILED status.

    The primary discriminator is the *step* identifier (``<order>_<Tool>``,
    e.g. ``005_XTB``), not the bare tool name. In a chain of the same tool
    (two Panda steps, two XTB scorers), an upstream-failure row carries the
    upstream step's id (``002_XTB``); comparing it against the current step's
    id (``007_XTB``) correctly treats it as upstream-propagated and excuses
    it, instead of mistaking it for a local failure of this step. The bare
    tool name is also treated as local for back-compatibility with manifests
    written before producers stamped the step id.

    Returns an empty list if no missing table is declared, the file does
    not exist yet, or it lacks an ``id`` column.
    """
    tables = expected_outputs.get('tables')
    if not isinstance(tables, dict):
        return []
    missing_info = tables.get('missing')
    if not missing_info:
        return []
    if isinstance(missing_info, dict):
        missing_path = missing_info.get('path')
    else:
        missing_path = str(missing_info)
    if not missing_path or not os.path.exists(missing_path):
        return []
    try:
        import pandas as pd
        df = pd.read_csv(missing_path)
    except Exception as e:
        print(f"Warning: Could not read {os.path.basename(missing_path)}: {e}")
        return []
    if 'id' not in df.columns or 'removed_by' not in df.columns or 'kind' not in df.columns:
        return []
    local = {step_id, tool_name}
    excused = df[(~df['removed_by'].astype(str).isin(local))
                 | (df['kind'].astype(str) == 'filter')]
    return excused['id'].astype(str).tolist()


def _filter_expected_missing(missing_pairs: List[Tuple[Optional[str], str]],
                             expected_missing_ids: List[str]) -> Tuple[List[str], List[str]]:
    """
    Split missing (owner_id, path) pairs into (unexpected, expected) paths
    given the output ids that are excused.

    Matching is on the path's *owner id* — the id carried through from the
    ``<id>`` substitution — never on the filename. A template can put ``<id>``
    anywhere (``rank1_<id>.pdb``, ``<id>/model.pdb``); the owner id is the same
    regardless, so excusal does not depend on the id being the filename stem.
    A path with no owner id (legacy/shared-file/string entries) is never
    excused.

    A path is "expected-missing" when its owner id is an excused id itself OR
    *descends from* one — i.e. an excused id appears in the owner id's
    suffix-base ancestor chain (``map_table_ids_to_ids``). The descends-from
    test covers a fan-out whose count is unknown at config time:
    ``_remap_missing_to_output_ids`` can only pre-expand a *deterministic*
    fan-out (``prot+lig2_<1..2>``) into the declared output ids; a lazy fan-out
    (``Panda_5[_rank<...>]``) collapses to its prefix ``Panda_5``, so the owner
    id of ``Panda_5_rank001.pdb`` would never exact-match. The ancestor chain
    excuses it, because a child of an upstream-filtered parent legitimately
    cannot exist. Ancestor matching is parent/child only — it walks the owner
    id's own suffix-strip chain, never the sibling tier — so it does not
    over-match a different design sharing only a top-level base
    (``Panda_6_rank001`` is NOT excused by ``Panda_5``). This reuses the suffix
    machinery ``get_mapped_ids`` applies for its child/parent tiers; there is
    one id-matching path.
    """
    if not expected_missing_ids:
        return [p for (_oid, p) in missing_pairs], []
    from biopipelines.id_map_utils import map_table_ids_to_ids
    expected_set = set(expected_missing_ids)
    unexpected = []
    expected = []
    for owner_id, path in missing_pairs:
        if owner_id is None:
            unexpected.append(path)
            continue
        ancestors = map_table_ids_to_ids(owner_id, {"*": "*_<S>"})
        is_expected = any(a in expected_set for a in ancestors)
        (expected if is_expected else unexpected).append(path)
    return unexpected, expected


def _collect_declared_output_ids(expected_outputs: Dict[str, Any]) -> Tuple[List[str], List[str]]:
    """Every declared output id across the tool's streams plus their map_table
    paths (deterministic expansion only; lazy patterns contribute their prefix).

    The ids are the targets for remapping an excused id into this tool's output
    id space — a product tool declares ``prot+lig2`` here even though the
    upstream manifest only names ``lig2``. The map_table paths let the remap use
    provenance matching for renames whose link lives only in a stream map
    (e.g. Panda's ``Panda_1`` -> ``LID_001_1``), not just suffix/product shape.
    """
    ids: List[str] = []
    seen = set()
    map_tables: List[str] = []
    seen_tables = set()
    for category, value in expected_outputs.items():
        if category in _RESERVED_TOP_LEVEL_KEYS or not isinstance(value, dict):
            continue
        mt = value.get('map_table')
        if mt and mt not in seen_tables and os.path.exists(mt):
            seen_tables.add(mt)
            map_tables.append(mt)
        raw_ids = value.get('ids', [])
        if not raw_ids:
            continue
        expanded, _complete = id_patterns.try_expand_ids(raw_ids)
        for eid in expanded:
            if eid not in seen:
                seen.add(eid)
                ids.append(eid)
    return ids, map_tables


def _remap_missing_to_output_ids(expected_missing_ids: List[str],
                                 expected_outputs: Dict[str, Any]) -> List[str]:
    """Map input-axis excused ids into this tool's output id space.

    The missing manifest names ids in the UPSTREAM tool's id space (an input
    axis, e.g. ``lig2``); this tool's expected files are keyed by its OWN output
    ids (a product ``prot+lig2``, a ``_N`` multiplier ``prot_1``, a group key).
    ``get_mapped_ids`` already understands ``+`` products, ``_`` suffixes, and
    provenance, so it bridges the two id spaces. An excused id that matches no
    declared output id is kept as-is (covers 1:1 tools and exact-name matches).
    """
    if not expected_missing_ids:
        return expected_missing_ids
    output_ids, map_tables = _collect_declared_output_ids(expected_outputs)
    if not output_ids:
        return expected_missing_ids
    from biopipelines.id_map_utils import get_mapped_ids
    mapped = get_mapped_ids(expected_missing_ids, output_ids, unique=False,
                            map_table_paths=map_tables or None)
    out: List[str] = []
    seen = set()
    for mid in expected_missing_ids:
        matches = mapped.get(mid, [])
        # No output matched (e.g. an exact-name 1:1 tool, or a stale row) → keep
        # the original id so exact-basename excusal still works.
        for cand in (matches or [mid]):
            if cand not in seen:
                seen.add(cand)
                out.append(cand)
    return out


def _step_id(output_folder: str, tool_name: str) -> str:
    """The step identifier producers write into ``removed_by``.

    Matches ``os.path.basename(self.output_folder)`` used by tools (e.g.
    ``005_XTB``). Falls back to the bare tool name for folders without a
    numeric step prefix, so exact-name producers still match.
    """
    folder_name = os.path.basename(output_folder.rstrip(os.sep))
    if '_' in folder_name and folder_name.split('_')[0].isdigit():
        return folder_name
    return tool_name


def check_expected_outputs(expected_outputs: Dict[str, Any],
                           tool_name: str,
                           output_folder: str = "") -> Tuple[bool, Dict[str, List[str]]]:
    """
    Check that every declared output file exists.

    Iterates *all* stream-shaped top-level keys, not a hardcoded list. A
    stream-shaped entry is either a list of paths (legacy) or a dict with
    ``files`` / ``ids`` (new format) — both handled by ``extract_id_file_pairs``.
    This means streams like ``fasta``, ``msas``, ``images``, ``plots`` are
    checked the same way as the historically-privileged ``structures`` /
    ``compounds`` / ``sequences``.

    Missing-aware: files corresponding to IDs flagged in ``tables.missing``
    as either propagated upstream (``removed_by != tool_name``) or locally
    filtered (``kind == "filter"``) are *expected* to be absent and do not
    count as failures. Local failure rows (``removed_by == tool_name`` and
    ``kind == "failure"``) are NOT excused — their missing files still
    trigger a FAILED status.

    Args:
        expected_outputs: Dictionary with standardized output format
        tool_name: Name of the current tool (used to distinguish local
            failure rows from upstream/local-filter rows in missing.csv)

    Returns:
        Tuple of (all_exist: bool, missing_by_category: Dict[str, List[str]])
    """
    missing_by_category = {}
    all_exist = True

    step_id = _step_id(output_folder, tool_name)
    expected_missing_ids = _load_expected_missing_ids(expected_outputs, step_id, tool_name)
    expected_missing_ids = _remap_missing_to_output_ids(expected_missing_ids, expected_outputs)
    if expected_missing_ids:
        print(f"Found {len(expected_missing_ids)} IDs expected to be missing "
              f"(upstream propagated or local filter)")

    # Every non-reserved key whose value reduces to a file list is a stream.
    for category, value in expected_outputs.items():
        if category in _RESERVED_TOP_LEVEL_KEYS:
            continue
        pairs = extract_id_file_pairs(value)
        if not pairs:
            continue
        exists, missing = check_pairs_exist(pairs)
        if exists:
            continue
        unexpected, expected = _filter_expected_missing(missing, expected_missing_ids)
        if expected:
            print(f"Info: {len(expected)} {category} files missing as expected "
                  f"(filtered out upstream)")
        if unexpected:
            missing_by_category[category] = unexpected
            all_exist = False

    # Tables (handle both old list-of-paths and new dict-of-TableInfo formats).
    if 'tables' in expected_outputs:
        tables = expected_outputs['tables']
        table_files = []

        if isinstance(tables, dict):
            for name, info in tables.items():
                if isinstance(info, dict) and 'path' in info:
                    table_files.append(info['path'])
                else:
                    table_files.append(str(info))
        elif isinstance(tables, list):
            table_files = tables

        if table_files:
            exists, missing = check_files_exist(table_files)
            if not exists:
                missing_by_category['tables'] = missing
                all_exist = False

    return all_exist, missing_by_category

def create_status_file(output_folder: str, tool_name: str, status: str, details: Dict[str, Any] = None) -> str:
    """
    Create a status file indicating completion or failure.
    
    Args:
        output_folder: Tool's output folder
        tool_name: Name of the tool (e.g., "RFdiffusion", "ProteinMPNN")
        status: "COMPLETED" or "FAILED"
        details: Optional details about the status
        
    Returns:
        Path to created status file
    """
    # Get step number from output folder name (e.g., "1_RFdiffusion" -> "1")
    folder_name = os.path.basename(output_folder)
    if '_' in folder_name and folder_name.split('_')[0].isdigit():
        step_number = folder_name.split('_')[0]
        status_filename = f"{step_number}_{tool_name}_{status}"
    else:
        status_filename = f"{tool_name}_{status}"
    
    # Create status file in the parent directory (pipeline level)
    parent_dir = os.path.dirname(output_folder)
    status_file_path = os.path.join(parent_dir, status_filename)
    
    # Create the status file with optional details
    with open(status_file_path, 'w') as f:
        f.write(f"Status: {status}\n")
        f.write(f"Tool: {tool_name}\n")
        f.write(f"Output folder: {output_folder}\n")
        f.write(f"Timestamp: {__import__('datetime').datetime.now().isoformat()}\n")
        
        if details:
            f.write(f"\nDetails:\n")
            if 'missing_files' in details:
                f.write(f"Missing files:\n")
                for category, files in details['missing_files'].items():
                    f.write(f"  {category}:\n")
                    for file_path in files:
                        f.write(f"    - {file_path}\n")
            
            if 'error_message' in details:
                f.write(f"Error: {details['error_message']}\n")
    
    return status_file_path

def check_completion_status(output_folder: str, tool_name: str) -> Optional[str]:
    """
    Check if completion status file already exists.

    Args:
        output_folder: Tool's output folder
        tool_name: Name of the tool

    Returns:
        Status ("COMPLETED" or "FAILED") if exists, None otherwise
    """
    parent_dir = os.path.dirname(output_folder)
    folder_name = os.path.basename(output_folder)

    if '_' in folder_name and folder_name.split('_')[0].isdigit():
        step_number = folder_name.split('_')[0]
        completed_file = os.path.join(parent_dir, f"{step_number}_{tool_name}_COMPLETED")
        failed_file = os.path.join(parent_dir, f"{step_number}_{tool_name}_FAILED")
    else:
        completed_file = os.path.join(parent_dir, f"{tool_name}_COMPLETED")
        failed_file = os.path.join(parent_dir, f"{tool_name}_FAILED")

    if os.path.exists(completed_file):
        return "COMPLETED"
    elif os.path.exists(failed_file):
        return "FAILED"

    return None

def clean_old_status_files(output_folder: str, tool_name: str) -> None:
    """
    Remove old status files (both COMPLETED and FAILED) to allow fresh evaluation.
    
    Args:
        output_folder: Tool's output folder
        tool_name: Name of the tool
    """
    parent_dir = os.path.dirname(output_folder)
    folder_name = os.path.basename(output_folder)
    
    # Determine status file patterns
    if '_' in folder_name and folder_name.split('_')[0].isdigit():
        step_number = folder_name.split('_')[0]
        completed_file = os.path.join(parent_dir, f"{step_number}_{tool_name}_COMPLETED")
        failed_file = os.path.join(parent_dir, f"{step_number}_{tool_name}_FAILED")
    else:
        completed_file = os.path.join(parent_dir, f"{tool_name}_COMPLETED")
        failed_file = os.path.join(parent_dir, f"{tool_name}_FAILED")
    
    # Derive warning file path
    warning_file = completed_file.replace("_COMPLETED", "_WARNING")

    # Remove old status files if they exist
    for status_file in [completed_file, failed_file, warning_file]:
        if os.path.exists(status_file):
            try:
                os.remove(status_file)
                print(f"Removed old status file: {os.path.basename(status_file)}")
            except OSError as e:
                print(f"Warning: Could not remove {os.path.basename(status_file)}: {e}")

def main():
    parser = argparse.ArgumentParser(description="Check pipeline tool completion status")
    parser.add_argument("output_folder", help="Tool's output folder path")
    parser.add_argument("tool_name", help="Name of the tool (e.g., RFdiffusion, ProteinMPNN)")
    parser.add_argument("expected_outputs", help="JSON file or string with expected outputs")
    parser.add_argument("--check-only", action="store_true", 
                       help="Only check status, don't create status files")
    parser.add_argument("--force", action="store_true",
                       help="Force check even if status file exists")
    parser.add_argument("--job-name", 
                       help="Job name for filter manifest lookup (auto-detected if not provided)")
    
    args = parser.parse_args()
    
    # Auto-detect job name from output folder if not provided
    job_name = args.job_name
    if not job_name:
        # Try to extract from folder structure or use tool name as fallback
        folder_name = os.path.basename(args.output_folder)
        if '_' in folder_name and len(folder_name.split('_')) > 1:
            # Format like "1_ToolName" - use ToolName as job name
            job_name = '_'.join(folder_name.split('_')[1:])
        else:
            job_name = args.tool_name
    
    # Clean up old status files to allow fresh evaluation
    clean_old_status_files(args.output_folder, args.tool_name)
    
    # Check if we already have a status (unless forced)
    if not args.force:
        existing_status = check_completion_status(args.output_folder, args.tool_name)
        if existing_status == "COMPLETED":
            print(f"Tool {args.tool_name} already completed")
            sys.exit(0)
        elif existing_status == "FAILED":
            print(f"Tool {args.tool_name} previously failed")
            sys.exit(1)
    
    # Load expected outputs
    try:
        if os.path.isfile(args.expected_outputs):
            # Load from JSON file
            with open(args.expected_outputs, 'r') as f:
                expected_outputs = json.load(f)
        else:
            # Parse as JSON string
            expected_outputs = json.loads(args.expected_outputs)
    except (json.JSONDecodeError, FileNotFoundError) as e:
        print(f"Error parsing expected outputs: {e}", file=sys.stderr)
        sys.exit(2)

    # Unwrap envelope format if present (tool_name, tool_class, output_structure wrapper)
    if 'output_structure' in expected_outputs and 'tool_name' in expected_outputs:
        expected_outputs = expected_outputs['output_structure']

    success, missing_by_category = check_expected_outputs(expected_outputs, args.tool_name, args.output_folder)

    if success:
        print(f"Required outputs found for {args.tool_name}")
        if not args.check_only:
            status_file = create_status_file(args.output_folder, args.tool_name, "COMPLETED")
            print(f"Created completed status file: {os.path.basename(status_file)}")
        sys.exit(0)
    else:
        print(f"Missing outputs for {args.tool_name}:")
        for category, files in missing_by_category.items():
            print(f"  {category}:")
            for file_path in files:
                print(f"    - {file_path}")

        if not args.check_only:
            details = {"missing_files": missing_by_category}
            status_file = create_status_file(args.output_folder, args.tool_name, "FAILED", details)
            print(f"Created failure status file: {os.path.basename(status_file)}")
        sys.exit(1)

if __name__ == "__main__":
    main()