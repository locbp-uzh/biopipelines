#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Utility functions for ID mapping pattern matching.

This module provides a unified approach to parsing and applying ID mapping patterns
across different tools. The pattern format is {"*": "<pattern>"} where <pattern>
describes how to map structure IDs to table IDs.

Pattern syntax:
- "*" in both key and value represents the base ID
- "<N>" represents a numeric suffix (one or more digits)
- "<S>" represents any non-underscore segment (alphanumeric suffixes like "19A")
- Everything else is treated as literal text
- Patterns are applied RECURSIVELY until no more matches

Examples:
- {"*": "*"}           -> no mapping (identity)
- {"*": "*_<S>"}       -> strip ALL trailing "_segment" suffixes recursively (DEFAULT)
                        -> "protein_1_19A" -> "protein_1" -> "protein"
                        -> "rifampicin_1_2_3" -> "rifampicin"
- {"*": "*_<N>"}       -> strip only trailing "_123" numeric suffixes recursively
                        -> "protein_1_19A" -> NO MATCH (19A is not purely numeric)
- {"*": "*-<N>"}       -> strip ALL trailing "-123" suffixes recursively
- {"*": "*_seq_<N>"}   -> strip ALL trailing "_seq_123" patterns recursively

Note: The default pattern {"*": "*_<S>"} works for any kind of suffix (numeric
or alphanumeric). Use {"*": "*_<N>"} if you need to restrict to numeric-only suffixes.

Multi-axis combinatorics IDs use "+" as a separator (e.g., "prot1+lig1" for a
protein-ligand pair). When generating candidate IDs, "+" components are expanded
so that "prot1+lig1" also yields "prot1" and "lig1" as candidates. Each component
also undergoes recursive suffix stripping. This makes provenance columns unnecessary
for matching multi-axis IDs against their constituent input IDs.

Matching priority order (when using get_mapped_ids with unique=True):
    1. Exact match (source_id == target_id)
    2. Exact provenance match (via map_table provenance columns)
    3. For each candidate id in [source_id, *provenance_ids]:
       a. Child match (target derives from candidate, closest first)
       b. Parent match (candidate derives from target, closest first)
    4. For each candidate id in [source_id, *provenance_ids]:
       Sibling match (common ancestor, closest combined distance first)
"""

import re
from functools import lru_cache
from typing import Dict, List, Optional, Tuple, Union


def parse_id_map_pattern(id_map: Dict[str, str]) -> Optional[re.Pattern]:
    """
    Parse ID mapping pattern and generate regex for extracting base ID.

    Args:
        id_map: ID mapping dictionary (e.g., {"*": "*_<N>_<N>"})

    Returns:
        Compiled regex pattern that captures the base ID, or None if no mapping needed

    Examples:
        >>> pattern = parse_id_map_pattern({"*": "*_<N>"})
        >>> match = pattern.match("rifampicin_1")
        >>> match.group(1)
        'rifampicin'

        >>> pattern = parse_id_map_pattern({"*": "*_<S>"})
        >>> match = pattern.match("protein_19A")
        >>> match.group(1)
        'protein'

        >>> pattern = parse_id_map_pattern({"*": "*-seq-<N>"})
        >>> match = pattern.match("protein-seq-42")
        >>> match.group(1)
        'protein'
    """
    if not id_map or "*" not in id_map:
        return None

    pattern_str = id_map["*"]

    # Check if pattern is just "*" (identity mapping - no transformation)
    if pattern_str == "*":
        return None

    # Verify pattern starts with "*"
    if not pattern_str.startswith("*"):
        raise ValueError(f"ID map pattern must start with '*': {pattern_str}")

    # Extract the suffix pattern (everything after the first "*")
    suffix_pattern = pattern_str[1:]  # Remove leading "*"

    if not suffix_pattern:
        # Pattern is just "*" - no mapping
        return None

    # Build regex pattern by:
    # 1. Replace "<N>" and "<S>" placeholders with temporary markers
    # 2. Escape special regex characters in the literal parts
    # 3. Replace markers back with regex patterns
    # 4. Anchor to end of string

    # Use placeholders that won't appear in normal text
    DIGIT_PLACEHOLDER = '\x00DIGIT_PATTERN\x00'
    SEGMENT_PLACEHOLDER = '\x00SEGMENT_PATTERN\x00'

    # Replace <N> and <S> with placeholders before escaping
    temp_pattern = suffix_pattern.replace('<N>', DIGIT_PLACEHOLDER)
    temp_pattern = temp_pattern.replace('<S>', SEGMENT_PLACEHOLDER)

    # Escape regex special characters
    regex_suffix = re.escape(temp_pattern)

    # Replace placeholders with regex patterns
    regex_suffix = regex_suffix.replace(re.escape(DIGIT_PLACEHOLDER), r'\d+')
    regex_suffix = regex_suffix.replace(re.escape(SEGMENT_PLACEHOLDER), r'[^_]+')

    # Build full pattern: capture everything before the suffix
    full_pattern = f'^(.+){regex_suffix}$'

    try:
        compiled_pattern = re.compile(full_pattern)
        return compiled_pattern
    except re.error as e:
        raise ValueError(f"Invalid ID map pattern '{pattern_str}': {e}")


@lru_cache(maxsize=200_000)
def _strip_suffixes_recursive(
    start_id: str,
    delimiter: str,
    has_s: bool,
    literal_part: str,
) -> Tuple[str, ...]:
    """
    Recursively strip suffixes from an ID according to the pattern parameters.

    Returns a tuple of progressively stripped IDs (NOT including start_id itself).

    Pure function of its arguments; memoized because it is called many times on
    a small set of distinct ids during cross-tool ID matching.
    """
    results = []
    current_id = start_id

    while True:
        parts = current_id.split(delimiter)

        if len(parts) <= 1:
            break

        stripped = False

        if literal_part:
            if len(parts) >= 2:
                last_part = parts[-1]
                second_last = parts[-2] if len(parts) >= 2 else ""

                if has_s:
                    if last_part and second_last == literal_part:
                        parts = parts[:-2]
                        stripped = True
                else:
                    if last_part.isdigit() and second_last == literal_part:
                        parts = parts[:-2]
                        stripped = True
        else:
            if has_s:
                parts = parts[:-1]
                stripped = True
            else:
                if parts[-1].isdigit():
                    parts = parts[:-1]
                    stripped = True

        if not stripped or not parts:
            break

        current_id = delimiter.join(parts)
        if not results or current_id != results[-1]:
            results.append(current_id)

    return tuple(results)


# Separator used for multi-axis combinatorics IDs (e.g., "prot1+lig1")
MULTI_AXIS_SEPARATOR = '+'


def map_table_ids_to_ids(structure_id: str, id_map: Dict[str, str]) -> list:
    """
    Generate all candidate IDs for lookup by recursively applying id_map pattern.

    This function generates a list of IDs to try when looking up a structure ID in a table.
    It starts with the original ID and progressively strips suffixes according to the
    id_map pattern, generating all intermediate IDs in priority order (most specific first).

    For multi-axis combinatorics IDs (containing "+"), the function generates all
    contiguous sub-sequences of "+" components (from longest to shortest), plus
    suffix-stripped variants of each. This allows "prot1+lig1+lig2" to match against
    "prot1+lig1", "lig1+lig2", "prot1", "lig1", or "lig2".

    Args:
        structure_id: Structure ID (e.g., "RFDAA_Hit_Screen_007_1_1" or "prot1+lig1_2")
        id_map: ID mapping dictionary (e.g., {"*": "*_<N>"} or {"*": "*_<S>"})

    Returns:
        List of candidate IDs to try, from most specific to least specific.

    Examples:
        >>> map_table_ids_to_ids("rifampicin_1_2", {"*": "*_<N>"})
        ['rifampicin_1_2', 'rifampicin_1', 'rifampicin']

        >>> map_table_ids_to_ids("protein_1_19A", {"*": "*_<S>"})
        ['protein_1_19A', 'protein_1', 'protein']

        >>> map_table_ids_to_ids("prot1+lig1", {"*": "*_<S>"})
        ['prot1+lig1', 'prot1', 'lig1']

        >>> map_table_ids_to_ids("prot1+lig1_2", {"*": "*_<S>"})
        ['prot1+lig1_2', 'prot1+lig1', 'prot1', 'lig1_2', 'lig1']

        >>> map_table_ids_to_ids("protein-seq-42", {"*": "*-seq-<N>"})
        ['protein-seq-42', 'protein']

        >>> map_table_ids_to_ids("no_change", {"*": "*"})
        ['no_change']
    """
    if not id_map or "*" not in id_map:
        return [structure_id]
    # Delegate to a memoized core keyed on (id, pattern). map_table_ids_to_ids
    # is a pure function and is called many times on a small set of distinct ids
    # during cross-tool ID matching; caching keeps repeated merges cheap. Return
    # a copy so callers that mutate the result (e.g. .extend) don't corrupt the
    # cache entry.
    return list(_map_table_ids_to_ids_cached(structure_id, id_map["*"]))


@lru_cache(maxsize=200_000)
def _map_table_ids_to_ids_cached(structure_id: str, pattern_str: str) -> Tuple[str, ...]:
    # Check if pattern is just "*" (identity mapping)
    if pattern_str == "*":
        return (structure_id,)

    # Verify pattern starts with "*"
    if not pattern_str.startswith("*"):
        return (structure_id,)

    # Extract the suffix pattern
    suffix_pattern = pattern_str[1:]  # Remove leading "*"

    if not suffix_pattern:
        return (structure_id,)

    # Determine placeholder type: <S> (any segment) or <N> (digits only)
    has_s = '<S>' in suffix_pattern
    has_n = '<N>' in suffix_pattern

    if not has_s and not has_n:
        return (structure_id,)

    # Find the position of the placeholder
    placeholder = '<S>' if has_s else '<N>'
    p_pos = suffix_pattern.find(placeholder)

    # Extract delimiter
    if p_pos == 0:
        return (structure_id,)  # Pattern like "*<N>" doesn't make sense

    delimiter = suffix_pattern[p_pos - 1]
    literal_with_delim = suffix_pattern[:p_pos - 1] if p_pos > 1 else ""
    if literal_with_delim.startswith(delimiter):
        literal_part = literal_with_delim[1:]
    else:
        literal_part = literal_with_delim

    # Generate all intermediate IDs by progressively stripping suffixes
    candidates = [structure_id]

    if MULTI_AXIS_SEPARATOR in structure_id:
        # For multi-axis IDs, generate all contiguous sub-sequences of + components
        # with suffix-stripped variants, ordered by decreasing length (more specific first)
        parts = structure_id.split(MULTI_AXIS_SEPARATOR)
        n = len(parts)
        seen = {structure_id}

        # Generate sub-sequences from longest to shortest
        for length in range(n, 0, -1):
            for start in range(n - length + 1):
                sub_parts = parts[start:start + length]
                sub_id = MULTI_AXIS_SEPARATOR.join(sub_parts)
                if sub_id not in seen:
                    candidates.append(sub_id)
                    seen.add(sub_id)
                # Suffix-strip the last component to generate variants
                last = sub_parts[-1]
                for stripped_last in _strip_suffixes_recursive(last, delimiter, has_s, literal_part):
                    variant_parts = sub_parts[:-1] + [stripped_last]
                    variant = MULTI_AXIS_SEPARATOR.join(variant_parts)
                    if variant not in seen:
                        candidates.append(variant)
                        seen.add(variant)
    else:
        stripped = _strip_suffixes_recursive(structure_id, delimiter, has_s, literal_part)
        candidates.extend(stripped)

    return tuple(candidates)


# Per-process cache of provenance lookup tables, keyed by map_table CSV
# path. Each entry maps `source_id` (string) to its list of provenance
# identity values pulled from `<stream>.id` columns. Built lazily once
# per CSV in _get_provenance_ids; the underlying CSV is itself cached
# in biopipelines.biopipelines_io._table_cache, so an existing entry
# stays valid as long as that DataFrame is still cached.
_provenance_lookup_cache: Dict[str, Dict[str, List[str]]] = {}


def _build_provenance_lookup(df) -> Dict[str, List[str]]:
    """Build {source_id -> [provenance ids]} from a map_table DataFrame.

    Provenance columns are any '<stream>.id' columns other than 'id'
    itself. Each row's id column maps to the list of non-empty,
    non-self provenance values across those columns.
    """
    lookup: Dict[str, List[str]] = {}
    if "id" not in df.columns:
        return lookup
    prov_cols = [c for c in df.columns if c.endswith(".id") and c != "id"]
    if not prov_cols:
        return lookup
    ids = [str(v) for v in df["id"].tolist()]
    cols = {c: [str(v) for v in df[c].tolist()] for c in prov_cols}
    for i, sid in enumerate(ids):
        prov: List[str] = []
        for c in prov_cols:
            val = cols[c][i]
            if val and val != sid and val != "nan":
                prov.append(val)
        if prov:
            lookup[sid] = prov
    return lookup


def _get_provenance_ids(source_id, map_table_paths):
    """
    Get all provenance identity values for a source_id from map_tables.

    Reads provenance columns ({stream}.id) from map_table CSVs and returns
    all alternate identities for the given source_id, without filtering
    against any target set.

    Args:
        source_id: The ID to look up
        map_table_paths: List of map_table CSV paths to search

    Returns:
        List of provenance identity strings (may be empty)
    """
    if not map_table_paths:
        return []
    from biopipelines.biopipelines_io import _table_cache
    import pandas as pd
    import os
    provenance_ids: List[str] = []
    for path in map_table_paths:
        if not path or not os.path.exists(path):
            continue
        # Build / reuse the per-source-id lookup index for this CSV.
        if path not in _provenance_lookup_cache:
            if path in _table_cache:
                df = _table_cache[path]
            else:
                try:
                    df = pd.read_csv(path)
                    _table_cache[path] = df
                except Exception:
                    continue
            _provenance_lookup_cache[path] = _build_provenance_lookup(df)
        lookup = _provenance_lookup_cache[path]
        provenance_ids.extend(lookup.get(source_id, ()))
    # Deduplicate while preserving order (a few CSVs may carry the same id).
    seen = set()
    deduped: List[str] = []
    for pid in provenance_ids:
        if pid not in seen:
            seen.add(pid)
            deduped.append(pid)
    return deduped


@lru_cache(maxsize=4096)
def _build_target_index(target_ids: Tuple[str, ...], pattern_str: str):
    """Build the target-side lookup structures used by get_mapped_ids.

    Returns ``(target_set, target_bases_cache, base_to_targets)`` where:
      * ``target_set`` — set of target ids (for exact / membership checks).
      * ``target_bases_cache`` — {target_id: tuple of its suffix-base ids,
        most-specific first} (for sibling matching).
      * ``base_to_targets`` — inverted index {base_id: [(target_id, distance), …]}
        with insertion order preserved per base (historical tie-break:
        first-encountered target wins).

    Memoized on (target_ids, pattern) so that repeated calls against the
    same target list — e.g. resolving many source ids against one tool's
    output, or many per-merge lookups in a campaign — build the index once
    instead of once per call. The returned structures must be treated as
    read-only by callers (get_mapped_ids only reads them or copies before
    mutating).
    """
    id_map = {"*": pattern_str}
    target_set = set(target_ids)
    target_bases_cache: Dict[str, Tuple[str, ...]] = {}
    base_to_targets: Dict[str, List[Tuple[str, int]]] = {}
    for target_id in target_ids:
        bases = tuple(map_table_ids_to_ids(target_id, id_map))
        target_bases_cache[target_id] = bases
        for distance, base in enumerate(bases):
            base_to_targets.setdefault(base, []).append((target_id, distance))
    return target_set, target_bases_cache, base_to_targets


def get_mapped_ids(
    source_ids: List[str],
    target_ids: List[str],
    id_map: Dict[str, str] = None,
    unique: bool = True,
    map_table_paths: Optional[List[str]] = None
) -> Union[Dict[str, Optional[str]], Dict[str, List[str]]]:
    """
    Match source IDs to target IDs using id_map patterns and provenance.

    For each source_id, finds target_ids that match according to the id_map.
    Matching uses a priority-based strategy:
    1. Exact match: source_id == target_id (highest priority)
    2. Provenance match (via map_table provenance columns — most reliable signal)
    3. Child match: target derives from source (target = source + suffix)
    4. Parent match: source derives from target (source = target + suffix)
    5. Sibling match: common ancestor

    Args:
        source_ids: List of source IDs to match from
        target_ids: List of target IDs to match against
        id_map: ID mapping pattern (default: {"*": "*_<S>"})
                Supports recursive suffix matching
        unique: If True (default), return the single most specific match or None.
                If False, return list of all matches.
        map_table_paths: Optional list of map_table CSV paths for provenance-based
                         matching. When suffix-based strategies fail (e.g., Panda
                         auto-renamed IDs like Panda_1), provenance columns in
                         map_tables (e.g., structures.id) are used to trace back
                         to original source IDs.

    Returns:
        If unique=True: Dict mapping each source_id to best matching target_id (or None)
        If unique=False: Dict mapping each source_id to list of matching target_ids

    Priority order (when unique=True):
        1. Exact match (source_id == target_id)
        2. Exact provenance match (via map_table provenance columns)
        3. For each candidate in [source_id, *provenance_ids]:
           a. Child match (target derives from candidate, closest first)
           b. Parent match (candidate derives from target, closest first)
        4. For each candidate in [source_id, *provenance_ids]:
           Sibling match (common ancestor, closest combined distance first)

    Examples:
        >>> # Exact match
        >>> get_mapped_ids(["protein_1"], ["protein_1"])
        {'protein_1': 'protein_1'}

        >>> # Child match (target = source + _N)
        >>> get_mapped_ids(["protein_1"], ["protein_1_1", "protein_1_2"])
        {'protein_1': 'protein_1_1'}

        >>> # Parent match (source = target + _N)
        >>> get_mapped_ids(["protein_1_1"], ["protein_1"])
        {'protein_1_1': 'protein_1'}

        >>> # Sibling match (common ancestor)
        >>> get_mapped_ids(["protein_1_1"], ["protein_1_2"])
        {'protein_1_1': 'protein_1_2'}

        >>> # unique=False - returns all matches
        >>> get_mapped_ids(["protein_1"], ["protein_1_1", "protein_1_2"], unique=False)
        {'protein_1': ['protein_1_1', 'protein_1_2']}

        >>> # Provenance match (Panda auto-renamed IDs)
        >>> get_mapped_ids(["Panda_1"], ["LID_001_1"], map_table_paths=["/path/to/map.csv"])
        {'Panda_1': 'LID_001_1'}
    """
    if id_map is None:
        id_map = {"*": "*_<S>"}

    result = {}

    # Build (or reuse) the target-side index. This is keyed only on the
    # target id list and the pattern, so callers that invoke get_mapped_ids
    # many times against the *same* target list pay the O(N) index build only
    # once instead of once per call.
    target_set, target_bases_cache, base_to_targets = _build_target_index(
        tuple(target_ids), id_map.get("*", "*_<S>")
    )

    for source_id in source_ids:
        # Priority 1: Exact match
        if source_id in target_set:
            if unique:
                result[source_id] = source_id
            else:
                # Collect exact match + all other targets whose base list
                # contains the source id.
                matches = [source_id]
                seen_match = {source_id}
                for tid, _dist in base_to_targets.get(source_id, ()):
                    if tid not in seen_match:
                        matches.append(tid)
                        seen_match.add(tid)
                result[source_id] = matches
            continue

        # Priority 2: Provenance match via map_table columns
        if map_table_paths:
            from biopipelines.biopipelines_io import resolve_id_by_provenance
            matched = resolve_id_by_provenance(source_id, target_set, map_table_paths)
            if matched is not None:
                if unique:
                    result[source_id] = matched
                else:
                    result[source_id] = [matched]
                continue

        # Build candidate IDs: source_id + provenance identities
        candidate_ids = [source_id]
        if map_table_paths:
            for prov_id in _get_provenance_ids(source_id, map_table_paths):
                if prov_id not in candidate_ids:
                    candidate_ids.append(prov_id)

        # Priority 3: Child/Parent match for each candidate id
        # For each candidate, try child then parent before moving to next candidate
        found = False
        for cand_id in candidate_ids:
            # 3a: Target is a "child" of candidate (target = cand_id + suffix)
            # base_to_targets gives us all targets whose suffix-base list
            # contains cand_id, paired with the distance to it.
            child_matches = list(base_to_targets.get(cand_id, ()))

            if child_matches:
                child_matches.sort(key=lambda x: x[1])
                if unique:
                    result[source_id] = child_matches[0][0]
                else:
                    result[source_id] = [m[0] for m in child_matches]
                found = True
                break

            # 3b: Candidate is a "child" of target (cand_id = target + suffix)
            cand_bases = map_table_ids_to_ids(cand_id, id_map)
            parent_matches = []
            for i, base in enumerate(cand_bases[1:], start=1):
                if base in target_set:
                    parent_matches.append((base, i))

            if parent_matches:
                parent_matches.sort(key=lambda x: x[1])
                if unique:
                    result[source_id] = parent_matches[0][0]
                else:
                    result[source_id] = [m[0] for m in parent_matches]
                found = True
                break

        if found:
            continue

        # Priority 4: Sibling match for each candidate id
        # For each base reachable from cand_id, base_to_targets gives us
        # the targets that share it together with their target-side
        # distance — no per-target scan required.
        sibling_matches = []
        for cand_id in candidate_ids:
            cand_bases = map_table_ids_to_ids(cand_id, id_map)
            for source_dist, common_base in enumerate(cand_bases):
                shared = base_to_targets.get(common_base)
                if not shared:
                    continue
                for target_id, target_dist in shared:
                    combined_dist = source_dist + target_dist
                    sibling_matches.append((target_id, combined_dist, source_dist))

        if sibling_matches:
            # Sort by combined distance, then by source distance
            sibling_matches.sort(key=lambda x: (x[1], x[2]))
            if unique:
                result[source_id] = sibling_matches[0][0]
            else:
                # Remove duplicates while preserving order
                seen = set()
                unique_matches = []
                for m in sibling_matches:
                    if m[0] not in seen:
                        seen.add(m[0])
                        unique_matches.append(m[0])
                result[source_id] = unique_matches
            continue

        # No match found
        if unique:
            result[source_id] = None
        else:
            result[source_id] = []

    return result


def get_mapped_ids_grouped(
    source_ids: list,
    target_ids: list,
    id_map: Dict[str, str] = None
) -> Dict[str, Dict[str, list]]:
    """
    Match source IDs to target IDs, grouped by their common base.

    Similar to get_mapped_ids but returns results grouped by the common
    ancestor ID, which is useful when you need to know the relationship
    hierarchy.

    Args:
        source_ids: List of source IDs to match from
        target_ids: List of target IDs to match against
        id_map: ID mapping pattern (default: {"*": "*_<S>"})

    Returns:
        Dict with structure: {base_id: {"sources": [...], "targets": [...]}}

    Example:
        >>> get_mapped_ids_grouped(
        ...     ["protein_1", "protein_2"],
        ...     ["protein_1_1", "protein_1_2", "protein_2_1"],
        ...     {"*": "*_<S>"}
        ... )
        {
            'protein_1': {'sources': ['protein_1'], 'targets': ['protein_1_1', 'protein_1_2']},
            'protein_2': {'sources': ['protein_2'], 'targets': ['protein_2_1']}
        }
    """
    if id_map is None:
        id_map = {"*": "*_<S>"}

    # Build base-to-ids mapping for both source and target
    source_by_base = {}
    for sid in source_ids:
        bases = map_table_ids_to_ids(sid, id_map)
        root = bases[-1] if bases else sid  # Most stripped version
        if root not in source_by_base:
            source_by_base[root] = []
        source_by_base[root].append(sid)

    target_by_base = {}
    for tid in target_ids:
        bases = map_table_ids_to_ids(tid, id_map)
        root = bases[-1] if bases else tid
        if root not in target_by_base:
            target_by_base[root] = []
        target_by_base[root].append(tid)

    # Combine into result
    all_bases = set(source_by_base.keys()) | set(target_by_base.keys())
    result = {}
    for base in all_bases:
        result[base] = {
            "sources": source_by_base.get(base, []),
            "targets": target_by_base.get(base, [])
        }

    return result
