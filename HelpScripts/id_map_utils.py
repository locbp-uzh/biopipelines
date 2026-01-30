#!/usr/bin/env python3
"""
Utility functions for ID mapping pattern matching.

This module provides a unified approach to parsing and applying ID mapping patterns
across different tools. The pattern format is {"*": "<pattern>"} where <pattern>
describes how to map structure IDs to table IDs.

Pattern syntax:
- "*" in both key and value represents the base ID
- "<N>" represents a numeric suffix (one or more digits)
- Everything else is treated as literal text
- Patterns are applied RECURSIVELY until no more matches

Examples:
- {"*": "*"}           → no mapping (identity)
- {"*": "*_<N>"}       → strip ALL trailing "_123" suffixes recursively
                        → "rifampicin_1_2_3" → "rifampicin"
- {"*": "*-<N>"}       → strip ALL trailing "-123" suffixes recursively
- {"*": "*_seq_<N>"}   → strip ALL trailing "_seq_123" patterns recursively

Note: The default pattern {"*": "*_<N>"} works for any number of numeric suffixes,
so you don't need to specify {"*": "*_<N>_<N>"} or similar patterns.
"""

import re
from typing import Dict, List, Optional, Union


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

        >>> pattern = parse_id_map_pattern({"*": "*_<N>_<N>"})
        >>> match = pattern.match("rifampicin_1_2")
        >>> match.group(1)
        'rifampicin'

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
    # 1. Replace "<N>" placeholder with a temporary marker
    # 2. Escape special regex characters in the literal parts
    # 3. Replace marker back with "\d+" pattern
    # 4. Anchor to end of string

    # Use a placeholder that won't appear in normal text
    PLACEHOLDER = '\x00DIGIT_PATTERN\x00'

    # Replace <N> with placeholder before escaping
    temp_pattern = suffix_pattern.replace('<N>', PLACEHOLDER)

    # Escape regex special characters
    regex_suffix = re.escape(temp_pattern)

    # Replace placeholder with digit pattern
    regex_suffix = regex_suffix.replace(re.escape(PLACEHOLDER), r'\d+')

    # Build full pattern: capture everything before the suffix
    full_pattern = f'^(.+){regex_suffix}$'

    try:
        compiled_pattern = re.compile(full_pattern)
        return compiled_pattern
    except re.error as e:
        raise ValueError(f"Invalid ID map pattern '{pattern_str}': {e}")


def map_table_ids_to_ids(structure_id: str, id_map: Dict[str, str]) -> list:
    """
    Generate all candidate IDs for lookup by recursively applying id_map pattern.

    This function generates a list of IDs to try when looking up a structure ID in a table.
    It starts with the original ID and progressively strips suffixes according to the
    id_map pattern, generating all intermediate IDs in priority order (most specific first).

    Args:
        structure_id: Structure ID (e.g., "RFDAA_Hit_Screen_007_1_1")
        id_map: ID mapping dictionary (e.g., {"*": "*_<N>"})

    Returns:
        List of candidate IDs to try, from most specific to least specific
        (e.g., ["RFDAA_Hit_Screen_007_1_1", "RFDAA_Hit_Screen_007_1", "RFDAA_Hit_Screen_007", "RFDAA_Hit_Screen"])

    Examples:
        >>> map_table_ids_to_ids("rifampicin_1_2", {"*": "*_<N>"})
        ['rifampicin_1_2', 'rifampicin_1', 'rifampicin']

        >>> map_table_ids_to_ids("protein-seq-42", {"*": "*-seq-<N>"})
        ['protein-seq-42', 'protein']

        >>> map_table_ids_to_ids("no_change", {"*": "*"})
        ['no_change']
    """
    if not id_map or "*" not in id_map:
        return [structure_id]

    pattern_str = id_map["*"]

    # Check if pattern is just "*" (identity mapping)
    if pattern_str == "*":
        return [structure_id]

    # Verify pattern starts with "*"
    if not pattern_str.startswith("*"):
        return [structure_id]

    # Extract the suffix pattern
    suffix_pattern = pattern_str[1:]  # Remove leading "*"

    if not suffix_pattern:
        return [structure_id]

    # Find the first occurrence of <N>
    n_pos = suffix_pattern.find('<N>')
    if n_pos == -1:
        return [structure_id]

    # Extract delimiter
    if n_pos == 0:
        return [structure_id]  # Pattern like "*<N>" doesn't make sense

    delimiter = suffix_pattern[n_pos - 1]
    literal_with_delim = suffix_pattern[:n_pos - 1] if n_pos > 1 else ""
    if literal_with_delim.startswith(delimiter):
        literal_part = literal_with_delim[1:]
    else:
        literal_part = literal_with_delim

    # Generate all intermediate IDs by progressively stripping suffixes
    candidates = [structure_id]  # Start with original
    current_id = structure_id

    # Recursively strip suffixes
    while True:
        parts = current_id.split(delimiter)

        if len(parts) <= 1:
            break

        stripped = False

        if literal_part:
            # Pattern like "*-seq-<N>"
            # Check if we can strip literal + number at the end
            if len(parts) >= 2:
                last_part = parts[-1]
                second_last = parts[-2] if len(parts) >= 2 else ""

                if last_part.isdigit() and second_last == literal_part:
                    # Remove both literal and number parts
                    parts = parts[:-2]
                    stripped = True
        else:
            # Pattern like "*_<N>": strip trailing numeric part
            if parts[-1].isdigit():
                parts = parts[:-1]
                stripped = True

        if not stripped or not parts:
            break

        # Rejoin and add to candidates
        current_id = delimiter.join(parts)
        if current_id != candidates[-1]:  # Avoid duplicates
            candidates.append(current_id)

    return candidates


def get_mapped_ids(
    source_ids: List[str],
    target_ids: List[str],
    id_map: Dict[str, str] = None,
    unique: bool = True
) -> Union[Dict[str, Optional[str]], Dict[str, List[str]]]:
    """
    Match source IDs to target IDs using id_map patterns.

    For each source_id, finds target_ids that match according to the id_map.
    Matching uses a priority-based strategy:
    1. Exact match: source_id == target_id (highest priority)
    2. Target is derived from source: target_id = source_id + suffix
    3. Source is derived from target: source_id = target_id + suffix
    4. Common ancestor: both source and target derive from same base

    Args:
        source_ids: List of source IDs to match from
        target_ids: List of target IDs to match against
        id_map: ID mapping pattern (default: {"*": "*_<N>"})
                Supports recursive suffix matching
        unique: If True (default), return the single most specific match or None.
                If False, return list of all matches.

    Returns:
        If unique=True: Dict mapping each source_id to best matching target_id (or None)
        If unique=False: Dict mapping each source_id to list of matching target_ids

    Priority order (when unique=True):
        1. Exact match (source_id == target_id)
        2. Child match (target derives from source, closest first)
        3. Parent match (source derives from target, closest first)
        4. Sibling match (common ancestor, closest combined distance first)

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
    """
    if id_map is None:
        id_map = {"*": "*_<N>"}

    target_set = set(target_ids)
    result = {}

    # Pre-compute base forms for all targets (for sibling matching)
    target_bases_cache = {}
    for target_id in target_ids:
        target_bases_cache[target_id] = map_table_ids_to_ids(target_id, id_map)

    for source_id in source_ids:
        # Priority 1: Exact match
        if source_id in target_set:
            if unique:
                result[source_id] = source_id
            else:
                # Collect exact match + all other matches
                matches = [source_id]
                for target_id in target_ids:
                    if target_id == source_id:
                        continue
                    if source_id in target_bases_cache[target_id]:
                        matches.append(target_id)
                result[source_id] = matches
            continue

        # Priority 2: Target is a "child" of source (target = source + suffix)
        child_matches = []
        for target_id in target_ids:
            target_bases = target_bases_cache[target_id]
            if source_id in target_bases:
                distance = target_bases.index(source_id)
                child_matches.append((target_id, distance))

        if child_matches:
            child_matches.sort(key=lambda x: x[1])
            if unique:
                result[source_id] = child_matches[0][0]
            else:
                result[source_id] = [m[0] for m in child_matches]
            continue

        # Priority 3: Source is a "child" of target (source = target + suffix)
        source_bases = map_table_ids_to_ids(source_id, id_map)
        parent_matches = []
        for i, base in enumerate(source_bases[1:], start=1):
            if base in target_set:
                parent_matches.append((base, i))

        if parent_matches:
            parent_matches.sort(key=lambda x: x[1])
            if unique:
                result[source_id] = parent_matches[0][0]
            else:
                result[source_id] = [m[0] for m in parent_matches]
            continue

        # Priority 4: Sibling match (source and target share common ancestor)
        # Find targets that share any base with source
        source_bases_set = set(source_bases)
        sibling_matches = []
        for target_id in target_ids:
            target_bases = target_bases_cache[target_id]
            # Find common ancestors
            common = source_bases_set & set(target_bases)
            if common:
                # Find the most specific common ancestor (closest to both)
                for common_base in common:
                    source_dist = source_bases.index(common_base)
                    target_dist = target_bases.index(common_base)
                    # Combined distance: prefer matches where both are close to ancestor
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
        id_map: ID mapping pattern (default: {"*": "*_<N>"})

    Returns:
        Dict with structure: {base_id: {"sources": [...], "targets": [...]}}

    Example:
        >>> get_mapped_ids_grouped(
        ...     ["protein_1", "protein_2"],
        ...     ["protein_1_1", "protein_1_2", "protein_2_1"],
        ...     {"*": "*_<N>"}
        ... )
        {
            'protein_1': {'sources': ['protein_1'], 'targets': ['protein_1_1', 'protein_1_2']},
            'protein_2': {'sources': ['protein_2'], 'targets': ['protein_2_1']}
        }
    """
    if id_map is None:
        id_map = {"*": "*_<N>"}

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


if __name__ == "__main__":
    # Test cases for map_table_ids_to_ids - expecting lists of candidates
    test_cases = [
        ("rifampicin_1", {"*": "*_<N>"}, ["rifampicin_1", "rifampicin"]),
        ("rifampicin_1_2", {"*": "*_<N>"}, ["rifampicin_1_2", "rifampicin_1", "rifampicin"]),
        ("rifampicin_1_2_3", {"*": "*_<N>"}, ["rifampicin_1_2_3", "rifampicin_1_2", "rifampicin_1", "rifampicin"]),
        ("RFDAA_Hit_Screen_007_1_1", {"*": "*_<N>"}, ["RFDAA_Hit_Screen_007_1_1", "RFDAA_Hit_Screen_007_1", "RFDAA_Hit_Screen_007", "RFDAA_Hit_Screen"]),
        ("protein-seq-42", {"*": "*-seq-<N>"}, ["protein-seq-42", "protein"]),
        ("data_something_123", {"*": "*_something_<N>"}, ["data_something_123", "data"]),
        ("no_change", {"*": "*"}, ["no_change"]),
        ("no_match", {"*": "*_<N>"}, ["no_match"]),  # No numeric suffix
        ("complex_001_v2_3", {"*": "*_v2_<N>"}, ["complex_001_v2_3", "complex_001"]),
        # Edge cases
        ("test_1_2_foo", {"*": "*_<N>"}, ["test_1_2_foo"]),  # Non-numeric part at end
        ("test_foo_1", {"*": "*_<N>"}, ["test_foo_1", "test_foo"]),  # Non-numeric part in middle
    ]

    print("Running ID mapping tests...")
    all_passed = True

    for structure_id, id_map, expected in test_cases:
        result = map_table_ids_to_ids(structure_id, id_map)
        passed = result == expected
        all_passed = all_passed and passed

        status = "PASS" if passed else "FAIL"
        print(f"[{status}] {structure_id} with {id_map}")
        print(f"  Result:   {result}")
        print(f"  Expected: {expected}")

    if all_passed:
        print("\nAll map_table_ids_to_ids tests passed!")
    else:
        print("\nSome map_table_ids_to_ids tests failed!")
        exit(1)

    # Test cases for get_mapped_ids with unique=True (default)
    print("\n" + "="*60)
    print("Running get_mapped_ids tests (unique=True)...")

    get_mapped_unique_tests = [
        # (source_ids, target_ids, id_map, expected)
        # Case 1: Exact match - same IDs on both sides
        (
            ["protein_1", "protein_2"],
            ["protein_1", "protein_2"],
            {"*": "*_<N>"},
            {"protein_1": "protein_1", "protein_2": "protein_2"}
        ),
        # Case 2: Target is child of source (target = source + _N) - returns first/most specific
        (
            ["protein_1", "protein_2"],
            ["protein_1_1", "protein_1_2", "protein_2_1"],
            {"*": "*_<N>"},
            {"protein_1": "protein_1_1", "protein_2": "protein_2_1"}
        ),
        # Case 3: Source is child of target (source = target + _N)
        (
            ["protein_1_1", "protein_2_1"],
            ["protein_1", "protein_2"],
            {"*": "*_<N>"},
            {"protein_1_1": "protein_1", "protein_2_1": "protein_2"}
        ),
        # Case 4: Mixed - some exact, some parent-child
        (
            ["protein_1", "protein_2_1"],
            ["protein_1", "protein_2"],
            {"*": "*_<N>"},
            {"protein_1": "protein_1", "protein_2_1": "protein_2"}
        ),
        # Case 5: No matches - returns None
        (
            ["alpha_1", "alpha_2"],
            ["beta_1", "beta_2"],
            {"*": "*_<N>"},
            {"alpha_1": None, "alpha_2": None}
        ),
        # Case 6: Deep nesting - returns most specific parent match
        (
            ["protein_1_2_3"],
            ["protein_1", "protein_1_2"],
            {"*": "*_<N>"},
            {"protein_1_2_3": "protein_1_2"}
        ),
        # Case 7: Sibling match - both derive from same ancestor
        (
            ["protein_1_1"],
            ["protein_1_2"],
            {"*": "*_<N>"},
            {"protein_1_1": "protein_1_2"}
        ),
        # Case 8: Sibling match - prefer closer sibling
        (
            ["protein_1_1"],
            ["protein_1_2", "protein_2_1"],
            {"*": "*_<N>"},
            {"protein_1_1": "protein_1_2"}
        ),
        # Case 9: Deep sibling - common ancestor at different levels
        (
            ["protein_1_1_1"],
            ["protein_1_1_2", "protein_1_2"],
            {"*": "*_<N>"},
            {"protein_1_1_1": "protein_1_1_2"}  # Closer sibling (distance 1+1=2 vs 2+1=3)
        ),
    ]

    all_passed = True
    for source_ids, target_ids, id_map, expected in get_mapped_unique_tests:
        result = get_mapped_ids(source_ids, target_ids, id_map, unique=True)
        passed = result == expected
        all_passed = all_passed and passed

        status = "PASS" if passed else "FAIL"
        print(f"[{status}] sources={source_ids}, targets={target_ids}")
        if not passed:
            print(f"  Result:   {result}")
            print(f"  Expected: {expected}")

    if all_passed:
        print("\nAll get_mapped_ids (unique=True) tests passed!")
    else:
        print("\nSome get_mapped_ids (unique=True) tests failed!")
        exit(1)

    # Test cases for get_mapped_ids with unique=False
    print("\n" + "="*60)
    print("Running get_mapped_ids tests (unique=False)...")

    get_mapped_list_tests = [
        # (source_ids, target_ids, id_map, expected)
        # Case 1: Exact match - includes exact + children
        (
            ["protein_1"],
            ["protein_1", "protein_1_1", "protein_1_2"],
            {"*": "*_<N>"},
            {"protein_1": ["protein_1", "protein_1_1", "protein_1_2"]}
        ),
        # Case 2: Target is child of source (target = source + _N)
        (
            ["protein_1", "protein_2"],
            ["protein_1_1", "protein_1_2", "protein_2_1"],
            {"*": "*_<N>"},
            {"protein_1": ["protein_1_1", "protein_1_2"], "protein_2": ["protein_2_1"]}
        ),
        # Case 3: Source is child of target - returns all parent matches
        (
            ["protein_1_2_3"],
            ["protein_1", "protein_1_2"],
            {"*": "*_<N>"},
            {"protein_1_2_3": ["protein_1_2", "protein_1"]}
        ),
        # Case 4: No matches - returns empty list
        (
            ["alpha_1"],
            ["beta_1"],
            {"*": "*_<N>"},
            {"alpha_1": []}
        ),
        # Case 5: Sibling match - returns siblings sorted by distance
        (
            ["protein_1_1"],
            ["protein_1_2", "protein_1_3"],
            {"*": "*_<N>"},
            {"protein_1_1": ["protein_1_2", "protein_1_3"]}
        ),
    ]

    all_passed_list = True
    for source_ids, target_ids, id_map, expected in get_mapped_list_tests:
        result = get_mapped_ids(source_ids, target_ids, id_map, unique=False)
        passed = result == expected
        all_passed_list = all_passed_list and passed

        status = "PASS" if passed else "FAIL"
        print(f"[{status}] sources={source_ids}, targets={target_ids}")
        if not passed:
            print(f"  Result:   {result}")
            print(f"  Expected: {expected}")

    if all_passed_list:
        print("\nAll get_mapped_ids (unique=False) tests passed!")
    else:
        print("\nSome get_mapped_ids (unique=False) tests failed!")
        exit(1)

    print("\n" + "="*60)
    print("All tests passed!")
