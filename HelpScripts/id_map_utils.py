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
from typing import Dict, Optional


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


if __name__ == "__main__":
    # Test cases - now expecting lists of candidates
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
        print("\nAll tests passed!")
    else:
        print("\nSome tests failed!")
        exit(1)
