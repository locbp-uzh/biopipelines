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


def map_table_ids_to_ids(structure_id: str, id_map: Dict[str, str]) -> str:
    """
    Map structure ID to table ID using id_map pattern, applying recursively.

    This function extracts the base ID from a structure ID by repeatedly removing
    numeric suffixes specified in the id_map pattern. This allows a single pattern
    like {"*": "*_<N>"} to handle IDs that have been modified multiple times
    (e.g., rifampicin_1_2 → rifampicin).

    Args:
        structure_id: Structure ID (e.g., "rifampicin_1_2", "protein-seq-42")
        id_map: ID mapping dictionary (e.g., {"*": "*_<N>"})

    Returns:
        Mapped table ID (e.g., "rifampicin", "protein")
        If pattern doesn't match, returns structure_id unchanged

    Examples:
        >>> map_structure_id_to_table_id("rifampicin_1", {"*": "*_<N>"})
        'rifampicin'

        >>> map_structure_id_to_table_id("rifampicin_1_2", {"*": "*_<N>"})
        'rifampicin'

        >>> map_structure_id_to_table_id("rifampicin_1_2_3", {"*": "*_<N>"})
        'rifampicin'

        >>> map_structure_id_to_table_id("protein-seq-42", {"*": "*-seq-<N>"})
        'protein'

        >>> map_structure_id_to_table_id("no_change", {"*": "*"})
        'no_change'

        >>> map_structure_id_to_table_id("no_match", {"*": "*_<N>"})
        'no_match'
    """
    if not id_map or "*" not in id_map:
        return structure_id

    pattern_str = id_map["*"]

    # Check if pattern is just "*" (identity mapping)
    if pattern_str == "*":
        return structure_id

    # Verify pattern starts with "*"
    if not pattern_str.startswith("*"):
        return structure_id

    # Extract the suffix pattern (everything after the first "*")
    suffix_pattern = pattern_str[1:]  # Remove leading "*"

    if not suffix_pattern:
        return structure_id

    # Parse the pattern to extract delimiter and literal parts
    # e.g., "_<N>" → delimiter="_", no literal
    # e.g., "-seq-<N>" → delimiter="-", literal="seq"

    # Find the first occurrence of <N>
    n_pos = suffix_pattern.find('<N>')
    if n_pos == -1:
        # No <N> in pattern, use regex method
        pattern = parse_id_map_pattern(id_map)
        if pattern is None:
            return structure_id
        match = pattern.match(structure_id)
        if match:
            return match.group(1)
        return structure_id

    # Extract delimiter (character before <N> or at start)
    if n_pos == 0:
        return structure_id  # Pattern like "*<N>" doesn't make sense

    delimiter = suffix_pattern[n_pos - 1]
    # Extract literal part, excluding the leading delimiter if present
    literal_with_delim = suffix_pattern[:n_pos - 1] if n_pos > 1 else ""
    # Remove leading delimiter from literal if it exists
    if literal_with_delim.startswith(delimiter):
        literal_part = literal_with_delim[1:]
    else:
        literal_part = literal_with_delim

    # Split the ID by delimiter
    parts = structure_id.split(delimiter)

    if len(parts) <= 1:
        # No delimiter found, return as-is
        return structure_id

    # Remove trailing parts that match the pattern
    # Work from the end backwards to strip matching suffixes

    if literal_part:
        # Pattern like "*-seq-<N>" or "*_something_<N>"
        # We need to find and remove suffix that matches: delimiter + literal + delimiter + digits
        # e.g., "-seq-42" from "protein-seq-42"

        # Reconstruct what we're looking for at the end
        # The pattern after "*" is like "-seq-<N>"
        # So we're looking for parts ending with: literal part followed by digits

        # Work backwards: remove trailing numeric parts that follow the literal
        while len(parts) >= 2:
            last_part = parts[-1]
            second_last = parts[-2] if len(parts) >= 2 else ""

            # Check if last part is a number and second-to-last matches literal
            if last_part.isdigit() and second_last == literal_part:
                # Remove both parts
                parts = parts[:-2]
            elif last_part == literal_part and len(parts) >= 2:
                # Literal is at the end without a number, shouldn't match
                break
            else:
                break

    else:
        # Pattern like "*_<N>": strip all trailing numeric parts
        while len(parts) > 0 and parts[-1].isdigit():
            parts = parts[:-1]

    # Rejoin with delimiter
    if not parts:
        return structure_id

    return delimiter.join(parts)


if __name__ == "__main__":
    # Test cases
    test_cases = [
        ("rifampicin_1", {"*": "*_<N>"}, "rifampicin"),
        ("rifampicin_1_2", {"*": "*_<N>"}, "rifampicin"),  # Now works with single pattern
        ("rifampicin_1_2_3", {"*": "*_<N>"}, "rifampicin"),  # Multiple suffixes
        ("rifampicin_1_2", {"*": "*_<N>_<N>"}, "rifampicin"),  # Old style still works
        ("protein-seq-42", {"*": "*-seq-<N>"}, "protein"),
        ("data_something_123", {"*": "*_something_<N>"}, "data"),
        ("no_change", {"*": "*"}, "no_change"),
        ("no_match", {"*": "*_<N>"}, "no_match"),
        ("complex_001_v2_3", {"*": "*_v2_<N>"}, "complex_001"),
        # Edge cases
        ("test_1_2_foo", {"*": "*_<N>"}, "test_1_2_foo"),  # Non-numeric part at end
        ("test_foo_1", {"*": "*_<N>"}, "test_foo"),  # Non-numeric part in middle
    ]

    print("Running ID mapping tests...")
    all_passed = True

    for structure_id, id_map, expected in test_cases:
        result = map_table_ids_to_ids(structure_id, id_map)
        passed = result == expected
        all_passed = all_passed and passed

        status = "PASS" if passed else "FAIL"
        print(f"[{status}] {structure_id} with {id_map} -> {result} (expected: {expected})")

    if all_passed:
        print("\nAll tests passed!")
    else:
        print("\nSome tests failed!")
        exit(1)
