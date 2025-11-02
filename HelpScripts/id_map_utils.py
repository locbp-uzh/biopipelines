#!/usr/bin/env python3
"""
Utility functions for ID mapping pattern matching.

This module provides a unified approach to parsing and applying ID mapping patterns
across different tools. The pattern format is {"*": "<pattern>"} where <pattern>
describes how to map structure IDs to datasheet IDs.

Pattern syntax:
- "*" in both key and value represents the base ID
- "<N>" represents a numeric suffix (one or more digits)
- Everything else is treated as literal text

Examples:
- {"*": "*"}           → no mapping (identity)
- {"*": "*_<N>"}       → strip "_123" suffix
- {"*": "*-<N>"}       → strip "-123" suffix
- {"*": "*_<N>_<N>"}   → strip "_123_456" suffix
- {"*": "*_seq_<N>"}   → strip "_seq_123" suffix
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


def map_structure_id_to_datasheet_id(structure_id: str, id_map: Dict[str, str]) -> str:
    """
    Map structure ID to datasheet ID using id_map pattern.

    This function extracts the base ID from a structure ID by removing the suffix
    specified in the id_map pattern.

    Args:
        structure_id: Structure ID (e.g., "rifampicin_1_2", "protein-seq-42")
        id_map: ID mapping dictionary (e.g., {"*": "*_<N>_<N>"})

    Returns:
        Mapped datasheet ID (e.g., "rifampicin", "protein")
        If pattern doesn't match, returns structure_id unchanged

    Examples:
        >>> map_structure_id_to_datasheet_id("rifampicin_1", {"*": "*_<N>"})
        'rifampicin'

        >>> map_structure_id_to_datasheet_id("rifampicin_1_2", {"*": "*_<N>_<N>"})
        'rifampicin'

        >>> map_structure_id_to_datasheet_id("protein-seq-42", {"*": "*-seq-<N>"})
        'protein'

        >>> map_structure_id_to_datasheet_id("no_change", {"*": "*"})
        'no_change'

        >>> map_structure_id_to_datasheet_id("no_match", {"*": "*_<N>"})
        'no_match'
    """
    # Parse the pattern
    pattern = parse_id_map_pattern(id_map)

    # If no pattern (identity mapping), return as-is
    if pattern is None:
        return structure_id

    # Try to match and extract base ID
    match = pattern.match(structure_id)
    if match:
        return match.group(1)

    # No match - return structure_id unchanged
    return structure_id


if __name__ == "__main__":
    # Test cases
    test_cases = [
        ("rifampicin_1", {"*": "*_<N>"}, "rifampicin"),
        ("rifampicin_1_2", {"*": "*_<N>_<N>"}, "rifampicin"),
        ("rifampicin_1_2", {"*": "*_<N>"}, "rifampicin_1"),
        ("protein-seq-42", {"*": "*-seq-<N>"}, "protein"),
        ("data_something_123", {"*": "*_something_<N>"}, "data"),
        ("no_change", {"*": "*"}, "no_change"),
        ("no_match", {"*": "*_<N>"}, "no_match"),
        ("complex_001_v2_3", {"*": "*_v2_<N>"}, "complex_001"),
    ]

    print("Running ID mapping tests...")
    all_passed = True

    for structure_id, id_map, expected in test_cases:
        result = map_structure_id_to_datasheet_id(structure_id, id_map)
        passed = result == expected
        all_passed = all_passed and passed

        status = "PASS" if passed else "FAIL"
        print(f"[{status}] {structure_id} with {id_map} -> {result} (expected: {expected})")

    if all_passed:
        print("\nAll tests passed!")
    else:
        print("\nSome tests failed!")
        exit(1)
