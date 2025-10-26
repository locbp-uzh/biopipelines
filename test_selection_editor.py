"""
Test script for SelectionEditor tool.

Tests the selection modification logic with various operations.
"""

import sys
import os

# Add paths for imports
sys.path.insert(0, os.path.join(os.getcwd(), 'HelpScripts'))
sys.path.insert(0, os.path.join(os.getcwd(), 'PipelineScripts'))

from pipe_selection_editor import (
    parse_pymol_selection,
    format_pymol_selection,
    merge_overlapping_ranges,
    expand_selection,
    shrink_selection,
    shift_selection,
    invert_selection
)


def test_parsing():
    """Test PyMOL selection parsing."""
    print("Testing PyMOL selection parsing...")

    # Test simple range
    assert parse_pymol_selection("3-45") == [(3, 45)]
    print("  [OK] Simple range: '3-45'")

    # Test multiple ranges
    assert parse_pymol_selection("3-45+58-60") == [(3, 45), (58, 60)]
    print("  [OK] Multiple ranges: '3-45+58-60'")

    # Test single residue
    assert parse_pymol_selection("10") == [(10, 10)]
    print("  [OK] Single residue: '10'")

    # Test mixed
    assert parse_pymol_selection("1-4+6-10+15") == [(1, 4), (6, 10), (15, 15)]
    print("  [OK] Mixed format: '1-4+6-10+15'")

    # Test empty
    assert parse_pymol_selection("") == []
    print("  [OK] Empty selection: ''")

    print("[OK] All parsing tests passed!\n")


def test_formatting():
    """Test PyMOL selection formatting."""
    print("Testing PyMOL selection formatting...")

    # Test simple range
    assert format_pymol_selection([(3, 45)]) == "3-45"
    print("  [OK] Simple range")

    # Test multiple ranges
    assert format_pymol_selection([(3, 45), (58, 60)]) == "3-45+58-60"
    print("  [OK] Multiple ranges")

    # Test single residue
    assert format_pymol_selection([(10, 10)]) == "10"
    print("  [OK] Single residue")

    # Test empty
    assert format_pymol_selection([]) == ""
    print("  [OK] Empty selection")

    print("[OK] All formatting tests passed!\n")


def test_merging():
    """Test range merging."""
    print("Testing range merging...")

    # Test overlapping ranges
    ranges = [(1, 5), (4, 10)]
    merged = merge_overlapping_ranges(ranges)
    assert merged == [(1, 10)]
    print("  [OK] Overlapping ranges: (1,5) + (4,10) = (1,10)")

    # Test adjacent ranges
    ranges = [(1, 5), (6, 10)]
    merged = merge_overlapping_ranges(ranges)
    assert merged == [(1, 10)]
    print("  [OK] Adjacent ranges: (1,5) + (6,10) = (1,10)")

    # Test separate ranges
    ranges = [(1, 5), (10, 15)]
    merged = merge_overlapping_ranges(ranges)
    assert merged == [(1, 5), (10, 15)]
    print("  [OK] Separate ranges: (1,5) + (10,15) = unchanged")

    # Test complex merging
    ranges = [(1, 4), (6, 10), (8, 15), (20, 25)]
    merged = merge_overlapping_ranges(ranges)
    assert merged == [(1, 4), (6, 15), (20, 25)]
    print("  [OK] Complex merging: (1,4) + (6,10) + (8,15) + (20,25) = (1,4) + (6,15) + (20,25)")

    print("[OK] All merging tests passed!\n")


def test_expand():
    """Test expansion operation."""
    print("Testing expansion operation...")

    # Create mock valid residues (1-100)
    valid_residues = set(range(1, 101))

    # Test simple expansion
    ranges = [(10, 20)]
    expanded = expand_selection(ranges, 2, valid_residues)
    assert expanded == [(8, 22)]
    print("  [OK] Simple expansion: (10,20) + expand=2 = (8,22)")

    # Test expansion with merging
    ranges = [(10, 15), (17, 20)]
    expanded = expand_selection(ranges, 1, valid_residues)
    assert expanded == [(9, 21)]
    print("  [OK] Expansion with merging: (10,15) + (17,20) + expand=1 = (9,21)")

    # Test expansion at boundary
    ranges = [(3, 10)]
    expanded = expand_selection(ranges, 5, valid_residues)
    assert expanded == [(1, 15)]  # Can't go below 1
    print("  [OK] Expansion at boundary: (3,10) + expand=5 = (1,15)")

    print("[OK] All expansion tests passed!\n")


def test_shrink():
    """Test shrinking operation."""
    print("Testing shrinking operation...")

    # Create mock valid residues (1-100)
    valid_residues = set(range(1, 101))

    # Test simple shrinking
    ranges = [(10, 20)]
    shrunk = shrink_selection(ranges, 2, valid_residues)
    assert shrunk == [(12, 18)]
    print("  [OK] Simple shrinking: (10,20) + shrink=2 = (12,18)")

    # Test shrinking that removes range
    ranges = [(10, 12)]
    shrunk = shrink_selection(ranges, 2, valid_residues)
    assert shrunk == []
    print("  [OK] Shrinking removes range: (10,12) + shrink=2 = []")

    # Test shrinking multiple ranges
    ranges = [(10, 20), (30, 40)]
    shrunk = shrink_selection(ranges, 1, valid_residues)
    assert shrunk == [(11, 19), (31, 39)]
    print("  [OK] Shrinking multiple: (10,20) + (30,40) + shrink=1 = (11,19) + (31,39)")

    print("[OK] All shrinking tests passed!\n")


def test_shift():
    """Test shifting operation."""
    print("Testing shifting operation...")

    # Create mock valid residues (1-100)
    valid_residues = set(range(1, 101))

    # Test positive shift
    ranges = [(10, 20)]
    shifted = shift_selection(ranges, 5, valid_residues)
    assert shifted == [(15, 25)]
    print("  [OK] Positive shift: (10,20) + shift=5 = (15,25)")

    # Test negative shift
    ranges = [(20, 30)]
    shifted = shift_selection(ranges, -5, valid_residues)
    assert shifted == [(15, 25)]
    print("  [OK] Negative shift: (20,30) + shift=-5 = (15,25)")

    # Test multiple ranges
    ranges = [(10, 15), (20, 25)]
    shifted = shift_selection(ranges, 3, valid_residues)
    assert shifted == [(13, 18), (23, 28)]
    print("  [OK] Shifting multiple: (10,15) + (20,25) + shift=3 = (13,18) + (23,28)")

    print("[OK] All shifting tests passed!\n")


def test_invert():
    """Test inversion operation."""
    print("Testing inversion operation...")

    # Create small set of valid residues for easier testing
    valid_residues = set(range(1, 21))  # 1-20

    # Test simple inversion
    ranges = [(5, 10)]
    inverted = invert_selection(ranges, valid_residues)
    assert inverted == [(1, 4), (11, 20)]
    print("  [OK] Simple inversion: (5,10) inverted = (1,4) + (11,20)")

    # Test inversion of multiple ranges
    ranges = [(5, 7), (10, 12)]
    inverted = invert_selection(ranges, valid_residues)
    assert inverted == [(1, 4), (8, 9), (13, 20)]
    print("  [OK] Multiple range inversion: (5,7) + (10,12) inverted = (1,4) + (8,9) + (13,20)")

    # Test inversion of everything (should give empty)
    ranges = [(1, 20)]
    inverted = invert_selection(ranges, valid_residues)
    assert inverted == []
    print("  [OK] Full inversion: (1,20) inverted = []")

    print("[OK] All inversion tests passed!\n")


def test_combined_operations():
    """Test combined operations."""
    print("Testing combined operations...")

    valid_residues = set(range(1, 101))

    # Expand then merge example from requirements
    # "1-4+6-10" with expand=1 should give "1-11"
    ranges = [(1, 4), (6, 10)]
    expanded = expand_selection(ranges, 1, valid_residues)
    result = format_pymol_selection(expanded)
    assert result == "1-11"
    print(f"  [OK] Example from requirements: '1-4+6-10' + expand=1 = '{result}'")

    # Another example: "3-45" with expand=2 where residues 1,2 don't exist
    valid_residues_with_gap = set(range(3, 101))
    ranges = [(3, 45)]
    expanded = expand_selection(ranges, 2, valid_residues_with_gap)
    result = format_pymol_selection(expanded)
    assert result == "3-47"
    print(f"  [OK] Expansion with missing residues: '3-45' + expand=2 (no 1,2) = '{result}'")

    print("[OK] All combined operation tests passed!\n")


def run_all_tests():
    """Run all tests."""
    print("="*60)
    print("SelectionEditor Test Suite")
    print("="*60 + "\n")

    try:
        test_parsing()
        test_formatting()
        test_merging()
        test_expand()
        test_shrink()
        test_shift()
        test_invert()
        test_combined_operations()

        print("="*60)
        print("[OK] ALL TESTS PASSED!")
        print("="*60)
        return True

    except AssertionError as e:
        print(f"\n[FAIL] TEST FAILED: {e}")
        import traceback
        traceback.print_exc()
        return False
    except Exception as e:
        print(f"\n[FAIL] ERROR: {e}")
        import traceback
        traceback.print_exc()
        return False


if __name__ == "__main__":
    success = run_all_tests()
    sys.exit(0 if success else 1)

