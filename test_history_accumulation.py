#!/usr/bin/env python3
"""
Test script for SequenceMetricAnalysis history accumulation.

Tests the weighted aggregation logic with synthetic data.
"""

import pandas as pd
import numpy as np
import sys
import os

# Add HelpScripts to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'HelpScripts'))

from pipe_sequence_metric_analysis import accumulate_with_history

class MockLogger:
    def info(self, msg):
        print(f"INFO: {msg}")

def test_history_accumulation():
    """Test that history accumulation works correctly."""

    print("Testing history accumulation...")

    # Create current cycle data (10 observations of mutation A10G with affinity=5.0)
    current_df = pd.DataFrame([{
        'position': 10,
        'wt_aa': 'A',
        'mut_aa': 'G',
        'count': 10,
        'affinity_mean': 5.0,
        'affinity_std': 0.5,
        'affinity_min': 4.0,
        'affinity_max': 6.0
    }])

    # Create historical data (20 observations of same mutation with affinity=7.0)
    history_df = pd.DataFrame([{
        'position': 10,
        'wt_aa': 'A',
        'mut_aa': 'G',
        'count': 20,
        'affinity_mean': 7.0,
        'affinity_std': 1.0,
        'affinity_min': 5.0,
        'affinity_max': 9.0
    }])

    # Expected results:
    # Total count: 30
    # Weighted mean: (5.0*10 + 7.0*20) / 30 = 190/30 = 6.333...
    # Min: min(4.0, 5.0) = 4.0
    # Max: max(6.0, 9.0) = 9.0

    logger = MockLogger()
    result = accumulate_with_history(current_df, history_df, ['affinity'], logger)

    print("\nResults:")
    print(result)

    # Verify results
    assert len(result) == 1, f"Expected 1 row, got {len(result)}"

    row = result.iloc[0]
    assert row['count'] == 30, f"Expected count=30, got {row['count']}"

    expected_mean = (5.0 * 10 + 7.0 * 20) / 30
    assert abs(row['affinity_mean'] - expected_mean) < 0.001, \
        f"Expected mean={expected_mean:.3f}, got {row['affinity_mean']:.3f}"

    assert row['affinity_min'] == 4.0, f"Expected min=4.0, got {row['affinity_min']}"
    assert row['affinity_max'] == 9.0, f"Expected max=9.0, got {row['affinity_max']}"

    print(f"\n[OK] Count: {row['count']} (expected 30)")
    print(f"[OK] Mean: {row['affinity_mean']:.3f} (expected {expected_mean:.3f})")
    print(f"[OK] Std: {row['affinity_std']:.3f}")
    print(f"[OK] Min: {row['affinity_min']} (expected 4.0)")
    print(f"[OK] Max: {row['affinity_max']} (expected 9.0)")

    print("\n[OK] Test PASSED: History accumulation works correctly!")
    return True


def test_new_mutations():
    """Test that new mutations in current cycle are preserved."""

    print("\n" + "="*60)
    print("Testing new mutations (only in current cycle)...")

    # Current cycle has A10G
    current_df = pd.DataFrame([{
        'position': 10,
        'wt_aa': 'A',
        'mut_aa': 'G',
        'count': 5,
        'affinity_mean': 3.0,
        'affinity_std': 0.3,
        'affinity_min': 2.5,
        'affinity_max': 3.5
    }])

    # History has A10D (different mutation)
    history_df = pd.DataFrame([{
        'position': 10,
        'wt_aa': 'A',
        'mut_aa': 'D',
        'count': 10,
        'affinity_mean': 5.0,
        'affinity_std': 0.5,
        'affinity_min': 4.0,
        'affinity_max': 6.0
    }])

    logger = MockLogger()
    result = accumulate_with_history(current_df, history_df, ['affinity'], logger)

    print("\nResults:")
    print(result)

    # Should have 2 mutations
    assert len(result) == 2, f"Expected 2 rows, got {len(result)}"

    # Check A10G (only in current)
    aog = result[(result['mut_aa'] == 'G')]
    assert len(aog) == 1
    assert aog.iloc[0]['count'] == 5
    assert aog.iloc[0]['affinity_mean'] == 3.0

    # Check A10D (only in history)
    aod = result[(result['mut_aa'] == 'D')]
    assert len(aod) == 1
    assert aod.iloc[0]['count'] == 10
    assert aod.iloc[0]['affinity_mean'] == 5.0

    print("\n[OK] Test PASSED: New mutations preserved correctly!")
    return True


if __name__ == '__main__':
    try:
        test_history_accumulation()
        test_new_mutations()
        print("\n" + "="*60)
        print("All tests PASSED!")
        print("="*60)
    except AssertionError as e:
        print(f"\n[FAIL] Test FAILED: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"\n[ERROR] Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
