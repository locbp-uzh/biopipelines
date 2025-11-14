#!/usr/bin/env python3
"""
Sequence-metric analysis runtime script for BioPipelines.

Analyzes correlations between sequence mutations and metrics, calculates
statistical measures, and generates scored mutation tables for data-driven
sequence optimization.
"""

import argparse
import json
import pandas as pd
import numpy as np
from pathlib import Path
import logging
from typing import Dict, List, Tuple, Optional
import sys

# Standard amino acids
AMINO_ACIDS = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']


def setup_logging():
    """Set up logging configuration."""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    return logging.getLogger(__name__)


def load_reference_sequence(ref_input: str) -> str:
    """
    Load reference sequence from string or file.

    Args:
        ref_input: Either direct sequence or @path/to/sequences.csv

    Returns:
        Reference sequence string
    """
    if ref_input.startswith('@'):
        # It's a file path
        csv_path = ref_input[1:]  # Remove @ prefix
        df = pd.read_csv(csv_path)

        # Find sequence column
        seq_col = None
        for col in df.columns:
            if col.lower() in ['sequence', 'seq']:
                seq_col = col
                break

        if seq_col is None:
            raise ValueError(f"Could not find sequence column in {csv_path}")

        # Return first sequence
        return df[seq_col].iloc[0]
    else:
        # Direct sequence string
        return ref_input


def call_mutations(sequence: str, reference: str) -> List[Tuple[int, str, str]]:
    """
    Call mutations between sequence and reference.

    Args:
        sequence: Query sequence
        reference: Reference sequence

    Returns:
        List of (position, ref_aa, mut_aa) tuples (1-indexed positions)
    """
    mutations = []

    min_len = min(len(sequence), len(reference))

    for i in range(min_len):
        ref_aa = reference[i]
        mut_aa = sequence[i]

        if ref_aa != mut_aa:
            mutations.append((i + 1, ref_aa, mut_aa))  # 1-indexed

    return mutations


def accumulate_with_history(current_df: pd.DataFrame,
                           history_df: pd.DataFrame,
                           metric_columns: List[str],
                           logger) -> pd.DataFrame:
    """
    Accumulate current statistics with historical statistics using weighted aggregation.

    Uses proper formulas for combining means and variances from independent samples.

    Args:
        current_df: Current cycle aggregated statistics
        history_df: Historical aggregated statistics (from previous cycles)
        metric_columns: List of metric column names
        logger: Logger instance

    Returns:
        Combined DataFrame with accumulated statistics
    """
    # Merge current and history on mutation key
    merge_keys = ['position', 'wt_aa', 'mut_aa']

    # Separate mutations: only in current, only in history, in both
    merged = pd.merge(current_df, history_df, on=merge_keys, how='outer', suffixes=('_curr', '_hist'))

    accumulated_rows = []

    for _, row in merged.iterrows():
        position = row['position']
        wt_aa = row['wt_aa']
        mut_aa = row['mut_aa']

        # Get counts (use 0 if NaN)
        n_curr = row.get('count_curr', 0)
        n_hist = row.get('count_hist', 0)
        n_curr = 0 if pd.isna(n_curr) else n_curr
        n_hist = 0 if pd.isna(n_hist) else n_hist

        # Total count
        n_total = n_curr + n_hist

        accumulated_row = {
            'position': position,
            'wt_aa': wt_aa,
            'mut_aa': mut_aa,
            'count': n_total
        }

        # For each metric, combine statistics
        for metric in metric_columns:
            mean_curr_col = f"{metric}_mean_curr"
            mean_hist_col = f"{metric}_mean_hist"
            std_curr_col = f"{metric}_std_curr"
            std_hist_col = f"{metric}_std_hist"
            min_curr_col = f"{metric}_min_curr"
            min_hist_col = f"{metric}_min_hist"
            max_curr_col = f"{metric}_max_curr"
            max_hist_col = f"{metric}_max_hist"

            # Get values (handle NaN)
            mean_curr = row.get(mean_curr_col, np.nan)
            mean_hist = row.get(mean_hist_col, np.nan)
            std_curr = row.get(std_curr_col, np.nan)
            std_hist = row.get(std_hist_col, np.nan)
            min_curr = row.get(min_curr_col, np.nan)
            min_hist = row.get(min_hist_col, np.nan)
            max_curr = row.get(max_curr_col, np.nan)
            max_hist = row.get(max_hist_col, np.nan)

            # Case 1: Only in current cycle
            if n_hist == 0 or pd.isna(mean_hist):
                accumulated_row[f"{metric}_mean"] = mean_curr
                accumulated_row[f"{metric}_std"] = std_curr if not pd.isna(std_curr) else 0.0
                accumulated_row[f"{metric}_min"] = min_curr
                accumulated_row[f"{metric}_max"] = max_curr

            # Case 2: Only in history
            elif n_curr == 0 or pd.isna(mean_curr):
                accumulated_row[f"{metric}_mean"] = mean_hist
                accumulated_row[f"{metric}_std"] = std_hist if not pd.isna(std_hist) else 0.0
                accumulated_row[f"{metric}_min"] = min_hist
                accumulated_row[f"{metric}_max"] = max_hist

            # Case 3: In both - use weighted aggregation
            else:
                # Weighted mean
                combined_mean = (mean_curr * n_curr + mean_hist * n_hist) / n_total

                # Combined variance using parallel variance formula
                # Var = E[X²] - E[X]²
                # E[X²] = Var + Mean²
                var_curr = std_curr ** 2 if not pd.isna(std_curr) else 0.0
                var_hist = std_hist ** 2 if not pd.isna(std_hist) else 0.0

                # E[X²] for each group
                ex2_curr = var_curr + mean_curr ** 2
                ex2_hist = var_hist + mean_hist ** 2

                # Combined E[X²]
                combined_ex2 = (ex2_curr * n_curr + ex2_hist * n_hist) / n_total

                # Combined variance
                combined_var = combined_ex2 - combined_mean ** 2
                combined_std = np.sqrt(max(0, combined_var))  # Ensure non-negative

                # Combined min/max
                combined_min = np.nanmin([min_curr, min_hist])
                combined_max = np.nanmax([max_curr, max_hist])

                accumulated_row[f"{metric}_mean"] = combined_mean
                accumulated_row[f"{metric}_std"] = combined_std
                accumulated_row[f"{metric}_min"] = combined_min
                accumulated_row[f"{metric}_max"] = combined_max

        accumulated_rows.append(accumulated_row)

    result_df = pd.DataFrame(accumulated_rows)

    logger.info(f"Accumulated {len(result_df)} mutations (current: {len(current_df)}, historical: {len(history_df)})")

    return result_df


def analyze_sequence_metrics(sequences_df: pd.DataFrame,
                            metrics_df: pd.DataFrame,
                            reference_seq: str,
                            metric_columns: List[str],
                            history_df: Optional[pd.DataFrame] = None) -> pd.DataFrame:
    """
    Analyze correlations between mutations and metrics.

    Args:
        sequences_df: DataFrame with id and sequence columns
        metrics_df: DataFrame with id and metric columns
        reference_seq: Reference sequence string
        metric_columns: List of metric column names to analyze
        history_df: Previous mutation_statistics to accumulate with

    Returns:
        DataFrame with mutation statistics
    """
    logger = logging.getLogger(__name__)

    # Merge sequences with metrics on id
    merged_df = pd.merge(sequences_df, metrics_df, on='id', how='inner')

    logger.info(f"Analyzing {len(merged_df)} sequences with metrics")

    # Find sequence column
    seq_col = None
    for col in sequences_df.columns:
        if col.lower() in ['sequence', 'seq']:
            seq_col = col
            break

    if seq_col is None:
        raise ValueError("Could not find sequence column in sequences dataframe")

    # Initialize statistics storage
    mutation_stats = []

    # Process each sequence
    for _, row in merged_df.iterrows():
        sequence = row[seq_col]
        mutations = call_mutations(sequence, reference_seq)

        # Add each mutation with its metrics
        for pos, ref_aa, mut_aa in mutations:
            mutation_record = {
                'position': pos,
                'wt_aa': ref_aa,
                'mut_aa': mut_aa
            }

            # Add metric values
            for metric in metric_columns:
                if metric in row:
                    mutation_record[metric] = row[metric]

            mutation_stats.append(mutation_record)

    # Convert to DataFrame
    stats_df = pd.DataFrame(mutation_stats)

    if len(stats_df) == 0:
        logger.warning("No mutations found in sequences")
        return pd.DataFrame()

    # Group by position and mutation, calculate statistics for current cycle
    group_cols = ['position', 'wt_aa', 'mut_aa']

    aggregation_dict = {'mut_aa': 'count'}  # Count observations

    # Add aggregations for each metric
    for metric in metric_columns:
        if metric in stats_df.columns:
            aggregation_dict[metric] = ['mean', 'std', 'min', 'max']

    # Perform aggregation on current cycle
    grouped = stats_df.groupby(group_cols).agg(aggregation_dict)

    # Flatten multi-level columns
    grouped.columns = ['_'.join(col).strip('_') if col[1] else col[0] for col in grouped.columns.values]

    # Rename count column
    grouped = grouped.rename(columns={'mut_aa_count': 'count'})

    # Reset index to get position, wt_aa, mut_aa as columns
    grouped = grouped.reset_index()

    # Fill NaN std values with 0 (happens when count=1)
    for metric in metric_columns:
        std_col = f"{metric}_std"
        if std_col in grouped.columns:
            grouped[std_col] = grouped[std_col].fillna(0.0)

    # Accumulate with history using weighted aggregation
    if history_df is not None and len(history_df) > 0:
        logger.info(f"Accumulating with {len(history_df)} historical mutation records")
        grouped = accumulate_with_history(grouped, history_df, metric_columns, logger)

    logger.info(f"Found {len(grouped)} unique mutations across {grouped['position'].nunique()} positions")

    return grouped


def calculate_delta_scores(stats_df: pd.DataFrame,
                           primary_metric: str,
                           mode: str,
                           min_observations: int) -> pd.DataFrame:
    """
    Calculate delta scores (improvement vs wild-type).

    Args:
        stats_df: Mutation statistics DataFrame
        primary_metric: Metric to use for scoring
        mode: "minimize" or "maximize"
        min_observations: Minimum count to assign non-zero scores

    Returns:
        DataFrame in mutation_deltas format (position, wt_aa, A, C, D, ...)
    """
    logger = logging.getLogger(__name__)

    # Calculate wild-type mean for each position
    wt_means = {}
    for pos in stats_df['position'].unique():
        pos_data = stats_df[stats_df['position'] == pos]
        wt_aa = pos_data['wt_aa'].iloc[0]

        # Check if wild-type was observed
        wt_obs = pos_data[pos_data['mut_aa'] == wt_aa]
        if len(wt_obs) > 0:
            wt_means[pos] = wt_obs[f"{primary_metric}_mean"].iloc[0]
        else:
            # No wild-type observations - use overall mean at this position
            wt_means[pos] = pos_data[f"{primary_metric}_mean"].mean()

    # Build delta scores table
    delta_rows = []

    for pos in sorted(stats_df['position'].unique()):
        pos_data = stats_df[stats_df['position'] == pos]
        wt_aa = pos_data['wt_aa'].iloc[0]
        wt_mean = wt_means[pos]

        delta_row = {'position': pos, 'wt_aa': wt_aa}

        for aa in AMINO_ACIDS:
            # Find mutation data
            mut_data = pos_data[pos_data['mut_aa'] == aa]

            if len(mut_data) > 0:
                count = mut_data['count'].iloc[0]
                mut_mean = mut_data[f"{primary_metric}_mean"].iloc[0]

                # Calculate delta based on mode
                if mode == "minimize":
                    # Lower is better: delta = wt_mean - mut_mean
                    # Positive delta = improvement
                    raw_delta = wt_mean - mut_mean
                else:  # maximize
                    # Higher is better: delta = mut_mean - wt_mean
                    # Positive delta = improvement
                    raw_delta = mut_mean - wt_mean

                # Apply minimum observations threshold
                if count >= min_observations:
                    delta_row[aa] = raw_delta
                else:
                    # Not enough observations - neutral score
                    delta_row[aa] = 0.0
            else:
                # Mutation not observed - neutral score
                delta_row[aa] = 0.0

        delta_rows.append(delta_row)

    delta_df = pd.DataFrame(delta_rows)

    logger.info(f"Calculated delta scores for {len(delta_df)} positions")

    return delta_df


def calculate_zscore_scores(stats_df: pd.DataFrame,
                            primary_metric: str,
                            mode: str,
                            min_observations: int) -> pd.DataFrame:
    """
    Calculate z-score scores (standardized values).

    Args:
        stats_df: Mutation statistics DataFrame
        primary_metric: Metric to use for scoring
        mode: "minimize" or "maximize"
        min_observations: Minimum count to assign non-zero scores

    Returns:
        DataFrame in mutation_zscores format (position, wt_aa, A, C, D, ...)
    """
    logger = logging.getLogger(__name__)

    # Calculate global mean and std for the metric
    global_mean = stats_df[f"{primary_metric}_mean"].mean()
    global_std = stats_df[f"{primary_metric}_mean"].std()

    if global_std == 0:
        logger.warning("Global std is 0 - all mutations have same metric value")
        global_std = 1.0  # Avoid division by zero

    # Build z-score table
    zscore_rows = []

    for pos in sorted(stats_df['position'].unique()):
        pos_data = stats_df[stats_df['position'] == pos]
        wt_aa = pos_data['wt_aa'].iloc[0]

        zscore_row = {'position': pos, 'wt_aa': wt_aa}

        for aa in AMINO_ACIDS:
            # Find mutation data
            mut_data = pos_data[pos_data['mut_aa'] == aa]

            if len(mut_data) > 0:
                count = mut_data['count'].iloc[0]
                mut_mean = mut_data[f"{primary_metric}_mean"].iloc[0]

                # Calculate z-score
                zscore = (mut_mean - global_mean) / global_std

                # Adjust sign based on mode
                if mode == "minimize":
                    # Lower values are better, so negate z-score
                    # Positive z-score = improvement
                    zscore = -zscore
                # For maximize, positive z-score already means improvement

                # Apply minimum observations threshold
                if count >= min_observations:
                    zscore_row[aa] = zscore
                else:
                    zscore_row[aa] = 0.0
            else:
                # Mutation not observed
                zscore_row[aa] = 0.0

        zscore_rows.append(zscore_row)

    zscore_df = pd.DataFrame(zscore_rows)

    logger.info(f"Calculated z-scores for {len(zscore_df)} positions")

    return zscore_df


def calculate_top_mutations(stats_df: pd.DataFrame,
                            delta_df: pd.DataFrame,
                            zscore_df: pd.DataFrame,
                            primary_metric: str,
                            metric_columns: List[str]) -> pd.DataFrame:
    """
    Calculate best mutation at each position.

    Args:
        stats_df: Mutation statistics DataFrame
        delta_df: Delta scores DataFrame
        zscore_df: Z-score DataFrame
        primary_metric: Primary metric used for scoring
        metric_columns: All metric columns

    Returns:
        DataFrame with top mutations per position
    """
    top_rows = []

    for pos in sorted(stats_df['position'].unique()):
        pos_data = stats_df[stats_df['position'] == pos]
        wt_aa = pos_data['wt_aa'].iloc[0]

        # Find mutation with highest delta score
        delta_row = delta_df[delta_df['position'] == pos].iloc[0]
        zscore_row = zscore_df[zscore_df['position'] == pos].iloc[0]

        best_aa = None
        best_delta = -np.inf
        best_zscore = 0.0
        best_count = 0

        for aa in AMINO_ACIDS:
            delta = delta_row[aa]
            if delta > best_delta:
                best_delta = delta
                best_aa = aa
                best_zscore = zscore_row[aa]

                # Get count from stats
                mut_data = pos_data[pos_data['mut_aa'] == aa]
                if len(mut_data) > 0:
                    best_count = mut_data['count'].iloc[0]

        # Get metric means for best mutation
        top_row = {
            'position': pos,
            'wt_aa': wt_aa,
            'best_mutation': best_aa if best_aa else wt_aa,
            'count': best_count,
            'delta_score': best_delta,
            'zscore': best_zscore
        }

        # Add metric means
        if best_aa:
            mut_data = pos_data[pos_data['mut_aa'] == best_aa]
            if len(mut_data) > 0:
                for metric in metric_columns:
                    mean_col = f"{metric}_mean"
                    if mean_col in mut_data.columns:
                        top_row[f"{metric}_mean"] = mut_data[mean_col].iloc[0]

        top_rows.append(top_row)

    return pd.DataFrame(top_rows)


def calculate_coverage(stats_df: pd.DataFrame) -> pd.DataFrame:
    """
    Calculate position exploration statistics.

    Args:
        stats_df: Mutation statistics DataFrame

    Returns:
        DataFrame with coverage statistics per position
    """
    coverage_rows = []

    for pos in sorted(stats_df['position'].unique()):
        pos_data = stats_df[stats_df['position'] == pos]
        wt_aa = pos_data['wt_aa'].iloc[0]

        n_observations = pos_data['count'].sum()
        n_mutations_tested = len(pos_data)
        n_mutations_possible = 20  # 20 amino acids

        coverage_fraction = n_mutations_tested / n_mutations_possible
        max_count = pos_data['count'].max()
        min_count = pos_data['count'].min()
        mean_count = pos_data['count'].mean()

        coverage_rows.append({
            'position': pos,
            'wt_aa': wt_aa,
            'n_observations': n_observations,
            'n_mutations_tested': n_mutations_tested,
            'coverage_fraction': coverage_fraction,
            'max_count': max_count,
            'min_count': min_count,
            'mean_count': mean_count
        })

    return pd.DataFrame(coverage_rows)


def main():
    """Main execution function."""
    parser = argparse.ArgumentParser(description='Analyze sequence-metric correlations')
    parser.add_argument('--config', required=True, help='Configuration JSON file')

    args = parser.parse_args()
    logger = setup_logging()

    try:
        # Load configuration
        with open(args.config, 'r') as f:
            config = json.load(f)

        logger.info(f"Loading sequences from: {config['sequences_path']}")
        logger.info(f"Loading metrics from: {config['metrics_path']}")

        # Load data
        sequences_df = pd.read_csv(config['sequences_path'])
        metrics_df = pd.read_csv(config['metrics_path'])

        # Load reference sequence
        reference_seq = load_reference_sequence(config['reference_sequence'])
        logger.info(f"Reference sequence length: {len(reference_seq)}")

        # Load history if provided
        history_df = None
        if config['history_path']:
            logger.info(f"Loading history from: {config['history_path']}")
            history_df = pd.read_csv(config['history_path'])

        logger.info(f"Loaded {len(sequences_df)} sequences and {len(metrics_df)} metric records")

        # Analyze mutations
        logger.info("Analyzing sequence-metric correlations...")
        stats_df = analyze_sequence_metrics(
            sequences_df,
            metrics_df,
            reference_seq,
            config['metric_columns'],
            history_df
        )

        if len(stats_df) == 0:
            logger.error("No mutations found - cannot proceed")
            sys.exit(1)

        # Calculate delta scores
        logger.info("Calculating delta scores...")
        delta_df = calculate_delta_scores(
            stats_df,
            config['primary_metric'],
            config['mode'],
            config['min_observations']
        )

        # Calculate z-scores
        logger.info("Calculating z-scores...")
        zscore_df = calculate_zscore_scores(
            stats_df,
            config['primary_metric'],
            config['mode'],
            config['min_observations']
        )

        # Calculate top mutations
        logger.info("Identifying top mutations per position...")
        top_df = calculate_top_mutations(
            stats_df,
            delta_df,
            zscore_df,
            config['primary_metric'],
            config['metric_columns']
        )

        # Calculate coverage
        logger.info("Calculating coverage statistics...")
        coverage_df = calculate_coverage(stats_df)

        # Save results
        logger.info(f"Saving mutation statistics to: {config['mutation_statistics_output']}")
        stats_df.to_csv(config['mutation_statistics_output'], index=False)

        logger.info(f"Saving mutation deltas to: {config['mutation_deltas_output']}")
        delta_df.to_csv(config['mutation_deltas_output'], index=False)

        logger.info(f"Saving mutation z-scores to: {config['mutation_zscores_output']}")
        zscore_df.to_csv(config['mutation_zscores_output'], index=False)

        logger.info(f"Saving top mutations to: {config['top_mutations_output']}")
        top_df.to_csv(config['top_mutations_output'], index=False)

        logger.info(f"Saving coverage to: {config['coverage_output']}")
        coverage_df.to_csv(config['coverage_output'], index=False)

        # Print summary
        primary_metric = config['primary_metric']
        primary_metric_mean_col = f"{primary_metric}_mean"

        logger.info("Analysis complete!")
        logger.info(f"  Total unique mutations: {len(stats_df)}")
        logger.info(f"  Positions analyzed: {stats_df['position'].nunique()}")
        logger.info(f"  Total observations: {stats_df['count'].sum()}")
        logger.info(f"  Average observations per mutation: {stats_df['count'].mean():.2f}")
        logger.info(f"  Primary metric ({primary_metric}):")
        logger.info(f"    Mean: {stats_df[primary_metric_mean_col].mean():.3f}")
        logger.info(f"    Std: {stats_df[primary_metric_mean_col].std():.3f}")

    except Exception as e:
        logger.error(f"Error during sequence-metric analysis: {str(e)}")
        import traceback
        logger.error(traceback.format_exc())
        sys.exit(1)


if __name__ == "__main__":
    main()
