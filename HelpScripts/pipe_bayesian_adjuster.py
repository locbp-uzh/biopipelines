#!/usr/bin/env python3
# Copyright (c) 2026 Gianluca Quargnali @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Bayesian frequency adjustment runtime script for BioPipelines.

Applies Bayesian log-odds updates to mutation frequency tables using correlation signals.

Formula:
    p(i,aa|c) = σ(σ⁻¹(p₀(i,aa)) + γ·c(i,aa))

where:
    - p₀(i,aa) = prior probability from MutationProfiler frequencies
    - c(i,aa) = correlation signal from SequenceMetricCorrelation
    - γ = strength hyperparameter
    - σ(x) = 1/(1+e⁻ˣ) = sigmoid function
    - σ⁻¹(p) = log(p/(1-p)) = logit function
"""

import argparse
import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import logomaker
from pathlib import Path
import logging
from typing import Dict, List, Tuple
import sys

# Standard amino acids
AMINO_ACIDS = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

# Small epsilon to avoid log(0) and division by zero
EPSILON = 1e-9


def parse_positions_selection(selection_str):
    """
    Parse PyMOL-style position selection string.

    Args:
        selection_str: Selection string like "141+143+145+147-149+151-152"
                      '+' separates individual positions or ranges
                      '-' indicates a range (e.g., "147-149" means 147, 148, 149)

    Returns:
        Sorted list of positions (1-indexed), or None if selection_str is None/empty
    """
    if not selection_str or not str(selection_str).strip():
        return None

    positions = []
    parts = str(selection_str).split('+')

    for part in parts:
        part = part.strip()
        if not part:
            continue

        if '-' in part and not part.startswith('-'):
            # Range like "147-149"
            range_parts = part.split('-')
            if len(range_parts) == 2:
                try:
                    start = int(range_parts[0])
                    end = int(range_parts[1])
                    positions.extend(range(start, end + 1))
                except ValueError:
                    pass  # Skip malformed ranges
        else:
            # Single position
            try:
                positions.append(int(part))
            except ValueError:
                pass  # Skip malformed positions

    return sorted(set(positions)) if positions else None


def setup_logging():
    """Set up logging configuration."""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    return logging.getLogger(__name__)


def sigmoid(x):
    """Sigmoid function: σ(x) = 1/(1+e⁻ˣ)"""
    return 1.0 / (1.0 + np.exp(-x))


def logit(p):
    """Logit function (inverse sigmoid): σ⁻¹(p) = log(p/(1-p))"""
    # Clip to avoid log(0) or log(negative)
    p_clipped = np.clip(p, EPSILON, 1 - EPSILON)
    return np.log(p_clipped / (1 - p_clipped))


def bayesian_update(prior_prob: float, correlation: float, gamma: float) -> float:
    """
    Apply Bayesian log-odds update to a single probability.

    Args:
        prior_prob: Prior probability p₀(i,aa)
        correlation: Correlation signal c(i,aa)
        gamma: Strength hyperparameter

    Returns:
        Updated probability p(i,aa|c)
    """
    # Handle edge cases
    if prior_prob <= 0:
        prior_prob = EPSILON
    elif prior_prob >= 1:
        prior_prob = 1 - EPSILON

    # Apply Bayesian update
    logit_prior = logit(prior_prob)
    logit_updated = logit_prior + gamma * correlation
    updated_prob = sigmoid(logit_updated)

    return updated_prob


def load_and_merge_tables(frequencies_path: str, correlations_path: str,
                          logger) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Load frequency and correlation tables.

    Args:
        frequencies_path: Path to frequencies CSV
        correlations_path: Path to correlations CSV
        logger: Logger instance

    Returns:
        Tuple of (frequencies_df, correlations_df)
    """
    logger.info(f"Loading frequencies from: {frequencies_path}")
    frequencies_df = pd.read_csv(frequencies_path)

    logger.info(f"Loading correlations from: {correlations_path}")
    correlations_df = pd.read_csv(correlations_path)

    # Validate columns
    required_cols = ['position', 'original'] + AMINO_ACIDS
    for col in required_cols:
        if col not in frequencies_df.columns:
            raise ValueError(f"Missing column '{col}' in frequencies table")
        if col.replace('original', 'wt_aa') not in correlations_df.columns:
            raise ValueError(f"Missing column '{col}' in correlations table")

    logger.info(f"Loaded {len(frequencies_df)} positions from frequencies table")
    logger.info(f"Loaded {len(correlations_df)} positions from correlations table")

    return frequencies_df, correlations_df


def apply_pseudocounts(frequencies_df: pd.DataFrame, correlations_df: pd.DataFrame,
                       pseudocount: float, mode: str, logger) -> pd.DataFrame:
    """
    Apply pseudocounts only to amino acids with good correlation signals and renormalize.

    Only amino acids where c(i,aa) is nonzero and good (positive for max mode, negative for min mode)
    receive pseudocounts. This prevents all amino acids from having nonzero probability after adjustment.

    Args:
        frequencies_df: DataFrame with prior frequencies
        correlations_df: DataFrame with correlation signals
        pseudocount: Pseudocount to add to amino acids with good correlations
        mode: "min" or "max" (determines what "good" correlation means)
        logger: Logger instance

    Returns:
        DataFrame with selective pseudocounts applied and renormalized
    """
    if pseudocount == 0:
        logger.info("Pseudocount is 0, skipping pseudocount application")
        return frequencies_df.copy()

    logger.info(f"Applying pseudocount: {pseudocount} (mode: {mode})")

    # Sign factor determines what "good" means
    # min mode: negative correlations are good (lower metric is better)
    # max mode: positive correlations are good (higher metric is better)
    sign_factor = -1.0 if mode == "min" else 1.0

    pseudocount_df = frequencies_df.copy()
    total_pseudocounts_added = 0

    for idx, row in pseudocount_df.iterrows():
        position = row['position']

        # Find matching correlation row
        corr_row = correlations_df[correlations_df['position'] == position]
        if corr_row.empty:
            logger.warning(f"No correlation data for position {position}, skipping pseudocount")
            continue

        corr_row = corr_row.iloc[0]

        # Calculate original sum across all amino acids
        original_sum = sum(row[aa] for aa in AMINO_ACIDS)

        # Add pseudocount only to amino acids with good correlations
        for aa in AMINO_ACIDS:
            correlation = corr_row[aa] * sign_factor
            # Good correlation: positive after sign correction, and nonzero
            if correlation > EPSILON:
                pseudocount_df.at[idx, aa] = row[aa] + pseudocount
                total_pseudocounts_added += 1
            # else: keep original frequency (no pseudocount added)

        # Calculate new sum after adding pseudocounts
        new_sum = sum(pseudocount_df.loc[idx, aa] for aa in AMINO_ACIDS)

        # Renormalize to preserve original sum
        if new_sum > 0 and original_sum > 0:
            normalization_factor = original_sum / new_sum
            for aa in AMINO_ACIDS:
                pseudocount_df.at[idx, aa] = pseudocount_df.at[idx, aa] * normalization_factor

    logger.info(f"Pseudocounts applied selectively to {total_pseudocounts_added} amino acids with good correlations")
    return pseudocount_df


def apply_bayesian_adjustment(frequencies_df: pd.DataFrame,
                              correlations_df: pd.DataFrame,
                              mode: str,
                              gamma: float,
                              pseudocount: float,
                              logger) -> Tuple[pd.DataFrame, List[Dict]]:
    """
    Apply Bayesian adjustment to all positions and amino acids.

    Args:
        frequencies_df: DataFrame with prior frequencies
        correlations_df: DataFrame with correlation signals
        mode: "min" or "max" (determines sign of correlations)
        gamma: Strength hyperparameter
        pseudocount: Pseudocount to add only to amino acids with good correlations
        logger: Logger instance

    Returns:
        Tuple of (adjusted_probabilities_df, adjustment_log)
    """
    logger.info("Applying Bayesian adjustments...")

    # If gamma is 0, return original frequencies unchanged
    if gamma == 0:
        logger.info("Gamma is 0, returning original frequencies without adjustment")
        adjusted_df = frequencies_df.copy()

        # Create empty adjustment log
        adjustment_log = []
        for _, row in adjusted_df.iterrows():
            position = row['position']
            original_aa = row['original']
            for aa in AMINO_ACIDS:
                adjustment_log.append({
                    'position': position,
                    'wt_aa': original_aa,
                    'aa': aa,
                    'prior_freq': row[aa],
                    'correlation': 0.0,
                    'adjusted_prob': row[aa],
                    'change': 0.0
                })

        return adjusted_df, adjustment_log

    # Apply pseudocounts selectively to amino acids with good correlations
    frequencies_df = apply_pseudocounts(frequencies_df, correlations_df, pseudocount, mode, logger)

    adjusted_data = []
    adjustment_log = []

    # Sign correction based on mode
    sign_factor = -1.0 if mode == "min" else 1.0
    logger.info(f"Mode: {mode}, sign factor: {sign_factor}")

    for _, freq_row in frequencies_df.iterrows():
        position = freq_row['position']
        original_aa = freq_row['original']

        # Find matching correlation row
        corr_row = correlations_df[correlations_df['position'] == position]
        if corr_row.empty:
            logger.warning(f"No correlation data for position {position}, skipping")
            continue

        corr_row = corr_row.iloc[0]

        # Verify original amino acid matches
        if corr_row.get('wt_aa', corr_row.get('original')) != original_aa:
            logger.warning(f"Original AA mismatch at position {position}: "
                          f"freq={original_aa}, corr={corr_row.get('wt_aa')}")

        # Adjust each amino acid
        adjusted_row = {'position': position, 'original': original_aa}

        for aa in AMINO_ACIDS:
            prior_freq = freq_row[aa]
            correlation = corr_row[aa] * sign_factor  # Apply sign correction

            # Apply Bayesian update
            adjusted_prob = bayesian_update(prior_freq, correlation, gamma)

            adjusted_row[aa] = adjusted_prob

            # Log the adjustment
            adjustment_log.append({
                'position': position,
                'wt_aa': original_aa,
                'aa': aa,
                'prior_freq': prior_freq,
                'correlation': correlation,
                'adjusted_prob': adjusted_prob,
                'change': adjusted_prob - prior_freq
            })

        adjusted_data.append(adjusted_row)

    adjusted_df = pd.DataFrame(adjusted_data)

    logger.info(f"Adjusted {len(adjusted_df)} positions")

    return adjusted_df, adjustment_log


def normalize_absolute(adjusted_df: pd.DataFrame, logger) -> pd.DataFrame:
    """
    Normalize adjusted probabilities as absolute probabilities.

    This makes them comparable to MutationProfiler's absolute_frequencies.
    Each value is normalized by the sum across all amino acids at that position.

    Args:
        adjusted_df: DataFrame with adjusted probabilities
        logger: Logger instance

    Returns:
        DataFrame with absolute probabilities
    """
    logger.info("Normalizing as absolute probabilities...")

    absolute_df = adjusted_df.copy()

    for idx, row in absolute_df.iterrows():
        position = row['position']
        original_aa = row['original']

        # Calculate sum across all amino acids
        total = sum(row[aa] for aa in AMINO_ACIDS)

        if total > 0:
            # Normalize each amino acid by total
            for aa in AMINO_ACIDS:
                absolute_df.at[idx, aa] = row[aa] / total
        else:
            logger.warning(f"Position {position}: total probability is zero, skipping normalization")

    return absolute_df


def normalize_relative(adjusted_df: pd.DataFrame, logger) -> pd.DataFrame:
    """
    Normalize adjusted probabilities as relative probabilities.

    This makes them comparable to MutationProfiler's relative_frequencies.
    - Original amino acid is set to 0.0
    - Other amino acids are normalized to sum to 1.0

    Args:
        adjusted_df: DataFrame with adjusted probabilities
        logger: Logger instance

    Returns:
        DataFrame with relative probabilities
    """
    logger.info("Normalizing as relative probabilities...")

    relative_df = adjusted_df.copy()

    for idx, row in relative_df.iterrows():
        position = row['position']
        original_aa = row['original']

        # Set original amino acid to 0.0
        relative_df.at[idx, original_aa] = 0.0

        # Calculate sum across non-original amino acids
        total = sum(row[aa] for aa in AMINO_ACIDS if aa != original_aa)

        if total > 0:
            # Normalize non-original amino acids
            for aa in AMINO_ACIDS:
                if aa != original_aa:
                    relative_df.at[idx, aa] = row[aa] / total
        else:
            logger.warning(f"Position {position}: total non-original probability is zero")

    return relative_df


def create_logo_plot(prob_df: pd.DataFrame, title: str, ylabel: str,
                    svg_path: str, png_path: str, logger, positions_filter=None, pseudocount=None):
    """
    Create sequence logo visualization for probabilities.

    Follows the same pattern as MutationProfiler:
    - Filters positions with significant mutations
    - Excludes original amino acid from visualization
    - Uses logomaker with chemistry color scheme

    Args:
        prob_df: DataFrame with probabilities
        title: Plot title
        ylabel: Y-axis label
        svg_path: Output SVG file path
        png_path: Output PNG file path
        logger: Logger instance
        positions_filter: List of positions to include (1-indexed). If None, includes positions where any amino acid has p(i,aa) > threshold.
        pseudocount: Threshold for significance. If None, uses 0.01 (1%).
    """
    logger.info(f"Creating logo plot: {title}")

    # Determine threshold: pseudocount if given, otherwise 1%
    threshold = pseudocount if pseudocount is not None else 0.01

    # Filter positions based on positions_filter or significance threshold
    positions_to_plot = []

    for i, row in prob_df.iterrows():
        position = int(row['position'])
        original_aa = row['original']

        # Apply position filter if specified
        if positions_filter is not None:
            if position not in positions_filter:
                continue
            # Include all specified positions
        else:
            # Default: only include positions where at least one non-original amino acid
            # has probability > threshold (this includes resurrected mutations)
            has_significant = any(row[aa] > threshold for aa in AMINO_ACIDS if aa != original_aa)
            if not has_significant:
                continue

        positions_to_plot.append(row)

    if not positions_to_plot:
        # Create empty plot
        fig, ax = plt.subplots(figsize=(10, 4))
        ax.text(0.5, 0.5, 'No significant probabilities found\n(all values below threshold)',
                ha='center', va='center', transform=ax.transAxes, fontsize=16)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis('off')
        plt.tight_layout()
        plt.savefig(svg_path, format='svg', dpi=300, bbox_inches='tight')
        plt.savefig(png_path, format='png', dpi=300, bbox_inches='tight')
        plt.close()
        logger.warning("No positions with significant probabilities for logo plot")
        return

    # Prepare data for logomaker
    logo_data = []
    x_labels = []

    for row in positions_to_plot:
        position = row['position']
        original_aa = row['original']

        # Create row excluding original amino acid
        aa_row = {}
        for aa in AMINO_ACIDS:
            if aa != original_aa:
                aa_row[aa] = row[aa]
            else:
                aa_row[aa] = 0.0  # Exclude original from visualization

        logo_data.append(aa_row)
        x_labels.append(f"{original_aa}{int(position)}")

    # Create figure
    num_positions = len(positions_to_plot)
    fig, ax = plt.subplots(figsize=(max(12, num_positions * 0.8), 6))

    # Convert to DataFrame
    logo_df = pd.DataFrame(logo_data)
    logo_df.index = range(len(logo_df))

    # Create logo
    logo = logomaker.Logo(logo_df, ax=ax, color_scheme='chemistry')

    # Customize plot
    ax.set_ylabel(ylabel, fontsize=14)
    ax.set_title(title, fontsize=16, fontweight='bold')
    ax.set_xlabel('Position', fontsize=14)

    # Set x-ticks
    ax.set_xticks(range(len(x_labels)))
    ax.set_xticklabels(x_labels, rotation=45, ha='right')

    # Clean up appearance
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    plt.savefig(svg_path, format='svg', dpi=300, bbox_inches='tight')
    plt.savefig(png_path, format='png', dpi=300, bbox_inches='tight')
    plt.close()

    logger.info(f"Saved logo plot to {svg_path} and {png_path}")


def main():
    """Main execution function."""
    parser = argparse.ArgumentParser(description='Apply Bayesian frequency adjustment')
    parser.add_argument('--config', required=True, help='Configuration JSON file')

    args = parser.parse_args()
    logger = setup_logging()

    try:
        # Load configuration
        with open(args.config, 'r') as f:
            config = json.load(f)

        logger.info("Starting Bayesian frequency adjustment...")

        # Load tables
        frequencies_df, correlations_df = load_and_merge_tables(
            config['frequencies_path'],
            config['correlations_path'],
            logger
        )

        # Apply Bayesian adjustment
        adjusted_df, adjustment_log = apply_bayesian_adjustment(
            frequencies_df,
            correlations_df,
            config['mode'],
            config['gamma'],
            config['pseudocount'],
            logger
        )

        # Save adjusted probabilities (raw Bayesian output)
        logger.info(f"Saving adjusted probabilities to: {config['adjusted_probabilities_output']}")
        adjusted_df.to_csv(config['adjusted_probabilities_output'], index=False)

        # Normalize as absolute probabilities
        absolute_df = normalize_absolute(adjusted_df, logger)
        logger.info(f"Saving absolute probabilities to: {config['absolute_probabilities_output']}")
        absolute_df.to_csv(config['absolute_probabilities_output'], index=False)

        # Normalize as relative probabilities
        relative_df = normalize_relative(adjusted_df, logger)
        logger.info(f"Saving relative probabilities to: {config['relative_probabilities_output']}")
        relative_df.to_csv(config['relative_probabilities_output'], index=False)

        # Save adjustment log
        logger.info(f"Saving adjustment log to: {config['adjustment_log_output']}")
        adjustment_log_df = pd.DataFrame(adjustment_log)
        adjustment_log_df.to_csv(config['adjustment_log_output'], index=False)

        # Parse positions filter if provided
        positions_filter = parse_positions_selection(config.get('positions'))
        if positions_filter:
            logger.info(f"Using positions filter: {positions_filter} ({len(positions_filter)} positions)")

        # Get pseudocount for threshold
        pseudocount = config.get('pseudocount')

        # Create logo plots
        create_logo_plot(
            adjusted_df,
            'Bayesian Adjusted Probabilities',
            'Adjusted Probability',
            config['adjusted_logo_svg'],
            config['adjusted_logo_png'],
            logger,
            positions_filter=positions_filter,
            pseudocount=pseudocount
        )

        create_logo_plot(
            absolute_df,
            'Absolute Probabilities',
            'Absolute Probability',
            config['absolute_logo_svg'],
            config['absolute_logo_png'],
            logger,
            positions_filter=positions_filter,
            pseudocount=pseudocount
        )

        create_logo_plot(
            relative_df,
            'Relative Probabilities',
            'Relative Probability',
            config['relative_logo_svg'],
            config['relative_logo_png'],
            logger,
            positions_filter=positions_filter,
            pseudocount=pseudocount
        )

        # Print summary
        logger.info("Adjustment complete!")
        logger.info(f"  Positions adjusted: {len(adjusted_df)}")
        logger.info(f"  Mode: {config['mode']}")
        logger.info(f"  Gamma: {config['gamma']}")
        logger.info(f"  Pseudocount: {config['pseudocount']}")

        # Calculate and log some statistics
        total_change = adjustment_log_df['change'].abs().sum()
        max_change = adjustment_log_df['change'].abs().max()
        logger.info(f"  Total absolute change: {total_change:.3f}")
        logger.info(f"  Maximum absolute change: {max_change:.3f}")

    except Exception as e:
        logger.error(f"Error during Bayesian adjustment: {str(e)}")
        import traceback
        logger.error(traceback.format_exc())
        sys.exit(1)


if __name__ == "__main__":
    main()
