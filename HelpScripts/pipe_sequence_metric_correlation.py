#!/usr/bin/env python3
"""
Sequence-metric correlation runtime script for BioPipelines.

Computes correlation signals c(i) and c(i,aa) where:

c(i) = (m̄_i-mut - m̄_i-wt) / sqrt(s²_i-mut + s²_i-wt)

where:
  - m̄_i-mut = average metric over sequences where position i is mutated
  - m̄_i-wt  = average metric over sequences where position i is NOT mutated
  - s²_i-mut = unbiased variance of metric for mutated group (ddof=1)
  - s²_i-wt  = unbiased variance of metric for wild-type group (ddof=1)

c(i,aa) = (m̄_i,aa - m̄_i,¬aa) / sqrt(s²_i,aa + s²_i,¬aa)

where:
  - m̄_i,aa   = average metric over sequences with amino acid 'aa' at position i
  - m̄_i,¬aa  = average metric over sequences with other amino acids at position i
  - s²_i,aa  = unbiased variance for sequences with 'aa' at position i
  - s²_i,¬aa = unbiased variance for sequences without 'aa' at position i

When n <= 1 for either group, variance cannot be reliably computed and correlation is set to 0.
"""

import argparse
import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import logomaker
from pathlib import Path
import logging
from typing import Dict, List, Tuple, Optional
import sys

# Standard amino acids
AMINO_ACIDS = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']


def parse_positions_selection(selection_str: Optional[str]) -> Optional[List[int]]:
    """
    Parse PyMOL-style position selection string.

    Args:
        selection_str: Selection string like "141+143+145+147-149+151-152"
                      '+' separates individual positions or ranges
                      '-' indicates a range (e.g., "147-149" means 147, 148, 149)

    Returns:
        Sorted list of positions (1-indexed), or None if selection_str is None/empty
    """
    if not selection_str or not selection_str.strip():
        return None

    positions = []
    parts = selection_str.split('+')

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


def load_and_merge_data(mutants_paths: List[str], data_paths: List[str],
                       metric: str, logger) -> pd.DataFrame:
    """
    Load and merge all mutants and data tables.

    Args:
        mutants_paths: List of paths to mutants tables
        data_paths: List of paths to data tables
        metric: Metric column name
        logger: Logger instance

    Returns:
        Merged DataFrame with id, sequence, and metric columns
    """
    all_data = []

    for mutants_path, data_path in zip(mutants_paths, data_paths):
        logger.info(f"Loading mutants from: {mutants_path}")
        logger.info(f"Loading data from: {data_path}")

        mutants_df = pd.read_csv(mutants_path)
        data_df = pd.read_csv(data_path)

        # Find sequence column in mutants
        seq_col = None
        for col in mutants_df.columns:
            if col.lower() in ['sequence', 'seq']:
                seq_col = col
                break

        if seq_col is None:
            raise ValueError(f"Could not find sequence column in {mutants_path}")

        # Check metric column exists in data
        if metric not in data_df.columns:
            raise ValueError(f"Metric column '{metric}' not found in {data_path}. Available columns: {list(data_df.columns)}")

        # Merge on id
        merged = pd.merge(mutants_df[['id', seq_col]], data_df[['id', metric]], on='id', how='inner')
        merged = merged.rename(columns={seq_col: 'sequence'})

        all_data.append(merged)

    # Concatenate all data
    combined_df = pd.concat(all_data, ignore_index=True)

    logger.info(f"Loaded {len(combined_df)} total sequences with metric values")

    return combined_df


def compute_correlation_1d(merged_df: pd.DataFrame, reference_seq: str,
                          metric: str, logger) -> pd.DataFrame:
    """
    Compute 1D correlation signal c(i) for each position.

    c(i) = (mean_mutated - mean_wt) / sqrt(var_mutated + var_wt)

    Args:
        merged_df: DataFrame with id, sequence, and metric columns
        reference_seq: Reference sequence
        metric: Metric column name
        logger: Logger instance

    Returns:
        DataFrame with correlation_1d results
    """
    logger.info("Computing 1D correlation signals...")

    correlation_1d_data = []
    seq_length = len(reference_seq)

    for pos in range(seq_length):
        position_label = pos + 1  # 1-indexed
        ref_aa = reference_seq[pos]

        # Split sequences into mutated vs wild-type at this position
        mutated_values = []
        wt_values = []

        for _, row in merged_df.iterrows():
            sequence = row['sequence']
            metric_value = row[metric]

            if pd.isna(metric_value):
                continue

            if pos < len(sequence):
                current_aa = sequence[pos]
                if current_aa != ref_aa:
                    mutated_values.append(metric_value)
                else:
                    wt_values.append(metric_value)

        n_mutated = len(mutated_values)
        n_wt = len(wt_values)

        # Compute correlation only if both groups have observations
        if n_mutated > 0 and n_wt > 0:
            mean_mutated = np.mean(mutated_values)
            mean_wt = np.mean(wt_values)

            # Compute unbiased variance (ddof=1)
            # If n <= 1, var returns NaN which we'll handle
            var_mutated = np.var(mutated_values, ddof=1) if n_mutated > 1 else 0.0
            var_wt = np.var(wt_values, ddof=1) if n_wt > 1 else 0.0

            # Handle cases where variance cannot be computed reliably
            # When either group has n <= 1, we can't compute unbiased variance
            if (n_mutated <= 1 or n_wt <= 1):
                correlation = 0.0
            elif (var_mutated + var_wt) == 0:
                # Both variances are zero (all values identical in both groups)
                correlation = 0.0
            else:
                # Compute effect size: (difference in means) / sqrt(sum of variances)
                # This is the standardized mean difference accounting for variance in both groups
                correlation = (mean_mutated - mean_wt) / np.sqrt(var_mutated + var_wt)
        else:
            mean_mutated = 0.0
            mean_wt = 0.0
            var_mutated = 0.0
            var_wt = 0.0
            correlation = 0.0

        correlation_1d_data.append({
            'position': position_label,
            'wt_aa': ref_aa,
            'correlation': correlation,
            'mean_mutated': mean_mutated,
            'mean_wt': mean_wt,
            'var_mutated': var_mutated,
            'var_wt': var_wt,
            'n_mutated': n_mutated,
            'n_wt': n_wt
        })

    correlation_1d_df = pd.DataFrame(correlation_1d_data)

    logger.info(f"Computed 1D correlations for {len(correlation_1d_df)} positions")

    return correlation_1d_df


def compute_correlation_2d(merged_df: pd.DataFrame, reference_seq: str,
                          metric: str, logger) -> pd.DataFrame:
    """
    Compute 2D correlation signal c(i,aa) for each position and amino acid.

    c(i,aa) = (mean_aa - mean_non_aa) / sqrt(var_aa + var_non_aa)

    Args:
        merged_df: DataFrame with id, sequence, and metric columns
        reference_seq: Reference sequence
        metric: Metric column name
        logger: Logger instance

    Returns:
        DataFrame with correlation_2d results (position, wt_aa, A, C, D, ...)
    """
    logger.info("Computing 2D correlation signals...")

    correlation_2d_data = []
    seq_length = len(reference_seq)

    for pos in range(seq_length):
        position_label = pos + 1  # 1-indexed
        ref_aa = reference_seq[pos]

        # Collect all metric values grouped by amino acid at this position
        # This includes both mutated and wild-type sequences
        position_aa_metrics = {}

        for _, row in merged_df.iterrows():
            sequence = row['sequence']
            metric_value = row[metric]

            if pd.isna(metric_value):
                continue

            if pos < len(sequence):
                current_aa = sequence[pos]

                # Group metric values by the amino acid at this position
                if current_aa not in position_aa_metrics:
                    position_aa_metrics[current_aa] = []
                position_aa_metrics[current_aa].append(metric_value)

        # Compute correlation for each amino acid
        row_data = {'position': position_label, 'wt_aa': ref_aa}

        for aa in AMINO_ACIDS:
            # Get values for sequences with this aa at position i vs all other aa
            aa_values = position_aa_metrics.get(aa, [])
            non_aa_values = []
            for other_aa, values in position_aa_metrics.items():
                if other_aa != aa:
                    non_aa_values.extend(values)

            n_aa = len(aa_values)
            n_non_aa = len(non_aa_values)

            if n_aa > 0 and n_non_aa > 0:
                mean_aa = np.mean(aa_values)
                mean_non_aa = np.mean(non_aa_values)

                # Compute unbiased variance (ddof=1)
                var_aa = np.var(aa_values, ddof=1) if n_aa > 1 else 0.0
                var_non_aa = np.var(non_aa_values, ddof=1) if n_non_aa > 1 else 0.0

                # Handle cases where variance cannot be computed reliably
                # When either group has n <= 1, we can't compute unbiased variance
                if (n_aa <= 1 or n_non_aa <= 1):
                    correlation = 0.0
                elif (var_aa + var_non_aa) == 0:
                    # Both variances are zero (all values identical in both groups)
                    correlation = 0.0
                else:
                    # Compute effect size: (difference in means) / sqrt(sum of variances)
                    # Positive correlation: this amino acid improves the metric
                    # Negative correlation: this amino acid worsens the metric
                    correlation = (mean_aa - mean_non_aa) / np.sqrt(var_aa + var_non_aa)
            else:
                correlation = 0.0

            row_data[aa] = correlation

        correlation_2d_data.append(row_data)

    correlation_2d_df = pd.DataFrame(correlation_2d_data)

    logger.info(f"Computed 2D correlations for {len(correlation_2d_df)} positions")

    return correlation_2d_df


def create_correlation_logo(correlation_2d_df: pd.DataFrame, correlation_1d_df: pd.DataFrame,
                           svg_path: str, png_path: str, logger, positions_filter: Optional[List[int]] = None):
    """
    Create sequence logo visualization for correlations.

    Positive correlations (c > 0) are stacked above the x-axis.
    Negative correlations (c < 0) are stacked below the x-axis.
    Height is proportional to |c(i,aa)|.

    Args:
        correlation_2d_df: DataFrame with 2D correlations
        correlation_1d_df: DataFrame with 1D correlations
        svg_path: Output SVG file path
        png_path: Output PNG file path
        logger: Logger instance
        positions_filter: List of positions to include (1-indexed). If None, includes all positions in the sequence (both mutated >1% and positions with correlations).
    """
    logger.info("Creating correlation logo plot...")

    # Filter positions based on positions_filter or mutation/correlation
    positions_to_plot = []

    for i, row in correlation_2d_df.iterrows():
        position = int(row['position'])
        corr_1d_row = correlation_1d_df.iloc[i]

        # Apply position filter if specified
        if positions_filter is not None:
            if position not in positions_filter:
                continue
        else:
            # Default: include positions that either have mutation or correlation
            has_correlation = any(abs(row[aa]) > 1e-6 for aa in AMINO_ACIDS)
            was_mutated = corr_1d_row['n_mutated'] > 0

            if not (has_correlation or was_mutated):
                continue

        positions_to_plot.append((row, corr_1d_row))

    if not positions_to_plot:
        # Create empty plot
        fig, ax = plt.subplots(figsize=(10, 4))
        ax.text(0.5, 0.5, 'No correlations found\n(all c(i,aa) = 0 at all positions)',
                ha='center', va='center', transform=ax.transAxes, fontsize=16)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis('off')
        plt.tight_layout()
        plt.savefig(svg_path, format='svg', dpi=300, bbox_inches='tight')
        plt.savefig(png_path, format='png', dpi=300, bbox_inches='tight')
        plt.close()
        logger.warning("All correlations are zero. Possible causes: insufficient data (n≤1), no metric variation, or all sequences identical.")
        return

    # Prepare data for logomaker
    # Create a single dataframe with all correlation values (both positive and negative)
    logo_data = []
    x_labels = []

    for row in positions_to_plot:
        position = row['position']
        wt_aa = row['wt_aa']

        aa_row = {}
        for aa in AMINO_ACIDS:
            aa_row[aa] = row[aa]  # Keep both positive and negative values

        logo_data.append(aa_row)
        x_labels.append(f"{wt_aa}{position}")

    # Create figure
    num_positions = len(positions_to_plot)
    fig, ax = plt.subplots(figsize=(max(12, num_positions * 0.8), 6))

    # Convert to DataFrame
    logo_df = pd.DataFrame(logo_data)
    logo_df.index = range(len(logo_df))

    # Create logo with both positive and negative values
    # logomaker handles positive and negative values automatically
    # flip_below=False keeps negative letters upright (not flipped upside down)
    logo = logomaker.Logo(logo_df, ax=ax, color_scheme='chemistry',
                         baseline_width=0, show_spines=False, flip_below=False)

    # Customize plot
    ax.set_ylabel('Correlation Signal c(i,aa)', fontsize=14)
    ax.set_title('Mutation-Metric Correlation Logo', fontsize=16, fontweight='bold')
    ax.set_xlabel('Position', fontsize=14)

    # Set x-ticks
    ax.set_xticks(range(len(x_labels)))
    ax.set_xticklabels(x_labels, rotation=45, ha='right')

    # Add horizontal line at y=0
    ax.axhline(y=0, color='black', linestyle='-', linewidth=1.5)

    # Clean up appearance
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    plt.savefig(svg_path, format='svg', dpi=300, bbox_inches='tight')
    plt.savefig(png_path, format='png', dpi=300, bbox_inches='tight')
    plt.close()

    logger.info(f"Saved correlation logo to {svg_path} and {png_path}")


def main():
    """Main execution function."""
    parser = argparse.ArgumentParser(description='Compute sequence-metric correlations')
    parser.add_argument('--config', required=True, help='Configuration JSON file')

    args = parser.parse_args()
    logger = setup_logging()

    try:
        # Load configuration
        with open(args.config, 'r') as f:
            config = json.load(f)

        logger.info("Starting correlation analysis...")

        # Load reference sequence
        reference_seq = load_reference_sequence(config['original_sequence'])
        logger.info(f"Reference sequence length: {len(reference_seq)}")

        # Load and merge all data
        merged_df = load_and_merge_data(
            config['mutants_paths'],
            config['data_paths'],
            config['metric'],
            logger
        )

        # Compute 1D correlations
        correlation_1d_df = compute_correlation_1d(merged_df, reference_seq, config['metric'], logger)

        # Compute 2D correlations
        correlation_2d_df = compute_correlation_2d(merged_df, reference_seq, config['metric'], logger)

        # Save results
        logger.info(f"Saving 1D correlations to: {config['correlation_1d_output']}")
        correlation_1d_df.to_csv(config['correlation_1d_output'], index=False)

        logger.info(f"Saving 2D correlations to: {config['correlation_2d_output']}")
        correlation_2d_df.to_csv(config['correlation_2d_output'], index=False)

        # Create correlation logo
        # Parse positions filter if provided
        positions_filter = parse_positions_selection(config.get('positions'))
        if positions_filter:
            logger.info(f"Using positions filter: {positions_filter} ({len(positions_filter)} positions)")

        create_correlation_logo(
            correlation_2d_df,
            correlation_1d_df,
            config['logo_svg_output'],
            config['logo_png_output'],
            logger,
            positions_filter=positions_filter
        )

        # Print summary
        logger.info("Analysis complete!")
        logger.info(f"  Positions analyzed: {len(correlation_1d_df)}")
        logger.info(f"  Mean absolute correlation (1D): {correlation_1d_df['correlation'].abs().mean():.3f}")
        logger.info(f"  Max absolute correlation (1D): {correlation_1d_df['correlation'].abs().max():.3f}")

    except Exception as e:
        logger.error(f"Error during correlation analysis: {str(e)}")
        import traceback
        logger.error(traceback.format_exc())
        sys.exit(1)


if __name__ == "__main__":
    main()
