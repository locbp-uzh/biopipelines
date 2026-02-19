#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Mutation profiler runtime script for BioPipelines.

Analyzes mutation patterns between original and mutant sequences, calculates
frequency statistics, and generates sequence logo visualizations.
"""

import argparse
import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib import colors
import seaborn as sns
from pathlib import Path
import logging
from typing import Dict, List, Tuple, Optional
import sys
import logomaker

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

# Color scheme for amino acids (based on chemical properties)
AMINO_ACID_COLORS = {
    # Hydrophobic (blue shades)
    'A': '#1f77b4', 'V': '#aec7e8', 'I': '#1f77b4', 'L': '#aec7e8', 
    'M': '#1f77b4', 'F': '#aec7e8', 'W': '#1f77b4', 'P': '#aec7e8',
    # Polar (green shades)
    'S': '#2ca02c', 'T': '#98df8a', 'Y': '#2ca02c', 'N': '#98df8a', 'Q': '#2ca02c',
    # Charged positive (red shades)
    'K': '#d62728', 'R': '#ff7f7f', 'H': '#d62728',
    # Charged negative (orange shades)
    'D': '#ff7f0e', 'E': '#ffbb78',
    # Special cases
    'C': '#9467bd', 'G': '#c5b0d5'  # Cysteine (purple), Glycine (light purple)
}


def setup_logging():
    """Set up logging configuration."""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    return logging.getLogger(__name__)


def load_sequences_from_csv(csv_path: str) -> List[Tuple[str, str]]:
    """
    Load sequences from CSV file.
    
    Args:
        csv_path: Path to CSV file with id and sequence columns
        
    Returns:
        List of (id, sequence) tuples
    """
    df = pd.read_csv(csv_path)
    
    # Handle different possible column names
    id_col = None
    seq_col = None
    
    for col in df.columns:
        col_lower = col.lower()
        if col_lower in ['id', 'name', 'identifier']:
            id_col = col
        elif col_lower in ['sequence', 'seq']:
            seq_col = col
    
    if id_col is None or seq_col is None:
        raise ValueError(f"Could not find id and sequence columns in {csv_path}. "
                        f"Available columns: {list(df.columns)}")
    
    sequences = [(row[id_col], row[seq_col]) for _, row in df.iterrows()]
    return sequences


def align_sequences(sequences: List[str]) -> List[str]:
    """
    Align sequences to same length by padding with gaps.
    For now, assumes sequences are already aligned or same length.
    
    Args:
        sequences: List of protein sequences
        
    Returns:
        List of aligned sequences
    """
    max_length = max(len(seq) for seq in sequences)
    aligned = []
    
    for seq in sequences:
        if len(seq) < max_length:
            # Pad with gaps at the end
            aligned_seq = seq + '-' * (max_length - len(seq))
        else:
            aligned_seq = seq
        aligned.append(aligned_seq)
    
    return aligned


def analyze_mutations(original_seqs: List[str], mutant_seqs: List[str], 
                     include_original: bool = True) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Analyze mutations between original and mutant sequences.
    
    Args:
        original_seqs: List of original sequences
        mutant_seqs: List of mutant sequences
        include_original: Whether to include original sequences in analysis
        
    Returns:
        Tuple of (profile_df, mutations_df, absolute_freq_df, relative_freq_df)
    """
    # Get reference sequence (first original sequence)
    if not original_seqs:
        raise ValueError("No original sequences provided")
    
    reference_seq = original_seqs[0]
    seq_length = len(reference_seq)
    
    # Combine sequences for analysis
    all_sequences = mutant_seqs.copy()
    if include_original:
        all_sequences.extend(original_seqs)
    
    total_sequences = len(all_sequences)
    
    # Initialize data structures
    profile_data = []
    mutations_data = []
    absolute_freq_data = []
    relative_freq_data = []
    
    # Analyze each position
    for pos in range(seq_length):
        position_label = pos + 1  # 1-indexed position
        reference_aa = reference_seq[pos]
        
        # Count amino acids at this position
        aa_counts = {aa: 0 for aa in AMINO_ACIDS}
        aa_counts['-'] = 0  # Gap character
        
        mutations_count = 0
        
        for seq in all_sequences:
            if pos < len(seq):
                current_aa = seq[pos]
                if current_aa != reference_aa:
                    mutations_count += 1
                
                if current_aa in aa_counts:
                    aa_counts[current_aa] += 1
                elif current_aa != reference_aa:  # Only count non-reference mutations
                    # Handle unknown amino acids as mutations
                    mutations_count += 1
        
        # Calculate frequencies
        mutations_frequency = mutations_count / total_sequences if total_sequences > 0 else 0
        
        # Profile data - simplified structure
        profile_data.append({
            'position': position_label,
            'original': reference_aa,
            'count': mutations_count,
            'frequency': mutations_frequency
        })
        
        # Raw mutations data (raw counts)
        mutations_row = {'position': position_label, 'original': reference_aa}
        for aa in AMINO_ACIDS:
            mutations_row[aa] = aa_counts[aa]
        mutations_data.append(mutations_row)
        
        # Absolute frequencies (relative to total sequences)
        abs_freq_row = {'position': position_label, 'original': reference_aa}
        for aa in AMINO_ACIDS:
            abs_freq_row[aa] = aa_counts[aa] / total_sequences if total_sequences > 0 else 0
        absolute_freq_data.append(abs_freq_row)
        
        # Relative frequencies (relative to mutations at this position)
        rel_freq_row = {'position': position_label, 'original': reference_aa}
        for aa in AMINO_ACIDS:
            if mutations_count > 0:
                # Only count mutations (not the reference amino acid)
                if aa != reference_aa:
                    rel_freq_row[aa] = aa_counts[aa] / mutations_count
                else:
                    rel_freq_row[aa] = 0.0
            else:
                rel_freq_row[aa] = 0.0
        relative_freq_data.append(rel_freq_row)
    
    # Create DataFrames
    profile_df = pd.DataFrame(profile_data)
    mutations_df = pd.DataFrame(mutations_data)
    absolute_freq_df = pd.DataFrame(absolute_freq_data)
    relative_freq_df = pd.DataFrame(relative_freq_data)
    
    return profile_df, mutations_df, absolute_freq_df, relative_freq_df


def create_sequence_logo(relative_freq_df: pd.DataFrame, absolute_freq_df: pd.DataFrame, profile_df: pd.DataFrame,
                        relative_svg: str, relative_png: str, absolute_svg: str, absolute_png: str,
                        counts_svg: str, counts_png: str, positions_filter: Optional[List[int]] = None):
    """
    Create sequence logo visualizations showing mutation patterns.

    Args:
        relative_freq_df: DataFrame with relative frequencies per position
        absolute_freq_df: DataFrame with absolute frequencies per position
        profile_df: DataFrame with mutation counts per position
        relative_svg: Output SVG file path for relative frequencies
        relative_png: Output PNG file path for relative frequencies
        absolute_svg: Output SVG file path for absolute frequencies
        absolute_png: Output PNG file path for absolute frequencies
        counts_svg: Output SVG file path for counts
        counts_png: Output PNG file path for counts
        positions_filter: List of positions to include (1-indexed). If None, includes only mutated positions (>1% frequency).
    """
    # Filter positions based on positions_filter or mutation threshold
    mutation_positions = []
    for i, row in relative_freq_df.iterrows():
        position = int(row['position'])

        # Apply position filter if specified
        if positions_filter is not None:
            if position not in positions_filter:
                continue
            # Include all specified positions, even if no mutations
        else:
            # Default: only show positions with >1% mutation frequency
            total_freq = sum(row[aa] for aa in AMINO_ACIDS)
            if total_freq <= 0.01:
                continue

        abs_row = absolute_freq_df.iloc[i]
        count_row = profile_df.iloc[i]
        mutation_positions.append((row, abs_row, count_row))
    
    if not mutation_positions:
        # Create empty plots for all three types
        for svg_path, png_path, title in [(relative_svg, relative_png, 'No significant mutations found (Relative)'),
                                          (absolute_svg, absolute_png, 'No significant mutations found (Absolute)'), 
                                          (counts_svg, counts_png, 'No significant mutations found (Counts)')]:
            fig, ax = plt.subplots(figsize=(10, 4))
            ax.text(0.5, 0.5, title, ha='center', va='center', transform=ax.transAxes, fontsize=16)
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1)
            ax.axis('off')
            plt.tight_layout()
            plt.savefig(svg_path, format='svg', dpi=300, bbox_inches='tight')
            plt.savefig(png_path, format='png', dpi=300, bbox_inches='tight')
            plt.close()
        return
    
    # Create three separate plots for different visualizations
    num_positions = len(mutation_positions)
    
    # Plot 1: Relative frequencies
    fig, ax = plt.subplots(figsize=(max(12, num_positions * 0.8), 6))
    _create_logo_plot(ax, mutation_positions, 'relative', 'Relative Mutation Frequency', 'Relative Frequency Sequence Logo')
    plt.tight_layout()
    plt.savefig(relative_svg, format='svg', dpi=300, bbox_inches='tight')
    plt.savefig(relative_png, format='png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Plot 2: Absolute frequencies  
    fig, ax = plt.subplots(figsize=(max(12, num_positions * 0.8), 6))
    _create_logo_plot(ax, mutation_positions, 'absolute', 'Absolute Mutation Frequency', 'Absolute Frequency Sequence Logo')
    plt.tight_layout()
    plt.savefig(absolute_svg, format='svg', dpi=300, bbox_inches='tight')
    plt.savefig(absolute_png, format='png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Plot 3: Counts
    fig, ax = plt.subplots(figsize=(max(12, num_positions * 0.8), 6))
    _create_logo_plot(ax, mutation_positions, 'count', 'Mutation Count', 'Mutation Count Sequence Logo')
    plt.tight_layout()
    plt.savefig(counts_svg, format='svg', dpi=300, bbox_inches='tight')
    plt.savefig(counts_png, format='png', dpi=300, bbox_inches='tight')
    plt.close()


def _create_logo_plot(ax, mutation_positions, data_type, ylabel, title):
    """Create a single sequence logo plot using logomaker."""
    if not mutation_positions:
        ax.text(0.5, 0.5, 'No significant mutations found', 
                ha='center', va='center', transform=ax.transAxes, fontsize=16)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        return
    
    # Prepare data for logomaker
    logo_data = []
    x_labels = []
    
    for rel_data, abs_data, count_data in mutation_positions:
        position = rel_data['position']
        original_aa = rel_data.get('original', abs_data.get('original', ''))
        
        # Choose data source and exclude original amino acid
        row_data = {'pos': position}
        
        if data_type == 'relative':
            # For relative: only mutations (exclude original)
            for aa in AMINO_ACIDS:
                if aa in rel_data and aa != original_aa:
                    row_data[aa] = rel_data[aa]
                else:
                    row_data[aa] = 0.0
        elif data_type == 'absolute':
            # For absolute: mutation frequencies (exclude original)  
            for aa in AMINO_ACIDS:
                if aa in abs_data and aa != original_aa:
                    row_data[aa] = abs_data[aa]
                else:
                    row_data[aa] = 0.0
        else:  # count
            # For count: raw mutation counts (exclude original)
            # Calculate total sequences from profile data
            if count_data['frequency'] > 0:
                total_seqs = count_data['count'] / count_data['frequency']
            else:
                total_seqs = 1
            for aa in AMINO_ACIDS:
                if aa in abs_data and aa != original_aa:
                    row_data[aa] = abs_data[aa] * total_seqs
                else:
                    row_data[aa] = 0.0
        
        logo_data.append(row_data)
        x_labels.append(f"{original_aa}{position}")
    
    # Convert to DataFrame for logomaker with sequential indexing
    logo_df = pd.DataFrame(logo_data)
    # Use sequential indexing (0, 1, 2, ...) instead of actual positions
    logo_df.index = range(len(logo_df))
    logo_df = logo_df.drop('pos', axis=1)  # Remove the position column
    
    # Create logo using logomaker
    logo = logomaker.Logo(logo_df, ax=ax, color_scheme='chemistry')
    
    # Customize the plot
    ax.set_ylabel(ylabel, fontsize=14)
    ax.set_title(title, fontsize=16, fontweight='bold')
    ax.set_xlabel('Position', fontsize=14)
    
    # Set x-ticks to sequential positions but label with original AA + actual position
    ax.set_xticks(range(len(x_labels)))
    ax.set_xticklabels(x_labels, rotation=45, ha='right')
    
    # Clean up appearance
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


def main():
    """Main execution function."""
    parser = argparse.ArgumentParser(description='Analyze mutation patterns between sequences')
    parser.add_argument('--config', required=True, help='Configuration JSON file')
    
    args = parser.parse_args()
    logger = setup_logging()
    
    try:
        # Load configuration
        with open(args.config, 'r') as f:
            config = json.load(f)
        
        logger.info(f"Loading original sequences from: {config['original_sequences']}")
        logger.info(f"Loading mutant sequences from: {config['mutants_sequences']}")
        
        # Load sequences
        original_seqs_data = load_sequences_from_csv(config['original_sequences'])
        mutant_seqs_data = load_sequences_from_csv(config['mutants_sequences'])
        
        # Extract just the sequences
        original_seqs = [seq for _, seq in original_seqs_data]
        mutant_seqs = [seq for _, seq in mutant_seqs_data]
        
        logger.info(f"Loaded {len(original_seqs)} original sequences and {len(mutant_seqs)} mutant sequences")
        
        # Align sequences (ensure same length)
        all_seqs = original_seqs + mutant_seqs
        aligned_seqs = align_sequences(all_seqs)
        aligned_original = aligned_seqs[:len(original_seqs)]
        aligned_mutants = aligned_seqs[len(original_seqs):]
        
        logger.info("Analyzing mutation patterns...")
        
        # Analyze mutations
        profile_df, mutations_df, absolute_freq_df, relative_freq_df = analyze_mutations(
            aligned_original, aligned_mutants, config['include_original']
        )
        
        logger.info(f"Found mutations at {len(profile_df)} positions")
        
        # Save results
        logger.info(f"Saving profile data to: {config['profile_output']}")
        profile_df.to_csv(config['profile_output'], index=False)
        
        logger.info(f"Saving mutations data to: {config['mutations_output']}")
        mutations_df.to_csv(config['mutations_output'], index=False)
        
        logger.info(f"Saving absolute frequencies to: {config['absolute_frequencies_output']}")
        absolute_freq_df.to_csv(config['absolute_frequencies_output'], index=False)
        
        logger.info(f"Saving relative frequencies to: {config['relative_frequencies_output']}")
        relative_freq_df.to_csv(config['relative_frequencies_output'], index=False)
        
        # Create sequence logo
        logger.info("Creating sequence logo visualizations...")
        # Parse positions filter if provided
        positions_filter = parse_positions_selection(config.get('positions'))
        if positions_filter:
            logger.info(f"Using positions filter: {positions_filter} ({len(positions_filter)} positions)")

        create_sequence_logo(relative_freq_df, absolute_freq_df, profile_df,
                           config['sequence_logo_relative_svg'], config['sequence_logo_relative_png'],
                           config['sequence_logo_absolute_svg'], config['sequence_logo_absolute_png'],
                           config['sequence_logo_counts_svg'], config['sequence_logo_counts_png'],
                           positions_filter=positions_filter)
        logger.info(f"Sequence logos saved to:")
        logger.info(f"  Relative frequencies: {config['sequence_logo_relative_svg']} and {config['sequence_logo_relative_png']}")
        logger.info(f"  Absolute frequencies: {config['sequence_logo_absolute_svg']} and {config['sequence_logo_absolute_png']}")
        logger.info(f"  Counts: {config['sequence_logo_counts_svg']} and {config['sequence_logo_counts_png']}")
        
        # Print summary statistics
        total_mutations = profile_df['count'].sum()
        avg_mutations_per_position = profile_df['count'].mean()
        positions_with_mutations = (profile_df['count'] > 0).sum()
        
        logger.info(f"Analysis complete!")
        logger.info(f"  Total mutations: {total_mutations}")
        logger.info(f"  Average mutations per position: {avg_mutations_per_position:.2f}")
        logger.info(f"  Positions with mutations: {positions_with_mutations}")
        
    except Exception as e:
        logger.error(f"Error during mutation analysis: {str(e)}")
        import traceback
        logger.error(traceback.format_exc())
        sys.exit(1)


if __name__ == "__main__":
    main()