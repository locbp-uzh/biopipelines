#!/usr/bin/env python3
"""
Mutation composer runtime script for BioPipelines.

Generates new sequences based on mutation frequency profiles using various
composition strategies.
"""

import argparse
import json
import pandas as pd
import numpy as np
from pathlib import Path
import logging
from typing import Dict, List, Tuple, Optional
import sys
import random

# Standard amino acids
AMINO_ACIDS = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']


def setup_logging():
    """Set up logging configuration."""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    return logging.getLogger(__name__)


def load_frequency_data(csv_path: str) -> pd.DataFrame:
    """
    Load mutation frequency data from CSV file.

    Args:
        csv_path: Path to CSV file with frequency data

    Returns:
        DataFrame with frequency data
    """
    df = pd.read_csv(csv_path)

    # Validate that we have expected columns
    if 'position' not in df.columns:
        raise ValueError(f"Missing 'position' column in {csv_path}")

    # Check for amino acid columns
    missing_aas = [aa for aa in AMINO_ACIDS if aa not in df.columns]
    if missing_aas:
        logging.warning(f"Missing amino acid columns: {missing_aas}")

    return df


def load_and_combine_frequency_data(csv_paths: List[str], strategy: str = "average") -> pd.DataFrame:
    """
    Load and combine multiple frequency datasheets using specified strategy.

    Args:
        csv_paths: List of paths to CSV files with frequency data
        strategy: Combination strategy ("average", "maximum", "stack", "round_robin")

    Returns:
        Combined DataFrame with frequency data
    """
    logger = logging.getLogger(__name__)

    if len(csv_paths) == 1:
        # Single datasheet - just load it
        return load_frequency_data(csv_paths[0])

    logger.info(f"Loading {len(csv_paths)} frequency datasheets for combination")

    # Load all datasheets
    dfs = []
    for i, csv_path in enumerate(csv_paths):
        df = load_frequency_data(csv_path)
        logger.info(f"Loaded datasheet {i+1}: {csv_path} with {len(df)} positions")
        dfs.append(df)

    # Validate compatibility (same positions and reference sequence)
    base_df = dfs[0]
    base_positions = set(base_df['position'].tolist())

    for i, df in enumerate(dfs[1:], 1):
        positions = set(df['position'].tolist())
        if positions != base_positions:
            raise ValueError(f"Position mismatch between datasheet 0 and {i}: {base_positions.symmetric_difference(positions)}")

        # Check reference sequence compatibility if 'original' column exists
        if 'original' in base_df.columns and 'original' in df.columns:
            base_orig = base_df.set_index('position')['original']
            df_orig = df.set_index('position')['original']
            mismatches = (base_orig != df_orig).sum()
            if mismatches > 0:
                raise ValueError(f"Reference sequence mismatch between datasheet 0 and {i}: {mismatches} positions differ")

    logger.info("All datasheets are compatible - proceeding with combination")

    if strategy in ["average", "maximum"]:
        # Combine frequencies by averaging or taking maximum
        combined_df = base_df.copy()

        for aa in AMINO_ACIDS:
            if aa not in combined_df.columns:
                continue

            # Collect values for this amino acid across all datasheets
            aa_values = []
            for df in dfs:
                if aa in df.columns:
                    aa_values.append(df[aa].values)
                else:
                    aa_values.append(np.zeros(len(df)))

            # Combine based on strategy
            aa_array = np.array(aa_values)  # Shape: (n_datasheets, n_positions)

            if strategy == "average":
                combined_df[aa] = np.mean(aa_array, axis=0)
            elif strategy == "maximum":
                combined_df[aa] = np.max(aa_array, axis=0)

        logger.info(f"Combined frequencies using {strategy} strategy")
        return combined_df

    else:
        # For "stack" and "round_robin", we'll handle combination during sequence generation
        # Return the list of datasheets as a special combined format
        combined_df = base_df.copy()
        combined_df.attrs['strategy'] = strategy
        combined_df.attrs['source_dfs'] = dfs
        logger.info(f"Prepared datasheets for {strategy} strategy")
        return combined_df


def convert_absolute_to_relative_frequencies(freq_df: pd.DataFrame) -> pd.DataFrame:
    """
    Convert absolute frequencies to relative frequencies for hotspot-focused mode.

    For each position, converts absolute frequencies (mutation_count/total_sequences)
    to relative frequencies (mutation_count/mutations_at_position).

    Args:
        freq_df: DataFrame with absolute frequencies

    Returns:
        DataFrame with relative frequencies (mutations only, original AA = 0)
    """
    logger = logging.getLogger(__name__)
    logger.info("Converting absolute frequencies to relative frequencies for hotspot mode")

    relative_df = freq_df.copy()

    for idx, row in freq_df.iterrows():
        # Get original amino acid
        if 'original' in row:
            original_aa = row['original']
        else:
            # Fallback: find AA with highest frequency (likely original)
            original_aa = None
            max_freq = 0
            for aa in AMINO_ACIDS:
                if aa in row and row[aa] > max_freq:
                    max_freq = row[aa]
                    original_aa = aa

        if original_aa is None:
            continue

        # Calculate total mutations at this position (exclude original)
        mutation_count = sum(row[aa] for aa in AMINO_ACIDS if aa in row and aa != original_aa)

        # Convert to relative frequencies
        for aa in AMINO_ACIDS:
            if aa in row:
                if aa == original_aa:
                    # Set original AA to 0 for relative frequencies
                    relative_df.at[idx, aa] = 0.0
                elif mutation_count > 0:
                    # Convert to relative frequency
                    relative_df.at[idx, aa] = row[aa] / mutation_count
                else:
                    # No mutations at this position
                    relative_df.at[idx, aa] = 0.0

    return relative_df


def get_reference_sequence(freq_df: pd.DataFrame) -> str:
    """
    Reconstruct reference sequence from frequency data using the 'original' column.
    
    Args:
        freq_df: DataFrame with frequency data including 'original' column
        
    Returns:
        Reference sequence string
    """
    reference_seq = []
    
    for _, row in freq_df.iterrows():
        # Use the 'original' column if available
        if 'original' in row:
            reference_aa = row['original']
        else:
            # Fallback: find amino acid with zero frequency (original) or most common if all have mutations
            min_freq = float('inf')
            reference_aa = 'X'  # Default unknown
            
            for aa in AMINO_ACIDS:
                if aa in row:
                    freq = row[aa]
                    if freq == 0:
                        reference_aa = aa
                        break
                    elif freq < min_freq:
                        min_freq = freq
                        reference_aa = aa
        
        reference_seq.append(reference_aa)
    
    return ''.join(reference_seq)


def generate_single_point_mutations(freq_df: pd.DataFrame, num_sequences: int,
                                   min_frequency: float, prefix: str = "", combination_strategy: str = "average") -> List[Tuple[str, str, str, List[int]]]:
    """
    Generate single-point mutations ordered by frequency.
    Each sequence contains exactly one mutation from the original sequence.

    Args:
        freq_df: DataFrame with frequency data (possibly combined)
        num_sequences: Number of sequences to generate
        min_frequency: Minimum frequency threshold
        prefix: Prefix for sequence IDs
        combination_strategy: Strategy used for multi-datasheet combination

    Returns:
        List of (sequence_id, sequence, mutations, mutation_positions) tuples
    """
    logger = logging.getLogger(__name__)

    # Handle special combination strategies
    if combination_strategy in ["stack", "round_robin"] and hasattr(freq_df, 'attrs') and 'source_dfs' in freq_df.attrs:
        return generate_multi_datasheet_single_point_mutations(
            freq_df.attrs['source_dfs'], num_sequences, min_frequency, prefix, combination_strategy
        )

    # Standard single-datasheet or pre-combined datasheet processing
    reference_seq = get_reference_sequence(freq_df)

    # Collect all mutations above threshold
    mutations = []
    for _, row in freq_df.iterrows():
        position = int(row['position']) - 1  # Convert to 0-indexed
        
        # Get original amino acid from the 'original' column or reference sequence
        if 'original' in row:
            original_aa = row['original']
        else:
            original_aa = reference_seq[position]
        
        # Find all mutations at this position
        for aa in AMINO_ACIDS:
            if aa in row and aa != original_aa:
                freq = row[aa]
                if freq >= min_frequency:
                    mutations.append((position, original_aa, aa, freq))
    
    # Sort by frequency (highest first)
    mutations.sort(key=lambda x: x[3], reverse=True)
    
    # Generate sequences - each with exactly one mutation from the original
    sequences = []
    for i in range(min(num_sequences, len(mutations))):
        pos, orig_aa, mut_aa, freq = mutations[i]
        
        # Start with the original reference sequence
        seq_list = list(reference_seq)
        # Apply only this single mutation
        seq_list[pos] = mut_aa
        mutated_seq = ''.join(seq_list)
        
        sequence_id = f"{prefix}_{i+1:03d}" if prefix else f"mut_{i+1:03d}"
        mutation_desc = f"{orig_aa}{pos+1}{mut_aa}"
        mutation_positions = [pos + 1]  # 1-indexed for output
        
        sequences.append((sequence_id, mutated_seq, mutation_desc, mutation_positions))
    
    return sequences


def generate_multi_datasheet_single_point_mutations(source_dfs: List[pd.DataFrame], num_sequences: int,
                                                   min_frequency: float, prefix: str, strategy: str) -> List[Tuple[str, str, str, List[int]]]:
    """
    Generate single-point mutations using multiple datasheets with stack or round_robin strategies.

    Args:
        source_dfs: List of source DataFrames with frequency data
        num_sequences: Number of sequences to generate
        min_frequency: Minimum frequency threshold
        prefix: Prefix for sequence IDs
        strategy: Either "stack" or "round_robin"

    Returns:
        List of (sequence_id, sequence, mutations, mutation_positions) tuples
    """
    logger = logging.getLogger(__name__)
    logger.info(f"Generating multi-datasheet single-point mutations using {strategy} strategy")

    # Get reference sequence (should be same across all datasheets)
    reference_seq = get_reference_sequence(source_dfs[0])

    if strategy == "stack":
        # Each sequence gets one mutation from each datasheet
        logger.info(f"Stack strategy: each sequence will have {len(source_dfs)} mutations")

        # Collect top mutations from each datasheet
        datasheet_mutations = []
        for df_idx, freq_df in enumerate(source_dfs):
            mutations = []
            for _, row in freq_df.iterrows():
                position = int(row['position']) - 1  # Convert to 0-indexed

                # Get original amino acid
                if 'original' in row:
                    original_aa = row['original']
                else:
                    original_aa = reference_seq[position]

                # Find all mutations at this position
                for aa in AMINO_ACIDS:
                    if aa in row and aa != original_aa:
                        freq = row[aa]
                        if freq >= min_frequency:
                            mutations.append((position, original_aa, aa, freq))

            # Sort by frequency (highest first)
            mutations.sort(key=lambda x: x[3], reverse=True)
            datasheet_mutations.append(mutations)
            logger.info(f"Datasheet {df_idx+1}: found {len(mutations)} valid mutations")

        # Generate sequences with one mutation from each datasheet
        sequences = []
        for seq_idx in range(num_sequences):
            seq_list = list(reference_seq)
            mutations_made = []
            mutation_positions = []
            used_positions = set()  # Track positions already mutated

            # Take mutations from each datasheet, handling conflicts
            for df_idx, mutations in enumerate(datasheet_mutations):
                selected_mutation = None

                # Try to find a non-conflicting mutation from this datasheet
                for mutation_idx in range(seq_idx, len(mutations)):
                    pos, orig_aa, mut_aa, freq = mutations[mutation_idx]

                    if pos not in used_positions:
                        # No conflict - use this mutation
                        selected_mutation = (pos, orig_aa, mut_aa, freq)
                        break
                    else:
                        # Conflict detected - log and try next mutation
                        existing_mut = None
                        for made_mut in mutations_made:
                            if f"{pos+1}" in made_mut:
                                existing_mut = made_mut
                                break
                        logger.info(f"Conflict at position {pos+1}: datasheet {df_idx+1} wants {orig_aa}->{mut_aa}, "
                                  f"but position already has {existing_mut}. Trying next best...")

                # Apply the selected mutation if found
                if selected_mutation:
                    pos, orig_aa, mut_aa, freq = selected_mutation
                    seq_list[pos] = mut_aa
                    mutations_made.append(f"{orig_aa}{pos+1}{mut_aa}")
                    mutation_positions.append(pos + 1)
                    used_positions.add(pos)
                    logger.debug(f"Applied mutation from datasheet {df_idx+1}: {orig_aa}{pos+1}{mut_aa}")
                else:
                    logger.warning(f"No non-conflicting mutation found for datasheet {df_idx+1} in sequence {seq_idx+1}")

            # Create sequence only if we have mutations
            if mutations_made:
                mutated_seq = ''.join(seq_list)
                sequence_id = f"{prefix}_{seq_idx+1:03d}" if prefix else f"stack_{seq_idx+1:03d}"
                mutations_desc = ",".join(mutations_made)
                sequences.append((sequence_id, mutated_seq, mutations_desc, mutation_positions))
                logger.info(f"Generated sequence {seq_idx+1} with {len(mutations_made)} mutations: {mutations_desc}")
            else:
                logger.warning(f"No mutations could be applied for sequence {seq_idx+1} - all positions conflicted")

        logger.info(f"Generated {len(sequences)} stacked sequences")
        return sequences

    elif strategy == "round_robin":
        # Alternate between datasheets for single mutations
        logger.info("Round-robin strategy: alternating between datasheets for single mutations")

        # Collect all mutations from all datasheets with source info
        all_mutations = []
        for df_idx, freq_df in enumerate(source_dfs):
            for _, row in freq_df.iterrows():
                position = int(row['position']) - 1

                # Get original amino acid
                if 'original' in row:
                    original_aa = row['original']
                else:
                    original_aa = reference_seq[position]

                # Find all mutations at this position
                for aa in AMINO_ACIDS:
                    if aa in row and aa != original_aa:
                        freq = row[aa]
                        if freq >= min_frequency:
                            all_mutations.append((position, original_aa, aa, freq, df_idx))

        # Sort by frequency (highest first)
        all_mutations.sort(key=lambda x: x[3], reverse=True)

        # Group by datasheet and maintain order
        datasheet_mutations = [[] for _ in source_dfs]
        for mutation in all_mutations:
            pos, orig_aa, mut_aa, freq, df_idx = mutation
            datasheet_mutations[df_idx].append((pos, orig_aa, mut_aa, freq))

        # Generate sequences by round-robin through datasheets
        sequences = []
        for seq_idx in range(num_sequences):
            df_idx = seq_idx % len(source_dfs)  # Round-robin datasheet selection
            mutation_idx = seq_idx // len(source_dfs)  # Which mutation within that datasheet

            if mutation_idx < len(datasheet_mutations[df_idx]):
                pos, orig_aa, mut_aa, freq = datasheet_mutations[df_idx][mutation_idx]

                # Apply single mutation
                seq_list = list(reference_seq)
                seq_list[pos] = mut_aa
                mutated_seq = ''.join(seq_list)

                sequence_id = f"{prefix}_{seq_idx+1:03d}" if prefix else f"rr_{seq_idx+1:03d}"
                mutations_desc = f"{orig_aa}{pos+1}{mut_aa}"
                mutation_positions = [pos + 1]

                sequences.append((sequence_id, mutated_seq, mutations_desc, mutation_positions))

        logger.info(f"Generated {len(sequences)} round-robin sequences")
        return sequences

    else:
        raise ValueError(f"Unknown multi-datasheet strategy: {strategy}")


def generate_multi_datasheet_top_mutations(source_dfs: List[pd.DataFrame], num_sequences: int,
                                         min_frequency: float, prefix: str, strategy: str) -> List[Tuple[str, str, str, List[int]]]:
    """
    Generate sequences using top mutations from multiple datasheets with stack or round_robin strategies.

    Args:
        source_dfs: List of source DataFrames with frequency data
        num_sequences: Number of sequences to generate
        min_frequency: Minimum frequency threshold
        prefix: Prefix for sequence IDs
        strategy: Either "stack" or "round_robin"

    Returns:
        List of (sequence_id, sequence, mutations, mutation_positions) tuples
    """
    logger = logging.getLogger(__name__)
    logger.info(f"Generating multi-datasheet top mutations using {strategy} strategy")

    # Get reference sequence (should be same across all datasheets)
    reference_seq = get_reference_sequence(source_dfs[0])

    if strategy == "stack":
        # Each sequence gets top mutations from each datasheet stacked together
        logger.info(f"Stack strategy: each sequence will combine top mutations from {len(source_dfs)} datasheets")

        # For each datasheet, get top mutations per position
        datasheet_top_mutations = []
        for df_idx, freq_df in enumerate(source_dfs):
            position_mutations = {}
            for _, row in freq_df.iterrows():
                position = int(row['position']) - 1
                if 'original' in row:
                    original_aa = row['original']
                else:
                    original_aa = reference_seq[position]

                aa_freqs = []
                for aa in AMINO_ACIDS:
                    if aa in row and aa != original_aa:
                        freq = row[aa]
                        if freq >= min_frequency:
                            aa_freqs.append((aa, freq))

                # Sort by frequency and take top 3
                aa_freqs.sort(key=lambda x: x[1], reverse=True)
                position_mutations[position] = aa_freqs[:3]

            datasheet_top_mutations.append(position_mutations)
            logger.info(f"Datasheet {df_idx+1}: found top mutations for {len(position_mutations)} positions")

        # Generate sequences by stacking mutations from all datasheets
        sequences = []
        for seq_idx in range(num_sequences):
            seq_list = list(reference_seq)
            mutations_made = []
            mutation_positions = []
            used_positions = set()

            # Apply top mutations from each datasheet
            for df_idx, position_mutations in enumerate(datasheet_top_mutations):
                for pos, mut_list in position_mutations.items():
                    if pos not in used_positions and mut_list:
                        # Use the top mutation from this datasheet for this position
                        mut_idx = seq_idx % len(mut_list)
                        selected_aa, freq = mut_list[mut_idx]

                        seq_list[pos] = selected_aa
                        original_aa = reference_seq[pos]
                        mutations_made.append(f"{original_aa}{pos+1}{selected_aa}")
                        mutation_positions.append(pos + 1)
                        used_positions.add(pos)
                        break  # Only one mutation per datasheet per sequence

            if mutations_made:
                mutated_seq = ''.join(seq_list)
                sequence_id = f"{prefix}_{seq_idx+1:03d}" if prefix else f"top_stack_{seq_idx+1:03d}"
                mutations_desc = ",".join(mutations_made)
                sequences.append((sequence_id, mutated_seq, mutations_desc, mutation_positions))

        logger.info(f"Generated {len(sequences)} stacked top-mutation sequences")
        return sequences

    elif strategy == "round_robin":
        # Alternate between datasheets, taking top mutations from each
        logger.info("Round-robin strategy: alternating between datasheets for top mutations")

        # Collect top mutations from all datasheets with datasheet info
        all_top_mutations = []
        for df_idx, freq_df in enumerate(source_dfs):
            for _, row in freq_df.iterrows():
                position = int(row['position']) - 1
                if 'original' in row:
                    original_aa = row['original']
                else:
                    original_aa = reference_seq[position]

                aa_freqs = []
                for aa in AMINO_ACIDS:
                    if aa in row and aa != original_aa:
                        freq = row[aa]
                        if freq >= min_frequency:
                            aa_freqs.append((aa, freq))

                # Sort by frequency and take top mutation
                aa_freqs.sort(key=lambda x: x[1], reverse=True)
                if aa_freqs:
                    top_aa, top_freq = aa_freqs[0]
                    all_top_mutations.append((position, original_aa, top_aa, top_freq, df_idx))

        # Sort by frequency and group by datasheet
        all_top_mutations.sort(key=lambda x: x[3], reverse=True)
        datasheet_mutations = [[] for _ in source_dfs]
        for mutation in all_top_mutations:
            pos, orig_aa, mut_aa, freq, df_idx = mutation
            datasheet_mutations[df_idx].append((pos, orig_aa, mut_aa, freq))

        # Generate sequences by round-robin through datasheets
        sequences = []
        for seq_idx in range(num_sequences):
            df_idx = seq_idx % len(source_dfs)
            mutation_idx = seq_idx // len(source_dfs)

            if mutation_idx < len(datasheet_mutations[df_idx]):
                pos, orig_aa, mut_aa, freq = datasheet_mutations[df_idx][mutation_idx]

                # Apply single mutation
                seq_list = list(reference_seq)
                seq_list[pos] = mut_aa
                mutated_seq = ''.join(seq_list)
                sequence_id = f"{prefix}_{seq_idx+1:03d}" if prefix else f"top_rr_{seq_idx+1:03d}"
                mutations_desc = f"{orig_aa}{pos+1}{mut_aa}"
                mutation_positions = [pos + 1]
                sequences.append((sequence_id, mutated_seq, mutations_desc, mutation_positions))

        logger.info(f"Generated {len(sequences)} round-robin top-mutation sequences")
        return sequences

    else:
        raise ValueError(f"Unknown multi-datasheet strategy: {strategy}")


def generate_multi_datasheet_weighted_random(source_dfs: List[pd.DataFrame], num_sequences: int,
                                           min_frequency: float, max_mutations: Optional[int], random_seed: Optional[int],
                                           prefix: str, strategy: str) -> List[Tuple[str, str, str, List[int]]]:
    """
    Generate weighted random sequences using multiple datasheets with stack or round_robin strategies.

    Args:
        source_dfs: List of source DataFrames with frequency data
        num_sequences: Number of sequences to generate
        min_frequency: Minimum frequency threshold
        max_mutations: Maximum mutations per sequence
        random_seed: Random seed for reproducibility
        prefix: Prefix for sequence IDs
        strategy: Either "stack" or "round_robin"

    Returns:
        List of (sequence_id, sequence, mutations, mutation_positions) tuples
    """
    logger = logging.getLogger(__name__)
    logger.info(f"Generating multi-datasheet weighted random mutations using {strategy} strategy")

    if random_seed is not None:
        random.seed(random_seed)
        np.random.seed(random_seed)

    # Get reference sequence (should be same across all datasheets)
    reference_seq = get_reference_sequence(source_dfs[0])

    if strategy == "stack":
        # Each sequence gets weighted random mutations from each datasheet combined
        logger.info(f"Stack strategy: combining weighted mutations from {len(source_dfs)} datasheets")

        sequences = []
        for seq_idx in range(num_sequences):
            seq_list = list(reference_seq)
            mutations_made = []
            mutation_positions = []
            used_positions = set()

            # Apply weighted random mutations from each datasheet
            for df_idx, freq_df in enumerate(source_dfs):
                # Collect weighted mutations from this datasheet
                weighted_mutations = []
                for _, row in freq_df.iterrows():
                    position = int(row['position']) - 1
                    if position in used_positions:
                        continue  # Skip positions already mutated

                    if 'original' in row:
                        original_aa = row['original']
                    else:
                        original_aa = reference_seq[position]

                    for aa in AMINO_ACIDS:
                        if aa in row and aa != original_aa:
                            freq = row[aa]
                            if freq >= min_frequency:
                                weighted_mutations.append((position, original_aa, aa, freq))

                # Apply weighted random selection from this datasheet
                if weighted_mutations:
                    # Calculate weights and select randomly
                    positions, orig_aas, mut_aas, weights = zip(*weighted_mutations)
                    total_weight = sum(weights)

                    if total_weight > 0:
                        # Select one mutation per datasheet
                        probs = [w/total_weight for w in weights]
                        selected_idx = np.random.choice(len(weighted_mutations), p=probs)
                        pos, orig_aa, mut_aa, freq = weighted_mutations[selected_idx]

                        seq_list[pos] = mut_aa
                        mutations_made.append(f"{orig_aa}{pos+1}{mut_aa}")
                        mutation_positions.append(pos + 1)
                        used_positions.add(pos)

            if mutations_made:
                mutated_seq = ''.join(seq_list)
                sequence_id = f"{prefix}_{seq_idx+1:03d}" if prefix else f"wr_stack_{seq_idx+1:03d}"
                mutations_desc = ",".join(mutations_made)
                sequences.append((sequence_id, mutated_seq, mutations_desc, mutation_positions))

        logger.info(f"Generated {len(sequences)} stacked weighted random sequences")
        return sequences

    elif strategy == "round_robin":
        # Alternate between datasheets for weighted random mutations
        logger.info("Round-robin strategy: alternating between datasheets for weighted random mutations")

        # Collect all weighted mutations from all datasheets
        datasheet_mutations = []
        for df_idx, freq_df in enumerate(source_dfs):
            weighted_mutations = []
            for _, row in freq_df.iterrows():
                position = int(row['position']) - 1
                if 'original' in row:
                    original_aa = row['original']
                else:
                    original_aa = reference_seq[position]

                for aa in AMINO_ACIDS:
                    if aa in row and aa != original_aa:
                        freq = row[aa]
                        if freq >= min_frequency:
                            weighted_mutations.append((position, original_aa, aa, freq))

            datasheet_mutations.append(weighted_mutations)
            logger.info(f"Datasheet {df_idx+1}: found {len(weighted_mutations)} weighted mutations")

        # Generate sequences by round-robin through datasheets
        sequences = []
        for seq_idx in range(num_sequences):
            df_idx = seq_idx % len(source_dfs)
            weighted_mutations = datasheet_mutations[df_idx]

            if weighted_mutations:
                # Apply weighted random selection from selected datasheet
                positions, orig_aas, mut_aas, weights = zip(*weighted_mutations)
                total_weight = sum(weights)

                if total_weight > 0:
                    probs = [w/total_weight for w in weights]
                    selected_idx = np.random.choice(len(weighted_mutations), p=probs)
                    pos, orig_aa, mut_aa, freq = weighted_mutations[selected_idx]

                    # Apply single mutation
                    seq_list = list(reference_seq)
                    seq_list[pos] = mut_aa
                    mutated_seq = ''.join(seq_list)
                    sequence_id = f"{prefix}_{seq_idx+1:03d}" if prefix else f"wr_rr_{seq_idx+1:03d}"
                    mutations_desc = f"{orig_aa}{pos+1}{mut_aa}"
                    mutation_positions = [pos + 1]
                    sequences.append((sequence_id, mutated_seq, mutations_desc, mutation_positions))

        logger.info(f"Generated {len(sequences)} round-robin weighted random sequences")
        return sequences

    else:
        raise ValueError(f"Unknown multi-datasheet strategy: {strategy}")


def generate_weighted_random_mutations(freq_df: pd.DataFrame, num_sequences: int,
                                     min_frequency: float, max_mutations: Optional[int],
                                     random_seed: Optional[int], prefix: str = "", combination_strategy: str = "average") -> List[Tuple[str, str, str, List[int]]]:
    """
    Generate random sequences weighted by mutation frequencies.
    
    Args:
        freq_df: DataFrame with frequency data
        num_sequences: Number of sequences to generate
        min_frequency: Minimum frequency threshold
        max_mutations: Maximum mutations per sequence
        random_seed: Random seed for reproducibility
        
    Returns:
        List of (sequence_id, sequence, mutations, mutation_positions) tuples
    """
    logger = logging.getLogger(__name__)

    # Handle multi-datasheet combination strategies
    if combination_strategy in ["stack", "round_robin"] and hasattr(freq_df, 'attrs') and 'source_dfs' in freq_df.attrs:
        return generate_multi_datasheet_weighted_random(
            freq_df.attrs['source_dfs'], num_sequences, min_frequency, max_mutations, random_seed, prefix, combination_strategy
        )

    if random_seed is not None:
        random.seed(random_seed)
        np.random.seed(random_seed)

    reference_seq = get_reference_sequence(freq_df)
    sequences = []
    
    for seq_idx in range(num_sequences):
        seq_list = list(reference_seq)
        mutations_made = []
        mutation_positions = []
        
        # Determine number of mutations for this sequence
        if max_mutations is None:
            max_muts_this_seq = len(reference_seq)
        else:
            max_muts_this_seq = max_mutations

        num_mutations = random.randint(1, min(max_muts_this_seq, 5))  # Cap at 5 for reasonable sequences

        # Step 1: Calculate absolute mutation frequencies per position for position selection
        position_absolute_freqs = []
        position_data = []

        for row_idx, row in freq_df.iterrows():
            pos = int(row['position']) - 1  # Convert to 0-indexed
            if pos >= len(reference_seq):
                continue

            # Get original amino acid
            if 'original' in row:
                original_aa = row['original']
            else:
                original_aa = reference_seq[pos]

            # Calculate total absolute frequency for mutations at this position
            absolute_freq = 0.0
            valid_mutations = []

            for aa in AMINO_ACIDS:
                if aa in row and aa != original_aa:
                    freq = row[aa]
                    if freq >= min_frequency:
                        absolute_freq += freq
                        valid_mutations.append((aa, freq))

            # Only include positions that have valid mutations
            if valid_mutations:
                position_absolute_freqs.append(absolute_freq)
                position_data.append((pos, original_aa, valid_mutations, absolute_freq))

        # If no valid positions found, skip this sequence
        if not position_data:
            continue

        # Step 2: Select positions to mutate based on absolute frequencies
        total_absolute_freq = sum(position_absolute_freqs)
        if total_absolute_freq == 0:
            continue

        # Normalize position probabilities
        position_probs = [freq / total_absolute_freq for freq in position_absolute_freqs]

        # Select positions without replacement using weighted sampling
        selected_positions = []
        remaining_positions = list(range(len(position_data)))
        remaining_probs = position_probs.copy()

        for _ in range(min(num_mutations, len(position_data))):
            if not remaining_positions:
                break

            # Normalize remaining probabilities
            total_remaining = sum(remaining_probs)
            if total_remaining == 0:
                break
            normalized_probs = [p / total_remaining for p in remaining_probs]

            # Select position
            selected_idx = np.random.choice(remaining_positions, p=normalized_probs)
            selected_positions.append(selected_idx)

            # Remove from remaining
            list_idx = remaining_positions.index(selected_idx)
            remaining_positions.pop(list_idx)
            remaining_probs.pop(list_idx)

        # Step 3: For each selected position, choose amino acid based on relative frequencies
        for pos_idx in selected_positions:
            pos, original_aa, valid_mutations, _ = position_data[pos_idx]

            # Calculate relative frequencies for amino acid selection
            total_mut_freq = sum(freq for aa, freq in valid_mutations)
            if total_mut_freq == 0:
                continue

            aa_choices = [aa for aa, freq in valid_mutations]
            aa_probs = [freq / total_mut_freq for aa, freq in valid_mutations]

            # Select amino acid based on relative probability
            selected_aa = np.random.choice(aa_choices, p=aa_probs)

            # Apply mutation
            seq_list[pos] = selected_aa
            mutations_made.append(f"{original_aa}{pos+1}{selected_aa}")
            mutation_positions.append(pos + 1)
        
        # Create final sequence
        mutated_seq = ''.join(seq_list)
        sequence_id = f"{prefix}_{seq_idx+1:03d}" if prefix else f"rand_{seq_idx+1:03d}"
        mutations_desc = ",".join(mutations_made) if mutations_made else "none"
        
        sequences.append((sequence_id, mutated_seq, mutations_desc, mutation_positions))
    
    return sequences


def generate_multi_datasheet_hotspot_focused(source_dfs: List[pd.DataFrame], num_sequences: int,
                                           min_frequency: float, prefix: str, hotspot_count: int,
                                           strategy: str) -> List[Tuple[str, str, str, List[int]]]:
    """
    Generate hotspot-focused sequences using multiple datasheets with stack or round_robin strategies.

    Args:
        source_dfs: List of source DataFrames with frequency data
        num_sequences: Number of sequences to generate
        min_frequency: Minimum frequency threshold
        prefix: Prefix for sequence IDs
        hotspot_count: Number of top hotspot positions to consider
        strategy: Either "stack" or "round_robin"

    Returns:
        List of (sequence_id, sequence, mutations, mutation_positions) tuples
    """
    logger = logging.getLogger(__name__)
    logger.info(f"Generating multi-datasheet hotspot-focused mutations using {strategy} strategy")

    # Get reference sequence (should be same across all datasheets)
    reference_seq = get_reference_sequence(source_dfs[0])

    if strategy == "stack":
        # Identify hotspots from each datasheet and stack mutations from them
        logger.info(f"Stack strategy: combining hotspots from {len(source_dfs)} datasheets")

        # For each datasheet, identify hotspots
        datasheet_hotspots = []
        for df_idx, freq_df in enumerate(source_dfs):
            # Convert to relative frequencies if needed for hotspot detection
            working_df = freq_df.copy()
            is_absolute_freq = False

            # Check if this is absolute frequencies and convert to relative
            for _, row in freq_df.head(20).iterrows():
                if 'original' in row:
                    orig_aa = row['original']
                    if orig_aa in row and row[orig_aa] > 0.8:  # High frequency for original indicates absolute
                        is_absolute_freq = True
                        break

            if is_absolute_freq:
                logger.info(f"Datasheet {df_idx+1}: Converting absolute to relative frequencies for hotspot detection")
                working_df = convert_absolute_to_relative_frequencies(working_df)

            # Calculate mutation rates per position
            position_mutation_rates = {}
            for _, row in working_df.iterrows():
                position = int(row['position']) - 1
                if 'original' in row:
                    original_aa = row['original']
                else:
                    original_aa = reference_seq[position]

                mutation_rate = 0.0
                for aa in AMINO_ACIDS:
                    if aa in row and aa != original_aa:
                        freq = row[aa]
                        if freq >= min_frequency:
                            mutation_rate += freq

                position_mutation_rates[position] = mutation_rate

            # Get top hotspots
            sorted_positions = sorted(position_mutation_rates.items(), key=lambda x: x[1], reverse=True)
            hotspot_positions = [pos for pos, rate in sorted_positions[:hotspot_count] if rate > 0]

            # Collect mutations from hotspot positions
            hotspot_mutations = []
            for _, row in working_df.iterrows():
                position = int(row['position']) - 1
                if position in hotspot_positions:
                    if 'original' in row:
                        original_aa = row['original']
                    else:
                        original_aa = reference_seq[position]

                    for aa in AMINO_ACIDS:
                        if aa in row and aa != original_aa:
                            freq = row[aa]
                            if freq >= min_frequency:
                                hotspot_mutations.append((position, original_aa, aa, freq))

            # Sort by frequency
            hotspot_mutations.sort(key=lambda x: x[3], reverse=True)
            datasheet_hotspots.append(hotspot_mutations)
            logger.info(f"Datasheet {df_idx+1}: found {len(hotspot_mutations)} hotspot mutations from {len(hotspot_positions)} positions")

        # Generate sequences by stacking hotspot mutations
        sequences = []
        for seq_idx in range(num_sequences):
            seq_list = list(reference_seq)
            mutations_made = []
            mutation_positions = []
            used_positions = set()

            # Apply mutations from each datasheet's hotspots
            for df_idx, hotspot_mutations in enumerate(datasheet_hotspots):
                if hotspot_mutations:
                    # Find a mutation that doesn't conflict
                    for mutation_idx in range(seq_idx, len(hotspot_mutations)):
                        pos, orig_aa, mut_aa, freq = hotspot_mutations[mutation_idx % len(hotspot_mutations)]
                        if pos not in used_positions:
                            seq_list[pos] = mut_aa
                            mutations_made.append(f"{orig_aa}{pos+1}{mut_aa}")
                            mutation_positions.append(pos + 1)
                            used_positions.add(pos)
                            break

            if mutations_made:
                mutated_seq = ''.join(seq_list)
                sequence_id = f"{prefix}_{seq_idx+1:03d}" if prefix else f"hs_stack_{seq_idx+1:03d}"
                mutations_desc = ",".join(mutations_made)
                sequences.append((sequence_id, mutated_seq, mutations_desc, mutation_positions))

        logger.info(f"Generated {len(sequences)} stacked hotspot-focused sequences")
        return sequences

    elif strategy == "round_robin":
        # Alternate between datasheets for hotspot mutations
        logger.info("Round-robin strategy: alternating between datasheets for hotspot mutations")

        # Collect hotspot mutations from all datasheets
        datasheet_hotspots = []
        for df_idx, freq_df in enumerate(source_dfs):
            # Convert to relative frequencies if needed
            working_df = freq_df.copy()
            is_absolute_freq = False

            for _, row in freq_df.head(20).iterrows():
                if 'original' in row:
                    orig_aa = row['original']
                    if orig_aa in row and row[orig_aa] > 0.8:
                        is_absolute_freq = True
                        break

            if is_absolute_freq:
                working_df = convert_absolute_to_relative_frequencies(working_df)

            # Calculate mutation rates and identify hotspots
            position_mutation_rates = {}
            for _, row in working_df.iterrows():
                position = int(row['position']) - 1
                if 'original' in row:
                    original_aa = row['original']
                else:
                    original_aa = reference_seq[position]

                mutation_rate = sum(row[aa] for aa in AMINO_ACIDS
                                  if aa in row and aa != original_aa and row[aa] >= min_frequency)
                position_mutation_rates[position] = mutation_rate

            # Get top hotspots and their mutations
            sorted_positions = sorted(position_mutation_rates.items(), key=lambda x: x[1], reverse=True)
            hotspot_positions = [pos for pos, rate in sorted_positions[:hotspot_count] if rate > 0]

            hotspot_mutations = []
            for _, row in working_df.iterrows():
                position = int(row['position']) - 1
                if position in hotspot_positions:
                    if 'original' in row:
                        original_aa = row['original']
                    else:
                        original_aa = reference_seq[position]

                    for aa in AMINO_ACIDS:
                        if aa in row and aa != original_aa:
                            freq = row[aa]
                            if freq >= min_frequency:
                                hotspot_mutations.append((position, original_aa, aa, freq))

            hotspot_mutations.sort(key=lambda x: x[3], reverse=True)
            datasheet_hotspots.append(hotspot_mutations)

        # Generate sequences by round-robin
        sequences = []
        for seq_idx in range(num_sequences):
            df_idx = seq_idx % len(source_dfs)
            hotspot_mutations = datasheet_hotspots[df_idx]
            mutation_idx = seq_idx // len(source_dfs)

            if mutation_idx < len(hotspot_mutations):
                pos, orig_aa, mut_aa, freq = hotspot_mutations[mutation_idx]

                # Apply single mutation
                seq_list = list(reference_seq)
                seq_list[pos] = mut_aa
                mutated_seq = ''.join(seq_list)
                sequence_id = f"{prefix}_{seq_idx+1:03d}" if prefix else f"hs_rr_{seq_idx+1:03d}"
                mutations_desc = f"{orig_aa}{pos+1}{mut_aa}"
                mutation_positions = [pos + 1]
                sequences.append((sequence_id, mutated_seq, mutations_desc, mutation_positions))

        logger.info(f"Generated {len(sequences)} round-robin hotspot-focused sequences")
        return sequences

    else:
        raise ValueError(f"Unknown multi-datasheet strategy: {strategy}")


def generate_hotspot_focused_mutations(freq_df: pd.DataFrame, num_sequences: int,
                                     min_frequency: float, prefix: str = "", hotspot_count: int = 10, combination_strategy: str = "average") -> List[Tuple[str, str, str, List[int]]]:
    """
    Generate mutations focused on hotspot positions with highest mutation rates.

    Args:
        freq_df: DataFrame with frequency data
        num_sequences: Number of sequences to generate
        min_frequency: Minimum frequency threshold
        prefix: Prefix for sequence IDs
        hotspot_count: Number of top hotspot positions to consider

    Returns:
        List of (sequence_id, sequence, mutations, mutation_positions) tuples
    """
    logger = logging.getLogger(__name__)

    # Handle multi-datasheet combination strategies
    if combination_strategy in ["stack", "round_robin"] and hasattr(freq_df, 'attrs') and 'source_dfs' in freq_df.attrs:
        return generate_multi_datasheet_hotspot_focused(
            freq_df.attrs['source_dfs'], num_sequences, min_frequency, prefix, hotspot_count, combination_strategy
        )

    # DEBUG: Log input data structure
    logger.info(f"DEBUG: freq_df shape: {freq_df.shape}")
    logger.info(f"DEBUG: freq_df columns: {list(freq_df.columns)}")
    logger.info(f"DEBUG: First few rows of freq_df:")
    logger.info(f"{freq_df.head()}")
    logger.info(f"DEBUG: min_frequency threshold: {min_frequency}")
    logger.info(f"DEBUG: hotspot_count: {hotspot_count}")

    # Detect if this is absolute frequencies and convert to relative for hotspot detection
    working_df = freq_df.copy()
    is_absolute_freq = False

    # Check if any original amino acid has frequency close to 1.0 (indicating absolute frequencies)
    for _, row in freq_df.head(20).iterrows():  # Check first 20 rows
        if 'original' in row:
            original_aa = row['original']
            if original_aa in row and row[original_aa] > 0.9:  # Original AA freq > 0.9 suggests absolute
                is_absolute_freq = True
                break

    if is_absolute_freq:
        logger.info("DEBUG: Detected absolute frequencies, converting to relative frequencies for hotspot detection")
        working_df = convert_absolute_to_relative_frequencies(freq_df)
        logger.info("DEBUG: Conversion complete")
    else:
        logger.info("DEBUG: Using frequencies as-is (appear to be relative frequencies)")

    reference_seq = get_reference_sequence(freq_df)
    logger.info(f"DEBUG: Reference sequence: {reference_seq[:50]}...")  # First 50 chars

    # Calculate mutation frequency per position (exclude original AA)
    position_scores = []
    for _, row in working_df.iterrows():
        position = int(row['position']) - 1

        # Get original amino acid
        if 'original' in row:
            original_aa = row['original']
        else:
            original_aa = reference_seq[position]

        # Sum mutation frequencies (exclude original)
        mutation_freq = sum(row[aa] for aa in AMINO_ACIDS if aa in row and aa != original_aa)
        position_scores.append((position, mutation_freq))

    # DEBUG: Log position scores
    logger.info(f"DEBUG: Top 20 position scores (mutation frequency): {sorted(position_scores, key=lambda x: x[1], reverse=True)[:20]}")

    # Sort by mutation frequency (highest first)
    position_scores.sort(key=lambda x: x[1], reverse=True)

    sequences = []
    hotspot_positions = [pos for pos, score in position_scores[:hotspot_count] if score >= min_frequency]

    # DEBUG: Log hotspot positions
    logger.info(f"DEBUG: Selected hotspot positions (top {hotspot_count}): {hotspot_positions}")
    logger.info(f"DEBUG: Number of hotspot positions: {len(hotspot_positions)}")

    if not hotspot_positions:
        logger.error(f"DEBUG: No hotspot positions found with min_frequency {min_frequency} from top {hotspot_count}")
        logger.error(f"DEBUG: Best position scores: {position_scores[:5]}")
        raise ValueError(f"No hotspot positions found above minimum frequency threshold {min_frequency}")
    
    for seq_idx in range(num_sequences):
        seq_list = list(reference_seq)
        mutations_made = []
        mutation_positions = []

        # Strategy for diversity: vary hotspot positions and mutations per sequence
        if seq_idx == 0:
            # First sequence: use top 3 hotspots with best mutations
            num_hotspots_to_use = min(len(hotspot_positions), 3)
            selected_hotspots = hotspot_positions[:num_hotspots_to_use]
            mutation_strategy = "best"
        elif seq_idx < len(hotspot_positions):
            # Next sequences: use different combinations of hotspots
            num_hotspots_to_use = min(len(hotspot_positions), 3)
            # Rotate starting position to get different combinations
            start_idx = seq_idx % len(hotspot_positions)
            selected_hotspots = hotspot_positions[start_idx:start_idx + num_hotspots_to_use]
            if len(selected_hotspots) < num_hotspots_to_use:
                # Wrap around if needed
                selected_hotspots.extend(hotspot_positions[:num_hotspots_to_use - len(selected_hotspots)])
            mutation_strategy = "best"
        else:
            # Later sequences: use top hotspots but vary mutation selection
            num_hotspots_to_use = min(len(hotspot_positions), 3)
            selected_hotspots = hotspot_positions[:num_hotspots_to_use]
            mutation_strategy = "second_best" if seq_idx % 2 == 0 else "weighted_random"

        logger.info(f"DEBUG: Sequence {seq_idx+1} - Using hotspots: {selected_hotspots}, strategy: {mutation_strategy}")

        for pos in selected_hotspots:
            row_idx = pos
            if row_idx >= len(freq_df):
                logger.warning(f"DEBUG: Position {pos} exceeds dataframe length {len(freq_df)}")
                continue

            # Use original frequencies (not converted ones) for mutation selection
            row = freq_df.iloc[row_idx]

            # Use original column if available
            if 'original' in row:
                original_aa = row['original']
            else:
                original_aa = reference_seq[pos]

            logger.info(f"DEBUG: Position {pos+1}, original AA: {original_aa}")

            # Collect all valid mutations at this position
            valid_mutations = []
            for aa in AMINO_ACIDS:
                if aa in row and aa != original_aa:
                    freq = row[aa]
                    if freq >= min_frequency:
                        valid_mutations.append((aa, freq))

            # Sort by frequency (highest first)
            valid_mutations.sort(key=lambda x: x[1], reverse=True)

            selected_aa = original_aa
            if valid_mutations:
                if mutation_strategy == "best":
                    # Use best mutation
                    selected_aa, selected_freq = valid_mutations[0]
                elif mutation_strategy == "second_best" and len(valid_mutations) > 1:
                    # Use second-best mutation
                    selected_aa, selected_freq = valid_mutations[1]
                elif mutation_strategy == "weighted_random":
                    # Weighted random selection
                    import random
                    import numpy as np
                    mutations_list = [aa for aa, freq in valid_mutations]
                    weights = [freq for aa, freq in valid_mutations]
                    # Normalize weights
                    total_weight = sum(weights)
                    if total_weight > 0:
                        weights = [w / total_weight for w in weights]
                        selected_aa = np.random.choice(mutations_list, p=weights)
                        selected_freq = dict(valid_mutations)[selected_aa]
                    else:
                        selected_aa, selected_freq = valid_mutations[0]
                else:
                    # Fallback to best
                    selected_aa, selected_freq = valid_mutations[0]
            else:
                selected_freq = 0

            # DEBUG: Log all amino acid frequencies at this position
            aa_freqs_debug = {}
            for aa in AMINO_ACIDS:
                if aa in row:
                    aa_freqs_debug[aa] = row[aa]

            logger.info(f"DEBUG: Position {pos+1} frequencies: {aa_freqs_debug}")
            logger.info(f"DEBUG: Position {pos+1}, selected mutation: {original_aa}->{selected_aa}, freq: {selected_freq}")

            # Apply mutation if found
            if selected_aa != original_aa:
                seq_list[pos] = selected_aa
                mutations_made.append(f"{original_aa}{pos+1}{selected_aa}")
                mutation_positions.append(pos + 1)
            else:
                logger.warning(f"DEBUG: No valid mutation found at position {pos+1}")

        # Create final sequence
        mutated_seq = ''.join(seq_list)
        sequence_id = f"{prefix}_{seq_idx+1:03d}" if prefix else f"hotspot_{seq_idx+1:03d}"
        mutations_desc = ",".join(mutations_made) if mutations_made else "none"

        logger.info(f"DEBUG: Generated sequence {seq_idx+1}: {mutations_desc}")
        sequences.append((sequence_id, mutated_seq, mutations_desc, mutation_positions))

    # Final validation - ensure we generated actual mutations
    mutations_found = any(mutations_desc != "none" for _, _, mutations_desc, _ in sequences)
    if not mutations_found:
        logger.error("DEBUG: No mutations were generated in any sequence!")
        logger.error(f"DEBUG: All sequences have mutations='none'")
        raise ValueError("Failed to generate any mutations. All sequences are identical to reference.")

    return sequences


def generate_top_mutations(freq_df: pd.DataFrame, num_sequences: int,
                          min_frequency: float, prefix: str = "", combination_strategy: str = "average") -> List[Tuple[str, str, str, List[int]]]:
    """
    Generate sequences using top N most frequent mutations per position.
    
    Args:
        freq_df: DataFrame with frequency data
        num_sequences: Number of sequences to generate
        min_frequency: Minimum frequency threshold
        
    Returns:
        List of (sequence_id, sequence, mutations, mutation_positions) tuples
    """
    logger = logging.getLogger(__name__)

    # Handle multi-datasheet combination strategies
    if combination_strategy in ["stack", "round_robin"] and hasattr(freq_df, 'attrs') and 'source_dfs' in freq_df.attrs:
        return generate_multi_datasheet_top_mutations(
            freq_df.attrs['source_dfs'], num_sequences, min_frequency, prefix, combination_strategy
        )

    reference_seq = get_reference_sequence(freq_df)
    sequences = []

    # For each position, get top mutations
    position_mutations = {}
    for _, row in freq_df.iterrows():
        position = int(row['position']) - 1
        # Use original column if available
        if 'original' in row:
            original_aa = row['original']
        else:
            original_aa = reference_seq[position]
        
        aa_freqs = []
        for aa in AMINO_ACIDS:
            if aa in row and aa != original_aa:
                freq = row[aa]
                if freq >= min_frequency:
                    aa_freqs.append((aa, freq))
        
        # Sort by frequency and take top 3
        aa_freqs.sort(key=lambda x: x[1], reverse=True)
        position_mutations[position] = aa_freqs[:3]
    
    # Generate sequences by combining top mutations
    for seq_idx in range(num_sequences):
        seq_list = list(reference_seq)
        mutations_made = []
        mutation_positions = []
        
        # Apply top mutation at each position (cycling through top choices)
        for pos, mut_list in position_mutations.items():
            if mut_list:
                # Cycle through top mutations
                mut_idx = seq_idx % len(mut_list)
                selected_aa, freq = mut_list[mut_idx]
                
                # Get original AA for this position from the data
                orig_row = None
                for _, row in freq_df.iterrows():
                    if int(row['position']) - 1 == pos:
                        orig_row = row
                        break
                
                if orig_row is not None and 'original' in orig_row:
                    original_aa = orig_row['original']
                else:
                    original_aa = reference_seq[pos]
                    
                seq_list[pos] = selected_aa
                mutations_made.append(f"{original_aa}{pos+1}{selected_aa}")
                mutation_positions.append(pos + 1)
        
        # Create final sequence
        mutated_seq = ''.join(seq_list)
        sequence_id = f"{prefix}_{seq_idx+1:03d}" if prefix else f"top_{seq_idx+1:03d}"
        mutations_desc = ",".join(mutations_made) if mutations_made else "none"
        
        sequences.append((sequence_id, mutated_seq, mutations_desc, mutation_positions))
    
    return sequences


def save_sequences(sequences: List[Tuple[str, str, str, List[int]]], csv_path: str, fasta_path: str):
    """
    Save sequences to CSV and FASTA files.
    
    Args:
        sequences: List of (sequence_id, sequence, mutations, mutation_positions) tuples
        csv_path: Output CSV file path
        fasta_path: Output FASTA file path
    """
    # Save CSV
    csv_data = []
    for seq_id, sequence, mutations, mut_positions in sequences:
        csv_data.append({
            'id': seq_id,
            'sequence': sequence,
            'mutations': mutations,
            'mutation_positions': ','.join(map(str, mut_positions))
        })
    
    df = pd.DataFrame(csv_data)
    df.to_csv(csv_path, index=False)
    
    # Save FASTA
    with open(fasta_path, 'w') as f:
        for seq_id, sequence, mutations, mut_positions in sequences:
            f.write(f">{seq_id} | {mutations}\n")
            # Write sequence in 80-character lines
            for i in range(0, len(sequence), 80):
                f.write(sequence[i:i+80] + '\n')


def main():
    """Main execution function."""
    parser = argparse.ArgumentParser(description='Generate sequences from mutation profiles')
    parser.add_argument('--config', required=True, help='Configuration JSON file')
    
    args = parser.parse_args()
    logger = setup_logging()
    
    try:
        # Load configuration
        with open(args.config, 'r') as f:
            config = json.load(f)
        
        # Load frequency data (single or multiple datasheets)
        frequency_datasheets = config['frequency_datasheets']
        combination_strategy = config.get('combination_strategy', 'average')

        if isinstance(frequency_datasheets, str):
            # Legacy single datasheet support
            frequency_datasheets = [frequency_datasheets]

        logger.info(f"Loading frequency data from {len(frequency_datasheets)} datasheet(s)")

        # Load and combine frequency data
        freq_df = load_and_combine_frequency_data(frequency_datasheets, combination_strategy)
        logger.info(f"Loaded/combined frequency data for {len(freq_df)} positions")
        
        # Generate sequences based on mode
        mode = config['mode']
        num_sequences = config['num_sequences']
        min_frequency = config['min_frequency']
        prefix = config.get('prefix', '')
        
        logger.info(f"Generating {num_sequences} sequences using mode: {mode}")
        if prefix:
            logger.info(f"Using prefix: {prefix}")
        
        if mode == "single_point":
            sequences = generate_single_point_mutations(freq_df, num_sequences, min_frequency, prefix, combination_strategy)
        elif mode == "weighted_random":
            sequences = generate_weighted_random_mutations(
                freq_df, num_sequences, min_frequency,
                config.get('max_mutations'), config.get('random_seed'), prefix, combination_strategy
            )
        elif mode == "hotspot_focused":
            hotspot_count = config.get('hotspot_count', 10)
            sequences = generate_hotspot_focused_mutations(freq_df, num_sequences, min_frequency, prefix, hotspot_count, combination_strategy)
        elif mode == "top_mutations":
            sequences = generate_top_mutations(freq_df, num_sequences, min_frequency, prefix, combination_strategy)
        else:
            raise ValueError(f"Unknown mode: {mode}")
        
        logger.info(f"Generated {len(sequences)} sequences")
        
        # Save sequences
        logger.info(f"Saving sequences to: {config['sequences_output']}")
        save_sequences(sequences, config['sequences_output'], config['sequences_fasta'])
        
        # Print summary and validate results
        if sequences:
            avg_mutations = np.mean([len(mut_pos) for _, _, _, mut_pos in sequences])
            logger.info(f"Average mutations per sequence: {avg_mutations:.1f}")

            sample_seq = sequences[0]
            logger.info(f"Example sequence: {sample_seq[0]} with mutations: {sample_seq[2]}")

            # DEBUG: Check if any mutations were actually generated
            mutations_found = any(mutations_desc != "none" for _, _, mutations_desc, _ in sequences)
            if not mutations_found:
                logger.error("CRITICAL ERROR: No mutations were generated!")
                logger.error("All generated sequences are identical to the reference sequence.")
                logger.error("This indicates a problem with the mutation frequency data or thresholds.")
                raise ValueError("Failed to generate any mutations. Check input data and parameters.")
        else:
            logger.error("CRITICAL ERROR: No sequences were generated!")
            raise ValueError("No sequences were generated")

        logger.info("Sequence generation complete!")
        
    except Exception as e:
        logger.error(f"Error during sequence generation: {str(e)}")
        import traceback
        logger.error(traceback.format_exc())
        sys.exit(1)


if __name__ == "__main__":
    main()