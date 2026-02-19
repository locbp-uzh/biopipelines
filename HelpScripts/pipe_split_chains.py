#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Runtime helper script for splitting concatenated single-chain sequences into multi-chain sequences.

This script takes sequences where multiple protein chains are concatenated into
a single sequence and splits them into separate chains based on specified
split positions.
"""

import os
import sys
import argparse
import json
import pandas as pd
from typing import List, Optional


def split_sequences(config_data: dict) -> None:
    """
    Split concatenated sequences into multiple chains.

    Args:
        config_data: Configuration dictionary with split parameters
    """
    input_csv = config_data['input_csv']
    output_csv = config_data['output_csv']
    split_positions = config_data['split_positions']
    chain_names = config_data.get('chain_names')

    print(f"Splitting sequences into {len(split_positions) + 1} chains")
    print(f"Input: {input_csv}")
    print(f"Split positions: {split_positions}")
    print(f"Output: {output_csv}")

    # Check input file exists
    if not os.path.exists(input_csv):
        raise FileNotFoundError(f"Input CSV not found: {input_csv}")

    # Read input sequences
    df = pd.read_csv(input_csv)

    # Validate input has required columns
    if 'id' not in df.columns or 'sequence' not in df.columns:
        raise ValueError("Input CSV must have 'id' and 'sequence' columns")

    print(f"Loaded {len(df)} sequences")

    # Process each sequence
    output_rows = []
    for idx, row in df.iterrows():
        seq_id = row['id']
        sequence = row['sequence']

        # Validate sequence length against split positions
        max_pos = max(split_positions)
        if max_pos >= len(sequence):
            print(f"Warning: Sequence {seq_id} (length {len(sequence)}) is shorter than split position {max_pos}, skipping")
            continue

        # Split sequence at specified positions
        chains = []
        start = 0
        for pos in split_positions:
            chains.append(sequence[start:pos])
            start = pos
        chains.append(sequence[start:])  # Last chain

        # Create output rows for each chain
        for chain_idx, chain_seq in enumerate(chains):
            if chain_names:
                chain_id = f"{seq_id}_{chain_names[chain_idx]}"
            else:
                chain_id = f"{seq_id}_{chain_idx + 1}"

            output_rows.append({
                'id': chain_id,
                'sequence': chain_seq,
                'source_id': seq_id,
                'complex_id': seq_id,  # Alias for Boltz2 multi-chain grouping
                'chain_index': chain_idx + 1,
                'chain_name': chain_names[chain_idx] if chain_names else str(chain_idx + 1),
                'chain_length': len(chain_seq)
            })

    # Create output dataframe
    output_df = pd.DataFrame(output_rows)

    # Create output directory if needed
    os.makedirs(os.path.dirname(output_csv), exist_ok=True)

    # Write output
    output_df.to_csv(output_csv, index=False)

    print(f"\nSplit {len(df)} sequences into {len(output_df)} chain sequences")
    print(f"Output written to: {output_csv}")

    # Show sample of results
    if not output_df.empty:
        print("\nFirst few results:")
        print(output_df.head(6).to_string(index=False))

        if len(output_df) > 6:
            print(f"... and {len(output_df) - 6} more rows")


def main():
    parser = argparse.ArgumentParser(description='Split concatenated sequences into multiple chains')
    parser.add_argument('--config', required=True, help='JSON config file with split parameters')

    args = parser.parse_args()

    # Load configuration
    if not os.path.exists(args.config):
        print(f"Error: Config file not found: {args.config}")
        sys.exit(1)

    try:
        with open(args.config, 'r') as f:
            config_data = json.load(f)
    except Exception as e:
        print(f"Error loading config: {e}")
        sys.exit(1)

    # Validate required parameters
    required_params = ['input_csv', 'output_csv', 'split_positions']
    for param in required_params:
        if param not in config_data:
            print(f"Error: Missing required parameter: {param}")
            sys.exit(1)

    try:
        split_sequences(config_data)
        print("\nSplitting completed successfully!")

    except Exception as e:
        print(f"Error splitting sequences: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
