#!/usr/bin/env python3
"""
Runtime helper script for Sequence tool.

Creates sequence CSV and FASTA files from raw sequence strings.
"""

import os
import sys
import argparse
import json
import pandas as pd
from typing import Dict, List, Any


def create_sequence_files(config_data: Dict[str, Any]) -> int:
    """
    Create sequence CSV and FASTA files from configuration.

    Args:
        config_data: Configuration dictionary with sequence parameters

    Returns:
        0 on success, 1 on failure
    """
    custom_ids = config_data['custom_ids']
    sequences = config_data['sequences']
    types = config_data['types']
    output_folder = config_data['output_folder']
    sequences_csv = config_data['sequences_csv']
    sequences_fasta = config_data['sequences_fasta']

    print(f"Creating sequence files for {len(sequences)} sequences")
    print(f"Output folder: {output_folder}")

    # Create output directory
    os.makedirs(output_folder, exist_ok=True)

    # Create CSV file
    data = []
    for seq_id, seq, seq_type in zip(custom_ids, sequences, types):
        # Clean sequence (remove whitespace)
        clean_seq = ''.join(seq.split())
        data.append({
            'id': seq_id,
            'sequence': clean_seq,
            'type': seq_type,
            'length': len(clean_seq)
        })

    df = pd.DataFrame(data)
    df.to_csv(sequences_csv, index=False)
    print(f"Created CSV: {sequences_csv}")

    # Create FASTA file
    with open(sequences_fasta, 'w') as f:
        for seq_id, seq, seq_type in zip(custom_ids, sequences, types):
            # Clean sequence (remove whitespace)
            clean_seq = ''.join(seq.split())
            # Write FASTA entry
            f.write(f">{seq_id}\n")
            # Write sequence in lines of 80 characters
            for i in range(0, len(clean_seq), 80):
                f.write(clean_seq[i:i+80] + "\n")

    print(f"Created FASTA: {sequences_fasta}")

    # Summary
    print(f"\n=== SEQUENCE SUMMARY ===")
    print(f"Total sequences: {len(sequences)}")
    for seq_id, seq, seq_type in zip(custom_ids, sequences, types):
        clean_seq = ''.join(seq.split())
        seq_preview = clean_seq[:30] + "..." if len(clean_seq) > 30 else clean_seq
        print(f"  {seq_id}: {seq_preview} ({seq_type}, {len(clean_seq)} residues)")

    return 0


def main():
    parser = argparse.ArgumentParser(description='Create sequence files from raw sequences')
    parser.add_argument('--config', required=True, help='JSON config file with sequence parameters')

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
    required_params = ['custom_ids', 'sequences', 'types', 'output_folder',
                       'sequences_csv', 'sequences_fasta']
    for param in required_params:
        if param not in config_data:
            print(f"Error: Missing required parameter: {param}")
            sys.exit(1)

    try:
        result = create_sequence_files(config_data)
        if result == 0:
            print("\nSequence files created successfully")
        sys.exit(result)

    except Exception as e:
        print(f"Error creating sequence files: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
