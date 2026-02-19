#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Mutagenesis helper script for generating amino acid substitutions.

This script generates systematic amino acid substitutions at specified positions
using various class-based strategies for mutagenesis studies.

Usage:
    python pipe_mutagenesis.py --sequences MAP_TABLE_CSV --position POS
           --mode MODE [--mutate-to AAS] --include-original [true|false]
           --exclude EXCLUDE --output OUTPUT_CSV

Arguments:
    --sequences: Path to map_table CSV with id and sequence columns
    --position: Position for mutagenesis (1-indexed)
    --mode: Mutagenesis mode (specific, saturation, hydrophobic, charged, etc.)
    --mutate-to: Target amino acid(s) for specific mode (e.g., "A" or "AV")
    --include-original: Whether to include original amino acid (true/false)
    --exclude: Amino acids to exclude (single string, e.g., "CEFGL")
    --output: Output CSV file path

Output:
    CSV file with columns: id, sequence, mutations, mutation_positions, original_aa, new_aa
"""

import sys
import os
import pandas as pd
import argparse
from typing import List, Dict


# Amino acid classifications
AMINO_ACID_CLASSES = {
    "saturation": "ACDEFGHIKLMNPQRSTVWY",  # All 20 amino acids
    "hydrophobic": "AFILMVWY",
    "hydrophilic": "DEHKNQRST",
    "charged": "DEHKR",
    "polar": "CDEHKNQRSTY",
    "nonpolar": "AFGILMPVW",
    "aromatic": "FHWY",
    "aliphatic": "AGILV",
    "positive": "HKR",
    "negative": "DE"
}


def load_sequences_from_map_table(file_path: str) -> List[Dict]:
    """
    Load sequences from a map_table CSV (id-value DataStream).

    Supports chaining: if the CSV has 'mutations' and 'mutation_positions'
    columns (from a previous Mutagenesis), these are preserved and accumulated.

    Args:
        file_path: Path to map_table CSV with id and sequence columns

    Returns:
        List of dicts with keys: id, sequence, prior_mutations, prior_positions
    """
    df = pd.read_csv(file_path)
    if len(df) == 0:
        raise ValueError("Empty sequences file")

    results = []
    for _, row in df.iterrows():
        seq_id = str(row['id'])
        sequence = str(row['sequence']).upper()

        # Get prior mutations/positions if chaining from a previous Mutagenesis
        prior_mutations = ""
        prior_positions = ""
        if 'mutations' in df.columns and pd.notna(row['mutations']):
            prior_mutations = str(row['mutations'])
        if 'mutation_positions' in df.columns and pd.notna(row['mutation_positions']):
            prior_positions = str(row['mutation_positions'])

        results.append({
            'id': seq_id,
            'sequence': sequence,
            'prior_mutations': prior_mutations,
            'prior_positions': prior_positions
        })

    return results


def generate_mutants(sequence: str, sequence_id: str, position: int, mode: str,
                    include_original: bool, exclude: str,
                    mutate_to: str = "",
                    prior_mutations: str = "",
                    prior_positions: str = "") -> List[Dict]:
    """
    Generate mutant sequences based on specified parameters.

    Args:
        sequence: Original protein sequence
        sequence_id: Identifier for the sequence
        position: Position for mutagenesis (1-indexed)
        mode: Mutagenesis mode
        include_original: Whether to include original amino acid
        exclude: Amino acids to exclude
        mutate_to: Target amino acid(s) for "specific" mode
        prior_mutations: Comma-separated prior mutations (from chaining)
        prior_positions: Plus-separated prior positions in PyMOL format (from chaining)

    Returns:
        List of dictionaries with mutant information
    """
    # Validate inputs
    if position <= 0 or position > len(sequence):
        raise ValueError(f"Position {position} is out of range for sequence of length {len(sequence)}")

    if mode == "specific":
        if not mutate_to:
            raise ValueError("mutate_to is required when mode is 'specific'")
        amino_acids = mutate_to.upper()
    elif mode in AMINO_ACID_CLASSES:
        amino_acids = AMINO_ACID_CLASSES[mode]
    else:
        raise ValueError(f"Unknown mode: {mode}. Available modes: ['specific'] + {list(AMINO_ACID_CLASSES.keys())}")

    # Remove excluded amino acids
    if exclude:
        amino_acids = ''.join(aa for aa in amino_acids if aa not in exclude.upper())

    # Get original amino acid at position
    original_aa = sequence[position - 1]

    # Remove original amino acid unless include_original=True
    if not include_original and original_aa in amino_acids:
        amino_acids = amino_acids.replace(original_aa, '')

    if not amino_acids:
        raise ValueError("No amino acids left after applying exclusions and original filtering")

    # Generate mutants
    mutants = []
    for new_aa in amino_acids:
        # Create mutant sequence
        mutant_sequence = sequence[:position-1] + new_aa + sequence[position:]

        # Generate mutant ID
        mutant_id = f"{sequence_id}_{position}{new_aa}"

        # Create mutation description (standard notation: OrigPosNew)
        this_mutation = f"{original_aa}{position}{new_aa}"

        # Accumulate mutations and positions from prior rounds
        if prior_mutations:
            mutations = f"{prior_mutations},{this_mutation}"
        else:
            mutations = this_mutation

        if prior_positions:
            mutation_positions = f"{prior_positions}+{position}"
        else:
            mutation_positions = str(position)

        mutants.append({
            'id': mutant_id,
            'original.id': sequence_id,
            'sequence': mutant_sequence,
            'mutations': mutations,
            'mutation_positions': mutation_positions,
            'original_aa': original_aa,
            'new_aa': new_aa
        })

    return mutants


def main():
    parser = argparse.ArgumentParser(description='Generate mutants')
    parser.add_argument('--sequences', required=True,
                       help='Path to map_table CSV with id and sequence columns')
    parser.add_argument('--position', type=int, required=True,
                       help='Position for mutagenesis (1-indexed)')
    parser.add_argument('--mode', required=True,
                       help='Mutagenesis mode')
    parser.add_argument('--mutate-to', default='',
                       help='Target amino acid(s) for specific mode (e.g., "A" or "AV")')
    parser.add_argument('--include-original', type=str, choices=['true', 'false'], default='false',
                       help='Whether to include original amino acid')
    parser.add_argument('--exclude', default='',
                       help='Amino acids to exclude (single string)')
    parser.add_argument('--output', required=True,
                       help='Output CSV file path')

    args = parser.parse_args()

    try:
        include_original = args.include_original.lower() == 'true'

        # Load all sequences from map_table
        sequences = load_sequences_from_map_table(args.sequences)

        # Generate mutants for all input sequences
        all_mutants = []
        for seq_info in sequences:
            sequence = seq_info['sequence']
            sequence_id = seq_info['id']

            print(f"Processing sequence: {sequence_id}")
            print(f"  Sequence length: {len(sequence)}")
            print(f"  Position: {args.position}")
            print(f"  Original amino acid: {sequence[args.position-1]}")
            if seq_info['prior_mutations']:
                print(f"  Prior mutations: {seq_info['prior_mutations']}")

            mutants = generate_mutants(
                sequence=sequence,
                sequence_id=sequence_id,
                position=args.position,
                mode=args.mode,
                include_original=include_original,
                exclude=args.exclude,
                mutate_to=args.mutate_to,
                prior_mutations=seq_info['prior_mutations'],
                prior_positions=seq_info['prior_positions']
            )

            all_mutants.extend(mutants)
            print(f"  Generated {len(mutants)} mutants")

        print(f"\nTotal mutants generated: {len(all_mutants)}")

        # Create output dataframe
        df = pd.DataFrame(all_mutants)

        # Ensure output directory exists
        os.makedirs(os.path.dirname(args.output), exist_ok=True)

        # Save to CSV
        df.to_csv(args.output, index=False)

        print(f"Mutants saved to: {args.output}")

        # Print summary
        print(f"\nSummary:")
        print(f"  Mode: {args.mode}")
        print(f"  Position: {args.position}")
        print(f"  Input sequences: {len(sequences)}")
        print(f"  Total mutants: {len(all_mutants)}")
        if all_mutants:
            print(f"  Example mutations: {', '.join([m['mutations'] for m in all_mutants[:5]])}")
            if len(all_mutants) > 5:
                print(f"  ... and {len(all_mutants) - 5} more")

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
