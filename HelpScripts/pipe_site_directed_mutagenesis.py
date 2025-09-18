#!/usr/bin/env python3
"""
Site-Directed Mutagenesis helper script for generating amino acid substitutions.

This script generates systematic amino acid substitutions at specified positions
using various class-based strategies for mutagenesis studies.

Usage:
    python pipe_site_directed_mutagenesis.py --sequence-source [direct|file]
           --sequence [SEQUENCE|PATH] --sequence-id ID --position POS
           --mode MODE --include-original [true|false] --exclude EXCLUDE
           --output OUTPUT_CSV

Arguments:
    --sequence-source: Source type ('direct' for sequence string, 'file' for CSV file)
    --sequence: Protein sequence string or path to sequences CSV file
    --sequence-id: Identifier for the sequence
    --position: Position for mutagenesis (1-indexed)
    --mode: Mutagenesis mode (saturation, hydrophobic, charged, etc.)
    --include-original: Whether to include original amino acid (true/false)
    --exclude: Amino acids to exclude (single string, e.g., "CEFGL")
    --output: Output CSV file path

Output:
    CSV file with columns: id, sequence, mutation, position, original_aa, new_aa
"""

import sys
import os
import pandas as pd
import argparse
from typing import List, Tuple, Dict


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


def load_sequence_from_file(file_path: str) -> Tuple[str, str]:
    """
    Load sequence from CSV file (typically from previous pipeline tool).

    Args:
        file_path: Path to CSV file with sequences

    Returns:
        Tuple of (sequence_id, sequence)
    """
    try:
        df = pd.read_csv(file_path)
        if len(df) == 0:
            raise ValueError("Empty sequences file")

        # Use first sequence if multiple are present
        if len(df) > 1:
            print(f"Warning: Multiple sequences found, using first one")

        row = df.iloc[0]

        # Try different possible column names for sequence ID
        seq_id = None
        for col in ['id', 'sequence_id', 'name']:
            if col in row:
                seq_id = str(row[col])
                break

        if seq_id is None:
            seq_id = "sequence"
            print("Warning: No ID column found, using 'sequence' as default")

        # Try different possible column names for sequence
        sequence = None
        for col in ['sequence', 'seq', 'protein_sequence']:
            if col in row:
                sequence = str(row[col]).upper()
                break

        if sequence is None:
            raise ValueError("No sequence column found in file")

        return seq_id, sequence

    except Exception as e:
        raise ValueError(f"Error loading sequence from file {file_path}: {e}")


def generate_mutants(sequence: str, sequence_id: str, position: int, mode: str,
                    include_original: bool, exclude: str) -> List[Dict]:
    """
    Generate mutant sequences based on specified parameters.

    Args:
        sequence: Original protein sequence
        sequence_id: Identifier for the sequence
        position: Position for mutagenesis (1-indexed)
        mode: Mutagenesis mode
        include_original: Whether to include original amino acid
        exclude: Amino acids to exclude

    Returns:
        List of dictionaries with mutant information
    """
    # Validate inputs
    if position <= 0 or position > len(sequence):
        raise ValueError(f"Position {position} is out of range for sequence of length {len(sequence)}")

    if mode not in AMINO_ACID_CLASSES:
        raise ValueError(f"Unknown mode: {mode}. Available modes: {list(AMINO_ACID_CLASSES.keys())}")

    # Get amino acids for this mode
    amino_acids = AMINO_ACID_CLASSES[mode]

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
        mutant_id = f"{sequence_id}_{original_aa}{position}{new_aa}"

        # Create mutation description
        mutation = f"{original_aa}{position}{new_aa}"

        mutants.append({
            'id': mutant_id,
            'sequence': mutant_sequence,
            'mutation': mutation,
            'position': position,
            'original_aa': original_aa,
            'new_aa': new_aa
        })

    return mutants


def main():
    parser = argparse.ArgumentParser(description='Generate site-directed mutants')
    parser.add_argument('--sequence-source', required=True, choices=['direct', 'file'],
                       help='Source type: direct sequence string or file path')
    parser.add_argument('--sequence', required=True,
                       help='Protein sequence string or path to sequences file')
    parser.add_argument('--sequence-id', required=True,
                       help='Identifier for the sequence')
    parser.add_argument('--position', type=int, required=True,
                       help='Position for mutagenesis (1-indexed)')
    parser.add_argument('--mode', required=True,
                       help='Mutagenesis mode')
    parser.add_argument('--include-original', type=str, choices=['true', 'false'], default='false',
                       help='Whether to include original amino acid')
    parser.add_argument('--exclude', default='',
                       help='Amino acids to exclude (single string)')
    parser.add_argument('--output', required=True,
                       help='Output CSV file path')

    args = parser.parse_args()

    try:
        # Convert string boolean to actual boolean
        include_original = args.include_original.lower() == 'true'

        # Load sequence
        if args.sequence_source == 'direct':
            sequence = args.sequence.upper()
            sequence_id = args.sequence_id
        else:  # file
            sequence_id, sequence = load_sequence_from_file(args.sequence)
            # Override with provided sequence_id if it's different from file
            if args.sequence_id != "sequence":
                sequence_id = args.sequence_id

        print(f"Processing sequence: {sequence_id}")
        print(f"Sequence length: {len(sequence)}")
        print(f"Position: {args.position}")
        print(f"Original amino acid: {sequence[args.position-1]}")
        print(f"Mode: {args.mode}")
        print(f"Include original: {include_original}")
        print(f"Exclude: {args.exclude}")

        # Generate mutants
        mutants = generate_mutants(
            sequence=sequence,
            sequence_id=sequence_id,
            position=args.position,
            mode=args.mode,
            include_original=include_original,
            exclude=args.exclude
        )

        print(f"Generated {len(mutants)} mutants")

        # Create output dataframe
        df = pd.DataFrame(mutants)

        # Ensure output directory exists
        os.makedirs(os.path.dirname(args.output), exist_ok=True)

        # Save to CSV
        df.to_csv(args.output, index=False)

        print(f"Mutants saved to: {args.output}")

        # Print summary
        print("\nSummary:")
        print(f"  Original sequence: {sequence_id}")
        print(f"  Position: {args.position} ({sequence[args.position-1]})")
        print(f"  Mode: {args.mode}")
        print(f"  Mutants generated: {len(mutants)}")
        if mutants:
            print(f"  Example mutants: {', '.join([m['mutation'] for m in mutants[:5]])}")
            if len(mutants) > 5:
                print(f"  ... and {len(mutants) - 5} more")

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()