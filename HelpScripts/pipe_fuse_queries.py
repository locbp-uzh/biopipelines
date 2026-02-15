#!/usr/bin/env python3
# Copyright (c) 2026 Gianluca Quargnali @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Runtime helper script for Fuse tool.

Generates fusion sequences by combining multiple sequences with linkers of variable lengths.
Outputs CSV with sequence information including sequence/linker position selections in PyMOL format.
"""

import os
import sys
import argparse
import pandas as pd
from typing import List, Tuple, Dict, Any
from itertools import product

# Import pdb_parser for extracting sequences from PDB files
try:
    from pdb_parser import parse_pdb_file, get_protein_sequence
except ImportError:
    # Direct path for testing
    sys.path.insert(0, os.path.dirname(__file__))
    from pdb_parser import parse_pdb_file, get_protein_sequence


def parse_length_spec(spec: str) -> List[int]:
    """
    Parse length specification like '1-6' or '3+5-7'.

    Args:
        spec: Length specification string

    Returns:
        List of integer lengths
    """
    lengths = []
    if '+' in spec:
        for part in spec.split('+'):
            lengths.extend(parse_length_spec(part.strip()))
    elif '-' in spec and not spec.startswith('-'):
        if spec.count('-') == 1:
            start, end = map(int, spec.split('-'))
            lengths.extend(range(start, end + 1))
        else:
            lengths.append(int(spec))
    else:
        lengths.append(int(spec))
    return lengths


def extract_sequence(protein_input: str) -> str:
    """
    Extract sequence from protein input (either a sequence string or PDB file path).

    Args:
        protein_input: Either an amino acid sequence or path to a PDB file

    Returns:
        Amino acid sequence string
    """
    # Check if it looks like a file path (ends with .pdb or .cif, or contains path separators)
    if protein_input.endswith('.pdb') or protein_input.endswith('.cif'):
        if not os.path.exists(protein_input):
            raise FileNotFoundError(f"PDB file not found: {protein_input}")

        atoms = parse_pdb_file(protein_input)
        sequences = get_protein_sequence(atoms)

        if not sequences:
            raise ValueError(f"No protein sequence found in: {protein_input}")

        # Concatenate all chains
        full_sequence = ''.join(sequences[chain] for chain in sorted(sequences.keys()))
        return full_sequence

    # Otherwise treat as sequence string
    return protein_input


def generate_fusion_sequences(
    sequences: List[str],
    linker: str,
    linker_lengths: List[str],
    name_base: str
) -> List[Dict[str, Any]]:
    """
    Generate all fusion sequence combinations with sequence/linker position information.

    Args:
        sequences: List of sequences (or PDB file paths)
        linker: Linker sequence template
        linker_lengths: List of length range specifications for each junction
        name_base: Base name for output sequence IDs

    Returns:
        List of dicts with sequence info including S1, L1, S2, L2, S3, etc.
    """
    # Extract sequences (handles both raw sequences and PDB files)
    input_sequences = []

    for i, seq_input in enumerate(sequences):
        seq = extract_sequence(seq_input)
        input_sequences.append(seq)

    # Parse linker length ranges
    length_ranges = [parse_length_spec(spec) for spec in linker_lengths]

    # Generate all combinations
    results = []

    for length_combo in product(*length_ranges):
        # Build the fused sequence and track positions
        fused_sequence = ""
        positions = []  # List of (type, idx, start, end) where type is 'S' or 'L'

        for i, seq in enumerate(input_sequences):
            # Add sequence
            start = len(fused_sequence) + 1  # 1-indexed
            fused_sequence += seq
            end = len(fused_sequence)
            positions.append(('S', i + 1, start, end))

            # Add linker if not last sequence
            if i < len(input_sequences) - 1:
                linker_len = length_combo[i]
                linker_seq = linker[:linker_len]

                start = len(fused_sequence) + 1
                fused_sequence += linker_seq
                end = len(fused_sequence)
                positions.append(('L', i + 1, start, end))

        # Build result dict
        # Format lengths string
        lengths_str = "-".join(str(l) for l in length_combo)

        # Build sequence ID
        seq_id = f"{name_base}_{lengths_str}" if name_base else f"fused_{lengths_str}"

        result = {
            'id': seq_id,
            'sequence': fused_sequence,
            'lengths': lengths_str
        }

        # Add S1, L1, S2, L2, S3, etc. columns in PyMOL selection format
        for pos_type, idx, start, end in positions:
            col_name = f"{pos_type}{idx}"
            if start == end:
                result[col_name] = str(start)
            else:
                result[col_name] = f"{start}-{end}"

        results.append(result)

    return results


def write_fasta(sequences: List[Dict[str, Any]], fasta_path: str) -> None:
    """
    Write sequences to FASTA file.

    Args:
        sequences: List of sequence dicts with 'id' and 'sequence' keys
        fasta_path: Output FASTA file path
    """
    with open(fasta_path, 'w') as f:
        for seq_data in sequences:
            f.write(f">{seq_data['id']}\n")
            # Write sequence in 80-character lines
            seq = seq_data['sequence']
            for i in range(0, len(seq), 80):
                f.write(f"{seq[i:i+80]}\n")


def main():
    parser = argparse.ArgumentParser(
        description='Generate fusion sequences with linker combinations'
    )
    parser.add_argument('name', help='Base name for output sequence IDs')
    parser.add_argument('sequences', help='Semicolon-separated list of sequences or PDB file paths')
    parser.add_argument('linker', help='Linker sequence template')
    parser.add_argument('linker_lengths', help='Linker length specifications (e.g., "L2-4;2-4L")')
    parser.add_argument('output_csv', help='Output CSV file path')
    parser.add_argument('output_fasta', help='Output FASTA file path')

    args = parser.parse_args()

    # Parse sequences
    sequences = args.sequences.split(';')
    print(f"Processing {len(sequences)} sequences")
    for i, s in enumerate(sequences):
        if s.endswith('.pdb') or s.endswith('.cif'):
            print(f"  Sequence {i+1}: {os.path.basename(s)} (PDB file)")
        else:
            print(f"  Sequence {i+1}: {len(s)} residues")

    # Parse linker lengths - strip leading/trailing L markers
    linker_lengths_str = args.linker_lengths
    if linker_lengths_str.startswith('L'):
        linker_lengths_str = linker_lengths_str[1:]
    if linker_lengths_str.endswith('L'):
        linker_lengths_str = linker_lengths_str[:-1]

    linker_lengths = linker_lengths_str.split(';')
    print(f"Linker: {args.linker}")
    print(f"Linker length ranges: {linker_lengths}")

    # Validate
    expected_junctions = len(sequences) - 1
    if len(linker_lengths) != expected_junctions:
        raise ValueError(
            f"Expected {expected_junctions} linker length specs for {len(sequences)} sequences, "
            f"got {len(linker_lengths)}"
        )

    # Generate fusion sequences
    results = generate_fusion_sequences(
        sequences=sequences,
        linker=args.linker,
        linker_lengths=linker_lengths,
        name_base=args.name
    )

    print(f"\nGenerated {len(results)} fusion sequences")

    # Create output directory if needed
    output_dir = os.path.dirname(args.output_csv)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    # Determine column order dynamically based on number of sequences
    num_sequences = len(sequences)
    position_columns = []
    for i in range(1, num_sequences + 1):
        position_columns.append(f"S{i}")
        if i < num_sequences:
            position_columns.append(f"L{i}")

    # Write CSV with ordered columns
    base_columns = ['id', 'sequence', 'lengths']
    all_columns = base_columns + position_columns

    df = pd.DataFrame(results)
    # Reorder columns
    df = df[all_columns]
    df.to_csv(args.output_csv, index=False)
    print(f"Saved CSV: {args.output_csv}")

    # Write FASTA
    write_fasta(results, args.output_fasta)
    print(f"Saved FASTA: {args.output_fasta}")

    # Show sample output
    if results:
        sample = results[0]
        print(f"\nSample output (first sequence):")
        print(f"  ID: {sample['id']}")
        print(f"  Length: {len(sample['sequence'])} residues")
        for col in position_columns:
            if col in sample:
                print(f"  {col}: {sample[col]}")


if __name__ == "__main__":
    main()
