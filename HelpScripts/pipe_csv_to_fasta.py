# Copyright (c) 2026 Gianluca Quargnali @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.


"""
Convert between CSV and FASTA formats for sequence files.

CSV format expects columns: id, sequence (and optionally others)
FASTA format: >id followed by sequence on next line
"""
import argparse
import pandas as pd
import os

parser = argparse.ArgumentParser(description='Convert between CSV and FASTA sequence formats')
parser.add_argument('input_file', type=str, help='Input file path (CSV or FASTA)')
parser.add_argument('output_file', type=str, help='Output file path')
parser.add_argument('--reverse', action='store_true',
                    help='Convert FASTA to CSV (default: CSV to FASTA)')
parser.add_argument('--id-column', type=str, default='id',
                    help='Column name for sequence ID in CSV (default: id)')
parser.add_argument('--seq-column', type=str, default='sequence',
                    help='Column name for sequence in CSV (default: sequence)')

args = parser.parse_args()


def csv_to_fasta(csv_file, fasta_file, id_col='id', seq_col='sequence'):
    """
    Convert CSV file with sequences to FASTA format.

    Args:
        csv_file: Input CSV file path
        fasta_file: Output FASTA file path
        id_col: Column name for sequence IDs
        seq_col: Column name for sequences
    """
    df = pd.read_csv(csv_file)

    if id_col not in df.columns:
        raise ValueError(f"ID column '{id_col}' not found in CSV. Available columns: {list(df.columns)}")
    if seq_col not in df.columns:
        raise ValueError(f"Sequence column '{seq_col}' not found in CSV. Available columns: {list(df.columns)}")

    with open(fasta_file, 'w') as f:
        for idx, row in df.iterrows():
            seq_id = str(row[id_col])
            sequence = str(row[seq_col])
            f.write(f">{seq_id}\n{sequence}\n")

    print(f"Converted {len(df)} sequences from CSV to FASTA")


def fasta_to_csv(fasta_file, csv_file, id_col='id', seq_col='sequence'):
    """
    Convert FASTA file to CSV format.

    Args:
        fasta_file: Input FASTA file path
        csv_file: Output CSV file path
        id_col: Column name for sequence IDs
        seq_col: Column name for sequences
    """
    sequences = []

    with open(fasta_file, 'r') as f:
        current_id = None
        current_seq = []

        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith('>'):
                # Save previous sequence if exists
                if current_id is not None:
                    sequences.append({
                        id_col: current_id,
                        seq_col: ''.join(current_seq)
                    })

                # Start new sequence
                current_id = line[1:].split()[0]  # Take first word after >
                current_seq = []
            else:
                current_seq.append(line)

        # Save last sequence
        if current_id is not None:
            sequences.append({
                id_col: current_id,
                seq_col: ''.join(current_seq)
            })

    df = pd.DataFrame(sequences)
    df.to_csv(csv_file, index=False)
    print(f"Converted {len(df)} sequences from FASTA to CSV")


if __name__ == "__main__":
    if args.reverse:
        fasta_to_csv(args.input_file, args.output_file, args.id_column, args.seq_column)
    else:
        csv_to_fasta(args.input_file, args.output_file, args.id_column, args.seq_column)
