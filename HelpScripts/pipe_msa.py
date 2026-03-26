#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
MSA format conversion between CSV and A3M.

Reads an input MSA table CSV and converts each referenced MSA file to the target format.

CSV format: columns `key, sequence` (key=-1 for all rows) — Boltz2 public server convention.
A3M format: FASTA-like with `>header\\nsequence\\n` pairs, query sequence first.
"""

import argparse
import json
import os
import sys

import pandas as pd


def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Convert MSA files between CSV and A3M formats'
    )
    parser.add_argument(
        '--config', required=True,
        help='Path to MSA conversion config JSON'
    )
    return parser.parse_args()


def convert_csv_to_a3m(csv_file: str, output_a3m: str) -> str:
    """
    Convert a Boltz2-style CSV MSA file to A3M format.

    Writes the A3M header line (``#<num_seqs>\\t<query_length>``) followed by
    FASTA-like ``>index\\nsequence\\n`` pairs.  The first sequence is treated as
    the query.

    Args:
        csv_file: Input CSV with columns `key, sequence`
        output_a3m: Output A3M file path

    Returns:
        Query sequence (first sequence in the MSA)
    """
    df = pd.read_csv(csv_file)

    if 'sequence' not in df.columns:
        raise ValueError(f"CSV file missing 'sequence' column: {csv_file}")

    sequences = df['sequence'].tolist()
    if not sequences:
        raise ValueError(f"No sequences found in CSV file: {csv_file}")

    query_len = len(sequences[0])

    with open(output_a3m, 'w') as f:
        f.write(f"#{len(sequences)}\t{query_len}\n")
        for i, seq in enumerate(sequences):
            f.write(f">{100 + i}\n{seq}\n")

    return sequences[0]


def convert_a3m_to_csv(a3m_file: str, output_csv: str) -> str:
    """
    Convert an A3M file to Boltz2-style CSV format.

    Skips comment/header lines (starting with ``#``) that carry A3M metadata
    such as sequence count and column width.

    Args:
        a3m_file: Input A3M file
        output_csv: Output CSV file path

    Returns:
        Query sequence (first sequence in the A3M)
    """
    sequences = []
    current_seq = ""

    with open(a3m_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            if line.startswith('>'):
                if current_seq:
                    sequences.append(current_seq)
                    current_seq = ""
            else:
                current_seq += line

    if current_seq:
        sequences.append(current_seq)

    if not sequences:
        raise ValueError(f"No sequences found in A3M file: {a3m_file}")

    msa_data = [{'key': -1, 'sequence': seq} for seq in sequences]
    msa_df = pd.DataFrame(msa_data)
    msa_df.to_csv(output_csv, index=False)

    return sequences[0]


def main():
    args = parse_arguments()

    with open(args.config, 'r') as f:
        config = json.load(f)

    input_msa_table = config['input_msa_table']
    convert = config['convert']
    output_folder = config['output_folder']
    output_msas_csv = config['output_msas_csv']

    os.makedirs(output_folder, exist_ok=True)

    if not os.path.exists(input_msa_table):
        print(f"Error: Input MSA table not found: {input_msa_table}")
        sys.exit(1)

    df = pd.read_csv(input_msa_table)

    if 'msa_file' not in df.columns:
        print(f"Error: Input MSA table missing 'msa_file' column")
        sys.exit(1)

    ext = ".a3m" if convert == "a3m" else ".csv"
    output_rows = []

    for _, row in df.iterrows():
        msa_file = row['msa_file']
        row_id = row.get('id', '')
        seq_id = row.get('sequences.id', row_id)

        if not msa_file or not os.path.exists(msa_file):
            print(f"Warning: MSA file not found: {msa_file}")
            continue

        output_file = os.path.join(output_folder, f"{seq_id}{ext}")

        if convert == "a3m":
            query_seq = convert_csv_to_a3m(msa_file, output_file)
            print(f"Converted CSV -> A3M: {os.path.basename(msa_file)} -> {os.path.basename(output_file)}")
        else:
            query_seq = convert_a3m_to_csv(msa_file, output_file)
            print(f"Converted A3M -> CSV: {os.path.basename(msa_file)} -> {os.path.basename(output_file)}")

        output_rows.append({
            'id': f"{seq_id}_msa",
            'sequences.id': seq_id,
            'sequence': query_seq,
            'msa_file': output_file
        })

    output_df = pd.DataFrame(output_rows)
    output_df.to_csv(output_msas_csv, index=False)
    print(f"\nConverted {len(output_rows)} MSA files to {convert} format")
    print(f"Output table: {output_msas_csv}")


if __name__ == "__main__":
    main()
