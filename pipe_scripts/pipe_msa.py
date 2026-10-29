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

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from biopipelines.id_map_utils import prune_redundant_provenance_columns  # noqa: E402


def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Convert MSA files between CSV and A3M formats'
    )
    parser.add_argument(
        '--config', required=True,
        help='Path to MSA conversion config JSON'
    )
    return parser.parse_args()


def a3m_meta_path(csv_file: str) -> str:
    """Sidecar beside a CSV MSA holding its A3M `#` meta line and `>` headers.

    Boltz rejects a CSV whose columns are not exactly ('key', 'sequence'), so the
    identifiers cannot be a column and live next to the file instead.
    """
    return os.path.splitext(csv_file)[0] + ".a3m_meta"


def read_a3m_meta(csv_file: str):
    """Return (meta_line, headers) for a CSV MSA, or ("", []) when there is no sidecar."""
    path = a3m_meta_path(csv_file)
    if not os.path.exists(path):
        return "", []
    with open(path) as f:
        lines = f.read().split("\n")
    meta = lines[0].strip() if lines else ""
    headers = list(lines[1:])
    while headers and headers[-1] == "":
        headers.pop()
    return meta, headers


def write_a3m_meta(csv_file: str, meta: str, headers) -> None:
    """Persist the A3M meta line and per-row headers beside a CSV MSA."""
    if not meta and not any(headers):
        return
    with open(a3m_meta_path(csv_file), 'w') as f:
        f.write((meta or "") + "\n")
        for hdr in headers:
            f.write(hdr + "\n")


def convert_csv_to_a3m(csv_file: str, output_a3m: str) -> str:
    """
    Convert a Boltz2-style CSV MSA file to A3M format.

    Writes the A3M header line (``#<query_length>\\t<cardinality>``) followed by
    FASTA-like ``>index\\nsequence\\n`` pairs.  The first sequence is treated as
    the query.

    Field order matters: ColabFold splits that line on tab and reads field 1 as
    the per-chain query length(s) and field 2 as the chain cardinality, then
    treats anything with cardinality != 1 as a multi-chain complex.

    A CSV with an `.a3m_meta` sidecar round-trips its real UniRef identifiers; without
    one the headers are synthesised as indices, which ColabFold cannot pair chains on.

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

    meta, headers = read_a3m_meta(csv_file)
    headers = (headers + [""] * len(sequences))[:len(sequences)]
    if not meta:
        meta = f"#{len(sequences[0])}\t1"

    with open(output_a3m, 'w') as f:
        f.write(meta + "\n")
        for i, (seq, hdr) in enumerate(zip(sequences, headers)):
            f.write(f">{hdr or (100 + i)}\n{seq}\n")

    return sequences[0]


def convert_a3m_to_csv(a3m_file: str, output_csv: str) -> str:
    """
    Convert an A3M file to Boltz2-style CSV format.

    The ``#`` metadata line and every ``>`` header are kept in an `.a3m_meta` sidecar,
    so converting the CSV back to A3M restores the original identifiers.

    Args:
        a3m_file: Input A3M file
        output_csv: Output CSV file path

    Returns:
        Query sequence (first sequence in the A3M)
    """
    sequences = []
    headers = []
    meta = ""
    current_seq = ""

    with open(a3m_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('#'):
                if not meta:
                    meta = line
                continue
            if line.startswith('>'):
                if current_seq:
                    sequences.append(current_seq)
                    current_seq = ""
                headers.append(line[1:])
            else:
                current_seq += line

    if current_seq:
        sequences.append(current_seq)

    if not sequences:
        raise ValueError(f"No sequences found in A3M file: {a3m_file}")

    headers = (headers + [""] * len(sequences))[:len(sequences)]
    msa_data = [{'key': -1, 'sequence': seq} for seq in sequences]
    msa_df = pd.DataFrame(msa_data)
    msa_df.to_csv(output_csv, index=False)

    write_a3m_meta(output_csv, meta, headers)

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

    if 'file' not in df.columns:
        print(f"Error: Input MSA table missing 'file' column")
        sys.exit(1)

    ext = ".a3m" if convert == "a3m" else ".csv"
    output_rows = []

    for _, row in df.iterrows():
        msa_file = row['file']
        out_id = str(row.get('id', '') or '')
        provenance = row.get('sequences.id', '')
        provenance = '' if pd.isna(provenance) else str(provenance)

        if not msa_file or not os.path.exists(msa_file):
            print(f"Warning: MSA file not found: {msa_file}")
            continue

        output_file = os.path.join(output_folder, f"{out_id}{ext}")

        if convert == "a3m":
            query_seq = convert_csv_to_a3m(msa_file, output_file)
            print(f"Converted CSV -> A3M: {os.path.basename(msa_file)} -> {os.path.basename(output_file)}")
        else:
            query_seq = convert_a3m_to_csv(msa_file, output_file)
            print(f"Converted A3M -> CSV: {os.path.basename(msa_file)} -> {os.path.basename(output_file)}")

        output_rows.append({
            'id': out_id,
            'sequences.id': provenance,
            'sequence': query_seq,
            'file': output_file
        })

    output_df = prune_redundant_provenance_columns(
        pd.DataFrame(output_rows, columns=['id', 'sequences.id', 'sequence', 'file']))
    output_df.to_csv(output_msas_csv, index=False)
    print(f"\nConverted {len(output_rows)} MSA files to {convert} format")
    print(f"Output table: {output_msas_csv}")


if __name__ == "__main__":
    main()
