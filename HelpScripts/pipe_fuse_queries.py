#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Runtime helper script for Fuse tool.

Generates fusion sequences by combining multiple sequences with linkers of variable lengths.
Reads a JSON config file specifying per-slot DataStream info and linker parameters.
Outputs CSV with sequence information including provenance and position selections in PyMOL format.
"""

import os
import sys
import json
import argparse
import pandas as pd
from typing import List, Dict, Any
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


def load_sequences_from_slot(slot: Dict[str, Any]) -> Dict[str, str]:
    """
    Load sequences from a slot definition.

    A slot has: ids, map_table, files.
    - If map_table is available, load id→sequence from CSV
    - If files are available, extract sequences from PDB/CIF files

    Args:
        slot: Slot dict with keys: ids, map_table, files

    Returns:
        Dict mapping ID → sequence string
    """
    ids = slot["ids"]
    map_table = slot.get("map_table", "")
    files = slot.get("files", [])

    # Try map_table first (CSV with id, sequence columns)
    if map_table and os.path.exists(map_table):
        df = pd.read_csv(map_table)
        if 'id' in df.columns and 'sequence' in df.columns:
            id_to_seq = {}
            for _, row in df.iterrows():
                row_id = str(row['id'])
                if row_id in ids:
                    id_to_seq[row_id] = str(row['sequence']).upper()
            if id_to_seq:
                return id_to_seq

    # Try files (PDB/CIF)
    if files:
        id_to_seq = {}
        for i, file_path in enumerate(files):
            if i >= len(ids):
                break
            seq_id = ids[i]
            if file_path.endswith('.pdb') or file_path.endswith('.cif'):
                if os.path.exists(file_path):
                    atoms = parse_pdb_file(file_path)
                    sequences = get_protein_sequence(atoms)
                    if sequences:
                        full_sequence = ''.join(sequences[chain] for chain in sorted(sequences.keys()))
                        id_to_seq[seq_id] = full_sequence
            elif file_path.endswith('.csv') and os.path.exists(file_path):
                df = pd.read_csv(file_path)
                if 'sequence' in df.columns:
                    id_to_seq[seq_id] = str(df['sequence'].iloc[0]).upper()
        if id_to_seq:
            return id_to_seq

    raise ValueError(f"Could not load sequences for slot with ids: {ids}")


def generate_fusion_sequences(
    slots: List[Dict[str, Any]],
    linker: str,
    linker_lengths: List[str],
    name_base: str
) -> List[Dict[str, Any]]:
    """
    Generate all fusion sequence combinations from multi-ID slots.

    The output is the cartesian product of all slot options × linker length options.

    Args:
        slots: List of slot dicts with ids, map_table, files
        linker: Linker sequence template
        linker_lengths: List of length range specifications for each junction
        name_base: Base name (unused in ID, IDs come from slot IDs)

    Returns:
        List of dicts with sequence info including provenance and position columns
    """
    # Load sequences for each slot
    slot_sequences = []
    for i, slot in enumerate(slots):
        print(f"Loading sequences for slot {i+1}...")
        id_to_seq = load_sequences_from_slot(slot)
        print(f"  Loaded {len(id_to_seq)} sequences")
        slot_sequences.append(id_to_seq)

    # Get per-slot ID lists
    slot_id_lists = [slot["ids"] for slot in slots]

    # Parse linker length ranges
    junction_length_lists = []
    for spec in linker_lengths:
        lengths = parse_length_spec(spec)
        junction_length_lists.append(lengths)

    # Build cartesian product axes: [slot0_ids, junction0_lengths, slot1_ids, ...]
    axes = []
    for i, slot_ids in enumerate(slot_id_lists):
        axes.append(slot_ids)
        if i < len(junction_length_lists):
            axes.append(junction_length_lists[i])

    # Generate all combinations
    results = []
    num_slots = len(slots)

    for combo in product(*axes):
        # combo alternates: (slot0_id, junc0_len, slot1_id, junc1_len, slot2_id, ...)
        # Extract slot IDs and junction lengths
        slot_ids_in_combo = [combo[i * 2] for i in range(num_slots)]
        junction_lengths_in_combo = [combo[i * 2 + 1] for i in range(num_slots - 1)]

        # Build the fused sequence and track positions
        fused_sequence = ""
        positions = []  # List of (type, idx, start, end)

        for i, seq_id in enumerate(slot_ids_in_combo):
            seq = slot_sequences[i].get(seq_id, "")
            if not seq:
                print(f"Warning: No sequence found for {seq_id} in slot {i+1}, skipping combo")
                break

            # Add sequence
            start = len(fused_sequence) + 1  # 1-indexed
            fused_sequence += seq
            end = len(fused_sequence)
            positions.append(('S', i + 1, start, end))

            # Add linker if not last slot
            if i < num_slots - 1:
                linker_len = junction_lengths_in_combo[i]
                linker_seq = linker[:linker_len]

                start = len(fused_sequence) + 1
                fused_sequence += linker_seq
                end = len(fused_sequence)
                positions.append(('L', i + 1, start, end))
        else:
            # Build result dict
            seq_id = "_".join(str(part) for part in combo)

            # Lengths string from junction lengths
            lengths_str = "-".join(str(l) for l in junction_lengths_in_combo)

            result = {
                'id': seq_id,
                'sequence': fused_sequence,
                'lengths': lengths_str
            }

            # Add provenance columns
            for slot_idx in range(num_slots):
                result[f'sequences_{slot_idx+1}.id'] = slot_ids_in_combo[slot_idx]

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
    parser.add_argument('--config', required=True,
                       help='Path to JSON config file')

    args = parser.parse_args()

    # Load config
    with open(args.config, 'r') as f:
        config = json.load(f)

    name = config["name"]
    linker = config["linker"]
    linker_lengths = config["linker_lengths"]
    slots = config["slots"]
    output_csv = config["output_csv"]
    output_fasta = config["output_fasta"]

    num_slots = len(slots)
    print(f"Processing {num_slots} slots")
    for i, slot in enumerate(slots):
        print(f"  Slot {i+1}: {len(slot['ids'])} sequences")

    print(f"Linker: {linker}")
    print(f"Linker length ranges: {linker_lengths}")

    # Validate
    expected_junctions = num_slots - 1
    if len(linker_lengths) != expected_junctions:
        raise ValueError(
            f"Expected {expected_junctions} linker length specs for {num_slots} slots, "
            f"got {len(linker_lengths)}"
        )

    # Generate fusion sequences
    results = generate_fusion_sequences(
        slots=slots,
        linker=linker,
        linker_lengths=linker_lengths,
        name_base=name
    )

    print(f"\nGenerated {len(results)} fusion sequences")

    # Create output directory if needed
    output_dir = os.path.dirname(output_csv)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    # Determine column order dynamically
    provenance_columns = [f'sequences_{i+1}.id' for i in range(num_slots)]
    position_columns = []
    for i in range(1, num_slots + 1):
        position_columns.append(f"S{i}")
        if i < num_slots:
            position_columns.append(f"L{i}")

    base_columns = ['id', 'sequence', 'lengths']
    all_columns = base_columns + provenance_columns + position_columns

    # Write CSV with ordered columns
    df = pd.DataFrame(results)
    df = df[all_columns]
    df.to_csv(output_csv, index=False)
    print(f"Saved CSV: {output_csv}")

    # Write FASTA
    write_fasta(results, output_fasta)
    print(f"Saved FASTA: {output_fasta}")

    # Show sample output
    if results:
        sample = results[0]
        print(f"\nSample output (first sequence):")
        print(f"  ID: {sample['id']}")
        print(f"  Length: {len(sample['sequence'])} residues")
        for col in provenance_columns + position_columns:
            if col in sample:
                print(f"  {col}: {sample[col]}")


if __name__ == "__main__":
    main()
