#!/usr/bin/env python3
"""
Runtime helper script for StitchSequences tool.

This script stitches sequences from multiple sequence generation tools by combining
a base sequence with overlay sequences at specified positions.
"""

import os
import sys
import argparse
import json
import pandas as pd
from typing import Dict, List, Any, Optional, Tuple


def parse_position_string(position_str: str) -> List[int]:
    """
    Parse position string into list of residue numbers.

    Args:
        position_str: Position string like '10-20+30-40+145'

    Returns:
        List of residue numbers (1-indexed)
    """
    positions = []

    if not position_str or position_str.lower() in ['', 'none', 'null']:
        return positions

    # Split by '+' to get individual parts
    parts = position_str.split('+')

    for part in parts:
        part = part.strip()
        if '-' in part:
            # Range specification like '10-20'
            try:
                start, end = part.split('-')
                start_num = int(start.strip())
                end_num = int(end.strip())
                positions.extend(range(start_num, end_num + 1))
            except ValueError:
                print(f"Warning: Could not parse range: {part}")
        else:
            # Single residue number
            try:
                positions.append(int(part))
            except ValueError:
                print(f"Warning: Could not parse position: {part}")

    return sorted(list(set(positions)))  # Remove duplicates and sort


def load_positions_from_datasheet(datasheet_path: str, column_name: str) -> Dict[str, List[int]]:
    """
    Load position specifications from datasheet CSV file.

    Args:
        datasheet_path: Path to CSV file
        column_name: Column containing position specifications

    Returns:
        Dictionary mapping sequence IDs to position lists
    """
    if not os.path.exists(datasheet_path):
        raise FileNotFoundError(f"Datasheet file not found: {datasheet_path}")

    df = pd.read_csv(datasheet_path)
    if column_name not in df.columns:
        raise ValueError(f"Column '{column_name}' not found in datasheet. Available columns: {list(df.columns)}")

    # Assuming the first column contains IDs
    id_column = df.columns[0]
    positions_map = {}

    for _, row in df.iterrows():
        sequence_id = row[id_column]
        position_value = row[column_name]

        # Parse position specification
        if pd.isna(position_value) or position_value == '':
            positions = []
        else:
            positions = parse_position_string(str(position_value))

        positions_map[str(sequence_id)] = positions

    print(f"Loaded positions for {len(positions_map)} sequences from {datasheet_path}")

    return positions_map


def stitch_sequences(base_sequence: str, overlay_sequence: str, positions: List[int]) -> str:
    """
    Stitch overlay sequence into base sequence at specified positions.

    Args:
        base_sequence: Base sequence string
        overlay_sequence: Overlay sequence string
        positions: List of 1-indexed positions where overlay should be applied

    Returns:
        Stitched sequence string
    """
    if not positions:
        return base_sequence

    # Convert to lists for easy modification
    base_chars = list(base_sequence)
    overlay_chars = list(overlay_sequence)

    # Sort positions to handle them in order
    sorted_positions = sorted(positions)

    # Map positions to overlay sequence indices
    overlay_index = 0

    for pos in sorted_positions:
        # Convert 1-indexed to 0-indexed
        base_index = pos - 1

        # Check bounds
        if base_index < 0 or base_index >= len(base_chars):
            print(f"Warning: Position {pos} out of range for base sequence (length {len(base_chars)})")
            continue

        if overlay_index >= len(overlay_chars):
            print(f"Warning: Overlay sequence too short for position {pos}")
            break

        # Replace base character with overlay character
        base_chars[base_index] = overlay_chars[overlay_index]
        overlay_index += 1

    return ''.join(base_chars)


def load_sequences_from_csv(csv_path: str) -> Dict[str, str]:
    """
    Load sequences from CSV file.

    Args:
        csv_path: Path to sequences CSV file

    Returns:
        Dictionary mapping sequence IDs to sequences
    """
    if not os.path.exists(csv_path):
        raise FileNotFoundError(f"Sequences file not found: {csv_path}")

    df = pd.read_csv(csv_path)

    # Expect 'id' and 'sequence' columns
    if 'id' not in df.columns or 'sequence' not in df.columns:
        raise ValueError(f"Expected 'id' and 'sequence' columns in {csv_path}. Found: {list(df.columns)}")

    sequences = {}
    for _, row in df.iterrows():
        seq_id = str(row['id'])
        sequence = str(row['sequence'])
        sequences[seq_id] = sequence

    return sequences


def stitch_sequences_from_config(config_data: Dict[str, Any]) -> None:
    """
    Stitch sequences based on configuration.

    Args:
        config_data: Configuration dictionary with stitching parameters
    """
    sequence_sources = config_data['sequence_sources']
    position_specs = config_data['position_specs']
    output_csv = config_data['output_csv']

    print(f"Stitching sequences from {len(sequence_sources)} sources")

    # Load all sequences
    all_sequences = []
    all_sequence_ids = set()

    for source in sequence_sources:
        index = source['index']
        sequences_file = source['sequences_file']
        tool_name = source['tool_name']

        print(f"Loading sequences from {tool_name}: {sequences_file}")

        sequences = load_sequences_from_csv(sequences_file)
        all_sequences.append(sequences)
        all_sequence_ids.update(sequences.keys())

        print(f"  Loaded {len(sequences)} sequences")

    if not all_sequences:
        raise ValueError("No sequences loaded")

    # Get base sequences (first source)
    base_sequences = all_sequences[0]
    print(f"Using {len(base_sequences)} base sequences from {sequence_sources[0]['tool_name']}")

    # Process position specifications
    position_maps = []

    for spec in position_specs:
        index = spec['index']
        if spec['type'] == 'fixed':
            # Fixed position string
            fixed_positions = parse_position_string(spec['value'])
            # Create map for all sequence IDs
            position_map = {seq_id: fixed_positions for seq_id in all_sequence_ids}
            position_maps.append(position_map)
        elif spec['type'] == 'datasheet':
            # Load from datasheet
            datasheet_path = spec['datasheet_path']
            column_name = spec['column_name']
            position_map = load_positions_from_datasheet(datasheet_path, column_name)
            position_maps.append(position_map)
        else:
            # Empty positions (no overlay)
            position_map = {seq_id: [] for seq_id in all_sequence_ids}
            position_maps.append(position_map)

    # Perform stitching
    stitched_results = []

    for seq_id in base_sequences.keys():
        print(f"\nProcessing sequence: {seq_id}")

        # Start with base sequence
        current_sequence = base_sequences[seq_id]
        print(f"  Base sequence: {current_sequence[:50]}{'...' if len(current_sequence) > 50 else ''}")

        # Apply overlays from each source
        for overlay_index in range(1, len(all_sequences)):
            overlay_sequences = all_sequences[overlay_index]
            position_map = position_maps[overlay_index]
            source_name = sequence_sources[overlay_index]['tool_name']

            if seq_id not in overlay_sequences:
                print(f"  Warning: No overlay sequence from {source_name} for {seq_id}")
                continue

            if seq_id not in position_map:
                print(f"  Warning: No position specification for {seq_id}")
                continue

            overlay_sequence = overlay_sequences[seq_id]
            positions = position_map[seq_id]

            if positions:
                print(f"  Overlaying {source_name} at positions: {positions[:10]}{'...' if len(positions) > 10 else ''}")
                current_sequence = stitch_sequences(current_sequence, overlay_sequence, positions)
            else:
                print(f"  No overlay positions for {source_name}")

        print(f"  Final sequence: {current_sequence[:50]}{'...' if len(current_sequence) > 50 else ''}")

        stitched_results.append({
            'id': seq_id,
            'sequence': current_sequence
        })

    # Save results
    if stitched_results:
        df = pd.DataFrame(stitched_results)

        # Create output directory
        output_dir = os.path.dirname(output_csv)
        print(f"Creating output directory: {output_dir}")
        os.makedirs(output_dir, exist_ok=True)

        # Save results
        print(f"Writing stitched sequences to: {output_csv}")
        df.to_csv(output_csv, index=False)

        print(f"\nSequence stitching completed successfully!")
        print(f"Generated {len(stitched_results)} stitched sequences")
        print(f"Results saved to: {output_csv}")

        # Show results summary
        print(f"\nResults summary:")
        print(df)

        # Sequence length statistics
        lengths = [len(result['sequence']) for result in stitched_results]
        print(f"\nSequence length statistics:")
        print(f"  Min: {min(lengths)} residues")
        print(f"  Max: {max(lengths)} residues")
        print(f"  Mean: {sum(lengths)/len(lengths):.1f} residues")

    else:
        raise ValueError("No stitched sequences generated")


def main():
    parser = argparse.ArgumentParser(description='Stitch sequences from multiple sequence generation tools')
    parser.add_argument('--config', required=True, help='JSON config file with stitching parameters')

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
    required_params = ['sequence_sources', 'position_specs', 'output_csv']
    for param in required_params:
        if param not in config_data:
            print(f"Error: Missing required parameter: {param}")
            sys.exit(1)

    try:
        stitch_sequences_from_config(config_data)

    except Exception as e:
        print(f"Error stitching sequences: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()