#!/usr/bin/env python3
"""
Runtime helper script for StitchSequences tool.

This script stitches sequences by replacing specific regions of a template
sequence with alternative segments, generating all Cartesian product combinations.

Supports:
- Raw sequence strings as template or substitutions
- ToolOutput CSV files as template or substitutions
- Fixed position strings (e.g., "11-19")
- Dynamic positions from table columns
- Variable length segment substitutions
"""

import os
import sys
import argparse
import json
import re
import pandas as pd
from typing import Dict, List, Any, Tuple
from itertools import product


def parse_position_range(pos_range: str) -> List[Tuple[int, int]]:
    """
    Parse position range string into list of (start, end) tuples.

    Args:
        pos_range: Position string like '10-20' or '10-20+30-40'

    Returns:
        List of (start, end) tuples, 1-indexed inclusive
    """
    ranges = []
    parts = pos_range.split('+')

    for part in parts:
        part = part.strip()
        if '-' in part:
            start, end = map(int, part.split('-'))
            ranges.append((start, end))
        else:
            pos = int(part)
            ranges.append((pos, pos))

    return sorted(ranges, key=lambda x: x[0])


def substitute_segments(template: str, substitutions: Dict[str, str]) -> str:
    """
    Replace segments of template sequence at specified positions.

    Args:
        template: Base sequence string
        substitutions: Dict mapping position ranges to replacement sequences

    Returns:
        New sequence with substitutions applied
    """
    # Parse all position ranges and sort by start position (descending)
    replacements = []
    for pos_range, replacement in substitutions.items():
        ranges = parse_position_range(pos_range)
        overall_start = ranges[0][0]
        overall_end = ranges[-1][1]
        replacements.append((overall_start, overall_end, ranges, replacement))

    # Sort by start position descending (process from end first)
    replacements.sort(key=lambda x: x[0], reverse=True)

    result = template
    for overall_start, overall_end, ranges, replacement in replacements:
        start_idx = overall_start - 1
        end_idx = overall_end
        result = result[:start_idx] + replacement + result[end_idx:]

    return result


def load_sequences_from_csv(csv_path: str) -> Dict[str, str]:
    """Load sequences from CSV file."""
    if not os.path.exists(csv_path):
        raise FileNotFoundError(f"Sequences file not found: {csv_path}")

    df = pd.read_csv(csv_path)

    if 'id' not in df.columns or 'sequence' not in df.columns:
        raise ValueError(f"Expected 'id' and 'sequence' columns in {csv_path}. Found: {list(df.columns)}")

    sequences = {}
    for _, row in df.iterrows():
        seq_id = str(row['id'])
        sequence = str(row['sequence'])
        sequences[seq_id] = sequence

    return sequences


def load_positions_from_table(table_path: str, column_name: str) -> Dict[str, str]:
    """
    Load position specifications from table CSV file.

    Returns dict mapping sequence IDs to position strings.
    """
    if not os.path.exists(table_path):
        raise FileNotFoundError(f"Table file not found: {table_path}")

    df = pd.read_csv(table_path)
    if column_name not in df.columns:
        raise ValueError(f"Column '{column_name}' not found. Available: {list(df.columns)}")

    id_column = df.columns[0]
    positions_map = {}

    for _, row in df.iterrows():
        seq_id = str(row[id_column])
        pos_value = row[column_name]
        if pd.notna(pos_value) and pos_value != '':
            positions_map[seq_id] = str(pos_value)

    return positions_map


def load_sequences_from_info(info: Dict[str, Any]) -> Dict[str, str]:
    """Load sequences from source info dict."""
    source_type = info.get("type")

    if source_type == "raw":
        sequences = info.get("sequences", [])
        return {f"seq_{i+1}": seq for i, seq in enumerate(sequences)}

    elif source_type == "raw_list":
        sequences = info.get("sequences", [])
        return {f"seq_{i+1}": seq for i, seq in enumerate(sequences)}

    elif source_type == "tool_output":
        csv_path = info.get("sequences_file")
        if not csv_path:
            raise ValueError("tool_output source missing sequences_file")
        return load_sequences_from_csv(csv_path)

    else:
        raise ValueError(f"Unknown source type: {source_type}")


def extract_base_id(seq_id: str) -> str:
    """Extract base ID by stripping trailing _N suffix."""
    match = re.match(r'^(.+)_\d+$', seq_id)
    return match.group(1) if match else seq_id


def group_by_base_id(sequences: Dict[str, str]) -> Dict[str, Dict[str, str]]:
    """Group sequences by their base ID."""
    grouped = {}
    for seq_id, sequence in sequences.items():
        base_id = extract_base_id(seq_id)
        if base_id not in grouped:
            grouped[base_id] = {}
        grouped[base_id][seq_id] = sequence
    return grouped


def stitch_sequences_from_config(config_data: Dict[str, Any]) -> None:
    """Stitch sequences based on configuration."""
    template_info = config_data['template']
    substitutions_info = config_data.get('substitutions', {})
    id_map = config_data.get('id_map', {"*": "*_<N>"})
    output_csv = config_data['output_csv']

    print("Loading template sequences...")
    template_sequences = load_sequences_from_info(template_info)
    print(f"  Loaded {len(template_sequences)} template sequences")

    # Load substitution options and positions for each region
    print(f"\nLoading substitutions for {len(substitutions_info)} regions...")
    substitutions = []
    for key_str, sub_info in substitutions_info.items():
        pos_key_info = sub_info['position_key']
        seq_info = sub_info['sequences']

        # Load sequences for this substitution
        sequences = load_sequences_from_info(seq_info)
        print(f"  {key_str}: {len(sequences)} sequence options")

        substitutions.append({
            'key_str': key_str,
            'position_key': pos_key_info,
            'sequences': sequences
        })

    # Determine stitching strategy
    template_is_raw = template_info.get("type") in ("raw", "raw_list")
    all_positions_fixed = all(
        sub['position_key']['type'] == 'fixed'
        for sub in substitutions
    )
    all_subs_raw = all(
        sub_info['sequences'].get("type") in ("raw", "raw_list")
        for sub_info in substitutions_info.values()
    )

    if all_positions_fixed and all_subs_raw:
        # Raw substitutions with fixed positions - apply to all templates
        results = stitch_with_raw_substitutions(template_sequences, substitutions)
    else:
        # Complex case: match sequences by base ID
        results = stitch_matched_sequences(template_sequences, substitutions, id_map)

    # Save results
    if results:
        df = pd.DataFrame(results)
        output_dir = os.path.dirname(output_csv)
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)
        df.to_csv(output_csv, index=False)

        print(f"\nSequence stitching completed!")
        print(f"Generated {len(results)} stitched sequences")
        print(f"Results saved to: {output_csv}")

        lengths = [len(r['sequence']) for r in results]
        print(f"\nSequence lengths: min={min(lengths)}, max={max(lengths)}, mean={sum(lengths)/len(lengths):.1f}")
    else:
        raise ValueError("No stitched sequences generated")


def stitch_with_raw_substitutions(
    template_sequences: Dict[str, str],
    substitutions: List[Dict]
) -> List[Dict[str, Any]]:
    """Stitch templates with raw substitutions at fixed positions using Cartesian product."""
    results = []

    # Build option lists for Cartesian product (just the sequences, ignore their IDs)
    option_lists = []
    pos_ranges = []
    for sub in substitutions:
        pos_ranges.append(sub['position_key']['positions'])
        # Only use the sequences, not the IDs (raw sequences have generic IDs)
        option_lists.append(list(sub['sequences'].values()))

    for template_id, template_seq in template_sequences.items():
        if not option_lists:
            results.append({'id': template_id, 'sequence': template_seq})
            continue

        # Get base ID from template
        base_id = extract_base_id(template_id)

        combo_idx = 0
        for combo in product(*option_lists):
            combo_idx += 1
            subs = {}
            for i, opt_seq in enumerate(combo):
                subs[pos_ranges[i]] = opt_seq

            stitched = substitute_segments(template_seq, subs)
            output_id = f"{base_id}_{combo_idx}"
            results.append({'id': output_id, 'sequence': stitched})

    print(f"\nGenerated {len(results)} stitched sequences")
    return results


def stitch_matched_sequences(
    template_sequences: Dict[str, str],
    substitutions: List[Dict],
    id_map: Dict[str, str]
) -> List[Dict[str, Any]]:
    """Stitch sequences by matching base IDs, with support for table-based positions."""
    results = []

    # Group templates by base ID
    template_grouped = group_by_base_id(template_sequences)

    # Group substitution sequences by base ID and load table positions if needed
    sub_grouped = []
    position_maps = []
    for sub in substitutions:
        grouped = group_by_base_id(sub['sequences'])
        sub_grouped.append(grouped)

        # Load position map for table-based positions
        pos_key = sub['position_key']
        if pos_key['type'] == 'table':
            pos_map = load_positions_from_table(pos_key['table_path'], pos_key['column'])
            print(f"  Loaded positions from table for {len(pos_map)} structures")
        else:
            # Fixed positions - same for all
            pos_map = None
        position_maps.append((pos_key, pos_map))

    print(f"\nProcessing {len(template_grouped)} base structure IDs...")

    for base_id in sorted(template_grouped.keys()):
        template_seqs = template_grouped[base_id]
        print(f"\n  Base ID: {base_id} ({len(template_seqs)} template sequences)")

        # Check if all substitutions have sequences for this base ID
        sub_seqs_for_base = []
        skip = False
        for i, grouped in enumerate(sub_grouped):
            if base_id in grouped:
                sub_seqs_for_base.append(grouped[base_id])
            else:
                print(f"    Skipping: no sequences in substitution {i}")
                skip = True
                break

        if skip:
            continue

        # Get positions for this base ID
        positions_for_base = []
        for pos_key, pos_map in position_maps:
            if pos_key['type'] == 'fixed':
                positions_for_base.append(pos_key['positions'])
            else:
                # Table-based: look up by base_id
                if base_id in pos_map:
                    positions_for_base.append(pos_map[base_id])
                else:
                    print(f"    Skipping: no position data for {base_id}")
                    skip = True
                    break

        if skip:
            continue

        # Build option lists for Cartesian product
        option_lists = [list(template_seqs.items())]
        for sub_seqs in sub_seqs_for_base:
            option_lists.append(list(sub_seqs.items()))

        # Generate all combinations
        combo_idx = 0
        for combo in product(*option_lists):
            combo_idx += 1
            template_id, template_seq = combo[0]

            subs = {}
            for i, (opt_id, opt_seq) in enumerate(combo[1:]):
                subs[positions_for_base[i]] = opt_seq

            stitched = substitute_segments(template_seq, subs)
            output_id = f"{base_id}_{combo_idx}"
            results.append({'id': output_id, 'sequence': stitched})

        print(f"    Generated {combo_idx} combinations")

    return results


def main():
    parser = argparse.ArgumentParser(description='Stitch sequences with segment substitutions')
    parser.add_argument('--config', required=True, help='JSON config file')
    args = parser.parse_args()

    if not os.path.exists(args.config):
        print(f"Error: Config file not found: {args.config}")
        sys.exit(1)

    try:
        with open(args.config, 'r') as f:
            config_data = json.load(f)
    except Exception as e:
        print(f"Error loading config: {e}")
        sys.exit(1)

    if 'template' not in config_data:
        print("Error: Missing required parameter: template")
        sys.exit(1)

    if 'output_csv' not in config_data:
        print("Error: Missing required parameter: output_csv")
        sys.exit(1)

    try:
        stitch_sequences_from_config(config_data)
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
