#!/usr/bin/env python3
# Copyright (c) 2026 Gianluca Quargnali @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# This software is freely available for use, modification, and redistribution.
# If you use this software or any derivative of it in published work,
# you must cite the original author and this repository.
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

# Import ID mapping utilities
# Handle both direct execution and import from other directories
try:
    from id_map_utils import get_mapped_ids, map_table_ids_to_ids
except ImportError:
    from HelpScripts.id_map_utils import get_mapped_ids, map_table_ids_to_ids


def sele_to_list(s: str) -> List[int]:
    """
    Convert a PyMOL-style selection string to a list of 1-indexed positions.

    Args:
        s: Selection string like '10-20' or '10-20+30-40' or '5+10+15'

    Returns:
        List of 1-indexed positions
    """
    def contiguous_sele_to_list(pp):
        if '-' in pp:
            min_val, max_val = pp.split('-')
            return [ri for ri in range(int(min_val), int(max_val) + 1)]
        else:
            return [int(pp)]

    a = []
    if s == "":
        return a
    elif '+' in s:
        plus_parts = s.split('+')
        for pp in plus_parts:
            a.extend(contiguous_sele_to_list(pp))
    else:
        a.extend(contiguous_sele_to_list(s))
    return a


def sele_to_segments(s: str) -> List[Tuple[int, int]]:
    """
    Parse selection string into list of (start, end) tuples for contiguous segments.

    Args:
        s: Selection string like '10-20' or '6-7+9-10+17-18'

    Returns:
        List of (start, end) tuples (1-indexed, inclusive)
    """
    segments = []
    if s == "":
        return segments

    parts = s.split('+')
    for part in parts:
        part = part.strip()
        if '-' in part:
            start, end = map(int, part.split('-'))
            segments.append((start, end))
        else:
            pos = int(part)
            segments.append((pos, pos))

    return segments


def apply_substitutions(template: str, substitutions: Dict[str, str],
                        debug_info: Dict[str, Any] = None) -> str:
    """
    Apply position-to-position substitutions from equal-length sequences.

    For each position in the selection, copy the residue at that position from the
    substitution sequence to the template. Template and substitution must be same length.

    Args:
        template: Base sequence string
        substitutions: Dict mapping position ranges to full replacement sequences.
                       Positions are 1-indexed (PyMOL style).
        debug_info: Optional dict with IDs for debug output.

    Returns:
        New sequence with substitutions applied
    """
    result = list(template)
    template_len = len(template)

    if debug_info:
        print(f"\n{'='*80}")
        print(f"SUBSTITUTIONS: {debug_info.get('template_id', 'unknown')}")
        print(f"  Template ID: {debug_info.get('template_id', 'N/A')}")
        print(f"  Template length: {template_len}")
        print(f"{'='*80}")

    for sub_idx, (pos_range, replacement_seq) in enumerate(substitutions.items()):
        # Warn if lengths mismatch
        if len(replacement_seq) != template_len:
            position_table_key = debug_info.get(f'sub_{sub_idx}_position_table_key') if debug_info else None
            print(f"  WARNING: Substitution sequence length mismatch:")
            print(f"    Template length: {template_len}")
            print(f"    Substitution sequence length: {len(replacement_seq)}")
            print(f"    Template ID: {debug_info.get('template_id', 'N/A') if debug_info else 'N/A'}")
            print(f"    Substitution ID: {debug_info.get(f'sub_{sub_idx}_id', 'N/A') if debug_info else 'N/A'}")
            print(f"    Position table key: {position_table_key if position_table_key else 'N/A'}")
            print(f"    Proceeding with available positions...")

        # Get list of 1-indexed positions to substitute
        positions = sele_to_list(pos_range)

        # Build a mask showing what will be replaced (for debug output)
        mask = ['-'] * template_len

        sub_id = debug_info.get(f'sub_{sub_idx}_id', 'N/A') if debug_info else 'N/A'
        position_table_key = debug_info.get(f'sub_{sub_idx}_position_table_key') if debug_info else None
        print(f"\n  Substitution {sub_idx + 1}: '{pos_range}'")
        print(f"    Substitution source ID: {sub_id}")
        if position_table_key:
            print(f"    Position table key: {position_table_key}")
        print(f"    Positions to replace: {len(positions)} residues")

        # Apply: for each position, copy from replacement_seq to result
        replacement_len = len(replacement_seq)
        for pos in positions:
            idx = pos - 1  # Convert 1-indexed to 0-indexed
            if idx < 0 or idx >= template_len:
                print(f"    WARNING: Position {pos} out of range for template of length {template_len}, skipping")
                continue
            if idx >= replacement_len:
                print(f"    WARNING: Position {pos} out of range for substitution sequence of length {replacement_len}, skipping")
                continue
            mask[idx] = replacement_seq[idx]
            result[idx] = replacement_seq[idx]

        # Show full alignment
        if positions:
            template_str = template
            mask_str = ''.join(mask)

            # Create match line
            match_line = ''.join(
                '|' if mask[i] != '-' and template[i] == mask[i] else
                '*' if mask[i] != '-' else ' '
                for i in range(template_len)
            )

            # Create numbering lines (mark every 5th position)
            fives_line = []
            tens_line = []
            hundreds_line = []
            for i in range(template_len):
                pos = i + 1  # 1-indexed position
                if pos % 5 == 0:
                    fives_line.append(str(pos % 10))
                    tens_line.append(str((pos // 10) % 10) if pos >= 10 else ' ')
                    hundreds_line.append(str((pos // 100) % 10) if pos >= 100 else ' ')
                else:
                    fives_line.append(' ')
                    tens_line.append(' ')
                    hundreds_line.append(' ')

            print(f"\n    Alignment (1-{template_len}):")
            if template_len >= 100:
                print(f"              {''.join(hundreds_line)}")
            if template_len >= 10:
                print(f"              {''.join(tens_line)}")
            print(f"              {''.join(fives_line)}")
            print(f"    Template: {template_str}")
            print(f"              {match_line}")
            print(f"    Substitut:{mask_str}")

    print(f"\n  Result length after substitutions: {len(result)}")
    print(f"{'='*80}\n")

    return ''.join(result)


def apply_indels(template: str, indels: Dict[str, str],
                 debug_info: Dict[str, Any] = None) -> str:
    """
    Apply segment replacements (insertions/deletions) to a sequence.

    Each contiguous segment in the selection is replaced with the given sequence.
    For example, "6-7+9-10+17-18": "GP" replaces segments 6-7, 9-10, and 17-18 each with "GP".

    Indels are processed from end to start to preserve position indices.

    Args:
        template: Base sequence string
        indels: Dict mapping position ranges to replacement sequences.
                Each contiguous segment is replaced with the full replacement.
        debug_info: Optional dict with IDs for debug output.

    Returns:
        New sequence with indels applied
    """
    result = template

    if debug_info:
        print(f"\n{'='*80}")
        print(f"INDELS: {debug_info.get('template_id', 'unknown')}")
        print(f"  Template ID: {debug_info.get('template_id', 'N/A')}")
        print(f"  Initial length: {len(template)}")
        print(f"{'='*80}")

    for indel_idx, (pos_range, replacement) in enumerate(indels.items()):
        segments = sele_to_segments(pos_range)

        indel_id = debug_info.get(f'indel_{indel_idx}_id', 'N/A') if debug_info else 'N/A'
        position_table_key = debug_info.get(f'indel_{indel_idx}_position_table_key') if debug_info else None
        print(f"\n  Indel {indel_idx + 1}: '{pos_range}'")
        print(f"    Indel source ID: {indel_id}")
        if position_table_key:
            print(f"    Position table key: {position_table_key}")
        print(f"    Replacement sequence: '{replacement}' ({len(replacement)} chars)")
        print(f"    Segments to replace: {segments}")

        # Process segments from end to start to preserve indices
        for start, end in sorted(segments, reverse=True):
            start_idx = start - 1  # Convert to 0-indexed
            end_idx = end          # end is inclusive, so end_idx = end for slicing

            old_segment = result[start_idx:end_idx]
            result = result[:start_idx] + replacement + result[end_idx:]

            print(f"      Replaced [{start}-{end}] '{old_segment}' -> '{replacement}'")

    print(f"\n  Result length after indels: {len(result)}")
    print(f"{'='*80}\n")

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
        # Include empty strings - they result in empty position lists (no substitution)
        # This is different from missing the ID entirely (which would skip the sequence)
        if pd.notna(pos_value):
            positions_map[seq_id] = str(pos_value)
        else:
            positions_map[seq_id] = ''

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


def extract_base_id(seq_id: str, id_map: Dict[str, str] = None) -> str:
    """
    Extract the most stripped base ID from a sequence ID.

    Args:
        seq_id: The sequence ID to strip
        id_map: ID mapping pattern (default: {"*": "*_<N>"})

    Returns:
        The base ID after recursively stripping suffixes
    """
    if id_map is None:
        id_map = {"*": "*_<N>"}
    bases = map_table_ids_to_ids(seq_id, id_map)
    return bases[-1] if bases else seq_id


def group_sequences_by_template(
    sequences: Dict[str, str],
    template_ids: List[str],
    id_map: Dict[str, str]
) -> Dict[str, Dict[str, str]]:
    """
    Group sequences by matching them to template IDs using id_map pattern.

    Uses get_mapped_ids from id_map_utils for flexible matching that supports:
    - Exact match (source == target)
    - Child match (target = source + suffix)
    - Parent match (source = target + suffix)
    - Sibling match (common ancestor)

    Args:
        sequences: Dict of sequence ID -> sequence
        template_ids: List of template IDs to match against
        id_map: ID mapping pattern (e.g., {"*": "*_<N>"})

    Returns:
        Dict mapping template_id -> {matched_seq_id: sequence}
    """
    seq_ids = list(sequences.keys())

    # Use get_mapped_ids to find matches for each template
    # Here template_ids are sources, seq_ids are targets
    matches = get_mapped_ids(template_ids, seq_ids, id_map, unique=False)

    grouped = {}
    for template_id in template_ids:
        matched_seq_ids = matches.get(template_id, [])
        if matched_seq_ids:
            grouped[template_id] = {
                seq_id: sequences[seq_id]
                for seq_id in matched_seq_ids
            }

    return grouped


def stitch_sequences_from_config(config_data: Dict[str, Any]) -> None:
    """Stitch sequences based on configuration."""
    template_info = config_data['template']
    substitutions_info = config_data.get('substitutions', {})
    indels_info = config_data.get('indels', {})
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

    # Load indel options and positions for each region
    print(f"\nLoading indels for {len(indels_info)} regions...")
    indels = []
    for key_str, indel_info in indels_info.items():
        pos_key_info = indel_info['position_key']
        seq_info = indel_info['sequences']

        # Load sequences for this indel
        sequences = load_sequences_from_info(seq_info)
        print(f"  {key_str}: {len(sequences)} sequence options")

        indels.append({
            'key_str': key_str,
            'position_key': pos_key_info,
            'sequences': sequences
        })

    # Determine stitching strategy
    all_positions_fixed = all(
        sub['position_key']['type'] == 'fixed'
        for sub in substitutions
    ) if substitutions else True
    all_subs_raw = all(
        sub_info['sequences'].get("type") in ("raw", "raw_list")
        for sub_info in substitutions_info.values()
    ) if substitutions_info else True

    all_indel_positions_fixed = all(
        indel['position_key']['type'] == 'fixed'
        for indel in indels
    ) if indels else True
    all_indels_raw = all(
        indel_info['sequences'].get("type") in ("raw", "raw_list")
        for indel_info in indels_info.values()
    ) if indels_info else True

    if all_positions_fixed and all_subs_raw and all_indel_positions_fixed and all_indels_raw:
        # Raw substitutions/indels with fixed positions - apply to all templates
        results = stitch_with_raw_operations(template_sequences, substitutions, indels, id_map)
    else:
        # Complex case: match sequences by base ID
        results = stitch_matched_sequences(template_sequences, substitutions, indels, id_map)

    # Remove duplicates if requested
    remove_duplicates = config_data.get('remove_duplicates', True)
    if remove_duplicates and results:
        seen_sequences = set()
        unique_results = []
        duplicates_removed = 0
        for r in results:
            if r['sequence'] not in seen_sequences:
                seen_sequences.add(r['sequence'])
                unique_results.append(r)
            else:
                print(f"Skipped duplicate: {r['id']}")
                duplicates_removed += 1
        if duplicates_removed > 0:
            print(f"Removed {duplicates_removed} duplicate sequences ({len(results)} -> {len(unique_results)})")
        results = unique_results

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


def stitch_with_raw_operations(
    template_sequences: Dict[str, str],
    substitutions: List[Dict],
    indels: List[Dict],
    id_map: Dict[str, str]
) -> List[Dict[str, Any]]:
    """Stitch templates with raw substitutions and indels at fixed positions using Cartesian product."""
    results = []

    # Build option lists for substitutions
    sub_option_lists = []
    sub_option_ids = []
    sub_pos_ranges = []
    for sub in substitutions:
        sub_pos_ranges.append(sub['position_key']['positions'])
        seq_items = list(sub['sequences'].items())
        sub_option_ids.append([item[0] for item in seq_items])
        sub_option_lists.append([item[1] for item in seq_items])

    # Build option lists for indels
    indel_option_lists = []
    indel_option_ids = []
    indel_pos_ranges = []
    for indel in indels:
        indel_pos_ranges.append(indel['position_key']['positions'])
        seq_items = list(indel['sequences'].items())
        indel_option_ids.append([item[0] for item in seq_items])
        indel_option_lists.append([item[1] for item in seq_items])

    for template_id, template_seq in template_sequences.items():
        base_id = extract_base_id(template_id, id_map)

        # Generate index ranges for substitutions and indels
        sub_index_ranges = [range(1, len(opts) + 1) for opts in sub_option_lists] if sub_option_lists else [range(1, 2)]
        indel_index_ranges = [range(1, len(opts) + 1) for opts in indel_option_lists] if indel_option_lists else [range(1, 2)]

        # Combine all index ranges for Cartesian product
        all_index_ranges = sub_index_ranges + indel_index_ranges

        for index_combo in product(*all_index_ranges):
            sub_indices = index_combo[:len(sub_option_lists)] if sub_option_lists else ()
            indel_indices = index_combo[len(sub_option_lists):] if indel_option_lists else ()

            current_seq = template_seq
            debug_info = {'template_id': template_id}

            # Apply substitutions first (position-to-position, same length)
            if sub_option_lists:
                subs = {}
                for i, idx in enumerate(sub_indices):
                    subs[sub_pos_ranges[i]] = sub_option_lists[i][idx - 1]
                    debug_info[f'sub_{i}_id'] = sub_option_ids[i][idx - 1]
                current_seq = apply_substitutions(current_seq, subs, debug_info=debug_info)

            # Apply indels second (segment replacement, can change length)
            if indel_option_lists:
                indel_ops = {}
                for i, idx in enumerate(indel_indices):
                    indel_ops[indel_pos_ranges[i]] = indel_option_lists[i][idx - 1]
                    debug_info[f'indel_{i}_id'] = indel_option_ids[i][idx - 1]
                current_seq = apply_indels(current_seq, indel_ops, debug_info=debug_info)

            # Build output ID suffix
            all_indices = list(sub_indices) + list(indel_indices)
            if all_indices:
                suffix = "_".join(str(idx) for idx in all_indices)
                output_id = f"{base_id}_{suffix}"
            else:
                output_id = template_id

            results.append({'id': output_id, 'sequence': current_seq})

    print(f"\nGenerated {len(results)} stitched sequences")
    return results


def stitch_matched_sequences(
    template_sequences: Dict[str, str],
    substitutions: List[Dict],
    indels: List[Dict],
    id_map: Dict[str, str]
) -> List[Dict[str, Any]]:
    """Stitch sequences by matching base IDs, with support for table-based positions."""
    results = []

    # Get all template IDs for matching
    all_template_ids = list(template_sequences.keys())

    # Group substitution sequences by template ID using id_map and load table positions if needed
    sub_grouped = []
    sub_position_maps = []
    for sub in substitutions:
        grouped = group_sequences_by_template(sub['sequences'], all_template_ids, id_map)
        sub_grouped.append(grouped)

        pos_key = sub['position_key']
        if pos_key['type'] == 'table':
            pos_map = load_positions_from_table(pos_key['table_path'], pos_key['column'])
            print(f"  Loaded substitution positions from table for {len(pos_map)} structures")
        else:
            pos_map = None
        sub_position_maps.append((pos_key, pos_map))

    # Group indel sequences by template ID using id_map and load table positions if needed
    indel_grouped = []
    indel_position_maps = []
    for indel in indels:
        grouped = group_sequences_by_template(indel['sequences'], all_template_ids, id_map)
        indel_grouped.append(grouped)

        pos_key = indel['position_key']
        if pos_key['type'] == 'table':
            pos_map = load_positions_from_table(pos_key['table_path'], pos_key['column'])
            print(f"  Loaded indel positions from table for {len(pos_map)} structures")
        else:
            pos_map = None
        indel_position_maps.append((pos_key, pos_map))

    print(f"\nProcessing {len(template_sequences)} template sequences...")

    # Process each template sequence directly
    for template_id, template_seq in sorted(template_sequences.items()):
        template_combo_count = 0
        base_id = extract_base_id(template_id, id_map)
        print(f"\n  Template: {template_id}")

        # Check if all substitutions have sequences for this template ID
        sub_seqs_for_template = []
        skip = False
        for i, grouped in enumerate(sub_grouped):
            if template_id in grouped:
                sub_seqs_for_template.append(grouped[template_id])
            else:
                print(f"    Skipping: no sequences in substitution {i} for template '{template_id}'")
                skip = True
                break

        if skip:
            continue

        # Check if all indels have sequences for this template ID
        indel_seqs_for_template = []
        for i, grouped in enumerate(indel_grouped):
            if template_id in grouped:
                indel_seqs_for_template.append(grouped[template_id])
            else:
                print(f"    Skipping: no sequences in indel {i} for template '{template_id}'")
                skip = True
                break

        if skip:
            continue

        # Get substitution positions for this template using id_map matching
        sub_positions_for_template = []
        sub_position_table_keys = []
        for pos_key, pos_map in sub_position_maps:
            if pos_key['type'] == 'fixed':
                sub_positions_for_template.append(pos_key['positions'])
                sub_position_table_keys.append(None)
            else:
                # Use get_mapped_ids for flexible matching
                pos_map_ids = list(pos_map.keys())
                match_result = get_mapped_ids([template_id], pos_map_ids, id_map, unique=True)
                matched_key = match_result.get(template_id)
                if matched_key:
                    sub_positions_for_template.append(pos_map[matched_key])
                    sub_position_table_keys.append(matched_key)
                else:
                    print(f"    Skipping: no substitution position data matching '{template_id}'")
                    skip = True
                    break

        if skip:
            continue

        # Get indel positions for this template using id_map matching
        indel_positions_for_template = []
        indel_position_table_keys = []
        for pos_key, pos_map in indel_position_maps:
            if pos_key['type'] == 'fixed':
                indel_positions_for_template.append(pos_key['positions'])
                indel_position_table_keys.append(None)
            else:
                # Use get_mapped_ids for flexible matching
                pos_map_ids = list(pos_map.keys())
                match_result = get_mapped_ids([template_id], pos_map_ids, id_map, unique=True)
                matched_key = match_result.get(template_id)
                if matched_key:
                    indel_positions_for_template.append(pos_map[matched_key])
                    indel_position_table_keys.append(matched_key)
                else:
                    print(f"    Skipping: no indel position data matching '{template_id}'")
                    skip = True
                    break

        if skip:
            continue

        # Generate combinations from all substitution and indel options
        # Build option lists for substitutions
        sub_option_lists = []
        sub_option_ids = []
        for sub_seqs in sub_seqs_for_template:
            seq_items = list(sub_seqs.items())
            sub_option_ids.append([item[0] for item in seq_items])
            sub_option_lists.append([item[1] for item in seq_items])

        # Build option lists for indels
        indel_option_lists = []
        indel_option_ids = []
        for indel_seqs in indel_seqs_for_template:
            seq_items = list(indel_seqs.items())
            indel_option_ids.append([item[0] for item in seq_items])
            indel_option_lists.append([item[1] for item in seq_items])

        # Generate index ranges for Cartesian product
        sub_index_ranges = [range(len(opts)) for opts in sub_option_lists] if sub_option_lists else [range(1)]
        indel_index_ranges = [range(len(opts)) for opts in indel_option_lists] if indel_option_lists else [range(1)]
        all_index_ranges = sub_index_ranges + indel_index_ranges

        for index_combo in product(*all_index_ranges):
            sub_indices = index_combo[:len(sub_option_lists)] if sub_option_lists else ()
            indel_indices = index_combo[len(sub_option_lists):] if indel_option_lists else ()

            current_seq = template_seq
            debug_info = {'template_id': template_id}

            # Apply substitutions
            if sub_option_lists:
                subs = {}
                for i, idx in enumerate(sub_indices):
                    positions_str = sub_positions_for_template[i]
                    subs[positions_str] = sub_option_lists[i][idx]
                    debug_info[f'sub_{i}_id'] = sub_option_ids[i][idx]
                    debug_info[f'sub_{i}_position_table_key'] = sub_position_table_keys[i]
                current_seq = apply_substitutions(current_seq, subs, debug_info=debug_info)

            # Apply indels
            if indel_option_lists:
                indel_ops = {}
                for i, idx in enumerate(indel_indices):
                    positions_str = indel_positions_for_template[i]
                    indel_ops[positions_str] = indel_option_lists[i][idx]
                    debug_info[f'indel_{i}_id'] = indel_option_ids[i][idx]
                    debug_info[f'indel_{i}_position_table_key'] = indel_position_table_keys[i]
                current_seq = apply_indels(current_seq, indel_ops, debug_info=debug_info)

            # Build output ID
            all_indices = list(sub_indices) + list(indel_indices)
            if all_indices and (len(sub_option_lists) > 0 or len(indel_option_lists) > 0):
                # Use 1-based indices for output
                suffix = "_".join(str(idx + 1) for idx in all_indices)
                output_id = f"{template_id}_{suffix}"
            else:
                output_id = template_id

            results.append({'id': output_id, 'sequence': current_seq})
            template_combo_count += 1

        print(f"    Generated {template_combo_count} combinations")

    print(f"\nTotal: {len(results)} stitched sequences")
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
