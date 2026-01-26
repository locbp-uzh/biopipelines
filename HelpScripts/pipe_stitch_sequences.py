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
        List of (start, end) tuples (PDB residue numbers, inclusive)
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


def build_residue_mapping(pdb_path: str) -> Dict[int, int]:
    """
    Build a mapping from PDB residue numbers to sequence string indices.

    Args:
        pdb_path: Path to PDB file

    Returns:
        Dict mapping PDB residue number -> string index (0-based)
    """
    from pdb_parser import parse_pdb_file

    # Standard amino acid codes
    aa_codes = {
        'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
        'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'
    }

    atoms = parse_pdb_file(pdb_path)

    # Collect unique residues (chain, res_num) for standard amino acids
    residues = {}
    for atom in atoms:
        if atom.res_name in aa_codes:
            key = (atom.chain, atom.res_num)
            if key not in residues:
                residues[key] = atom.res_num

    # Sort by chain and residue number (same order as get_protein_sequence)
    sorted_residues = sorted(residues.keys(), key=lambda x: (x[0], x[1]))

    # Build mapping: PDB residue number -> string index
    # Note: This concatenates all chains, matching get_protein_sequence behavior
    mapping = {}
    for idx, (chain, res_num) in enumerate(sorted_residues):
        mapping[res_num] = idx

    return mapping


def substitute_segments(template: str, substitutions: Dict[str, str],
                        residue_mapping: Dict[int, int] = None,
                        debug: bool = False,
                        debug_info: Dict[str, Any] = None) -> str:
    """
    Replace segments of template sequence at specified positions.

    Args:
        template: Base sequence string
        substitutions: Dict mapping position ranges to replacement sequences
        residue_mapping: Optional dict mapping PDB residue numbers to string indices.
                        If provided, positions are interpreted as PDB residue numbers.
                        If None, positions are interpreted as 1-indexed sequence positions.
        debug: If True, print alignment visualization for each substitution.
        debug_info: Optional dict with IDs for debug output.

    Returns:
        New sequence with substitutions applied
    """
    # Convert template to list for character-by-character replacement
    result = list(template)
    template_len = len(template)

    if debug and debug_info:
        print(f"\n{'='*80}")
        print(f"STITCH DEBUG: {debug_info.get('template_id', 'unknown')}")
        print(f"  Template ID: {debug_info.get('template_id', 'N/A')}")
        print(f"  Template length: {template_len}")
        if residue_mapping:
            print(f"  Residue mapping: PDB {min(residue_mapping.keys())}-{max(residue_mapping.keys())} -> idx 0-{len(residue_mapping)-1}")
        print(f"{'='*80}")

    # Process each substitution
    for sub_idx, (pos_range, replacement) in enumerate(substitutions.items()):
        ranges = parse_position_range(pos_range)

        # Build list of all positions (PDB residue numbers) to replace
        positions = []
        for start, end in ranges:
            positions.extend(range(start, end + 1))

        # Map PDB residue numbers to string indices
        if residue_mapping is not None:
            # Use PDB residue numbers via mapping
            indices = []
            for pos in positions:
                if pos not in residue_mapping:
                    raise ValueError(f"PDB residue {pos} not found in structure. "
                                   f"Available residue numbers: {min(residue_mapping.keys())}-{max(residue_mapping.keys())}")
                indices.append(residue_mapping[pos])
        else:
            # 1-indexed positions -> 0-indexed string indices
            indices = [pos - 1 for pos in positions]

        # Build masked replacement sequence: same length as template, with '-' except at substitution positions
        masked_replacement = ['-'] * template_len

        # Validate lengths match before building masked sequence
        if len(indices) != len(replacement):
            raise ValueError(
                f"Length mismatch for substitution at '{pos_range}':\n"
                f"  Selection spans {len(indices)} positions (indices {indices[0]}-{indices[-1]})\n"
                f"  But replacement sequence has {len(replacement)} characters\n"
                f"  PDB positions: {positions[0]}-{positions[-1]}\n"
                f"  Template ID: {debug_info.get('template_id', 'N/A') if debug_info else 'N/A'}\n"
                f"  Substitution ID: {debug_info.get(f'sub_{sub_idx}_id', 'N/A') if debug_info else 'N/A'}"
            )

        # Place replacement characters at the correct template indices
        for i, idx in enumerate(indices):
            if idx < 0 or idx >= template_len:
                raise ValueError(f"Index {idx} out of range for template of length {template_len}")
            masked_replacement[idx] = replacement[i]

        if debug:
            sub_id = debug_info.get(f'sub_{sub_idx}_id', 'N/A') if debug_info else 'N/A'
            print(f"\n  Substitution {sub_idx + 1}: '{pos_range}'")
            print(f"    Substitution source ID: {sub_id}")
            print(f"    PDB positions: {positions[0]}-{positions[-1]} ({len(positions)} residues)")
            print(f"    Template indices: {indices[0]}-{indices[-1]} ({len(indices)} positions)")
            print(f"    Replacement length: {len(replacement)}")

            # Show alignment around the substitution region with context
            context = 10
            start_idx = max(0, indices[0] - context)
            end_idx = min(template_len, indices[-1] + context + 1)

            template_slice = template[start_idx:end_idx]
            masked_slice = ''.join(masked_replacement[start_idx:end_idx])

            # Create match line
            match_line = ''.join(
                '|' if masked_replacement[i] != '-' and template[i] == masked_replacement[i] else
                '*' if masked_replacement[i] != '-' else ' '
                for i in range(start_idx, end_idx)
            )

            print(f"\n    Alignment (positions {start_idx}-{end_idx-1}, 0-indexed):")
            print(f"    Template:    {template_slice}")
            print(f"    Match:       {match_line}")
            print(f"    Replacement: {masked_slice}")

        # Apply substitution: replace template positions with non-'-' characters from masked sequence
        for idx in indices:
            result[idx] = masked_replacement[idx]

    if debug:
        print(f"\n  Final stitched sequence length: {len(result)}")
        print(f"{'='*80}\n")

    return ''.join(result)


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


def load_residue_mapping_from_info(info: Dict[str, Any]) -> Dict[str, Dict[int, int]]:
    """
    Load residue mappings from template info if structure files are available.

    Args:
        info: Template info dictionary

    Returns:
        Dict mapping sequence ID to residue mapping (PDB res num -> string index).
        Empty dict if no structure files available.
    """
    if info.get("type") != "tool_output":
        return {}

    structure_files = info.get("structure_files", [])
    structure_ids = info.get("structure_ids", [])

    if not structure_files:
        return {}

    mappings = {}
    for i, struct_file in enumerate(structure_files):
        if not os.path.exists(struct_file):
            print(f"  Warning: Structure file not found: {struct_file}")
            continue

        # Get ID for this structure
        struct_id = structure_ids[i] if i < len(structure_ids) else f"struct_{i+1}"

        try:
            mapping = build_residue_mapping(struct_file)
            mappings[struct_id] = mapping
            first_res = min(mapping.keys())
            last_res = max(mapping.keys())
            print(f"  Built residue mapping for {struct_id}: residues {first_res}-{last_res} -> indices 0-{len(mapping)-1}")
        except Exception as e:
            print(f"  Warning: Could not build residue mapping for {struct_file}: {e}")

    return mappings


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
    debug = config_data.get('debug', False)

    print("Loading template sequences...")
    template_sequences = load_sequences_from_info(template_info)
    print(f"  Loaded {len(template_sequences)} template sequences")

    # Load residue mappings if structure files are available
    print("Loading residue mappings from structure files...")
    residue_mappings = load_residue_mapping_from_info(template_info)
    if residue_mappings:
        print(f"  Loaded residue mappings for {len(residue_mappings)} structures")
    else:
        print("  No structure files available, using 1-indexed positions")

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
        results = stitch_with_raw_substitutions(template_sequences, substitutions, residue_mappings, debug=debug)
    else:
        # Complex case: match sequences by base ID
        results = stitch_matched_sequences(template_sequences, substitutions, id_map, residue_mappings, debug=debug)

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
    substitutions: List[Dict],
    residue_mappings: Dict[str, Dict[int, int]] = None,
    debug: bool = False
) -> List[Dict[str, Any]]:
    """Stitch templates with raw substitutions at fixed positions using Cartesian product."""
    results = []
    residue_mappings = residue_mappings or {}

    # Build option lists for Cartesian product (just the sequences, ignore their IDs)
    option_lists = []
    option_ids = []
    pos_ranges = []
    for sub in substitutions:
        pos_ranges.append(sub['position_key']['positions'])
        # Store both sequences and their IDs for debug output
        seq_items = list(sub['sequences'].items())
        option_ids.append([item[0] for item in seq_items])
        option_lists.append([item[1] for item in seq_items])

    for template_id, template_seq in template_sequences.items():
        if not option_lists:
            results.append({'id': template_id, 'sequence': template_seq})
            continue

        # Get base ID from template
        base_id = extract_base_id(template_id)

        # Get residue mapping for this template (try exact match, then base_id)
        residue_mapping = residue_mappings.get(template_id) or residue_mappings.get(base_id)

        # Generate index ranges for each substitution (1-indexed)
        index_ranges = [range(1, len(opts) + 1) for opts in option_lists]

        for index_combo in product(*index_ranges):
            # Build substitutions using the indices
            subs = {}
            debug_info = {'template_id': template_id}
            for i, idx in enumerate(index_combo):
                subs[pos_ranges[i]] = option_lists[i][idx - 1]  # Convert to 0-indexed for list access
                debug_info[f'sub_{i}_id'] = option_ids[i][idx - 1]
                debug_info[f'sub_{i}_positions'] = pos_ranges[i]

            stitched = substitute_segments(template_seq, subs, residue_mapping, debug=debug, debug_info=debug_info)
            # Format: <basename>_<A>_<B>_... where A, B are indices per substitution
            suffix = "_".join(str(idx) for idx in index_combo)
            output_id = f"{base_id}_{suffix}"
            results.append({'id': output_id, 'sequence': stitched})

    print(f"\nGenerated {len(results)} stitched sequences")
    return results


def stitch_matched_sequences(
    template_sequences: Dict[str, str],
    substitutions: List[Dict],
    id_map: Dict[str, str],
    residue_mappings: Dict[str, Dict[int, int]] = None,
    debug: bool = False
) -> List[Dict[str, Any]]:
    """Stitch sequences by matching base IDs, with support for table-based positions."""
    results = []
    residue_mappings = residue_mappings or {}

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

        # Get residue mapping for this base ID
        residue_mapping = residue_mappings.get(base_id)

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

        # Match sequences by ID: extract suffix from template_id and find matching substitution sequences
        combo_count = 0
        for template_id, template_seq in template_seqs.items():
            # Extract the suffix from template_id (everything after base_id)
            # Example: "Test_1_1" with base_id "Test_1" -> suffix "_1"
            if template_id.startswith(base_id + "_"):
                suffix = template_id[len(base_id) + 1:]  # +1 to skip the underscore
            else:
                # No suffix, skip this template
                print(f"    Warning: template_id '{template_id}' doesn't match base_id '{base_id}'")
                continue

            # Try to find matching substitution sequences with the same suffix
            matched_subs = []
            matched_sub_ids = []
            all_found = True
            for sub_seqs in sub_seqs_for_base:
                # Look for sequence with ID = base_id + "_" + suffix
                target_id = f"{base_id}_{suffix}"
                if target_id in sub_seqs:
                    matched_subs.append(sub_seqs[target_id])
                    matched_sub_ids.append(target_id)
                else:
                    print(f"    Warning: no matching substitution sequence for '{target_id}'")
                    all_found = False
                    break

            if not all_found:
                continue

            # Apply substitutions
            subs = {}
            debug_info = {'template_id': template_id}
            for i, sub_seq in enumerate(matched_subs):
                subs[positions_for_base[i]] = sub_seq
                debug_info[f'sub_{i}_id'] = matched_sub_ids[i]
                debug_info[f'sub_{i}_positions'] = positions_for_base[i]

            stitched = substitute_segments(template_seq, subs, residue_mapping, debug=debug, debug_info=debug_info)
            output_id = template_id
            results.append({'id': output_id, 'sequence': stitched})
            combo_count += 1

        print(f"    Generated {combo_count} combinations")

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
