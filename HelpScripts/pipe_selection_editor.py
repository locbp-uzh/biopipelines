#!/usr/bin/env python3
# Copyright (c) 2026 Gianluca Quargnali @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# This software is freely available for use, modification, and redistribution.
# If you use this software or any derivative of it in published work,
# you must cite the original author and this repository.
"""
Structure-aware PyMOL selection editor for BioPipelines.

Modifies PyMOL-formatted selection strings (e.g., "3-45+58-60") with operations
like expand, shrink, shift, and invert, while validating against actual PDB residue numbers.

Usage:
    python pipe_selection_editor.py <config_json>

Config JSON should contain:
    - selection_table: Path to CSV with selection column
    - selection_column: Name of column containing selections
    - structures_json: Path to DataStream JSON file with structures
    - expand: Residues to add on each side (int)
    - shrink: Residues to remove from each side (int)
    - shift: Residues to shift by (+/-) (int)
    - invert: Whether to invert selection (bool)
    - output_csv: Path to output CSV file
"""

import sys
import os
import json
import pandas as pd
from typing import List, Tuple, Set

# Import unified I/O utilities and PDB parser
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from biopipelines_io import load_datastream, iterate_files
from pdb_parser import parse_pdb_file, Atom


def parse_pymol_selection(selection: str) -> List[Tuple[int, int]]:
    """
    Parse PyMOL selection string into list of (start, end) tuples.

    Args:
        selection: PyMOL selection string (e.g., "3-45+58-60" or "10+15+20-25")

    Returns:
        List of (start, end) tuples representing ranges
    """
    if not selection or selection.strip() == "":
        return []

    ranges = []
    parts = selection.split('+')

    for part in parts:
        part = part.strip()
        if not part:
            continue

        if '-' in part:
            # Range format: "3-45"
            range_parts = part.split('-')
            if len(range_parts) != 2:
                raise ValueError(f"Invalid range format: {part}")
            start = int(range_parts[0])
            end = int(range_parts[1])
            if start > end:
                raise ValueError(f"Invalid range {part}: start > end")
            ranges.append((start, end))
        else:
            # Single residue: "10"
            res_num = int(part)
            ranges.append((res_num, res_num))

    return ranges


def format_pymol_selection(ranges: List[Tuple[int, int]]) -> str:
    """
    Format list of (start, end) tuples back to PyMOL selection string.

    Args:
        ranges: List of (start, end) tuples

    Returns:
        PyMOL selection string (e.g., "3-45+58-60")
    """
    if not ranges:
        return ""

    parts = []
    for start, end in ranges:
        if start == end:
            parts.append(f"{start}")
        else:
            parts.append(f"{start}-{end}")

    return "+".join(parts)


def get_protein_residues(pdb_path: str) -> Set[int]:
    """
    Get set of valid protein residue numbers from PDB file.

    Args:
        pdb_path: Path to PDB file

    Returns:
        Set of residue numbers present in the structure
    """
    atoms = parse_pdb_file(pdb_path)

    # Standard amino acid names
    standard_residues = {
        'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
        'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'
    }

    # Extract unique residue numbers for protein residues
    residue_numbers = set()
    for atom in atoms:
        if atom.res_name in standard_residues:
            residue_numbers.add(atom.res_num)

    if not residue_numbers:
        raise ValueError(f"No protein residues found in {pdb_path}")

    return residue_numbers


def merge_overlapping_ranges(ranges: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """
    Merge overlapping or adjacent ranges.

    Args:
        ranges: List of (start, end) tuples

    Returns:
        Merged list of ranges
    """
    if not ranges:
        return []

    # Sort ranges by start position
    sorted_ranges = sorted(ranges)

    merged = [sorted_ranges[0]]
    for current_start, current_end in sorted_ranges[1:]:
        last_start, last_end = merged[-1]

        # Check if ranges overlap or are adjacent
        if current_start <= last_end + 1:
            # Merge by extending the last range
            merged[-1] = (last_start, max(last_end, current_end))
        else:
            # No overlap, add as new range
            merged.append((current_start, current_end))

    return merged


def expand_selection(ranges: List[Tuple[int, int]], expand_by: int,
                    valid_residues: Set[int]) -> List[Tuple[int, int]]:
    """
    Expand selection ranges by adding residues on each side.

    Args:
        ranges: List of (start, end) tuples
        expand_by: Number of residues to add on each side
        valid_residues: Set of valid residue numbers from PDB

    Returns:
        Expanded and merged ranges
    """
    if expand_by == 0:
        return ranges

    expanded = []
    for start, end in ranges:
        # Attempt to expand
        new_start = start - expand_by
        new_end = end + expand_by

        # Find actual valid start (first valid residue >= new_start)
        actual_start = start
        for res in sorted(valid_residues):
            if res >= new_start:
                actual_start = res
                break

        # Find actual valid end (last valid residue <= new_end)
        actual_end = end
        for res in sorted(valid_residues, reverse=True):
            if res <= new_end:
                actual_end = res
                break

        expanded.append((actual_start, actual_end))

    # Merge any overlapping ranges created by expansion
    return merge_overlapping_ranges(expanded)


def shrink_selection(ranges: List[Tuple[int, int]], shrink_by: int,
                    valid_residues: Set[int]) -> List[Tuple[int, int]]:
    """
    Shrink selection ranges by removing residues from each side.

    Args:
        ranges: List of (start, end) tuples
        shrink_by: Number of residues to remove from each side
        valid_residues: Set of valid residue numbers from PDB

    Returns:
        Shrunk ranges (empty ranges are removed)
    """
    if shrink_by == 0:
        return ranges

    shrunk = []
    for start, end in ranges:
        # Calculate new boundaries
        new_start = start + shrink_by
        new_end = end - shrink_by

        # Skip if range would become invalid
        if new_start > new_end:
            continue

        # Validate against PDB residues
        if new_start not in valid_residues or new_end not in valid_residues:
            # Find nearest valid residues
            valid_in_range = sorted([r for r in valid_residues if new_start <= r <= new_end])
            if not valid_in_range:
                continue
            new_start = valid_in_range[0]
            new_end = valid_in_range[-1]

        shrunk.append((new_start, new_end))

    return shrunk


def shift_selection(ranges: List[Tuple[int, int]], shift_by: int,
                   valid_residues: Set[int]) -> List[Tuple[int, int]]:
    """
    Shift selection ranges by a fixed offset.

    Args:
        ranges: List of (start, end) tuples
        shift_by: Number of residues to shift (+/-)
        valid_residues: Set of valid residue numbers from PDB

    Returns:
        Shifted ranges (only includes valid ranges)
    """
    if shift_by == 0:
        return ranges

    shifted = []
    for start, end in ranges:
        new_start = start + shift_by
        new_end = end + shift_by

        # Check if shifted range is valid
        if new_start not in valid_residues or new_end not in valid_residues:
            # Check if any part of the shifted range is valid
            range_residues = [r for r in valid_residues if new_start <= r <= new_end]
            if not range_residues:
                # Skip this range entirely if no valid residues
                continue
            # Adjust to valid boundaries
            new_start = min(range_residues)
            new_end = max(range_residues)

        shifted.append((new_start, new_end))

    return shifted


def invert_selection(ranges: List[Tuple[int, int]],
                    valid_residues: Set[int]) -> List[Tuple[int, int]]:
    """
    Invert selection to get complement (all residues NOT in selection).

    Args:
        ranges: List of (start, end) tuples representing current selection
        valid_residues: Set of valid residue numbers from PDB

    Returns:
        Inverted selection ranges
    """
    # Get all selected residues
    selected = set()
    for start, end in ranges:
        for res in range(start, end + 1):
            if res in valid_residues:
                selected.add(res)

    # Get complement
    inverted_residues = valid_residues - selected

    if not inverted_residues:
        return []

    # Convert back to ranges
    sorted_residues = sorted(inverted_residues)
    ranges = []
    start = sorted_residues[0]
    end = sorted_residues[0]

    for res in sorted_residues[1:]:
        if res == end + 1:
            end = res
        else:
            ranges.append((start, end))
            start = res
            end = res

    ranges.append((start, end))
    return ranges


def modify_selection(selection: str, pdb_path: str, expand: int = 0,
                    shrink: int = 0, shift: int = 0, invert: bool = False) -> str:
    """
    Modify a PyMOL selection string with structure-aware operations.

    Args:
        selection: Original PyMOL selection string
        pdb_path: Path to PDB structure file
        expand: Residues to add on each side
        shrink: Residues to remove from each side
        shift: Residues to shift by (+/-)
        invert: Whether to invert selection

    Returns:
        Modified selection string
    """
    # Parse selection
    ranges = parse_pymol_selection(selection)

    # Get valid residues from PDB
    valid_residues = get_protein_residues(pdb_path)

    # Apply operations in order
    if expand > 0:
        ranges = expand_selection(ranges, expand, valid_residues)

    if shrink > 0:
        ranges = shrink_selection(ranges, shrink, valid_residues)

    if shift != 0:
        ranges = shift_selection(ranges, shift, valid_residues)

    if invert:
        ranges = invert_selection(ranges, valid_residues)

    # Format back to string
    return format_pymol_selection(ranges)


def process_selections(config_path: str):
    """
    Process selections according to config file.

    Args:
        config_path: Path to JSON config file
    """
    # Load config
    with open(config_path, 'r') as f:
        config = json.load(f)

    selection_table = config['selection_table']
    selection_column = config['selection_column']
    structures_json = config['structures_json']
    expand = config.get('expand', 0)
    shrink = config.get('shrink', 0)
    shift = config.get('shift', 0)
    invert = config.get('invert', False)
    output_csv = config['output_csv']

    # Load selection table
    if not os.path.exists(selection_table):
        raise ValueError(f"Selection table not found: {selection_table}")

    df = pd.read_csv(selection_table)

    # Validate column exists
    if selection_column not in df.columns:
        raise ValueError(
            f"Column '{selection_column}' not found in table. "
            f"Available columns: {list(df.columns)}"
        )

    # Load structures DataStream using pipe_biopipelines_io
    structures_ds = load_datastream(structures_json)

    # Create structure ID to path mapping from DataStream
    struct_map = {}
    for struct_id, struct_path in iterate_files(structures_ds):
        struct_map[struct_id] = struct_path

    print(f"Loaded {len(struct_map)} structures from DataStream")
    print(f"Processing {len(df)} selections")
    print(f"Operations: expand={expand}, shrink={shrink}, shift={shift}, invert={invert}")

    # Process each row
    results = []
    for idx, row in df.iterrows():
        struct_id = row['id']
        original_selection = row[selection_column]

        # Find corresponding PDB
        if struct_id not in struct_map:
            print(f"Warning: No PDB found for structure ID '{struct_id}', skipping")
            continue

        pdb_path = struct_map[struct_id]

        if not os.path.exists(pdb_path):
            print(f"Warning: PDB file not found: {pdb_path}, skipping")
            continue

        try:
            # Modify selection
            modified_selection = modify_selection(
                original_selection,
                pdb_path,
                expand=expand,
                shrink=shrink,
                shift=shift,
                invert=invert
            )

            results.append({
                'id': struct_id,
                'pdb': pdb_path,
                selection_column: modified_selection,
                f'original_{selection_column}': original_selection
            })

            print(f"  {struct_id}: '{original_selection}' -> '{modified_selection}'")

        except Exception as e:
            print(f"Error processing {struct_id}: {e}")
            raise

    if not results:
        raise ValueError("No selections were processed successfully")

    # Save results
    result_df = pd.DataFrame(results)
    os.makedirs(os.path.dirname(output_csv), exist_ok=True)
    result_df.to_csv(output_csv, index=False)

    print(f"\nSelection modification completed!")
    print(f"Processed {len(results)} structures")
    print(f"Results saved to: {output_csv}")


def main():
    if len(sys.argv) != 2:
        print("Usage: python pipe_selection_editor.py <config_json>")
        sys.exit(1)

    config_path = sys.argv[1]

    if not os.path.exists(config_path):
        print(f"Error: Config file not found: {config_path}")
        sys.exit(1)

    try:
        process_selections(config_path)
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
