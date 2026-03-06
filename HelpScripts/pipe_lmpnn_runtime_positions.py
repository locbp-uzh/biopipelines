#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Position determination for LigandMPNN from tables or direct specifications.

This script processes a DataStream JSON file to create a JSON file with LigandMPNN position arguments.

Usage:
    python pipe_lmpnn_runtime_positions.py <structures_json> <input_source> <input_table> <fixed_positions> <designed_positions> <ligand> <design_within> <output_file>

Arguments:
    structures_json: Path to DataStream JSON file with input structures
    input_source: "table" (e.g. RFdiffusion table), "selection" (use direct positions), or "ligand" (ligand-based)
    input_table: Path to table file (or "-" if not using table)
    fixed_positions: PyMOL selection for fixed positions (or "-")
    designed_positions: PyMOL selection for designed positions (or "-")
    ligand: Ligand identifier
    design_within: Distance cutoff for ligand-based design
    output_file: Path to output JSON file

Output:
    JSON file mapping design_id -> {fixed_option, redesigned_option} for LigandMPNN
"""

import sys
import os
import json
import pandas as pd
import glob
from pathlib import Path

# Import unified I/O utilities
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from biopipelines.biopipelines_io import load_datastream, iterate_files, load_table, lookup_table_value
from biopipelines.pdb_parser import parse_pdb_file


def sele_to_list(sele_str):
    """Convert selection string to list of residue numbers."""
    # Handle None, empty, or "-"
    if not sele_str or sele_str == "-" or pd.isna(sele_str):
        return []

    residues = []
    parts = str(sele_str).split('+')

    for part in parts:
        if '-' in part and not part.startswith('-'):
            # Range like "10-20"
            start, end = map(int, part.split('-', 1))
            residues.extend(range(start, end + 1))
        else:
            # Single residue
            residues.append(int(part))

    return sorted(residues)


def process_table_source(input_table, design_entries):
    """Process input table to get fixed/designed positions.

    Args:
        input_table: Path to CSV table with fixed/designed columns
        design_entries: List of (design_id, pdb_file) tuples from DataStream
    """
    if not os.path.exists(input_table):
        raise FileNotFoundError(f"Table not found: {input_table}")

    df = pd.read_csv(input_table)
    positions_data = {}

    for design_id, pdb_file in design_entries:
        # Find matching row in table - try by id column first, then by pdb column
        matching_rows = df[df['id'] == design_id] if 'id' in df.columns else pd.DataFrame()
        if matching_rows.empty and 'pdb' in df.columns:
            pdb_name = os.path.basename(pdb_file)
            matching_rows = df[df['pdb'] == pdb_name]

        if not matching_rows.empty:
            row = matching_rows.iloc[0]
            fixed_positions = sele_to_list(row.get('fixed', ''))
            designed_positions = sele_to_list(row.get('designed', ''))

            positions_data[design_id] = {
                'fixed_positions': fixed_positions,
                'designed_positions': designed_positions,
                'pdb_file': pdb_file
            }
        else:
            print(f"Warning: No table entry found for {design_id}")
            positions_data[design_id] = {
                'fixed_positions': [],
                'designed_positions': [],
                'pdb_file': pdb_file
            }

    return positions_data


def resolve_table_reference(reference, design_ids):
    """
    Resolve table reference to per-design selections.

    Args:
        reference: Either a table reference like "TABLE_REFERENCE:path:column" or direct PyMOL selection
        design_ids: List of design IDs from DataStream

    Returns:
        Dictionary mapping design IDs to position lists
    """
    if not reference.startswith("TABLE_REFERENCE:"):
        # Direct PyMOL selection - same for all designs
        return {design_id: sele_to_list(reference) for design_id in design_ids}

    # Use pipe_biopipelines_io to load table and column
    table, column_name = load_table(reference)
    positions_per_design = {}

    for design_id in design_ids:
        try:
            selection_value = lookup_table_value(table, design_id, column_name)
            positions_per_design[design_id] = sele_to_list(selection_value)
        except KeyError:
            print(f"Warning: No table entry found for {design_id} in column {column_name}")
            positions_per_design[design_id] = []

    return positions_per_design


def get_protein_residues_from_pdb(pdb_path, chain="A"):
    """Get all protein residue numbers for a specific chain from PDB file."""
    standard_residues = {
        'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
        'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'
    }
    atoms = parse_pdb_file(pdb_path)
    residues = set()
    for atom in atoms:
        if atom.chain == chain and atom.res_name in standard_residues:
            residues.add(atom.res_num)
    return sorted(list(residues))


def process_selection_source(fixed_positions, designed_positions, design_entries):
    """Process direct PyMOL selections or table references for fixed/designed positions.

    Args:
        fixed_positions: Fixed positions string or TABLE_REFERENCE
        designed_positions: Designed positions string or TABLE_REFERENCE
        design_entries: List of (design_id, pdb_file) tuples from DataStream
    """
    design_ids = [entry[0] for entry in design_entries]
    fixed_per_design = resolve_table_reference(fixed_positions, design_ids)
    designed_per_design = resolve_table_reference(designed_positions, design_ids)

    positions_data = {}
    for design_id, pdb_file in design_entries:
        fixed = fixed_per_design.get(design_id, [])
        designed = designed_per_design.get(design_id, [])

        # When redesigned resolves to empty and fixed is also empty,
        # fix all residues so the tool outputs the original sequence.
        # This handles the case where a table reference (e.g. DistanceSelector "within")
        # resolved to empty at runtime.
        if not designed and not fixed:
            all_residues = get_protein_residues_from_pdb(pdb_file)
            fixed = all_residues
            print(f"No positions to redesign for {design_id}, fixing all residues")

        positions_data[design_id] = {
            'fixed_positions': fixed,
            'designed_positions': designed,
            'pdb_file': pdb_file
        }

    return positions_data


def process_ligand_source(ligand, design_within, design_entries):
    """Process ligand-based design - determine positions based on distance to ligand.

    Args:
        ligand: Ligand identifier
        design_within: Distance cutoff in Angstroms
        design_entries: List of (design_id, pdb_file) tuples from DataStream
    """
    # For now, implement a simple fallback - in practice this would need
    # structural analysis to find residues near the ligand
    print(f"Warning: Ligand-based position determination not fully implemented. Using design_within={design_within}Å around {ligand}")

    positions_data = {}
    for design_id, pdb_file in design_entries:
        # Fallback: assume all positions are designable for ligand-based approach
        # In a full implementation, this would analyze the PDB structure
        positions_data[design_id] = {
            'fixed_positions': [],
            'designed_positions': [],  # Would be populated by structural analysis
            'pdb_file': pdb_file
        }

    return positions_data


def write_positions_json(positions_data, output_file):
    """Write positions data as JSON for runtime lookup.

    Args:
        positions_data: Dict mapping design_id -> {fixed_positions, designed_positions, pdb_file}
        output_file: Path to output JSON file
    """
    result = {}
    for design_id, positions in positions_data.items():
        # Build fixed positions option for this design
        fixed_option = ""
        if positions['fixed_positions']:
            fixed_str = " ".join([f"A{pos}" for pos in positions['fixed_positions']])
            fixed_option = f'--fixed_residues "{fixed_str}"'

        # Build designed positions option for this design
        redesigned_option = ""
        if positions['designed_positions']:
            designed_str = " ".join([f"A{pos}" for pos in positions['designed_positions']])
            redesigned_option = f'--redesigned_residues "{designed_str}"'

        result[design_id] = {
            "fixed_option": fixed_option,
            "redesigned_option": redesigned_option
        }

    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    with open(output_file, 'w') as f:
        json.dump(result, f, indent=2)

    print(f"Created LigandMPNN positions JSON: {output_file}")


def main():
    if len(sys.argv) != 9:
        print("Usage: python pipe_lmpnn_runtime_positions.py <structures_json> <input_source> <input_table> <fixed_positions> <designed_positions> <ligand> <design_within> <output_file>")
        sys.exit(1)

    structures_json = sys.argv[1]
    input_source = sys.argv[2]
    input_table = sys.argv[3]
    fixed_positions = sys.argv[4]
    designed_positions = sys.argv[5]
    ligand = sys.argv[6]
    design_within = float(sys.argv[7])
    output_file = sys.argv[8]

    # Load DataStream and get (id, file) pairs
    structures_ds = load_datastream(structures_json)
    design_entries = list(iterate_files(structures_ds))  # List of (design_id, pdb_file) tuples

    if not design_entries:
        raise ValueError(f"No structures found in DataStream: {structures_json}")

    print(f"Processing {len(design_entries)} structures with input_source='{input_source}'")

    # Process based on input source
    if input_source == "table" and input_table != "-":
        positions_data = process_table_source(input_table, design_entries)
    elif input_source == "selection":
        positions_data = process_selection_source(fixed_positions, designed_positions, design_entries)
    elif input_source == "ligand":
        positions_data = process_ligand_source(ligand, design_within, design_entries)
    else:
        raise ValueError(f"Invalid input_source: {input_source}")

    # Write positions JSON for runtime lookup
    write_positions_json(positions_data, output_file)


if __name__ == "__main__":
    main()