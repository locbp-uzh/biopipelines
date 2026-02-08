#!/usr/bin/env python3
# Copyright (c) 2026 Gianluca Quargnali @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# This software is freely available for use, modification, and redistribution.
# If you use this software or any derivative of it in published work,
# you must cite the original author and this repository.
"""
Position determination for LigandMPNN from tables or direct specifications.

This script processes a DataStream JSON file to create a bash script with LigandMPNN position arguments.

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
    output_file: Path to output bash script file

Output:
    Bash script with FIXED_OPTION and REDESIGNED_OPTION variables for LigandMPNN
"""

import sys
import os
import json
import pandas as pd
import glob
from pathlib import Path

# Import unified I/O utilities
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from biopipelines_io import load_datastream, iterate_files, load_table, lookup_table_value


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
        reference: Either a table reference like "DATASHEET_REFERENCE:path:column" or direct PyMOL selection
        design_ids: List of design IDs from DataStream

    Returns:
        Dictionary mapping design IDs to position lists
    """
    if not reference.startswith("DATASHEET_REFERENCE:"):
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


def process_selection_source(fixed_positions, designed_positions, design_entries):
    """Process direct PyMOL selections or table references for fixed/designed positions.

    Args:
        fixed_positions: Fixed positions string or DATASHEET_REFERENCE
        designed_positions: Designed positions string or DATASHEET_REFERENCE
        design_entries: List of (design_id, pdb_file) tuples from DataStream
    """
    design_ids = [entry[0] for entry in design_entries]
    fixed_per_design = resolve_table_reference(fixed_positions, design_ids)
    designed_per_design = resolve_table_reference(designed_positions, design_ids)

    positions_data = {}
    for design_id, pdb_file in design_entries:
        positions_data[design_id] = {
            'fixed_positions': fixed_per_design.get(design_id, []),
            'designed_positions': designed_per_design.get(design_id, []),
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


def create_ligandmpnn_bash_script(positions_data, output_file):
    """Create bash script that modifies the main script with position arguments.

    Args:
        positions_data: Dict mapping design_id -> {fixed_positions, designed_positions, pdb_file}
        output_file: Path to output bash script
    """
    # Build bash script content that does placeholder replacement
    script_content = "#!/bin/bash\n"
    script_content += "# LigandMPNN position replacement script - Generated by pipe_lmpnn_runtime_positions.py\n\n"

    script_content += "# Get the script file to modify\n"
    script_content += "SCRIPT_FILE=\"$1\"\n"
    script_content += "if [ -z \"$SCRIPT_FILE\" ]; then\n"
    script_content += "    echo \"Usage: $0 <script_file>\"\n"
    script_content += "    exit 1\n"
    script_content += "fi\n\n"

    script_content += "echo \"Updating $SCRIPT_FILE with LigandMPNN position arguments\"\n\n"

    # For each design, create the replacement commands
    for design_id, positions in positions_data.items():
        pdb_file = positions.get('pdb_file', '')
        pdb_name = os.path.basename(pdb_file) if pdb_file else f"{design_id}.pdb"

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

        # Add sed commands to replace placeholders - escape special characters
        script_content += f"# Replace placeholders for {design_id} ({pdb_name})\n"

        # Escape the fixed_option for sed (escape forward slashes, quotes, etc.)
        fixed_escaped = fixed_option.replace('/', '\\/').replace('"', '\\"')
        redesigned_escaped = redesigned_option.replace('/', '\\/').replace('"', '\\"')

        script_content += f"sed -i 's|{design_id}_FIXED_OPTION_PLACEHOLDER|{fixed_escaped}|g' \"$SCRIPT_FILE\"\n"
        script_content += f"sed -i 's|{design_id}_REDESIGNED_OPTION_PLACEHOLDER|{redesigned_escaped}|g' \"$SCRIPT_FILE\"\n"
        script_content += f"echo \"Updated placeholders for {design_id}: fixed='{fixed_option}' redesigned='{redesigned_option}'\"\n\n"

    script_content += "echo \"All position placeholders updated successfully\"\n"

    # Write bash script
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    with open(output_file, 'w') as f:
        f.write(script_content)

    # Make executable
    os.chmod(output_file, 0o755)

    print(f"Created LigandMPNN replacement script: {output_file}")


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

    # Create LigandMPNN bash script with position variables
    create_ligandmpnn_bash_script(positions_data, output_file)


if __name__ == "__main__":
    main()