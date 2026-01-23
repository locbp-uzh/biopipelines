#!/usr/bin/env python3
"""
Position determination for LigandMPNN from tables or direct specifications.

This script processes a table containing columns id, fixed, designed, direct position specifications, or
ligand-based design to create a bash script with LigandMPNN position arguments.

Usage:
    python pipe_lmpnn_runtime_positions.py <structures_list_file> <input_source> <input_table> <fixed_positions> <designed_positions> <ligand> <design_within> <output_file>

Arguments:
    structures_list_file: File containing list of PDB file paths (one per line)
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


def process_table_source(input_table, pdb_files):
    """Process input table to get fixed/designed positions."""
    if not os.path.exists(input_table):
        raise FileNotFoundError(f"Table not found: {input_table}")
    
    df = pd.read_csv(input_table)
    positions_data = {}
    
    for pdb_file in pdb_files:
        pdb_name = os.path.basename(pdb_file)
        pdb_base = os.path.splitext(pdb_name)[0]
        
        # Find matching row in table
        matching_rows = df[df['pdb'] == pdb_name]
        if matching_rows.empty:
            # Try with base name
            matching_rows = df[df['id'] == pdb_base]
        
        if not matching_rows.empty:
            row = matching_rows.iloc[0]
            fixed_positions = sele_to_list(row.get('fixed', ''))
            designed_positions = sele_to_list(row.get('designed', ''))
            
            positions_data[pdb_file] = {
                'fixed_positions': fixed_positions,
                'designed_positions': designed_positions
            }
        else:
            print(f"Warning: No table entry found for {pdb_name}")
            positions_data[pdb_file] = {
                'fixed_positions': [],
                'designed_positions': []
            }
    
    return positions_data


def resolve_table_reference(reference, pdb_files):
    """
    Resolve table reference to per-PDB selections.
    
    Args:
        reference: Either a table reference like "DATASHEET_REFERENCE:path:column" or direct PyMOL selection
        pdb_files: List of PDB files to process
        
    Returns:
        Dictionary mapping PDB files to position lists
    """
    if not reference.startswith("DATASHEET_REFERENCE:"):
        # Direct PyMOL selection
        return {pdb_file: sele_to_list(reference) for pdb_file in pdb_files}
    
    # Parse table reference: DATASHEET_REFERENCE:path:column
    _, table_path, column_name = reference.split(":", 2)
    
    if not os.path.exists(table_path):
        raise FileNotFoundError(f"Table not found: {table_path}")
    
    df = pd.read_csv(table_path)
    positions_per_pdb = {}
    
    for pdb_file in pdb_files:
        pdb_name = os.path.basename(pdb_file)
        pdb_base = os.path.splitext(pdb_name)[0]
        
        # Find matching row in table
        matching_rows = df[df['pdb'] == pdb_name]
        if matching_rows.empty:
            # Try with base name
            matching_rows = df[df['id'] == pdb_base]
        
        if not matching_rows.empty:
            row = matching_rows.iloc[0]
            selection_value = row.get(column_name, '')
            positions_per_pdb[pdb_file] = sele_to_list(selection_value)
        else:
            print(f"Warning: No table entry found for {pdb_name} in column {column_name}")
            positions_per_pdb[pdb_file] = []
    
    return positions_per_pdb


def process_selection_source(fixed_positions, designed_positions, pdb_files):
    """Process direct PyMOL selections or table references for fixed/designed positions.""" 
    fixed_per_pdb = resolve_table_reference(fixed_positions, pdb_files)
    designed_per_pdb = resolve_table_reference(designed_positions, pdb_files)
    
    positions_data = {}
    for pdb_file in pdb_files:
        positions_data[pdb_file] = {
            'fixed_positions': fixed_per_pdb.get(pdb_file, []),
            'designed_positions': designed_per_pdb.get(pdb_file, [])
        }
    
    return positions_data


def process_ligand_source(ligand, design_within, pdb_files):
    """Process ligand-based design - determine positions based on distance to ligand."""
    # For now, implement a simple fallback - in practice this would need 
    # structural analysis to find residues near the ligand
    print(f"Warning: Ligand-based position determination not fully implemented. Using design_within={design_within}Å around {ligand}")
    
    positions_data = {}
    for pdb_file in pdb_files:
        # Fallback: assume all positions are designable for ligand-based approach
        # In a full implementation, this would analyze the PDB structure
        positions_data[pdb_file] = {
            'fixed_positions': [],
            'designed_positions': []  # Would be populated by structural analysis
        }
    
    return positions_data


def create_ligandmpnn_bash_script(positions_data, output_file):
    """Create bash script that modifies the main script with position arguments."""
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
    
    # For each PDB, create the replacement commands
    for pdb_file, positions in positions_data.items():
        pdb_name = os.path.basename(pdb_file)
        pdb_id = os.path.splitext(pdb_name)[0]
        
        # Build fixed positions option for this PDB
        fixed_option = ""
        if positions['fixed_positions']:
            fixed_str = " ".join([f"A{pos}" for pos in positions['fixed_positions']])
            fixed_option = f'--fixed_residues "{fixed_str}"'
        
        # Build designed positions option for this PDB
        redesigned_option = ""
        if positions['designed_positions']:
            designed_str = " ".join([f"A{pos}" for pos in positions['designed_positions']])
            redesigned_option = f'--redesigned_residues "{designed_str}"'
        
        # Add sed commands to replace placeholders - escape special characters
        script_content += f"# Replace placeholders for {pdb_name}\n"
        
        # Escape the fixed_option for sed (escape forward slashes, quotes, etc.)
        fixed_escaped = fixed_option.replace('/', '\\/').replace('"', '\\"')
        redesigned_escaped = redesigned_option.replace('/', '\\/').replace('"', '\\"')
        
        script_content += f"sed -i 's|{pdb_id}_FIXED_OPTION_PLACEHOLDER|{fixed_escaped}|g' \"$SCRIPT_FILE\"\n"
        script_content += f"sed -i 's|{pdb_id}_REDESIGNED_OPTION_PLACEHOLDER|{redesigned_escaped}|g' \"$SCRIPT_FILE\"\n"
        script_content += f"echo \"Updated placeholders for {pdb_name}: fixed='{fixed_option}' redesigned='{redesigned_option}'\"\n\n"
    
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
        print("Usage: python pipe_lmpnn_runtime_positions.py <structures_list_file> <input_source> <input_table> <fixed_positions> <designed_positions> <ligand> <design_within> <output_file>")
        sys.exit(1)

    structures_list_file = sys.argv[1]
    input_source = sys.argv[2]
    input_table = sys.argv[3]
    fixed_positions = sys.argv[4]
    designed_positions = sys.argv[5]
    ligand = sys.argv[6]
    design_within = float(sys.argv[7])
    output_file = sys.argv[8]

    # Read structure files from list file (one path per line)
    with open(structures_list_file, 'r') as f:
        pdb_files = [line.strip() for line in f if line.strip()]
    if not pdb_files:
        raise ValueError(f"No structure files found in: {structures_list_file}")
    
    print(f"Processing {len(pdb_files)} PDB files with input_source='{input_source}'")
    
    # Process based on input source
    if input_source == "table" and input_table != "-":
        positions_data = process_table_source(input_table, pdb_files)
    elif input_source == "selection":
        positions_data = process_selection_source(fixed_positions, designed_positions, pdb_files)
    elif input_source == "ligand":
        positions_data = process_ligand_source(ligand, design_within, pdb_files)
    else:
        raise ValueError(f"Invalid input_source: {input_source}")
    
    # Create LigandMPNN bash script with position variables
    create_ligandmpnn_bash_script(positions_data, output_file)


if __name__ == "__main__":
    main()