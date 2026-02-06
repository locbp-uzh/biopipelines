#Copyright © 2024 LOCBP @ University of Zürich
#Distributed under MIT license
import argparse
import os
import sys

# Import PDB parser and I/O utilities
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from pdb_parser import parse_pdb_file
from biopipelines_io import load_datastream, iterate_files, load_table, lookup_table_value

parser = argparse.ArgumentParser(description='Establish residues for which to generate a sequence with ProteinMPNN. Can read from RFdiffusion table, use direct selections, or apply pLDDT threshold')
parser.add_argument('structures_json', type=str, help="Path to DataStream JSON file with input structures")
parser.add_argument('input_source', type=str, help = "'table' for reading from table, 'selection' for direct input, 'plddt' for threshold")
parser.add_argument('input_table', type=str, help='Path to input table (e.g., tool_results.csv) or "-" if not used')
parser.add_argument('pLDDT_thr', type=float, help="Will consider residues with pLDDT > threshold as fixed ")
parser.add_argument('FIXED', type=str, help="PyMOL selection for fixed positions or '-'")
parser.add_argument('DESIGNED', type=str, help="PyMOL selection for designed positions or '-'")
parser.add_argument('FIXED_CHAIN', type=str, help="A or B or whatever")
parser.add_argument('fixed_jsonl_file', type=str, help="output file")
parser.add_argument('sele_csv_file', type=str, help="output file with selections of fixed and mobile parts")
# Parse the arguments
args = parser.parse_args()
structures_json=args.structures_json
input_source=args.input_source
input_table=args.input_table
pLDDT_thr=args.pLDDT_thr
FIXED=args.FIXED
DESIGNED=args.DESIGNED
FIXED_CHAIN=args.FIXED_CHAIN
fixed_jsonl_file=args.fixed_jsonl_file
sele_csv_file=args.sele_csv_file

#Converts a selection string into an array
def sele_to_list(s):
    """Convert selection string to list of residue numbers."""
    a = []
    if not s or s == "":
        return a

    # Convert to string to handle numeric types
    s = str(s)
    if s == "nan":
        return a

    # Handle both '+' separated ranges and space-separated legacy format
    if '+' in s:
        parts = s.split('+')
    else:
        # Legacy space-separated format
        parts = s.split()

    for part in parts:
        part = part.strip()
        if not part:
            continue

        if '-' in part and not part.startswith('-'):
            # Range like "10-15"
            range_parts = part.split('-')
            if len(range_parts) == 2:
                min_val, max_val = range_parts
                for ri in range(int(min_val), int(max_val) + 1):
                    a.append(ri)
            else:
                print(f"Warning: Malformed range '{part}', skipping")
        else:
            # Single residue
            try:
                a.append(int(part))
            except ValueError:
                print(f"Warning: Could not parse '{part}' as integer, skipping")

    return sorted(a)

def list_to_sele(a):
    s = ""
    i = 0
    while i < len(a):
        if i > 0: s += "+"
        s += f"{a[i]}"
        #represent consecutive indeces with a dash
        if i < len(a) - 1:
            if int(a[i])+1 == int(a[i+1]):
                s += "-"
                j = i + 2
                while j < len(a):
                    if int(a[j]) != int(a[j-1])+1: break
                    j += 1
                i = j - 1
                s += f"{a[i]}"
        i += 1
    return s

def get_protein_residues_from_pdb(pdb_path, chain):
    """
    Get all protein residue numbers for a specific chain from PDB file.

    Args:
        pdb_path: Path to PDB file
        chain: Chain identifier (e.g., 'A')

    Returns:
        Sorted list of residue numbers for the specified chain
    """
    # Standard amino acid names
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

def compute_complement(all_residues, redesigned_residues):
    """
    Compute complement: all_residues - redesigned_residues.

    Args:
        all_residues: List of all protein residue numbers
        redesigned_residues: List of residues to redesign

    Returns:
        Sorted list of residues NOT in redesigned_residues
    """
    redesigned_set = set(redesigned_residues)
    complement = [res for res in all_residues if res not in redesigned_set]
    return sorted(complement)

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
            if pd.notna(selection_value) and selection_value != '':
                positions_per_design[design_id] = sele_to_list(str(selection_value))
            else:
                positions_per_design[design_id] = []
        except KeyError:
            print(f"Warning: No table entry found for {design_id} in column {column_name}")
            positions_per_design[design_id] = []

    return positions_per_design

import os
import pandas as pd
import pymol
from pymol import cmd

# Load DataStream and get (id, file) pairs
structures_ds = load_datastream(structures_json)
design_entries = list(iterate_files(structures_ds))  # List of (design_id, pdb_file) tuples
design_ids = [entry[0] for entry in design_entries]
design_files = {entry[0]: entry[1] for entry in design_entries}  # Map id -> file path

fixed_dict = dict()
mobile_dict = dict()

if input_source == "table": #derive from table
    try:
        # Read input table
        df = pd.read_csv(input_table)
        print(f"Reading from table: {input_table}")

        for design_id in design_ids:
            # Find matching row - try by id column first, then by pdb column
            matching_rows = df[df['id'] == design_id] if 'id' in df.columns else pd.DataFrame()
            if matching_rows.empty and 'pdb' in df.columns:
                # Try matching by pdb filename
                matching_rows = df[df['pdb'].apply(lambda x: os.path.splitext(os.path.basename(str(x)))[0]) == design_id]

            if not matching_rows.empty:
                row = matching_rows.iloc[0]
                fixed_dict[design_id] = dict()
                mobile_dict[design_id] = dict()

                # Parse fixed and designed selections from table
                if 'fixed' in df.columns and pd.notna(row['fixed']) and row['fixed'] != '':
                    fixed_dict[design_id][FIXED_CHAIN] = sele_to_list(str(row['fixed']))
                else:
                    fixed_dict[design_id][FIXED_CHAIN] = []

                if 'designed' in df.columns and pd.notna(row['designed']) and row['designed'] != '':
                    mobile_dict[design_id][FIXED_CHAIN] = sele_to_list(str(row['designed']))
                else:
                    mobile_dict[design_id][FIXED_CHAIN] = []

                print(f"Design: {design_id}, Fixed: {fixed_dict[design_id][FIXED_CHAIN]}, Designed: {mobile_dict[design_id][FIXED_CHAIN]}")
            else:
                print(f"Warning: No table entry found for {design_id}")

    except Exception as e:
        print(f"Error reading table {input_table}: {e}")
        print("Falling back to pLDDT threshold method...")
        input_source = "plddt"

if input_source == "selection":
    # Use direct selections or table references
    FIXED = FIXED if FIXED != '-' else ''
    DESIGNED = DESIGNED if DESIGNED != '-' else ''

    # Resolve table references if present
    fixed_per_design = resolve_table_reference(FIXED, design_ids) if FIXED else {design_id: [] for design_id in design_ids}
    designed_per_design = resolve_table_reference(DESIGNED, design_ids) if DESIGNED else {design_id: [] for design_id in design_ids}

    for design_id in design_ids:
        fixed_dict[design_id] = dict()
        mobile_dict[design_id] = dict()

        # Store original mobile/designed positions for documentation
        mobile_dict[design_id][FIXED_CHAIN] = designed_per_design[design_id]

        # Compute what ProteinMPNN should keep fixed:
        # Option C: Union of explicit fixed + complement of redesigned
        pdb_path = design_files[design_id]

        # Start with explicit fixed positions
        final_fixed = list(fixed_per_design[design_id])

        # If redesigned positions are specified, add their complement to fixed
        if designed_per_design[design_id]:
            # Get all protein residues from PDB
            all_residues = get_protein_residues_from_pdb(pdb_path, FIXED_CHAIN)

            # Compute complement of redesigned positions
            complement = compute_complement(all_residues, designed_per_design[design_id])

            # Union: fixed + complement
            final_fixed = sorted(list(set(final_fixed + complement)))

            print(f"Design: {design_id}, Selection-based - Explicit Fixed: {list_to_sele(fixed_per_design[design_id]) if fixed_per_design[design_id] else ''}, Redesigned: {list_to_sele(designed_per_design[design_id])}, Final Fixed (to ProteinMPNN): {list_to_sele(final_fixed)}")
        else:
            # No redesigned specified, use fixed as-is
            print(f"Design: {design_id}, Selection-based - Fixed: {list_to_sele(final_fixed) if final_fixed else ''}, Redesigned: ")

        # Store final fixed positions (what ProteinMPNN will use)
        fixed_dict[design_id][FIXED_CHAIN] = final_fixed
        
elif input_source == "plddt" or input_source == "table":  # table fallback
    # Use pLDDT threshold method
    FIXED = FIXED if FIXED != '-' else ''
    if pLDDT_thr < 100:
        try:
            # Initialize PyMOL in headless mode (no GUI)
            pymol.pymol_argv = ['pymol', '-c']  # -q for quiet, -c for no GUI
            pymol.finish_launching()
        except Exception as e:
            print("Error while initializing pymol")
            print(str(e))

    for design_id in design_ids:
        if design_id not in fixed_dict:  # Skip if already processed from table
            fixed_dict[design_id] = dict()
            mobile_dict[design_id] = dict()

            if pLDDT_thr < 100:
                pdb_file = design_files[design_id]
                try:
                    cmd.load(pdb_file,"prot")
                    fixed_residues = []
                    mobile_residues = []
                    atom_iterator = cmd.get_model("prot and name CA")
                    parfixed = sele_to_list(FIXED)
                    for atom in atom_iterator.atom:
                        resi = int(atom.resi)
                        if atom.b < pLDDT_thr and not resi in parfixed:
                            if not resi in mobile_residues:
                                mobile_residues.append(int(atom.resi))
                        else:
                            if not resi in fixed_residues:
                                fixed_residues.append(int(atom.resi))
                    cmd.delete("prot")
                    fixed_dict[design_id][FIXED_CHAIN] = fixed_residues[:]
                    mobile_dict[design_id][FIXED_CHAIN] = mobile_residues[:]
                    print(f"Design: {design_id}, pLDDT-based - Fixed: {len(fixed_residues)}, Mobile: {len(mobile_residues)}")
                except Exception as e:
                    print("Error while calculating fixed positions")
                    print(str(e))
                    fixed_dict[design_id][FIXED_CHAIN] = sele_to_list(FIXED)
                    mobile_dict[design_id][FIXED_CHAIN] = []
            else:
                fixed_dict[design_id][FIXED_CHAIN] = sele_to_list(FIXED)
                mobile_dict[design_id][FIXED_CHAIN] = sele_to_list(DESIGNED)
                print(f"Design: {design_id}, Selection-based - Fixed: {FIXED}, Redesigned: {DESIGNED}")

with open(fixed_jsonl_file,"w") as jsonl_file:
    #Python converts dictionaries to string having keys inside '', json only recognises ""
    jsonl_file.write(str(fixed_dict).replace("\'","\""))

with open(sele_csv_file,"w") as csv_file:
    csv_file.write("id,fixed,mobile")
    for id in fixed_dict.keys():
        fixed = list_to_sele(fixed_dict[id][FIXED_CHAIN]) if FIXED_CHAIN in fixed_dict[id].keys() else ""
        mobile = list_to_sele(mobile_dict[id][FIXED_CHAIN]) if FIXED_CHAIN in fixed_dict[id].keys() else ""
        csv_file.write("\n")
        csv_file.write(f"{id},{fixed},{mobile}")
        
print(f"Fixed positions written to: {fixed_jsonl_file}")
print(f"Selections summary written to: {sele_csv_file}")
        
#Dictionary of fixed positions looks like this
#{"design_0": {"A": [1, 2, 3, 7, 8, 9, 22, 25, 33], "B": []}, "design_1": {"A": [], "B": []}}
#Input modes:
# - table: Read fixed/designed from input table (e.g., from RFdiffusion or other tools)
# - selection: Use FIXED and DESIGNED parameters directly
# - plddt: Use pLDDT threshold to determine positions

if pLDDT_thr < 100:    
    cmd.quit() #MUST ALWAYS BE AT THE END OF THE SCRIPT