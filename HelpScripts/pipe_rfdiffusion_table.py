# Copyright © 2024 LOCBP @ University of Zürich
# Distributed under MIT license
"""
Create RFdiffusion results datasheet with id, pdb, fixed, designed columns.

This script processes RFdiffusion output to create a standardized CSV datasheet
that tracks design information for pipeline integration. Parses RFdiffusion log
files the same way pmpnn_fixed_positions.py does.
"""

import argparse
import pandas as pd
import os

def list_to_sele(a):
    """Convert list of residue numbers to PyMOL selection string."""
    if not a:
        return ""
    
    s = ""
    i = 0
    while i < len(a):   
        if i > 0: s += "+"
        s += f"{a[i]}"   
        # Represent consecutive indices with a dash
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

def parse_rfd_log(rfd_log_file, design_character='-'):
    """Parse RFdiffusion log file to extract fixed and designed regions."""
    fixed_dict = {}
    mobile_dict = {}
    
    if not os.path.exists(rfd_log_file):
        raise FileNotFoundError(f"RFdiffusion log file not found: {rfd_log_file}")
    
    # Reopen and parse
    with open(rfd_log_file, "r") as rfdlog:
        current_design_id = None
        design_processed = False
        
        for line in rfdlog:
            line = line.strip()
            
            # Look for design start
            if "Making design" in line and '/' in line:
                design_id = os.path.basename(line.split(" ")[-1])

                current_design_id = design_id
                design_processed = False
                fixed_dict[design_id] = {}
                mobile_dict[design_id] = {}
            
            # Look for timestep sequence (only if we have a current design and haven't processed it yet)
            elif current_design_id and not design_processed and "Timestep" in line:
                seq = line.split(" ")[-1].replace('a','')  # Remove ligand atoms
                if seq:  # Only process if sequence is not empty
                    # Chain A is default (same as pmpnn_fixed_positions.py with FIXED_CHAIN='A')
                    fixed_dict[current_design_id]['A'] = [i+1 for i, x in enumerate(seq) if x != design_character]
                    mobile_dict[current_design_id]['A'] = [i+1 for i, x in enumerate(seq) if x == design_character]
                    design_processed = True  # Mark this design as processed
    
    if not fixed_dict:
        raise ValueError(f"No design information found in log file {rfd_log_file}. Check RFdiffusion output format.")
    
    return fixed_dict, mobile_dict

def main():
    parser = argparse.ArgumentParser(description='Create RFdiffusion results datasheet')
    parser.add_argument('output_folder', type=str, help="RFdiffusion output folder")
    parser.add_argument('rfd_log_file', type=str, help="Path to RFdiffusion log file")
    parser.add_argument('design_character', type=str, help="Design character: '-' for RFdiffusion, '?' for RFdiffusion-AllAtom")
    parser.add_argument('pipeline_name', type=str, help="Pipeline name for ID generation")
    parser.add_argument('num_designs', type=int, help="Number of designs generated")
    parser.add_argument('datasheet_path', type=str, help="Output CSV datasheet path")
    parser.add_argument('design_startnum', type=int, nargs='?', default=0, help="Starting number for design numbering (default: 0 for backward compatibility)")
    
    args = parser.parse_args()
    
    # Parse RFdiffusion log file to get fixed/designed regions
    fixed_dict, mobile_dict = parse_rfd_log(args.rfd_log_file, args.design_character)
    
    # Extract design information
    designs = []
    
    for i in range(args.num_designs):
        design_num = args.design_startnum + i
        pdb_file = f"{args.pipeline_name}_{design_num}.pdb"
        pdb_path = os.path.join(args.output_folder, pdb_file)
        design_id = f"{args.pipeline_name}_{design_num}"
        
        # Get fixed/designed from log file using design ID (no fallback - fail if not found)
        if design_id not in fixed_dict:
            raise ValueError(f"Design ID {design_id} not found in parsed log data. Available design IDs: {list(fixed_dict.keys())}")
        
        fixed_list = fixed_dict[design_id].get('A', [])
        mobile_list = mobile_dict[design_id].get('A', [])
        fixed_regions = list_to_sele(fixed_list)
        designed_regions = list_to_sele(mobile_list)
        
        
        # Check if PDB file exists
        pdb_exists = os.path.exists(pdb_path)
        
        design_info = {
            'id': design_id,
            'pdb': pdb_file,
            'fixed': fixed_regions,
            'designed': designed_regions,
            'exists': pdb_exists
        }
        
        designs.append(design_info)
    
    # Create DataFrame and save
    df = pd.DataFrame(designs)
    
    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(args.datasheet_path), exist_ok=True)
    
    # Save datasheet
    df.to_csv(args.datasheet_path, index=False)
    
    # Save datasheet silently
    existing_pdbs = sum(1 for d in designs if d['exists'])

if __name__ == "__main__":
    main()