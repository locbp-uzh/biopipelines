#!/usr/bin/env python3
"""
Helper script to build Boltz2 YAML configuration with MSA paths.

This script creates a Boltz2 configuration that specifies MSA file paths
instead of triggering MSA calculation, enabling MSA recycling between runs.
"""

import os
import sys
import argparse
import pandas as pd
import yaml

def build_config_with_msas_from_csv(sequences_csv: str, msa_table: str, 
                                   ligand_smiles: str = None, affinity: bool = False) -> str:
    """
    Build Boltz2 YAML configuration with MSA paths from CSV files.
    
    Args:
        sequences_csv: Path to CSV with id,sequence columns 
        msa_table: Path to CSV with id,sequence_id,sequence,msa_file columns
        ligand_smiles: Ligand SMILES string (optional)
        affinity: Whether to calculate binding affinity
        
    Returns:
        YAML configuration string
    """
    # Read sequences CSV
    try:
        sequences_df = pd.read_csv(sequences_csv)
        print(f"Read {len(sequences_df)} sequences from {sequences_csv}")
    except Exception as e:
        print(f"Error reading sequences CSV: {e}")
        return None
    
    # Read MSA table to create dual lookup: by sequence_id and by sequence
    msa_lookup_by_id = {}
    msa_lookup_by_sequence = {}
    if msa_table and os.path.exists(msa_table):
        try:
            msa_df = pd.read_csv(msa_table)
            print(f"Read {len(msa_df)} MSA entries from {msa_table}")

            # Create dual lookup
            for _, row in msa_df.iterrows():
                sequence_id = row['sequence_id']
                sequence = row['sequence']
                msa_file = row['msa_file']
                if os.path.exists(msa_file):
                    msa_lookup_by_id[sequence_id] = msa_file
                    msa_lookup_by_sequence[sequence] = msa_file
                    print(f"MSA found for {sequence_id}: {os.path.basename(msa_file)}")
                else:
                    print(f"Warning: MSA file not found for {sequence_id}: {msa_file}")
        except Exception as e:
            print(f"Error reading MSA table {msa_table}: {e}")
    
    # Build configuration dictionary
    config = {'sequences': []}
    
    # Process each sequence
    for _, row in sequences_df.iterrows():
        seq_id = row['id']
        sequence = row['sequence']
        
        # Create protein entry
        protein_entry = {
            'protein': {
                'id': 'A',
                'sequence': sequence
            }
        }
        
        # Add MSA if available - priority: ID match -> sequence match -> error
        if msa_table and os.path.exists(msa_table):
            msa_file = None
            # Priority 1: Try exact sequence_id match
            if seq_id in msa_lookup_by_id:
                msa_file = msa_lookup_by_id[seq_id]
                print(f"Added MSA for sequence {seq_id} (matched by ID)")
            # Priority 2: Try sequence match
            elif sequence in msa_lookup_by_sequence:
                msa_file = msa_lookup_by_sequence[sequence]
                print(f"Added MSA for sequence {seq_id} (matched by sequence)")
            # Priority 3: Error - MSA was provided but cannot be matched
            else:
                raise ValueError(f"ERROR: MSA table provided but no MSA found for sequence {seq_id}. "
                               f"Tried matching by ID and by sequence. Cannot proceed with Boltz2.")

            protein_entry['protein']['msa'] = msa_file
        
        # Add ligand if provided
        if ligand_smiles:
            protein_entry['ligand'] = {
                'id': 'B',
                'smiles': ligand_smiles
            }
        
        config['sequences'].append(protein_entry)
    
    # Add affinity calculation if requested and ligand present
    if affinity and ligand_smiles:
        config['properties'] = [
            {
                'affinity': {
                    'binder': 'B'
                }
            }
        ]
    
    # Convert to YAML string
    yaml_str = yaml.dump(config, default_flow_style=False, sort_keys=False)
    return yaml_str

def build_config_with_msas(protein_sequence: str, ligand_smiles: str, 
                          msa_table: str, affinity: bool = False) -> str:
    """
    Build Boltz2 YAML configuration with MSA paths (legacy single sequence version).
    
    Args:
        protein_sequence: Protein amino acid sequence
        ligand_smiles: Ligand SMILES string
        msa_table: Path to CSV with MSA file mappings
        affinity: Whether to calculate binding affinity
        
    Returns:
        YAML configuration string
    """
    # Read MSA table to get MSA file path
    msa_file_path = None
    if msa_table and os.path.exists(msa_table):
        try:
            msa_df = pd.read_csv(msa_table)
            # Look for MSA file for protein sequence A (assuming single protein)
            protein_rows = msa_df[msa_df['sequence_id'] == 'A']
            if not protein_rows.empty:
                msa_file_path = protein_rows.iloc[0]['msa_file']
                print(f"Found MSA file for sequence A: {msa_file_path}")
            else:
                print(f"Warning: No MSA file found for sequence A in {msa_table}")
        except Exception as e:
            print(f"Error reading MSA table {msa_table}: {e}")
    
    # Build configuration dictionary
    config = {
        'sequences': [
            {
                'protein': {
                    'id': 'A',
                    'sequence': protein_sequence
                }
            }
        ]
    }
    
    # Add MSA path if found
    if msa_file_path and os.path.exists(msa_file_path):
        config['sequences'][0]['protein']['msa'] = msa_file_path
        print(f"Added MSA path to config: {msa_file_path}")
    else:
        print("Warning: MSA file not found or doesn't exist, Boltz2 will calculate MSAs")
    
    # Add ligand if provided
    if ligand_smiles:
        config['sequences'].append({
            'ligand': {
                'id': 'B',
                'smiles': ligand_smiles
            }
        })
    
    # Add affinity calculation if requested
    if affinity:
        config['properties'] = [
            {
                'affinity': {
                    'binder': 'B'
                }
            }
        ]
    
    # Convert to YAML string
    yaml_str = yaml.dump(config, default_flow_style=False, sort_keys=False)
    return yaml_str

def main():
    parser = argparse.ArgumentParser(description='Build Boltz2 config with MSA paths')
    
    # Support both single sequence and CSV-based approaches
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--protein-sequence', help='Single protein amino acid sequence')
    group.add_argument('--sequences-csv', help='CSV file with id,sequence columns')
    
    parser.add_argument('--ligand-smiles', help='Ligand SMILES string')
    parser.add_argument('--msa-table', help='Path to MSA table CSV')
    parser.add_argument('--output', required=True, help='Output YAML config file')
    parser.add_argument('--affinity', help='Enable affinity calculation (true/false)')
    
    args = parser.parse_args()
    
    # Parse affinity flag
    affinity_flag = args.affinity and args.affinity.lower() == 'true'
    
    # Build configuration
    try:
        if args.sequences_csv:
            # CSV-based approach for multiple sequences
            yaml_config = build_config_with_msas_from_csv(
                sequences_csv=args.sequences_csv,
                msa_table=args.msa_table or "",
                ligand_smiles=args.ligand_smiles or "",
                affinity=affinity_flag
            )
        else:
            # Single sequence approach (legacy)
            yaml_config = build_config_with_msas(
                protein_sequence=args.protein_sequence,
                ligand_smiles=args.ligand_smiles or "",
                msa_table=args.msa_table or "",
                affinity=affinity_flag
            )
        
        if yaml_config is None:
            print("Error: Failed to build configuration")
            sys.exit(1)
        
        # Write to output file
        with open(args.output, 'w') as f:
            f.write(yaml_config)
        
        print(f"Boltz2 configuration written to: {args.output}")
        print("Configuration preview:")
        print(yaml_config[:500] + "..." if len(yaml_config) > 500 else yaml_config)
        
    except Exception as e:
        print(f"Error building configuration: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()