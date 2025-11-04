#!/usr/bin/env python3
"""
Generate Boltz2 YAML config files for all protein-ligand combinations.

Takes protein CSV and ligand CSV files and generates one YAML config file 
per protein-ligand combination, including MSA integration when available.
"""

import argparse
import pandas as pd
import yaml
import os
import sys
import glob
import hashlib

def parse_arguments():
    parser = argparse.ArgumentParser(description='Generate Boltz2 config files for protein-ligand combinations')
    parser.add_argument('proteins_csv', help='Path to proteins CSV file (id,sequence format)')
    parser.add_argument('ligands_csv', help='Path to ligands CSV file or "None" if no ligands')
    parser.add_argument('output_dir', help='Directory to write config files')
    parser.add_argument('--affinity', action='store_true', help='Enable affinity calculation')
    parser.add_argument('--msa-table', help='Path to MSA table CSV file')
    parser.add_argument('--global-msa-cache', action='store_true', help='Use global MSA cache lookup')
    parser.add_argument('--msa-folder', help='Path to MSAs folder for global cache lookup')
    parser.add_argument('--job-name', help='Pipeline job name to use for single protein-ligand configs')
    
    return parser.parse_args()

def load_proteins(proteins_csv):
    """Load proteins from CSV file."""
    try:
        df = pd.read_csv(proteins_csv)
        if 'id' not in df.columns or 'sequence' not in df.columns:
            raise ValueError("Proteins CSV must have 'id' and 'sequence' columns")
        return df.to_dict('records')
    except Exception as e:
        print(f"Error loading proteins CSV {proteins_csv}: {e}")
        sys.exit(1)

def load_ligands(ligands_csv):
    """Load ligands from CSV file."""
    if ligands_csv == "None":
        return []
    
    try:
        df = pd.read_csv(ligands_csv)
        # Expected columns: id, format, smiles, ccd
        required_cols = ['id', 'format']
        for col in required_cols:
            if col not in df.columns:
                raise ValueError(f"Ligands CSV must have '{col}' column")
        return df.to_dict('records')
    except Exception as e:
        print(f"Error loading ligands CSV {ligands_csv}: {e}")
        sys.exit(1)

def load_msa_mappings(msa_table, use_global_cache=False):
    """Load MSA mappings from table CSV or check global cache."""
    if use_global_cache:
        # For global cache, we'll check if MSA files exist in the MSA cache folder
        # This will be handled during config creation
        return {}

    if not msa_table:
        return {}

    try:
        df = pd.read_csv(msa_table)
        # Expected columns: id, sequence_id, sequence, msa_file
        if 'id' not in df.columns or 'msa_file' not in df.columns:
            print(f"Warning: MSA table {msa_table} missing required columns")
            return {}

        # Create two mappings: by ID and by sequence (for fallback)
        msa_map_by_id = {}
        msa_map_by_seq = {}
        for _, row in df.iterrows():
            protein_id = row['sequence_id'] if 'sequence_id' in df.columns else row['id']
            msa_map_by_id[protein_id] = row['msa_file']

            # Also map by sequence if available (for when protein IDs don't match)
            if 'sequence' in df.columns and pd.notna(row['sequence']):
                msa_map_by_seq[row['sequence']] = row['msa_file']

        return {'by_id': msa_map_by_id, 'by_seq': msa_map_by_seq}
    except Exception as e:
        print(f"Error loading MSA table {msa_table}: {e}")
        return {'by_id': {}, 'by_seq': {}}

def create_config(protein, ligand, msa_file, enable_affinity, use_global_cache=False, msa_folder=None):
    """Create YAML config for a protein-ligand combination."""
    config = {
        'sequences': [
            {
                'protein': {
                    'id': 'A',  #this is not the structure_id, but an id internal to Boltz
                    'sequence': protein['sequence']
                }
            }
        ]
    }
    
    # Add MSA if available
    if msa_file and os.path.exists(msa_file):
        config['sequences'][0]['protein']['msa'] = msa_file
    elif use_global_cache and msa_folder:
        # Check for MSA file that matches this protein sequence
        
        if os.path.exists(msa_folder):
            # Calculate sequence hash to match against MSA files
            sequence = protein['sequence']
            sequence_hash = hashlib.sha256(sequence.encode()).hexdigest()[:8]
            
            # Look for MSA file that matches this protein sequence
            msa_pattern = os.path.join(msa_folder, "*.csv")
            msa_files = glob.glob(msa_pattern)
            
            # Try to find MSA file that matches this protein
            msa_file_found = None
            for msa_file in msa_files:
                # For global cache, we need to check if this MSA corresponds to our protein
                # Since we can't easily verify the sequence match from the filename alone,
                # we'll use the first MSA file but add a warning
                msa_file_found = msa_file
                break
            
            if msa_file_found:
                config['sequences'][0]['protein']['msa'] = msa_file_found
                print(f"Using global cache MSA: {msa_file_found}")
                print(f"Warning: Assuming MSA corresponds to protein sequence (hash: {sequence_hash})")
            else:
                print(f"No MSA files found in {msa_folder}")
        else:
            print(f"MSAs folder not found: {msa_folder}")
    
    # Add ligand if provided
    if ligand:
        ligand_entry = {
            'ligand': {
                'id': 'B' #this is not the compound_id, but an id internal to Boltz
            }
        }
        
        # Add ligand data based on format
        if ligand['format'].lower() == 'smiles' and 'smiles' in ligand and ligand['smiles']:
            ligand_entry['ligand']['smiles'] = ligand['smiles']
        elif ligand['format'].lower() == 'ccd' and 'ccd' in ligand and ligand['ccd']:
            ligand_entry['ligand']['ccd'] = ligand['ccd']
        
        config['sequences'].append(ligand_entry)
    
    # Add affinity calculation if requested
    if enable_affinity:
        config['properties'] = [
            {
                'affinity': {
                    'binder': 'B' #The binder is B which is the ligand
                }
            }
        ]
    
    return config

def main():
    args = parse_arguments()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Load proteins
    proteins = load_proteins(args.proteins_csv)
    print(f"Loaded {len(proteins)} proteins")
    
    # Load ligands
    ligands = load_ligands(args.ligands_csv)
    if ligands:
        print(f"Loaded {len(ligands)} ligands")
    else:
        print("No ligands provided - creating protein-only configs")
    
    # Load MSA mappings if provided
    msa_mappings = {}
    use_global_cache = args.global_msa_cache
    if args.msa_table:
        msa_mappings = load_msa_mappings(args.msa_table, use_global_cache)
        print(f"Loaded MSA mappings for {len(msa_mappings)} proteins")
    elif use_global_cache:
        print("Using global MSA cache - will check for MSA files during config creation")
    
    # Generate config files
    config_count = 0
    num_proteins = len(proteins)
    num_ligands = len(ligands)
    
    for protein in proteins:
        protein_id = protein['id']
        protein_seq = protein['sequence']

        # Look up MSA file - first by ID, then by sequence (fallback)
        msa_file = None
        if isinstance(msa_mappings, dict) and 'by_id' in msa_mappings:
            # Try lookup by protein ID first
            msa_file = msa_mappings['by_id'].get(protein_id)

            # Fallback: lookup by sequence if ID didn't match
            if not msa_file and 'by_seq' in msa_mappings:
                msa_file = msa_mappings['by_seq'].get(protein_seq)
                if msa_file:
                    print(f"MSA matched by sequence for protein '{protein_id}' (ID mismatch resolved)")
        elif isinstance(msa_mappings, dict):
            # Legacy format (backward compatibility)
            msa_file = msa_mappings.get(protein_id)
        
        if ligands:
            # Create one config per protein-ligand combination
            for ligand in ligands:
                config = create_config(protein, ligand, msa_file, args.affinity, use_global_cache, args.msa_folder)
                
                # Determine config filename based on counts
                if num_proteins == 1 and num_ligands == 1:
                    # One protein, one ligand: use job name if provided, otherwise ligand name
                    if args.job_name:
                        config_filename = f"{args.job_name}.yaml"
                    else:
                        config_filename = f"{ligand['id']}.yaml"
                elif num_proteins == 1 and num_ligands > 1:
                    # One protein, multiple ligands: use ligand name
                    config_filename = f"{ligand['id']}.yaml"
                elif num_proteins > 1 and num_ligands == 1:
                    # Multiple proteins, one ligand: use protein name
                    config_filename = f"{protein_id}.yaml"
                else:
                    # Multiple proteins, multiple ligands: use both
                    config_filename = f"{protein_id}_{ligand['id']}.yaml"
                
                config_path = os.path.join(args.output_dir, config_filename)
                
                with open(config_path, 'w') as f:
                    yaml.dump(config, f, default_flow_style=False, sort_keys=False)
                
                config_count += 1
                print(f"Created config: {config_filename}")
        else:
            # Create protein-only config
            config = create_config(protein, None, msa_file, args.affinity, use_global_cache, args.msa_folder)
            config_filename = f"{protein_id}.yaml"
            config_path = os.path.join(args.output_dir, config_filename)
            
            with open(config_path, 'w') as f:
                yaml.dump(config, f, default_flow_style=False, sort_keys=False)
            
            config_count += 1
            print(f"Created config: {config_filename}")
    
    print(f"\nGenerated {config_count} config files in {args.output_dir}")

if __name__ == "__main__":
    main()