#!/usr/bin/env python3

"""
Script to convert FASTA files to individual Boltz2 YAML configuration files.
Each sequence gets its own config file with proper chain IDs.

Usage:
    python pipe_boltz_fasta_to_configs.py FASTA_FILES CONFIG_DIR LIGAND_SMILES [AFFINITY]

Args:
    FASTA_FILES: Comma-separated list of FASTA file paths
    CONFIG_DIR: Directory where individual YAML config files will be created
    LIGAND_SMILES: SMILES string for ligand (optional, use "None" for protein-only)
    AFFINITY: Enable affinity calculation (optional, default: True)
"""

import argparse
import os
import yaml
from Bio import SeqIO

def main():
    parser = argparse.ArgumentParser(description='Convert FASTA files to Boltz2 YAML configs')
    parser.add_argument('fasta_files', help='Comma-separated list of FASTA file paths')
    parser.add_argument('config_dir', help='Output directory for config files')
    parser.add_argument('ligand_smiles', help='Ligand SMILES string or "None"')
    parser.add_argument('--affinity', action='store_true', default=True, help='Enable affinity calculation')
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.config_dir, exist_ok=True)
    
    # Parse FASTA file list
    fasta_file_list = [f.strip().strip("'\"") for f in args.fasta_files.split(',')]
    
    # Process each FASTA file
    configs_created = 0
    for fasta_file in fasta_file_list:
        try:
            if not os.path.exists(fasta_file):
                print(f"Warning: FASTA file not found: {fasta_file}")
                continue
                
            for record in SeqIO.parse(fasta_file, 'fasta'):
                try:
                    # Create YAML data for single sequence with proper chain ID
                    yaml_data = {'sequences': []}
                    
                    # Add protein with chain ID A
                    protein_entry = {
                        'protein': {
                            'id': 'A',  # Proper chain ID
                            'sequence': str(record.seq)
                        }
                    }
                    yaml_data['sequences'].append(protein_entry)
                    
                    # Add ligand if specified (with chain ID B)
                    if args.ligand_smiles and args.ligand_smiles.upper() != "NONE":
                        ligand_entry = {
                            'ligand': {
                                'id': 'B',  # Proper chain ID for ligand
                                'smiles': args.ligand_smiles
                            }
                        }
                        yaml_data['sequences'].append(ligand_entry)
                        
                        # Add affinity calculation if enabled
                        if args.affinity:
                            if 'properties' not in yaml_data:
                                yaml_data['properties'] = []
                            yaml_data['properties'].append({'affinity': {'binder': 'B'}})
                    
                    # Write individual config file named after record ID
                    config_file = os.path.join(args.config_dir, f"{record.id}.yaml")
                    with open(config_file, 'w') as f:
                        yaml.dump(yaml_data, f, default_flow_style=False)
                    
                    print(f"Created config: {config_file}")
                    configs_created += 1
                    
                except Exception as e:
                    print(f"Error processing sequence {record.id}: {e}")
                    continue
                    
        except Exception as e:
            print(f"Error reading FASTA file {fasta_file}: {e}")
            continue
    
    print(f"Successfully created {configs_created} config files in {args.config_dir}")
    return 0

if __name__ == "__main__":
    exit(main())