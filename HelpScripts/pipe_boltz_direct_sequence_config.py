#!/usr/bin/env python3
# Copyright (c) 2026 Gianluca Quargnali @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.


"""
Script to create a Boltz2 YAML configuration file from a direct protein sequence.

Usage:
    python pipe_boltz_direct_sequence_config.py SEQUENCE CONFIG_FILE LIGAND_SMILES [AFFINITY]

Args:
    SEQUENCE: Protein sequence string
    CONFIG_FILE: Output YAML config file path
    LIGAND_SMILES: SMILES string for ligand (optional, use "None" for protein-only)
    AFFINITY: Enable affinity calculation (optional, default: True)
"""

import argparse
import yaml
import os

def main():
    parser = argparse.ArgumentParser(description='Create Boltz2 YAML config from direct sequence')
    parser.add_argument('sequence', help='Protein sequence string')
    parser.add_argument('config_file', help='Output YAML config file path')
    parser.add_argument('ligand_smiles', help='Ligand SMILES string or "None"')
    parser.add_argument('--affinity', action='store_true', default=True, help='Enable affinity calculation')
    
    args = parser.parse_args()
    
    try:
        # Create YAML data for single sequence with proper chain ID
        yaml_data = {'sequences': []}
        
        # Add protein with chain ID A
        protein_entry = {
            'protein': {
                'id': 'A',  # Proper chain ID
                'sequence': args.sequence
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
        
        # Ensure output directory exists
        os.makedirs(os.path.dirname(args.config_file), exist_ok=True)
        
        # Write YAML config file
        with open(args.config_file, 'w') as f:
            yaml.dump(yaml_data, f, default_flow_style=False)
        
        print(f"Created config: {args.config_file}")
        return 0
        
    except Exception as e:
        print(f"Error creating config: {e}")
        return 1

if __name__ == "__main__":
    exit(main())