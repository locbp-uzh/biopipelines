#!/usr/bin/env python3
"""
Helper script to build Boltz2 YAML configuration with global MSA cache integration.

This script creates a Boltz2 configuration that uses cached MSA files from the global cache
when available, otherwise builds config without MSA paths (for calculation).
"""

import os
import sys
import argparse
import yaml

def build_config_with_global_cache(protein_sequence: str, ligand_smiles: str, 
                                  cached_msa_path: str, affinity: bool = False) -> str:
    """
    Build Boltz2 YAML configuration with global MSA cache integration.
    
    Args:
        protein_sequence: Protein amino acid sequence
        ligand_smiles: Ligand SMILES string
        cached_msa_path: Path to cached MSA file (may not exist)
        affinity: Whether to calculate binding affinity
        
    Returns:
        YAML configuration string
    """
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
    
    # Add MSA path if cached file exists
    if cached_msa_path and os.path.exists(cached_msa_path):
        config['sequences'][0]['protein']['msa'] = cached_msa_path
        print(f"Using cached MSA: {cached_msa_path}")
    else:
        print("No cached MSA found, Boltz2 will calculate MSAs")
    
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
    parser = argparse.ArgumentParser(description='Build Boltz2 config with global MSA cache')
    parser.add_argument('--protein-sequence', required=True, help='Protein amino acid sequence')
    parser.add_argument('--ligand-smiles', help='Ligand SMILES string')
    parser.add_argument('--cached-msa-path', help='Path to cached MSA file')
    parser.add_argument('--output', required=True, help='Output YAML config file')
    parser.add_argument('--affinity', help='Enable affinity calculation (true/false)')
    
    args = parser.parse_args()
    
    # Parse affinity flag
    affinity_flag = args.affinity and args.affinity.lower() == 'true'
    
    # Build configuration
    try:
        yaml_config = build_config_with_global_cache(
            protein_sequence=args.protein_sequence,
            ligand_smiles=args.ligand_smiles or "",
            cached_msa_path=args.cached_msa_path or "",
            affinity=affinity_flag
        )
        
        # Write to output file
        with open(args.output, 'w') as f:
            f.write(yaml_config)
        
        print(f"Boltz2 configuration written to: {args.output}")
        print("Configuration preview:")
        print(yaml_config)
        
    except Exception as e:
        print(f"Error building configuration: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()