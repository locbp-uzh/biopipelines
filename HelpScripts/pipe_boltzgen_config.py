#!/usr/bin/env python3
# Copyright (c) 2026 Gianluca Quargnali @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Generate BoltzGen YAML design specification at runtime.

Reads ligand information from compounds.csv (from Ligand tool) and generates
the appropriate YAML specification for BoltzGen binder design.

Supports both CCD codes and SMILES strings based on the format column.
"""

import argparse
import pandas as pd
import yaml
import os
import sys


def parse_arguments():
    parser = argparse.ArgumentParser(description='Generate BoltzGen design specification YAML')
    parser.add_argument('--ligand-csv', required=True, help='Path to ligand compounds.csv file')
    parser.add_argument('--output-yaml', required=True, help='Path to output YAML file')
    parser.add_argument('--binder-length-min', type=int, required=True, help='Minimum binder length')
    parser.add_argument('--binder-length-max', type=int, required=True, help='Maximum binder length')
    parser.add_argument('--target-structure', help='Path to target protein structure (for protein-anything mode)')
    parser.add_argument('--target-chain', default='A', help='Chain ID to include from target structure')
    parser.add_argument('--binding-region', help='Target binding region (e.g., A:9-140)')

    return parser.parse_args()


def load_ligand(ligand_csv):
    """
    Load ligand information from compounds.csv.

    Expected columns: id, format, code, lookup, source, ccd, cid, cas, smiles, name, formula, file_path

    Returns:
        dict with keys: id, format, smiles, ccd
    """
    try:
        df = pd.read_csv(ligand_csv)

        if len(df) == 0:
            raise ValueError("Ligand CSV is empty")

        # Take the first ligand
        row = df.iloc[0]

        ligand_info = {
            'id': row.get('id', 'ligand'),
            'format': row.get('format', 'smiles').lower(),
            'smiles': row.get('smiles', ''),
            'ccd': row.get('ccd', row.get('code', ''))  # Try 'ccd' then 'code'
        }

        # Validate we have the necessary data
        if ligand_info['format'] == 'smiles' and not ligand_info['smiles']:
            raise ValueError("Ligand format is 'smiles' but no SMILES string found")
        if ligand_info['format'] == 'ccd' and not ligand_info['ccd']:
            raise ValueError("Ligand format is 'ccd' but no CCD code found")

        return ligand_info

    except Exception as e:
        print(f"Error loading ligand CSV {ligand_csv}: {e}")
        sys.exit(1)


def create_design_spec(ligand_info, binder_min, binder_max, target_structure=None,
                       target_chain='A', binding_region=None):
    """
    Create BoltzGen design specification YAML.

    For protein-small_molecule protocol (ligand only, no target protein):
    entities:
      - protein:
          id: A
          sequence: 100..200
      - ligand:
          id: B
          smiles: "..." or ccd: ATP

    For protein-anything protocol (with target protein structure):
    entities:
      - file:
          path: /path/to/target.pdb
          include:
            - chain:
                id: A
      - protein:
          id: binder
          sequence: 100..200
      - ligand:  # optional
          ...
    """
    entities = []

    if target_structure:
        # Mode: protein binder design against a protein target
        # Add target structure
        target_entity = {
            'file': {
                'path': target_structure,
                'include': [{'chain': {'id': target_chain}}]
            }
        }
        entities.append(target_entity)

        # Add binder protein
        binder_entity = {
            'protein': {
                'id': 'binder',
                'sequence': f'{binder_min}..{binder_max}'
            }
        }
        entities.append(binder_entity)

    else:
        # Mode: protein-small_molecule (binder design against ligand)
        # Add binder protein first
        binder_entity = {
            'protein': {
                'id': 'A',
                'sequence': f'{binder_min}..{binder_max}'
            }
        }
        entities.append(binder_entity)

    # Add ligand entity
    ligand_entity = {
        'ligand': {
            'id': 'B' if not target_structure else 'LIG'
        }
    }

    # Add smiles or ccd based on format
    if ligand_info['format'] == 'smiles':
        ligand_entity['ligand']['smiles'] = ligand_info['smiles']
    elif ligand_info['format'] == 'ccd':
        ligand_entity['ligand']['ccd'] = ligand_info['ccd']
    else:
        # Default to smiles if available
        if ligand_info['smiles']:
            ligand_entity['ligand']['smiles'] = ligand_info['smiles']
        elif ligand_info['ccd']:
            ligand_entity['ligand']['ccd'] = ligand_info['ccd']
        else:
            raise ValueError(f"No valid ligand representation found (format: {ligand_info['format']})")

    entities.append(ligand_entity)

    # Build full specification
    spec = {'entities': entities}

    # Add bindings section if binding region is specified
    if binding_region and target_structure:
        if ':' in binding_region:
            chain, region = binding_region.split(':', 1)
            spec['bindings'] = [{
                'binder_chain': 'binder',
                'target_chains': [chain],
                'target_region': region
            }]

    return spec


def main():
    args = parse_arguments()

    print("="*60)
    print("BoltzGen Config Generation")
    print("="*60)

    # Load ligand information
    print(f"\nLoading ligand from: {args.ligand_csv}")
    ligand_info = load_ligand(args.ligand_csv)
    print(f"  ID: {ligand_info['id']}")
    print(f"  Format: {ligand_info['format']}")
    if ligand_info['format'] == 'smiles':
        smiles_preview = ligand_info['smiles'][:50] + '...' if len(ligand_info['smiles']) > 50 else ligand_info['smiles']
        print(f"  SMILES: {smiles_preview}")
    else:
        print(f"  CCD: {ligand_info['ccd']}")

    # Create design specification
    print(f"\nBinder length range: {args.binder_length_min}-{args.binder_length_max}")
    if args.target_structure:
        print(f"Target structure: {args.target_structure}")
        print(f"Target chain: {args.target_chain}")

    spec = create_design_spec(
        ligand_info,
        args.binder_length_min,
        args.binder_length_max,
        target_structure=args.target_structure,
        target_chain=args.target_chain,
        binding_region=args.binding_region
    )

    # Write YAML file
    output_dir = os.path.dirname(args.output_yaml)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    with open(args.output_yaml, 'w') as f:
        yaml.dump(spec, f, default_flow_style=False, sort_keys=False)

    print(f"\nGenerated design specification: {args.output_yaml}")
    print("\nContents:")
    print("-"*40)
    with open(args.output_yaml, 'r') as f:
        print(f.read())
    print("-"*40)

    print("\nBoltzGen config generation completed successfully")
    return 0


if __name__ == "__main__":
    sys.exit(main())
