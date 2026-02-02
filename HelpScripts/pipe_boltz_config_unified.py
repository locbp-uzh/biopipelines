#!/usr/bin/env python3
"""
Unified Boltz2 YAML config file generator.

Generates Boltz2 YAML configuration files based on combinatorics config,
supporting all input combinations (Bundle/Each for proteins and ligands).

Replaces:
- pipe_boltz_config.py (old library replacement)
- pipe_build_boltz_config_with_msas.py
- pipe_build_boltz_config_with_global_cache.py
- pipe_boltz_protein_ligand_configs.py
"""

import argparse
import json
import os
import sys
import csv
from typing import Dict, List, Any, Optional

import pandas as pd
import yaml


# Custom YAML representer for inline list formatting in constraints
class FlowList(list):
    """Wrapper class to force flow-style (inline) representation for lists in YAML"""
    pass


def flow_list_representer(dumper, data):
    """YAML representer that renders lists in flow style (inline format)"""
    return dumper.represent_sequence('tag:yaml.org,2002:seq', data, flow_style=True)


yaml.add_representer(FlowList, flow_list_representer)


def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Generate Boltz2 config files based on combinatorics config'
    )
    parser.add_argument(
        '--combinatorics-config', required=True,
        help='Path to combinatorics config JSON file'
    )
    parser.add_argument(
        '--output-dir', required=True,
        help='Directory to write YAML config files'
    )
    parser.add_argument(
        '--msa-table',
        help='Path to MSA table CSV file (columns: id, sequence_id, sequence, msa_file)'
    )
    parser.add_argument(
        '--affinity', action='store_true',
        help='Enable affinity calculation'
    )
    # Template parameters
    parser.add_argument('--template', help='Path to PDB template file')
    parser.add_argument('--template-chains', help='Comma-separated chain IDs for template')
    parser.add_argument('--template-force', action='store_true', help='Force template usage')
    parser.add_argument('--template-threshold', type=float, default=5.0, help='Template RMSD threshold')
    # Pocket constraint parameters
    parser.add_argument('--pocket-residues', help='Pocket residues as Python list string')
    parser.add_argument('--pocket-max-distance', type=float, default=7.0, help='Max distance for pocket')
    parser.add_argument('--pocket-force', action='store_true', help='Force pocket constraint')
    # Glycosylation parameters
    parser.add_argument('--glycosylation', help='JSON dict mapping chain IDs to Asn positions')
    # Covalent linkage parameters
    parser.add_argument('--covalent-linkage', help='JSON dict for covalent attachment')

    return parser.parse_args()


def load_combinatorics_config(config_path: str) -> Dict:
    """Load combinatorics configuration from JSON file."""
    with open(config_path, 'r') as f:
        return json.load(f)


def load_axis_data(axis_config: Dict) -> List[Dict]:
    """
    Load data from an axis's source CSV files.

    Args:
        axis_config: Dict with 'name', 'mode', 'sources' keys

    Returns:
        List of dicts, each representing a row from the CSVs
    """
    all_data = []
    sources = axis_config.get('sources', [])

    for source_path in sources:
        if not source_path or not os.path.exists(source_path):
            print(f"Warning: Source file not found: {source_path}")
            continue
        try:
            df = pd.read_csv(source_path)
            all_data.extend(df.to_dict('records'))
        except Exception as e:
            print(f"Error loading {source_path}: {e}")
            sys.exit(1)

    return all_data


def load_msa_mappings(msa_table: Optional[str]) -> Dict:
    """
    Load MSA mappings from table CSV.

    Returns dict with 'by_id' and 'by_seq' sub-dicts for flexible matching.
    """
    if not msa_table or not os.path.exists(msa_table):
        return {'by_id': {}, 'by_seq': {}}

    try:
        df = pd.read_csv(msa_table)
        msa_map_by_id = {}
        msa_map_by_seq = {}

        for _, row in df.iterrows():
            seq_id = row.get('sequence_id', row.get('id', ''))
            msa_file = row.get('msa_file', '')
            sequence = row.get('sequence', '')

            if seq_id and msa_file:
                msa_map_by_id[seq_id] = msa_file
            if sequence and msa_file:
                msa_map_by_seq[sequence] = msa_file

        return {'by_id': msa_map_by_id, 'by_seq': msa_map_by_seq}
    except Exception as e:
        print(f"Error loading MSA table: {e}")
        return {'by_id': {}, 'by_seq': {}}


def get_msa_file(protein_id: str, sequence: str, msa_mappings: Dict) -> Optional[str]:
    """Look up MSA file for a protein by ID or sequence."""
    # Try by ID first
    msa_file = msa_mappings['by_id'].get(protein_id)
    if msa_file and os.path.exists(msa_file):
        return msa_file

    # Try by sequence
    msa_file = msa_mappings['by_seq'].get(sequence)
    if msa_file and os.path.exists(msa_file):
        print(f"MSA matched by sequence for protein '{protein_id}'")
        return msa_file

    return None


def add_template_to_config(config: Dict, args) -> Dict:
    """Add template section to config if specified."""
    if not args.template:
        return config

    template_entry = {
        'pdb': args.template,
        'force': args.template_force,
        'threshold': args.template_threshold
    }
    if args.template_chains:
        template_entry['chain_id'] = [c.strip() for c in args.template_chains.split(',')]

    config['templates'] = [template_entry]
    return config


def add_pocket_to_config(config: Dict, args, ligand_chain_id: str) -> Dict:
    """Add pocket constraint to config if specified."""
    if not args.pocket_residues:
        return config

    import ast
    pocket_residues = ast.literal_eval(args.pocket_residues)

    constraint = {
        'pocket': {
            'binder': ligand_chain_id,
            'contacts': pocket_residues,
            'max_distance': args.pocket_max_distance,
            'force': args.pocket_force
        }
    }

    if 'constraints' not in config:
        config['constraints'] = []
    config['constraints'].append(constraint)

    return config


def add_glycosylation_to_config(config: Dict, glycosylation: Dict) -> Dict:
    """Add NAG ligands and bond constraints for N-glycosylation."""
    if not glycosylation:
        return config

    existing_chain_count = len(config['sequences'])

    for chain_id, positions in glycosylation.items():
        for pos in positions:
            nag_chain_id = chr(65 + existing_chain_count)

            config['sequences'].append({
                'ligand': {
                    'id': nag_chain_id,
                    'ccd': 'NAG'
                }
            })

            if 'constraints' not in config:
                config['constraints'] = []
            config['constraints'].append({
                'bond': {
                    'atom1': FlowList([chain_id, pos, 'ND2']),
                    'atom2': FlowList([nag_chain_id, 1, 'C1'])
                }
            })

            existing_chain_count += 1
            print(f"Added NAG ligand (chain '{nag_chain_id}') linked to {chain_id}:{pos}")

    return config


def add_covalent_linkage_to_config(config: Dict, covalent_linkage: Dict, ligand_chain_id: str) -> Dict:
    """Add bond constraint for covalent ligand attachment."""
    if not covalent_linkage:
        return config

    chain = covalent_linkage.get('chain')
    position = covalent_linkage.get('position')
    protein_atom = covalent_linkage.get('protein_atom', 'SG')
    ligand_atom = covalent_linkage.get('ligand_atom')
    ligand_residue = covalent_linkage.get('ligand_residue', 1)

    if not all([chain, position, ligand_atom]):
        print("Warning: covalent_linkage missing required fields")
        return config

    if 'constraints' not in config:
        config['constraints'] = []
    config['constraints'].append({
        'bond': {
            'atom1': FlowList([chain, position, protein_atom]),
            'atom2': FlowList([ligand_chain_id, ligand_residue, ligand_atom])
        }
    })

    print(f"Added covalent linkage: {chain}:{position} {protein_atom} -> {ligand_chain_id}:{ligand_residue} {ligand_atom}")
    return config


def build_protein_entry(protein: Dict, chain_id: str, msa_file: Optional[str]) -> Dict:
    """Build a protein entry for Boltz2 config."""
    entry = {
        'protein': {
            'id': chain_id,
            'sequence': protein['sequence']
        }
    }
    if msa_file:
        entry['protein']['msa'] = msa_file
    return entry


def build_ligand_entry(ligand: Dict, chain_id: str) -> Dict:
    """Build a ligand entry for Boltz2 config."""
    entry = {'ligand': {'id': chain_id}}

    fmt = ligand.get('format', '').lower()
    if fmt == 'smiles' and ligand.get('smiles'):
        entry['ligand']['smiles'] = ligand['smiles']
    elif fmt == 'ccd' and ligand.get('ccd'):
        entry['ligand']['ccd'] = ligand['ccd']

    return entry


def generate_config_id(
    proteins_mode: str, ligands_mode: str,
    protein_ids: List[str], ligand_ids: List[str],
    protein_idx: Optional[int] = None, ligand_idx: Optional[int] = None
) -> str:
    """
    Generate config ID based on combinatorics modes and indices.

    Rules:
    - bundle + bundle: "bundled_complex"
    - bundle + each: ligand_id
    - each + bundle: protein_id
    - each + each:
      - 1 protein: ligand_id
      - 1 ligand: protein_id
      - multiple both: protein_id_ligand_id
    """
    if proteins_mode == 'bundle' and ligands_mode == 'bundle':
        return 'bundled_complex'
    elif proteins_mode == 'bundle' and ligands_mode == 'each':
        return ligand_ids[ligand_idx]
    elif proteins_mode == 'each' and ligands_mode == 'bundle':
        return protein_ids[protein_idx]
    else:  # each x each
        if len(protein_ids) == 1:
            return ligand_ids[ligand_idx]
        elif len(ligand_ids) == 1:
            return protein_ids[protein_idx]
        else:
            return f"{protein_ids[protein_idx]}_{ligand_ids[ligand_idx]}"


def generate_configs(
    proteins: List[Dict], ligands: List[Dict],
    proteins_mode: str, ligands_mode: str,
    msa_mappings: Dict, args
) -> List[tuple]:
    """
    Generate all config dictionaries based on combinatorics modes.

    Returns list of (config_id, config_dict) tuples.
    """
    # Parse extra parameters
    glycosylation = json.loads(args.glycosylation) if args.glycosylation else None
    covalent_linkage = json.loads(args.covalent_linkage) if args.covalent_linkage else None

    protein_ids = [p['id'] for p in proteins] if proteins else []
    ligand_ids = [l['id'] for l in ligands] if ligands else []

    configs = []

    # No ligands axis - proteins only (apo mode)
    if ligands_mode is None:
        if proteins_mode == 'bundle':
            # Single config with all proteins bundled
            config = {'sequences': []}
            for idx, protein in enumerate(proteins):
                chain_id = chr(65 + idx)
                msa_file = get_msa_file(protein['id'], protein['sequence'], msa_mappings)
                config['sequences'].append(build_protein_entry(protein, chain_id, msa_file))
            config = add_template_to_config(config, args)
            config = add_glycosylation_to_config(config, glycosylation)
            configs.append(('bundled_complex', config))
        else:  # proteins_mode == 'each'
            # One config per protein
            for prot_idx, protein in enumerate(proteins):
                config = {'sequences': []}
                msa_file = get_msa_file(protein['id'], protein['sequence'], msa_mappings)
                config['sequences'].append(build_protein_entry(protein, 'A', msa_file))
                config = add_template_to_config(config, args)
                config = add_glycosylation_to_config(config, glycosylation)
                configs.append((protein_ids[prot_idx], config))

    elif proteins_mode == 'bundle' and ligands_mode == 'bundle':
        # Single config with all proteins and ligands
        config = {'sequences': []}

        for idx, protein in enumerate(proteins):
            chain_id = chr(65 + idx)
            msa_file = get_msa_file(protein['id'], protein['sequence'], msa_mappings)
            config['sequences'].append(build_protein_entry(protein, chain_id, msa_file))

        config = add_template_to_config(config, args)
        config = add_glycosylation_to_config(config, glycosylation)

        first_ligand_chain = None
        for ligand in ligands:
            chain_id = chr(65 + len(config['sequences']))
            if first_ligand_chain is None:
                first_ligand_chain = chain_id
            config['sequences'].append(build_ligand_entry(ligand, chain_id))

        if args.affinity and first_ligand_chain:
            config['properties'] = [{'affinity': {'binder': first_ligand_chain}}]

        config = add_covalent_linkage_to_config(config, covalent_linkage, first_ligand_chain or 'C')
        if first_ligand_chain:
            config = add_pocket_to_config(config, args, first_ligand_chain)

        configs.append(('bundled_complex', config))

    elif proteins_mode == 'bundle' and ligands_mode == 'each':
        # One config per ligand, all proteins bundled
        for lig_idx, ligand in enumerate(ligands):
            config = {'sequences': []}

            for idx, protein in enumerate(proteins):
                chain_id = chr(65 + idx)
                msa_file = get_msa_file(protein['id'], protein['sequence'], msa_mappings)
                config['sequences'].append(build_protein_entry(protein, chain_id, msa_file))

            config = add_template_to_config(config, args)
            config = add_glycosylation_to_config(config, glycosylation)

            ligand_chain_id = chr(65 + len(config['sequences']))
            config['sequences'].append(build_ligand_entry(ligand, ligand_chain_id))

            if args.affinity:
                config['properties'] = [{'affinity': {'binder': ligand_chain_id}}]

            config = add_covalent_linkage_to_config(config, covalent_linkage, ligand_chain_id)
            config = add_pocket_to_config(config, args, ligand_chain_id)

            config_id = generate_config_id(proteins_mode, ligands_mode, protein_ids, ligand_ids, ligand_idx=lig_idx)
            configs.append((config_id, config))

    elif proteins_mode == 'each' and ligands_mode == 'bundle':
        # One config per protein, all ligands bundled
        for prot_idx, protein in enumerate(proteins):
            config = {'sequences': []}

            msa_file = get_msa_file(protein['id'], protein['sequence'], msa_mappings)
            config['sequences'].append(build_protein_entry(protein, 'A', msa_file))

            config = add_template_to_config(config, args)
            config = add_glycosylation_to_config(config, glycosylation)

            first_ligand_chain = None
            for ligand in ligands:
                chain_id = chr(65 + len(config['sequences']))
                if first_ligand_chain is None:
                    first_ligand_chain = chain_id
                config['sequences'].append(build_ligand_entry(ligand, chain_id))

            if args.affinity and first_ligand_chain:
                config['properties'] = [{'affinity': {'binder': first_ligand_chain}}]

            config = add_covalent_linkage_to_config(config, covalent_linkage, first_ligand_chain or 'B')
            if first_ligand_chain:
                config = add_pocket_to_config(config, args, first_ligand_chain)

            config_id = generate_config_id(proteins_mode, ligands_mode, protein_ids, ligand_ids, protein_idx=prot_idx)
            configs.append((config_id, config))

    else:  # each x each (cartesian product)
        for prot_idx, protein in enumerate(proteins):
            ligands_to_process = ligands if ligands else [None]

            for lig_idx, ligand in enumerate(ligands_to_process):
                config = {'sequences': []}

                msa_file = get_msa_file(protein['id'], protein['sequence'], msa_mappings)
                config['sequences'].append(build_protein_entry(protein, 'A', msa_file))

                config = add_template_to_config(config, args)
                config = add_glycosylation_to_config(config, glycosylation)

                ligand_chain_id = None
                if ligand:
                    ligand_chain_id = chr(65 + len(config['sequences']))
                    config['sequences'].append(build_ligand_entry(ligand, ligand_chain_id))

                    if args.affinity:
                        config['properties'] = [{'affinity': {'binder': ligand_chain_id}}]

                    config = add_covalent_linkage_to_config(config, covalent_linkage, ligand_chain_id)
                    config = add_pocket_to_config(config, args, ligand_chain_id)

                if ligand:
                    config_id = generate_config_id(proteins_mode, ligands_mode, protein_ids, ligand_ids, prot_idx, lig_idx)
                else:
                    config_id = protein_ids[prot_idx]

                configs.append((config_id, config))

    return configs


def write_configs(configs: List[tuple], output_dir: str) -> List[str]:
    """Write config files and return list of generated IDs."""
    os.makedirs(output_dir, exist_ok=True)
    generated_ids = []

    for config_id, config in configs:
        config_path = os.path.join(output_dir, f"{config_id}.yaml")
        with open(config_path, 'w') as f:
            yaml.dump(config, f, default_flow_style=False, sort_keys=False)
        generated_ids.append(config_id)
        print(f"Created config: {config_id}.yaml")

    return generated_ids


def write_sequence_ids(ids: List[str], output_dir: str):
    """Write sequence_ids.csv file."""
    ids_file = os.path.join(os.path.dirname(output_dir), "sequence_ids.csv")
    with open(ids_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['id'])
        for seq_id in ids:
            writer.writerow([seq_id])
    print(f"Created sequence IDs file: {ids_file}")


def main():
    args = parse_arguments()

    # Load combinatorics config
    comb_config = load_combinatorics_config(args.combinatorics_config)
    axes = comb_config.get('axes', {})

    # Extract sequences and compounds axes
    if 'sequences' not in axes:
        print("Error: sequences axis is required")
        sys.exit(1)
    proteins_axis = axes['sequences']
    ligands_axis = axes.get('compounds')  # None if no ligands

    proteins_mode = proteins_axis.get('mode', 'each')
    ligands_mode = ligands_axis.get('mode', 'each') if ligands_axis else None

    print(f"Combinatorics modes: proteins={proteins_mode}, ligands={ligands_mode}")

    # Load data
    proteins = load_axis_data(proteins_axis)
    ligands = load_axis_data(ligands_axis) if ligands_axis else []

    print(f"Loaded {len(proteins)} proteins, {len(ligands)} ligands")

    # Load MSA mappings
    msa_mappings = load_msa_mappings(args.msa_table)
    if msa_mappings['by_id']:
        print(f"Loaded MSA mappings for {len(msa_mappings['by_id'])} proteins")

    # Generate configs
    configs = generate_configs(
        proteins, ligands,
        proteins_mode, ligands_mode,
        msa_mappings, args
    )

    # Write configs
    generated_ids = write_configs(configs, args.output_dir)

    # Write sequence_ids.csv
    write_sequence_ids(generated_ids, args.output_dir)

    print(f"\nGenerated {len(configs)} config files")


if __name__ == "__main__":
    main()
