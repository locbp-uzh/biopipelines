#!/usr/bin/env python3
# Copyright (c) 2026 Gianluca Quargnali @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

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

# Import shared ID prediction from biopipelines
_biopipelines_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'biopipelines')
sys.path.insert(0, _biopipelines_dir)
from combinatorics import predict_single_output_id, CombinatoricsConfig


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
    # Contact constraint parameters
    parser.add_argument('--contacts', help='JSON list of contact constraints')
    # Sequences output
    parser.add_argument('--sequences-csv', help='Path to write protein sequences CSV (id, sequence)')

    return parser.parse_args()


def load_combinatorics_config(config_path: str) -> Dict:
    """Load combinatorics configuration from JSON file."""
    with open(config_path, 'r') as f:
        return json.load(f)


def load_axis_data(axis_config: Dict) -> tuple:
    """
    Load data from an axis's source CSV files.

    Args:
        axis_config: Dict with 'name', 'mode', 'sources' keys
                    sources can be list of strings or list of dicts with 'path', 'iterate', and 'order'

    Returns:
        Tuple of (iterated_data, static_data, static_first):
        - iterated_data: List of dicts from sources with iterate=True
        - static_data: List of dicts from sources with iterate=False
        - static_first: True if static sources should be added before iterated (based on order)
    """
    iterated_data = []
    static_data = []
    sources = axis_config.get('sources', [])

    # Track order to determine if static should come first
    min_iterated_order = float('inf')
    min_static_order = float('inf')

    for source in sources:
        # Handle both old string format and new dict format
        if isinstance(source, dict):
            source_path = source.get('path')
            is_iterate = source.get('iterate', True)
            order = source.get('order', 0)
        else:
            source_path = source
            is_iterate = True
            order = 0

        if not source_path or not os.path.exists(source_path):
            print(f"Warning: Source file not found: {source_path}")
            continue
        try:
            df = pd.read_csv(source_path)
            records = df.to_dict('records')
            if is_iterate:
                iterated_data.extend(records)
                min_iterated_order = min(min_iterated_order, order)
            else:
                static_data.extend(records)
                min_static_order = min(min_static_order, order)
        except Exception as e:
            print(f"Error loading {source_path}: {e}")
            sys.exit(1)

    # Static comes first if its minimum order is less than iterated's minimum order
    static_first = min_static_order < min_iterated_order

    return iterated_data, static_data, static_first


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


def add_contacts_to_config(config: Dict, contacts: List[Dict]) -> Dict:
    """Add contact constraints to config."""
    if not contacts:
        return config

    if 'constraints' not in config:
        config['constraints'] = []

    # Determine which chain IDs are ligands (non-polymer)
    ligand_chain_ids = set()
    for seq in config.get('sequences', []):
        if 'ligand' in seq:
            ligand_chain_ids.add(seq['ligand']['id'])

    for contact in contacts:
        for key in ('token1', 'token2'):
            token = contact[key]
            chain_id = token[0]
            token_value = token[1]
            if chain_id in ligand_chain_ids and isinstance(token_value, int):
                raise ValueError(
                    f"Contact constraint {key}={token}: chain '{chain_id}' is a ligand, "
                    f"so the second element must be an atom name (string, e.g. 'C1'), "
                    f"not a residue index (integer). "
                    f"Check the CIF file of your ligand on RCSB to find standardized atom names, "
                    f"or run a prediction without constraints first and inspect the output."
                )

        entry = {
            'contact': {
                'token1': FlowList(contact['token1']),
                'token2': FlowList(contact['token2']),
            }
        }
        if 'max_distance' in contact:
            entry['contact']['max_distance'] = contact['max_distance']
        if 'force' in contact:
            entry['contact']['force'] = contact['force']
        config['constraints'].append(entry)

    print(f"Added {len(contacts)} contact constraint(s)")
    return config


def get_ligand_identity(ligand_entry: Dict) -> Optional[str]:
    """Get the identity string (smiles or ccd) of a ligand entry for duplicate detection."""
    lig = ligand_entry.get('ligand', {})
    if lig.get('smiles'):
        return f"smiles:{lig['smiles']}"
    if lig.get('ccd'):
        return f"ccd:{lig['ccd']}"
    return None


def has_duplicate_ligand(config: Dict, binder_chain_id: str) -> bool:
    """
    Check if the affinity binder ligand has duplicate copies in the config.

    Boltz does not allow affinity calculation for a ligand that appears
    multiple times in the same config. Returns True if duplicates are found.
    """
    # Find the binder ligand entry and its identity
    binder_identity = None
    for entry in config.get('sequences', []):
        if 'ligand' in entry and entry['ligand'].get('id') == binder_chain_id:
            binder_identity = get_ligand_identity(entry)
            break

    if binder_identity is None:
        return False

    # Count how many ligands share the same identity
    count = 0
    for entry in config.get('sequences', []):
        if 'ligand' in entry and get_ligand_identity(entry) == binder_identity:
            count += 1

    return count > 1


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


def build_dna_entry(dna: Dict, chain_id: str) -> Dict:
    """Build a DNA entry for Boltz2 config."""
    return {
        'dna': {
            'id': chain_id,
            'sequence': dna['sequence']
        }
    }


def build_rna_entry(rna: Dict, chain_id: str) -> Dict:
    """Build an RNA entry for Boltz2 config."""
    return {
        'rna': {
            'id': chain_id,
            'sequence': rna['sequence']
        }
    }


def build_ligand_entry(ligand: Dict, chain_id: str) -> Dict:
    """Build a ligand entry for Boltz2 config."""
    entry = {'ligand': {'id': chain_id}}

    fmt = ligand.get('format', '').lower()
    if fmt == 'smiles' and ligand.get('smiles'):
        entry['ligand']['smiles'] = ligand['smiles']
    elif fmt == 'ccd' and ligand.get('ccd'):
        entry['ligand']['ccd'] = ligand['ccd']

    return entry


# Entity types that support MSA lookup
MSA_ENTITY_TYPES = {'protein'}

# Entity types that are ligand-like (for affinity/pocket/covalent)
LIGAND_ENTITY_TYPES = {'ligand'}

# Backward-compatible mapping for old configs without entity_type
AXIS_NAME_TO_ENTITY_TYPE = {
    'proteins': 'protein',
    'sequences': 'protein',
    'dna': 'dna',
    'rna': 'rna',
    'ligands': 'ligand',
    'compounds': 'ligand',
}


def build_entry(entity_type: str, item: Dict, chain_id: str, msa_mappings: Dict) -> Dict:
    """Build a YAML sequence entry based on entity_type."""
    if entity_type == 'protein':
        msa_file = get_msa_file(item['id'], item.get('sequence', ''), msa_mappings)
        return build_protein_entry(item, chain_id, msa_file)
    elif entity_type == 'dna':
        return build_dna_entry(item, chain_id)
    elif entity_type == 'rna':
        return build_rna_entry(item, chain_id)
    elif entity_type == 'ligand':
        return build_ligand_entry(item, chain_id)
    else:
        raise ValueError(f"Unknown entity_type: {entity_type}")


def generate_configs(axis_data: Dict[str, Dict], msa_mappings: Dict, args) -> List[tuple]:
    """
    Generate all config dicts for N axes with independent Each/Bundle modes.

    For each axis:
    - 'each' mode: iterate over items (one per config)
    - 'bundle' mode without nested Each: all items static (appear in every config)
    - 'bundle' mode with nested Each: iterate over Each items, static items appear in every config

    The cartesian product is taken across all axes that have iteration.

    Args:
        axis_data: Dict mapping axis_name -> {entity_type, mode, iterated, static, static_first}
        msa_mappings: MSA file mappings for protein entries
        args: Command line arguments

    Returns:
        List of (config_id, config_dict) tuples.
    """
    import itertools

    glycosylation = json.loads(args.glycosylation) if args.glycosylation else None
    covalent_linkage = json.loads(args.covalent_linkage) if args.covalent_linkage else None
    contacts = json.loads(args.contacts) if args.contacts else None

    # Classify axes into iteration axes vs static-only axes
    axis_order = list(axis_data.keys())
    iteration_axes = []  # Axes contributing to the cartesian product
    static_only_axes = []  # Pure bundle axes (everything is static)

    for axis_name in axis_order:
        ad = axis_data[axis_name]
        entity_type = ad['entity_type']
        mode = ad['mode']
        iterated = ad['iterated']
        static = ad['static']
        static_first = ad['static_first']

        if mode == 'each':
            items_to_iterate = iterated if iterated else static
            static_items = static if iterated else []
            all_ids = [item['id'] for item in items_to_iterate]
            iteration_axes.append({
                'axis_name': axis_name, 'entity_type': entity_type,
                'items': items_to_iterate, 'static': static_items,
                'static_first': static_first, 'mode': mode, 'all_ids': all_ids,
            })
        elif mode == 'bundle':
            if iterated:
                # Bundle with nested Each: iterate over iterated items
                all_ids = [item['id'] for item in iterated]
                iteration_axes.append({
                    'axis_name': axis_name, 'entity_type': entity_type,
                    'items': iterated, 'static': static,
                    'static_first': static_first, 'mode': mode, 'all_ids': all_ids,
                })
            else:
                # Pure bundle: everything is static
                static_only_axes.append({
                    'axis_name': axis_name, 'entity_type': entity_type,
                    'items': static, 'static_first': static_first,
                })

    def next_chain_id(counter):
        cid = chr(65 + counter[0])
        counter[0] += 1
        return cid

    def add_axis_items_to_config(config, entity_type, items, chain_counter):
        """Add all items from an axis to the config, tracking first ligand chain."""
        first_ligand_chain = None
        for item in items:
            chain_id = next_chain_id(chain_counter)
            config['sequences'].append(build_entry(entity_type, item, chain_id, msa_mappings))
            if entity_type in LIGAND_ENTITY_TYPES and first_ligand_chain is None:
                first_ligand_chain = chain_id
        return first_ligand_chain

    def apply_decorations(config, first_ligand_chain):
        """Apply template, glycosylation, affinity, covalent, pocket, contacts."""
        config = add_template_to_config(config, args)
        config = add_glycosylation_to_config(config, glycosylation)
        if first_ligand_chain:
            if args.affinity:
                if has_duplicate_ligand(config, first_ligand_chain):
                    print(f"Warning: Affinity binder ligand (chain {first_ligand_chain}) has duplicates - skipping affinity")
                else:
                    config['properties'] = [{'affinity': {'binder': first_ligand_chain}}]
            config = add_covalent_linkage_to_config(config, covalent_linkage, first_ligand_chain)
            config = add_pocket_to_config(config, args, first_ligand_chain)
        config = add_contacts_to_config(config, contacts)
        return config

    # If no iteration axes, single config with everything bundled
    if not iteration_axes:
        config = {'sequences': []}
        chain_counter = [0]
        first_ligand_chain = None
        for sa in static_only_axes:
            flc = add_axis_items_to_config(config, sa['entity_type'], sa['items'], chain_counter)
            if first_ligand_chain is None:
                first_ligand_chain = flc
        config = apply_decorations(config, first_ligand_chain)
        return [('bundled_complex', config)]

    # Compute cartesian product of all iteration axes
    item_lists = [list(enumerate(ia['items'])) for ia in iteration_axes]
    product = list(itertools.product(*item_lists))

    configs = []
    for combo in product:
        config = {'sequences': []}
        chain_counter = [0]
        first_ligand_chain = None
        axis_selections = {}

        # Add entries from iteration axes
        for i, ia in enumerate(iteration_axes):
            idx, item = combo[i]
            axis_name = ia['axis_name']
            entity_type = ia['entity_type']
            static_items = ia['static']
            static_first = ia['static_first']

            axis_selections[axis_name] = (ia['mode'], ia['all_ids'], idx)

            if ia['mode'] == 'bundle' and static_first and static_items:
                # Static items first (for Bundle(cofactor, Each(library)) pattern)
                flc = add_axis_items_to_config(config, entity_type, static_items, chain_counter)
                if first_ligand_chain is None:
                    first_ligand_chain = flc
                # Then iterated item
                chain_id = next_chain_id(chain_counter)
                config['sequences'].append(build_entry(entity_type, item, chain_id, msa_mappings))
                if entity_type in LIGAND_ENTITY_TYPES and first_ligand_chain is None:
                    first_ligand_chain = chain_id
            else:
                # Iterated item first
                chain_id = next_chain_id(chain_counter)
                config['sequences'].append(build_entry(entity_type, item, chain_id, msa_mappings))
                if entity_type in LIGAND_ENTITY_TYPES and first_ligand_chain is None:
                    first_ligand_chain = chain_id
                # Then static items
                if static_items:
                    flc = add_axis_items_to_config(config, entity_type, static_items, chain_counter)
                    if first_ligand_chain is None:
                        first_ligand_chain = flc

        # Add entries from static-only (pure bundle) axes
        for sa in static_only_axes:
            axis_selections[sa['axis_name']] = ('bundle', [s['id'] for s in sa['items']], None)
            flc = add_axis_items_to_config(config, sa['entity_type'], sa['items'], chain_counter)
            if first_ligand_chain is None:
                first_ligand_chain = flc

        config = apply_decorations(config, first_ligand_chain)

        config_id = predict_single_output_id(
            bundled_name="bundled_complex",
            **axis_selections
        )
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


def write_protein_sequences(configs: List[tuple], output_path: str):
    """Write protein sequences CSV with id and sequence columns.

    Uses config_id as the sequence ID to ensure consistency with
    generated_ids used by the rest of the pipeline.
    """
    with open(output_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['id', 'sequence'])
        seen = set()
        for config_id, config in configs:
            if config_id and config_id not in seen:
                # Extract protein sequence from the first protein entry in the config
                seq = ''
                for entry in config.get('sequences', []):
                    if 'protein' in entry:
                        seq = entry['protein'].get('sequence', '')
                        break
                writer.writerow([config_id, seq])
                seen.add(config_id)
    print(f"Created protein sequences file: {output_path}")


def main():
    args = parse_arguments()

    # Load combinatorics config (includes pre-computed predictions and provenance)
    comb_config = load_combinatorics_config(args.combinatorics_config)
    stored_predicted_ids = comb_config.get('predicted_ids', [])
    stored_provenance = comb_config.get('provenance', {})
    axes = comb_config.get('axes', {})

    if not axes:
        print("Error: no axes found in combinatorics config")
        sys.exit(1)

    # Load all axes generically
    axis_data = {}
    for axis_name, axis_config in axes.items():
        # Resolve entity_type: use stored value, fall back to name-based mapping
        entity_type = axis_config.get('entity_type')
        if entity_type is None:
            entity_type = AXIS_NAME_TO_ENTITY_TYPE.get(axis_name)
        if entity_type is None:
            print(f"Error: cannot determine entity_type for axis '{axis_name}'. "
                  f"Known axis names: {list(AXIS_NAME_TO_ENTITY_TYPE.keys())}")
            sys.exit(1)

        mode = axis_config.get('mode', 'each')
        iterated, static, static_first = load_axis_data(axis_config)

        print(f"Axis '{axis_name}' (entity_type={entity_type}, mode={mode}): "
              f"{len(iterated)} iterated, {len(static)} static"
              f"{', static first' if static_first else ''}")

        axis_data[axis_name] = {
            'entity_type': entity_type,
            'mode': mode,
            'iterated': iterated,
            'static': static,
            'static_first': static_first,
        }

    # Load MSA mappings
    msa_mappings = load_msa_mappings(args.msa_table)
    if msa_mappings['by_id']:
        print(f"Loaded MSA mappings for {len(msa_mappings['by_id'])} proteins")

    # Generate configs
    configs = generate_configs(axis_data, msa_mappings, args)

    # Write configs
    generated_ids = write_configs(configs, args.output_dir)

    # Write sequence_ids.csv
    write_sequence_ids(generated_ids, args.output_dir)

    # Write protein sequences CSV with id and sequence columns
    if args.sequences_csv:
        write_protein_sequences(configs, args.sequences_csv)

    print(f"\nGenerated {len(configs)} config files")


if __name__ == "__main__":
    main()
