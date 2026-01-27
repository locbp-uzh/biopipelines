#!/usr/bin/env python3
"""
Generate Boltz2 YAML config files for all protein-ligand combinations.

Takes protein CSV and ligand CSV files and generates one YAML config file
per protein-ligand combination, including MSA integration when available.

Supports:
- Single protein + single/multiple ligands
- Multiple proteins (complexes) + single/multiple ligands
- SplitChains output with complex_id grouping
- Templates for structure guidance
- Pocket constraints for ligand binding
- Glycosylation (N-linked NAG)
- Covalent ligand attachment
"""

import argparse
import pandas as pd
import yaml
import os
import sys
import glob
import hashlib


# Custom YAML representer for inline list formatting in constraints
class FlowList(list):
    """Wrapper class to force flow-style (inline) representation for lists in YAML"""
    pass


def flow_list_representer(dumper, data):
    """YAML representer that renders lists in flow style (inline format)"""
    return dumper.represent_sequence('tag:yaml.org,2002:seq', data, flow_style=True)


# Register the custom representer
yaml.add_representer(FlowList, flow_list_representer)


def parse_arguments():
    parser = argparse.ArgumentParser(description='Generate Boltz2 config files for protein-ligand combinations')
    parser.add_argument('proteins_csv', help='Path to proteins CSV file (id,sequence format) or file containing list of CSV paths (one per line)')
    parser.add_argument('ligands_csv', help='Path to ligands CSV file or "None" if no ligands')
    parser.add_argument('output_dir', help='Directory to write config files')
    parser.add_argument('--affinity', action='store_true', help='Enable affinity calculation')
    parser.add_argument('--msa-table', help='Path to MSA table CSV file')
    parser.add_argument('--global-msa-cache', action='store_true', help='Use global MSA cache lookup')
    parser.add_argument('--msa-folder', help='Path to MSAs folder for global cache lookup')
    parser.add_argument('--job-name', help='Pipeline job name to use for single protein-ligand configs')
    # Template parameters
    parser.add_argument('--template', help='Path to PDB template file')
    parser.add_argument('--template-chains', help='Comma-separated list of chain IDs to apply template to')
    parser.add_argument('--template-force', action='store_true', help='Force template usage')
    parser.add_argument('--template-threshold', type=float, default=5.0, help='RMSD threshold for template')
    # Pocket constraint parameters
    parser.add_argument('--pocket-residues', help='Pocket residues as Python list string')
    parser.add_argument('--pocket-max-distance', type=float, default=7.0, help='Maximum distance for pocket constraint')
    parser.add_argument('--pocket-force', action='store_true', help='Force pocket constraint')
    # Glycosylation parameters
    parser.add_argument('--glycosylation', help='JSON string mapping chain IDs to list of Asn positions for N-glycosylation')
    # Covalent linkage parameters
    parser.add_argument('--covalent-linkage', help='JSON string with keys: chain, position, protein_atom, ligand_atom for covalent attachment')

    return parser.parse_args()


def load_proteins(proteins_csv):
    """Load proteins from CSV file or from a file containing list of CSV paths."""
    try:
        # Check if the file is a list file (contains paths, one per line)
        with open(proteins_csv, 'r') as f:
            first_line = f.readline().strip()

        # If first line looks like a file path (not a CSV header), treat as list file
        if first_line and (first_line.endswith('.csv') or first_line.endswith('.fasta') or first_line.endswith('.fa')):
            with open(proteins_csv, 'r') as f:
                csv_paths = [line.strip() for line in f if line.strip()]
            if not csv_paths:
                raise ValueError(f"No CSV paths found in list file: {proteins_csv}")
            # Concatenate all CSV files
            all_proteins = []
            for csv_path in csv_paths:
                df = pd.read_csv(csv_path)
                if 'id' not in df.columns or 'sequence' not in df.columns:
                    raise ValueError(f"Proteins CSV {csv_path} must have 'id' and 'sequence' columns")
                all_proteins.extend(df.to_dict('records'))
            return all_proteins
        else:
            # Direct CSV file
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


def add_glycosylation_to_config(config, glycosylation, num_existing_chains):
    """
    Add NAG ligands and bond constraints for N-glycosylation.

    Args:
        config: The config dict to modify
        glycosylation: Dict mapping chain IDs to list of Asn positions (e.g., {"A": [164]})
        num_existing_chains: Number of chains already in config (used to assign chain IDs)

    Returns:
        Updated config dict
    """
    if not glycosylation:
        return config

    # Start assigning NAG chain IDs after existing chains
    existing_chain_count = len(config['sequences'])

    for chain_id, positions in glycosylation.items():
        for pos in positions:
            # Use single-character chain ID (D, E, F, ...) for Boltz2 compatibility
            nag_chain_id = chr(65 + existing_chain_count)  # 65 is ASCII for 'A'

            # Add NAG ligand with CCD code
            config['sequences'].append({
                'ligand': {
                    'id': nag_chain_id,
                    'ccd': 'NAG'
                }
            })

            # Add bond constraint: Asn ND2 to NAG C1
            # Use FlowList to render as inline format: [A, 164, ND2]
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


def add_covalent_linkage_to_config(config, covalent_linkage, ligand_chain_id='C'):
    """
    Add bond constraint for covalent attachment of ligand to protein residue.

    Args:
        config: The config dict to modify
        covalent_linkage: Dict with keys:
                         - "chain" (protein chain ID)
                         - "position" (residue position)
                         - "protein_atom" (atom name in protein, e.g., 'SG' for Cys, 'OD1'/'OD2' for Asp)
                         - "ligand_atom" (atom name in ligand)
                         - "ligand_residue" (optional, defaults to 1)
                         If "protein_atom" is not specified, defaults to 'SG' (cysteine) for backwards compatibility
        ligand_chain_id: The chain ID of the ligand to attach (default: 'C')

    Returns:
        Updated config dict
    """
    if not covalent_linkage:
        return config

    chain_id = covalent_linkage.get('chain')
    position = covalent_linkage.get('position')
    protein_atom = covalent_linkage.get('protein_atom', 'SG')  # Default to SG for cysteine
    ligand_atom = covalent_linkage.get('ligand_atom')
    ligand_residue = covalent_linkage.get('ligand_residue', 1)  # Default to residue 1

    if not all([chain_id, position, ligand_atom]):
        print("Warning: covalent_linkage missing required fields (chain, position, ligand_atom)")
        return config

    # Add bond constraint: protein atom to ligand atom
    # Use FlowList to render as inline format: [A, 50, SG]
    if 'constraints' not in config:
        config['constraints'] = []
    config['constraints'].append({
        'bond': {
            'atom1': FlowList([chain_id, position, protein_atom]),
            'atom2': FlowList([ligand_chain_id, ligand_residue, ligand_atom])
        }
    })

    print(f"Added covalent linkage: {chain_id}:{position} {protein_atom} -> {ligand_chain_id}:{ligand_residue} {ligand_atom}")

    return config


def create_config(protein, ligand, msa_file, enable_affinity, use_global_cache=False, msa_folder=None,
                  template=None, template_chains=None, template_force=True, template_threshold=5.0,
                  pocket_residues=None, pocket_max_distance=7.0, pocket_force=True):
    """Create YAML config for a protein-ligand combination."""
    config = {
        'sequences': [
            {
                'protein': {
                    'id': 'A',  # Internal Boltz ID, not the structure_id
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
        # Use "C" as ligand chain ID for consistency with multi-protein case
        ligand_entry = {
            'ligand': {
                'id': 'C'  # Fixed single-character ID
            }
        }

        # Add ligand data based on format
        if ligand['format'].lower() == 'smiles' and 'smiles' in ligand and ligand['smiles']:
            ligand_entry['ligand']['smiles'] = ligand['smiles']
        elif ligand['format'].lower() == 'ccd' and 'ccd' in ligand and ligand['ccd']:
            ligand_entry['ligand']['ccd'] = ligand['ccd']

        config['sequences'].append(ligand_entry)

    # Add template if specified
    if template:
        template_entry = {
            'pdb': template,
            'force': template_force,
            'threshold': template_threshold
        }
        if template_chains:
            # Parse chain IDs from comma-separated string or list
            if isinstance(template_chains, str):
                chain_ids = [c.strip() for c in template_chains.split(',')]
            else:
                chain_ids = template_chains
            template_entry['chain_id'] = chain_ids

        config['templates'] = [template_entry]

    # Add pocket constraints if specified
    if pocket_residues and ligand:
        # Parse pocket residues if it's a string representation
        if isinstance(pocket_residues, str):
            import ast
            pocket_residues = ast.literal_eval(pocket_residues)

        constraint_entry = {
            'pocket': {
                'binder': 'C',  # Use fixed ligand chain ID
                'contacts': pocket_residues,
                'max_distance': pocket_max_distance,
                'force': pocket_force
            }
        }

        if 'constraints' not in config:
            config['constraints'] = []
        config['constraints'].append(constraint_entry)

    # Add affinity calculation if requested
    if enable_affinity and ligand:
        config['properties'] = [
            {
                'affinity': {
                    'binder': 'C'  # Use fixed ligand chain ID
                }
            }
        ]

    return config


def main():
    args = parse_arguments()

    # Parse glycosylation JSON if provided
    glycosylation = None
    if args.glycosylation:
        import json
        try:
            glycosylation = json.loads(args.glycosylation)
            print(f"Loaded glycosylation: {glycosylation}")
        except json.JSONDecodeError as e:
            print(f"Error parsing glycosylation JSON: {e}")
            sys.exit(1)

    # Parse covalent linkage JSON if provided
    covalent_linkage = None
    if args.covalent_linkage:
        import json
        try:
            covalent_linkage = json.loads(args.covalent_linkage)
            print(f"Loaded covalent_linkage: {covalent_linkage}")
        except json.JSONDecodeError as e:
            print(f"Error parsing covalent linkage JSON: {e}")
            sys.exit(1)

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
        if isinstance(msa_mappings, dict) and 'by_id' in msa_mappings:
            print(f"Loaded MSA mappings for {len(msa_mappings['by_id'])} proteins")
    elif use_global_cache:
        print("Using global MSA cache - will check for MSA files during config creation")

    # Generate config files
    config_count = 0
    num_proteins = len(proteins)
    num_ligands = len(ligands)

    # Check if proteins have complex_id column (from SplitChains tool)
    proteins_df = pd.DataFrame(proteins)
    has_complex_id = 'complex_id' in proteins_df.columns

    # Helper to look up MSA file for a protein
    def get_msa_file(protein_id, protein_seq):
        msa_file = None
        if isinstance(msa_mappings, dict) and 'by_id' in msa_mappings:
            msa_file = msa_mappings['by_id'].get(protein_id)
            if not msa_file and 'by_seq' in msa_mappings:
                msa_file = msa_mappings['by_seq'].get(protein_seq)
                if msa_file:
                    print(f"MSA matched by sequence for protein '{protein_id}'")
        elif isinstance(msa_mappings, dict):
            msa_file = msa_mappings.get(protein_id)
        return msa_file

    # Case 1: Multiple proteins WITH complex_id - group by complex_id (from SplitChains)
    if has_complex_id and num_proteins > 1:
        protein_groups = proteins_df.groupby('complex_id')

        for complex_id, group in protein_groups:
            group_proteins = group.to_dict('records')

            ligands_to_process = ligands if ligands else [None]
            for ligand in ligands_to_process:
                config = {'sequences': []}

                # Add all proteins in this complex
                for idx, protein in enumerate(group_proteins):
                    protein_id = protein['id']
                    msa_file = get_msa_file(protein_id, protein['sequence'])
                    chain_id = chr(65 + idx)  # A, B, C, ...

                    protein_entry = {
                        'protein': {
                            'id': chain_id,
                            'sequence': protein['sequence']
                        }
                    }

                    if msa_file and os.path.exists(msa_file):
                        protein_entry['protein']['msa'] = msa_file
                    elif use_global_cache and args.msa_folder and os.path.exists(args.msa_folder):
                        msa_pattern = os.path.join(args.msa_folder, "*.csv")
                        msa_files = glob.glob(msa_pattern)
                        if msa_files:
                            protein_entry['protein']['msa'] = msa_files[0]

                    config['sequences'].append(protein_entry)

                # Add template if specified
                if args.template:
                    template_entry = {
                        'pdb': args.template,
                        'force': args.template_force,
                        'threshold': args.template_threshold
                    }
                    if args.template_chains:
                        chain_ids = [c.strip() for c in args.template_chains.split(',')]
                        template_entry['chain_id'] = chain_ids
                    config['templates'] = [template_entry]

                # Add glycosylation BEFORE the ligand
                num_chains = len(group_proteins)
                config = add_glycosylation_to_config(config, glycosylation, num_chains)

                # Add ligand
                ligand_chain_id = None
                if ligand:
                    ligand_chain_id = chr(65 + len(config['sequences']))
                    ligand_entry = {'ligand': {'id': ligand_chain_id}}

                    if ligand['format'].lower() == 'smiles' and 'smiles' in ligand and ligand['smiles']:
                        ligand_entry['ligand']['smiles'] = ligand['smiles']
                    elif ligand['format'].lower() == 'ccd' and 'ccd' in ligand and ligand['ccd']:
                        ligand_entry['ligand']['ccd'] = ligand['ccd']

                    config['sequences'].append(ligand_entry)

                # Add pocket constraints
                if args.pocket_residues and ligand:
                    import ast
                    pocket_residues = ast.literal_eval(args.pocket_residues) if isinstance(args.pocket_residues, str) else args.pocket_residues
                    constraint_entry = {
                        'pocket': {
                            'binder': ligand_chain_id,
                            'contacts': pocket_residues,
                            'max_distance': args.pocket_max_distance,
                            'force': args.pocket_force
                        }
                    }
                    if 'constraints' not in config:
                        config['constraints'] = []
                    config['constraints'].append(constraint_entry)

                # Add covalent linkage
                if ligand:
                    config = add_covalent_linkage_to_config(config, covalent_linkage, ligand_chain_id)

                # Add affinity
                if args.affinity and ligand:
                    config['properties'] = [{'affinity': {'binder': ligand_chain_id}}]

                # Write config file
                if ligand and num_ligands > 1:
                    config_filename = f"{complex_id}_{ligand['id']}.yaml"
                else:
                    config_filename = f"{complex_id}.yaml"
                config_path = os.path.join(args.output_dir, config_filename)

                with open(config_path, 'w') as f:
                    yaml.dump(config, f, default_flow_style=False, sort_keys=False)

                config_count += 1
                print(f"Created config: {config_filename}")

    # Case 2: Multiple proteins WITHOUT complex_id + ligands - one config per ligand with ALL proteins
    elif num_proteins > 1 and num_ligands >= 1:
        for ligand in ligands:
            config = {'sequences': []}

            for idx, protein in enumerate(proteins):
                protein_id = protein['id']
                msa_file = get_msa_file(protein_id, protein['sequence'])
                chain_id = chr(65 + idx)

                protein_entry = {
                    'protein': {
                        'id': chain_id,
                        'sequence': protein['sequence']
                    }
                }

                if msa_file and os.path.exists(msa_file):
                    protein_entry['protein']['msa'] = msa_file
                elif use_global_cache and args.msa_folder and os.path.exists(args.msa_folder):
                    msa_pattern = os.path.join(args.msa_folder, "*.csv")
                    msa_files = glob.glob(msa_pattern)
                    if msa_files:
                        protein_entry['protein']['msa'] = msa_files[0]

                config['sequences'].append(protein_entry)

            # Add template
            if args.template:
                template_entry = {
                    'pdb': args.template,
                    'force': args.template_force,
                    'threshold': args.template_threshold
                }
                if args.template_chains:
                    chain_ids = [c.strip() for c in args.template_chains.split(',')]
                    template_entry['chain_id'] = chain_ids
                config['templates'] = [template_entry]

            # Add glycosylation
            num_chains = num_proteins
            config = add_glycosylation_to_config(config, glycosylation, num_chains)

            # Add ligand
            ligand_chain_id = chr(65 + len(config['sequences']))
            ligand_entry = {'ligand': {'id': ligand_chain_id}}

            if ligand['format'].lower() == 'smiles' and 'smiles' in ligand and ligand['smiles']:
                ligand_entry['ligand']['smiles'] = ligand['smiles']
            elif ligand['format'].lower() == 'ccd' and 'ccd' in ligand and ligand['ccd']:
                ligand_entry['ligand']['ccd'] = ligand['ccd']

            config['sequences'].append(ligand_entry)

            # Add pocket constraints
            if args.pocket_residues:
                import ast
                pocket_residues_parsed = ast.literal_eval(args.pocket_residues) if isinstance(args.pocket_residues, str) else args.pocket_residues
                constraint_entry = {
                    'pocket': {
                        'binder': ligand_chain_id,
                        'contacts': pocket_residues_parsed,
                        'max_distance': args.pocket_max_distance,
                        'force': args.pocket_force
                    }
                }
                if 'constraints' not in config:
                    config['constraints'] = []
                config['constraints'].append(constraint_entry)

            # Add covalent linkage
            config = add_covalent_linkage_to_config(config, covalent_linkage, ligand_chain_id)

            # Add affinity
            if args.affinity:
                config['properties'] = [{'affinity': {'binder': ligand_chain_id}}]

            config_filename = f"{ligand['id']}.yaml"
            config_path = os.path.join(args.output_dir, config_filename)

            with open(config_path, 'w') as f:
                yaml.dump(config, f, default_flow_style=False, sort_keys=False)

            config_count += 1
            print(f"Created config: {config_filename}")

    # Case 3: Multiple proteins without ligands - one config with all proteins
    elif not ligands and num_proteins > 1:
        config = {'sequences': []}

        for idx, protein in enumerate(proteins):
            protein_id = protein['id']
            msa_file = get_msa_file(protein_id, protein['sequence'])
            chain_id = chr(65 + idx)

            protein_entry = {
                'protein': {
                    'id': chain_id,
                    'sequence': protein['sequence']
                }
            }

            if msa_file and os.path.exists(msa_file):
                protein_entry['protein']['msa'] = msa_file
            elif use_global_cache and args.msa_folder and os.path.exists(args.msa_folder):
                msa_pattern = os.path.join(args.msa_folder, "*.csv")
                msa_files = glob.glob(msa_pattern)
                if msa_files:
                    protein_entry['protein']['msa'] = msa_files[0]

            config['sequences'].append(protein_entry)

        # Add template
        if args.template:
            template_entry = {
                'pdb': args.template,
                'force': args.template_force,
                'threshold': args.template_threshold
            }
            if args.template_chains:
                chain_ids = [c.strip() for c in args.template_chains.split(',')]
                template_entry['chain_id'] = chain_ids
            config['templates'] = [template_entry]

        config_filename = args.job_name + ".yaml" if args.job_name else "complex.yaml"
        config_path = os.path.join(args.output_dir, config_filename)

        with open(config_path, 'w') as f:
            yaml.dump(config, f, default_flow_style=False, sort_keys=False)

        config_count = 1
        print(f"Created multi-protein config: {config_filename}")

    # Case 4: Regular single protein case
    else:
        for protein in proteins:
            protein_id = protein['id']
            msa_file = get_msa_file(protein_id, protein['sequence'])

            if ligands:
                for ligand in ligands:
                    # Create config WITHOUT ligand first
                    config = create_config(protein, None, msa_file, args.affinity, use_global_cache, args.msa_folder,
                                         template=args.template, template_chains=args.template_chains,
                                         template_force=args.template_force, template_threshold=args.template_threshold,
                                         pocket_residues=args.pocket_residues, pocket_max_distance=args.pocket_max_distance,
                                         pocket_force=args.pocket_force)

                    # Add glycosylation BEFORE the ligand
                    num_chains = 1
                    config = add_glycosylation_to_config(config, glycosylation, num_chains)

                    # Add ligand with correct chain ID
                    ligand_chain_id = chr(65 + len(config['sequences']))
                    ligand_entry = {'ligand': {'id': ligand_chain_id}}

                    if ligand['format'].lower() == 'smiles' and 'smiles' in ligand and ligand['smiles']:
                        ligand_entry['ligand']['smiles'] = ligand['smiles']
                    elif ligand['format'].lower() == 'ccd' and 'ccd' in ligand and ligand['ccd']:
                        ligand_entry['ligand']['ccd'] = ligand['ccd']

                    config['sequences'].append(ligand_entry)

                    # Update affinity binder
                    if args.affinity:
                        config['properties'] = [{'affinity': {'binder': ligand_chain_id}}]

                    # Add pocket constraints with correct ligand chain ID
                    if args.pocket_residues:
                        import ast
                        pocket_residues_parsed = ast.literal_eval(args.pocket_residues) if isinstance(args.pocket_residues, str) else args.pocket_residues
                        constraint_entry = {
                            'pocket': {
                                'binder': ligand_chain_id,
                                'contacts': pocket_residues_parsed,
                                'max_distance': args.pocket_max_distance,
                                'force': args.pocket_force
                            }
                        }
                        if 'constraints' not in config:
                            config['constraints'] = []
                        config['constraints'].append(constraint_entry)

                    # Add covalent linkage
                    config = add_covalent_linkage_to_config(config, covalent_linkage, ligand_chain_id)

                    # Determine config filename
                    if num_proteins == 1 and num_ligands == 1:
                        if args.job_name:
                            config_filename = f"{args.job_name}.yaml"
                        else:
                            config_filename = f"{ligand['id']}.yaml"
                    elif num_proteins == 1 and num_ligands > 1:
                        config_filename = f"{ligand['id']}.yaml"
                    elif num_proteins > 1 and num_ligands == 1:
                        config_filename = f"{protein_id}.yaml"
                    else:
                        config_filename = f"{protein_id}_{ligand['id']}.yaml"

                    config_path = os.path.join(args.output_dir, config_filename)

                    with open(config_path, 'w') as f:
                        yaml.dump(config, f, default_flow_style=False, sort_keys=False)

                    config_count += 1
                    print(f"Created config: {config_filename}")
            else:
                # Protein-only config
                config = create_config(protein, None, msa_file, args.affinity, use_global_cache, args.msa_folder,
                                     template=args.template, template_chains=args.template_chains,
                                     template_force=args.template_force, template_threshold=args.template_threshold,
                                     pocket_residues=args.pocket_residues, pocket_max_distance=args.pocket_max_distance,
                                     pocket_force=args.pocket_force)

                # Add glycosylation
                num_chains = 1
                config = add_glycosylation_to_config(config, glycosylation, num_chains)

                config_filename = f"{protein_id}.yaml"
                config_path = os.path.join(args.output_dir, config_filename)

                with open(config_path, 'w') as f:
                    yaml.dump(config, f, default_flow_style=False, sort_keys=False)

                config_count += 1
                print(f"Created config: {config_filename}")

    print(f"\nGenerated {config_count} config files in {args.output_dir}")

    # Generate sequence_ids.csv for post-processing
    sequence_ids_file = os.path.join(os.path.dirname(args.output_dir), "sequence_ids.csv")

    sequence_ids = []
    if has_complex_id and num_proteins > 1:
        sequence_ids = proteins_df['complex_id'].unique().tolist()
    elif num_proteins > 1 and ligands:
        sequence_ids = [ligand['id'] for ligand in ligands]
    elif num_proteins == 1 and len(ligands) > 1:
        sequence_ids = [ligand['id'] for ligand in ligands]
    elif num_proteins > 1 and not ligands:
        sequence_ids = [protein['id'] for protein in proteins]
    elif num_proteins == 1 and len(ligands) == 1:
        if args.job_name:
            sequence_ids = [args.job_name]
        else:
            sequence_ids = [ligands[0]['id']]
    elif num_proteins == 1 and not ligands:
        sequence_ids = [proteins[0]['id']]

    if sequence_ids:
        import csv
        with open(sequence_ids_file, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['id'])
            for seq_id in sequence_ids:
                writer.writerow([seq_id])
        print(f"Created sequence IDs file: {sequence_ids_file}")
        print(f"Expected output structure IDs: {', '.join(sequence_ids)}")


if __name__ == "__main__":
    main()
