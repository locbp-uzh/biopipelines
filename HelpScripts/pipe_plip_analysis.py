#!/usr/bin/env python3
# Copyright (c) 2026 Gianluca Quargnali @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Runtime helper script for PLIP analysis tool.

This script processes PLIP raw outputs (XML, text files, PyMOL sessions) into
standardized CSV format for downstream analysis.
"""

import os
import sys
import argparse
import json
import pandas as pd
import xml.etree.ElementTree as ET
from typing import Dict, List, Any, Optional
from pathlib import Path

# Import unified I/O utilities
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from biopipelines_io import load_datastream, iterate_files


def parse_plip_xml(xml_file: str) -> List[Dict[str, Any]]:
    """
    Parse PLIP XML output file to extract interaction data.

    Args:
        xml_file: Path to PLIP XML output file

    Returns:
        List of interaction dictionaries
    """
    interactions = []

    try:
        tree = ET.parse(xml_file)
        root = tree.getroot()

        # Extract basic information
        pdb_id = root.get('pdbid', 'unknown')

        # Parse binding sites
        binding_sites = root.findall('.//bindingsite')

        for site in binding_sites:
            ligand_id = site.get('id', 'unknown')
            ligand_name = site.get('ligand_longname', 'unknown')

            # Parse different interaction types
            for interaction_type in ['hydrophobic_interaction', 'hydrogen_bond',
                                   'water_bridge', 'salt_bridge', 'pi_stacking',
                                   'pi_cation_interaction', 'halogen_bond', 'metal_complexation']:

                interactions_of_type = site.findall(f'.//{interaction_type}')

                for interaction in interactions_of_type:
                    interaction_data = {
                        'structure_id': pdb_id,
                        'ligand_id': ligand_id,
                        'ligand_name': ligand_name,
                        'interaction_type': interaction_type,
                        'residue_name': interaction.get('resname', ''),
                        'residue_number': interaction.get('resnr', ''),
                        'chain': interaction.get('reschain', ''),
                        'distance': float(interaction.get('dist', 0.0)) if interaction.get('dist') else None,
                        'angle': float(interaction.get('angle', 0.0)) if interaction.get('angle') else None,
                        'energy': float(interaction.get('energy', 0.0)) if interaction.get('energy') else None,
                        'ligand_atom': interaction.get('ligcoo', ''),
                        'protein_atom': interaction.get('protcoo', ''),
                        'ligand_idx': interaction.get('ligcarbonidx', ''),
                        'protein_idx': interaction.get('protcarbonidx', '')
                    }
                    interactions.append(interaction_data)

    except Exception as e:
        print(f"Error parsing XML file {xml_file}: {e}")

    return interactions


def parse_plip_txt(txt_file: str) -> Dict[str, Any]:
    """
    Parse PLIP text output file to extract summary information.

    Args:
        txt_file: Path to PLIP text output file

    Returns:
        Dictionary with summary information
    """
    summary = {}

    try:
        with open(txt_file, 'r') as f:
            content = f.read()

        # Extract basic statistics
        lines = content.split('\n')
        for line in lines:
            line = line.strip()
            if 'binding sites detected' in line.lower():
                summary['num_binding_sites'] = int(line.split()[0])
            elif 'interactions detected' in line.lower():
                summary['total_interactions'] = int(line.split()[0])

    except Exception as e:
        print(f"Error parsing text file {txt_file}: {e}")

    return summary


def process_plip_structure(structure_id: str, structure_path: str, raw_dir: str, ligand_filter: str = "") -> List[Dict[str, Any]]:
    """
    Process PLIP outputs for a single structure.

    Args:
        structure_id: ID of the structure from DataStream
        structure_path: Path to input structure file
        raw_dir: Directory containing PLIP raw outputs
        ligand_filter: Specific ligand ID to filter (empty = all ligands)

    Returns:
        List of interaction data
    """
    structure_output_dir = os.path.join(raw_dir, structure_id)

    interactions = []

    if not os.path.exists(structure_output_dir):
        print(f"Warning: No PLIP output directory found for {structure_id}")
        return interactions

    # Look for XML files (main source of interaction data)
    xml_files = []
    for file in os.listdir(structure_output_dir):
        if file.endswith('.xml'):
            xml_files.append(os.path.join(structure_output_dir, file))

    if not xml_files:
        print(f"Warning: No XML output files found for {structure_id}")
        return interactions

    # Process each XML file
    for xml_file in xml_files:
        xml_interactions = parse_plip_xml(xml_file)

        # Filter by ligand if specified
        if ligand_filter:
            xml_interactions = [i for i in xml_interactions if i['ligand_id'] == ligand_filter]

        interactions.extend(xml_interactions)

    # Add structure path information
    for interaction in interactions:
        interaction['structure_file'] = structure_path

    return interactions


def create_summary(all_interactions: List[Dict[str, Any]]) -> Dict[str, Any]:
    """
    Create summary statistics from all interactions.

    Args:
        all_interactions: List of all interaction data

    Returns:
        Summary statistics dictionary
    """
    if not all_interactions:
        return {"total_interactions": 0, "structures_analyzed": 0}

    df = pd.DataFrame(all_interactions)

    summary = {
        "total_interactions": len(all_interactions),
        "structures_analyzed": df['structure_id'].nunique(),
        "unique_ligands": df['ligand_id'].nunique(),
        "interaction_types": df['interaction_type'].value_counts().to_dict(),
        "interactions_per_structure": df.groupby('structure_id').size().to_dict(),
        "interactions_per_ligand": df.groupby('ligand_id').size().to_dict()
    }

    return summary


def create_summary_csv(all_interactions: List[Dict[str, Any]], structure_items: List[tuple], output_csv: str):
    """
    Create a summary CSV with aggregated interaction counts per structure.

    Args:
        all_interactions: List of all interaction data
        structure_items: List of (structure_id, structure_path) tuples from DataStream
        output_csv: Output CSV file path
    """
    # Map interaction types to column names
    type_mapping = {
        'hydrogen_bond': 'hbonds',
        'salt_bridge': 'saltbridges',
        'hydrophobic_interaction': 'hydrophobic',
        'pi_stacking': 'pistacking',
        'pi_cation_interaction': 'pication',
        'halogen_bond': 'halogen',
        'metal_complexation': 'metal',
        'water_bridge': 'waterbridges'
    }

    # Initialize results for all structures
    results = []
    for structure_id, structure_path in structure_items:
        results.append({
            'id': structure_id,
            'structure': structure_path,
            'hbonds': 0,
            'saltbridges': 0,
            'hydrophobic': 0,
            'pistacking': 0,
            'pication': 0,
            'halogen': 0,
            'metal': 0,
            'waterbridges': 0,
            'total_interactions': 0
        })

    # Create lookup dict for faster access
    results_dict = {r['id']: r for r in results}

    # Aggregate interaction counts
    for interaction in all_interactions:
        structure_id = interaction.get('structure_id', '')
        interaction_type = interaction.get('interaction_type', '')

        if structure_id in results_dict:
            column_name = type_mapping.get(interaction_type)
            if column_name and column_name in results_dict[structure_id]:
                results_dict[structure_id][column_name] += 1
            results_dict[structure_id]['total_interactions'] += 1

    # Create DataFrame and save
    df = pd.DataFrame(results)
    os.makedirs(os.path.dirname(output_csv), exist_ok=True)
    df.to_csv(output_csv, index=False)
    print(f"Saved summary CSV with {len(results)} structures to {output_csv}")


def copy_additional_outputs(raw_dir: str, processed_dir: str, structure_items: List[tuple]):
    """
    Copy additional PLIP outputs (PyMOL files, images) to processed directory.

    Args:
        raw_dir: Directory containing raw PLIP outputs
        processed_dir: Directory for processed outputs
        structure_items: List of (structure_id, structure_path) tuples from DataStream
    """
    os.makedirs(processed_dir, exist_ok=True)

    for structure_id, structure_path in structure_items:
        structure_output_dir = os.path.join(raw_dir, structure_id)

        if not os.path.exists(structure_output_dir):
            continue

        # Create structure-specific processed directory
        processed_structure_dir = os.path.join(processed_dir, structure_id)
        os.makedirs(processed_structure_dir, exist_ok=True)

        # Copy PyMOL session files
        for file in os.listdir(structure_output_dir):
            if file.endswith(('.pse', '.png', '.pml')):
                src = os.path.join(structure_output_dir, file)
                dst = os.path.join(processed_structure_dir, file)

                try:
                    import shutil
                    shutil.copy2(src, dst)
                    print(f"Copied {file} to processed outputs")
                except Exception as e:
                    print(f"Warning: Could not copy {file}: {e}")


def main():
    parser = argparse.ArgumentParser(description='Process PLIP outputs into standardized format')
    parser.add_argument('--structures', required=True, help='JSON file containing DataStream with structure files')
    parser.add_argument('--raw_dir', required=True, help='Directory containing PLIP raw outputs')
    parser.add_argument('--output_csv', required=True, help='Output CSV file for interactions')
    parser.add_argument('--summary_csv', required=True, help='Output CSV file for aggregated counts per structure')
    parser.add_argument('--summary_txt', required=True, help='Output text file for summary')
    parser.add_argument('--ligand', default='', help='Specific ligand ID to analyze (empty = all)')
    parser.add_argument('--processed_dir', required=True, help='Directory for processed outputs')

    args = parser.parse_args()

    # Load structures DataStream using pipe_biopipelines_io
    structures_ds = load_datastream(args.structures)

    # Build list of (struct_id, struct_path) tuples
    structure_items = list(iterate_files(structures_ds))

    if not structure_items:
        print(f"Error: No structures found in: {args.structures}")
        sys.exit(1)

    print(f"Processing PLIP outputs for {len(structure_items)} structures")

    # Process each structure
    all_interactions = []

    for struct_id, structure_path in structure_items:
        print(f"Processing structure: {struct_id}")
        interactions = process_plip_structure(struct_id, structure_path, args.raw_dir, args.ligand)
        all_interactions.extend(interactions)
        print(f"  Found {len(interactions)} interactions")

    # Create interactions CSV
    if all_interactions:
        df = pd.DataFrame(all_interactions)

        # Ensure output directory exists
        os.makedirs(os.path.dirname(args.output_csv), exist_ok=True)

        # Save to CSV
        df.to_csv(args.output_csv, index=False)
        print(f"Saved {len(all_interactions)} interactions to {args.output_csv}")
    else:
        print("Warning: No interactions found")
        # Create empty CSV with expected columns
        empty_df = pd.DataFrame(columns=[
            'structure_id', 'ligand_id', 'ligand_name', 'interaction_type',
            'residue_name', 'residue_number', 'chain', 'distance', 'angle', 'energy',
            'ligand_atom', 'protein_atom', 'ligand_idx', 'protein_idx', 'structure_file'
        ])
        os.makedirs(os.path.dirname(args.output_csv), exist_ok=True)
        empty_df.to_csv(args.output_csv, index=False)

    # Create summary CSV with aggregated counts per structure
    create_summary_csv(all_interactions, structure_items, args.summary_csv)

    # Create summary
    summary = create_summary(all_interactions)

    # Write summary text file
    os.makedirs(os.path.dirname(args.summary_txt), exist_ok=True)
    with open(args.summary_txt, 'w') as f:
        f.write("PLIP Analysis Summary\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Total interactions found: {summary.get('total_interactions', 0)}\n")
        f.write(f"Structures analyzed: {summary.get('structures_analyzed', 0)}\n")
        f.write(f"Unique ligands: {summary.get('unique_ligands', 0)}\n\n")

        if 'interaction_types' in summary:
            f.write("Interaction types:\n")
            for itype, count in summary['interaction_types'].items():
                f.write(f"  {itype}: {count}\n")
            f.write("\n")

        if 'interactions_per_structure' in summary:
            f.write("Interactions per structure:\n")
            for structure, count in summary['interactions_per_structure'].items():
                f.write(f"  {structure}: {count}\n")
            f.write("\n")

        if 'interactions_per_ligand' in summary:
            f.write("Interactions per ligand:\n")
            for ligand, count in summary['interactions_per_ligand'].items():
                f.write(f"  {ligand}: {count}\n")

    print(f"Saved summary to {args.summary_txt}")

    # Copy additional outputs (PyMOL files, images)
    copy_additional_outputs(args.raw_dir, args.processed_dir, structure_items)

    print("PLIP output processing completed successfully")


if __name__ == "__main__":
    main()