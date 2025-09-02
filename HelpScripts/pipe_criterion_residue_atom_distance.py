#!/usr/bin/env python3
"""
Runtime helper script for ResidueAtomDistance analysis.

This script analyzes protein structures to calculate distances between atoms and residues
and outputs results as CSV for downstream filtering and analysis.
"""

import os
import sys
import argparse
import json
import numpy as np
import pandas as pd
import glob
from typing import Dict, List, Any, Optional, Tuple

# Try to import structural biology libraries
try:
    import pymol
    from pymol import cmd
    PYMOL_AVAILABLE = True
except ImportError:
    PYMOL_AVAILABLE = False

try:
    import MDAnalysis as mda
    MDANALYSIS_AVAILABLE = True
except ImportError:
    MDANALYSIS_AVAILABLE = False


def parse_selection_string(selection: str) -> Dict[str, Any]:
    """
    Parse selection string into structured criteria.
    """
    if '.' in selection and ' in ' in selection:
        # Format: 'protein.D in TRGDTGH'
        parts = selection.split(' in ')
        if len(parts) == 2:
            entity_res, sequence = parts
            if '.' in entity_res:
                entity, residue_type = entity_res.split('.', 1)
                return {
                    "entity": entity.lower(),
                    "residue_type": residue_type,
                    "sequence_context": sequence,
                    "selection_type": "sequence_context"
                }
    elif '.' in selection:
        # Format: 'ligand.Cl' or 'protein.CA'
        parts = selection.split('.', 1)
        if len(parts) == 2:
            entity, atom_type = parts
            return {
                "entity": entity.lower(),
                "atom_type": atom_type,
                "selection_type": "entity_atom"
            }
    elif selection.lower().startswith('resname') or selection.lower().startswith('resid'):
        # MDAnalysis-style selection
        return {
            "raw_selection": selection,
            "selection_type": "mdanalysis"
        }
    
    # Default: treat as raw selection string
    return {
        "raw_selection": selection,
        "selection_type": "raw"
    }


def calculate_distance(structure_file: str, config: Dict[str, Any]) -> Optional[float]:
    """
    Calculate distance between atom and residue selections.
    """
    atom_selection = config.get('atom', '')
    residue_selection = config.get('residue', '')
    distance_metric = config.get('distance_metric', 'min')
    
    try:
        if MDANALYSIS_AVAILABLE:
            return calculate_distance_mdanalysis(structure_file, atom_selection, residue_selection, distance_metric)
        elif PYMOL_AVAILABLE:
            return calculate_distance_pymol(structure_file, atom_selection, residue_selection, distance_metric)
        else:
            print("Warning: Neither MDAnalysis nor PyMOL available - cannot calculate distances")
            return None
    
    except Exception as e:
        print(f"Error calculating distance for {structure_file}: {e}")
        return None


def calculate_distance_mdanalysis(structure_file: str, atom_sel: str, residue_sel: str, metric: str) -> Optional[float]:
    """Calculate distance using MDAnalysis."""
    u = mda.Universe(structure_file)
    
    # Parse selections
    atom_parsed = parse_selection_string(atom_sel)
    residue_parsed = parse_selection_string(residue_sel)
    
    # Build MDAnalysis selection strings
    atom_selection_str = build_mdanalysis_selection(atom_parsed)
    residue_selection_str = build_mdanalysis_selection(residue_parsed)
    
    if not atom_selection_str or not residue_selection_str:
        print(f"Could not build valid selections: atom='{atom_sel}', residue='{residue_sel}'")
        return None
    
    # Select atoms
    try:
        atom_group = u.select_atoms(atom_selection_str)
        residue_group = u.select_atoms(residue_selection_str)
    except Exception as e:
        print(f"Selection error: {e}")
        return None
    
    if len(atom_group) == 0 or len(residue_group) == 0:
        print(f"Empty selection: atoms={len(atom_group)}, residues={len(residue_group)}")
        return None
    
    # Calculate distances
    distances = []
    for atom in atom_group:
        for residue_atom in residue_group:
            dist = np.linalg.norm(atom.position - residue_atom.position)
            distances.append(dist)
    
    if not distances:
        return None
    
    # Apply distance metric
    if metric == 'min':
        return np.min(distances)
    elif metric == 'max':
        return np.max(distances)
    elif metric == 'mean':
        return np.mean(distances)
    elif metric == 'closest':
        return np.min(distances)  # Same as min
    else:
        return np.min(distances)


def build_mdanalysis_selection(parsed_sel: Dict[str, Any]) -> Optional[str]:
    """Build MDAnalysis selection string from parsed selection."""
    sel_type = parsed_sel.get('selection_type', 'raw')
    
    if sel_type == 'mdanalysis' or sel_type == 'raw':
        return parsed_sel.get('raw_selection', '')
    
    elif sel_type == 'entity_atom':
        entity = parsed_sel['entity']
        atom_type = parsed_sel['atom_type']
        
        if entity == 'ligand':
            # Assume ligand is HETATM
            if atom_type.upper() == 'ALL':
                return "not protein"
            else:
                return f"not protein and name {atom_type.upper()}"
        
        elif entity == 'protein':
            if atom_type.upper() == 'ALL':
                return "protein"
            else:
                return f"protein and name {atom_type.upper()}"
        
        else:
            # Try as resname
            if atom_type.upper() == 'ALL':
                return f"resname {entity.upper()}"
            else:
                return f"resname {entity.upper()} and name {atom_type.upper()}"
    
    elif sel_type == 'sequence_context':
        entity = parsed_sel['entity']
        residue_type = parsed_sel['residue_type']
        
        if entity == 'protein':
            return f"protein and resname {residue_type.upper()}"
        else:
            return f"resname {residue_type.upper()}"
    
    return None


def calculate_distance_pymol(structure_file: str, atom_sel: str, residue_sel: str, metric: str) -> Optional[float]:
    """Calculate distance using PyMOL (placeholder - would need full implementation)."""
    # This would require a full PyMOL implementation
    print("PyMOL distance calculation not yet implemented")
    return None


def main():
    parser = argparse.ArgumentParser(description='Generate distance analysis CSV')
    parser.add_argument('--input_dir', required=True, help='Directory containing structure files')
    parser.add_argument('--output_csv', required=True, help='Output CSV file for analysis results')
    parser.add_argument('--config', required=True, help='JSON config file with analysis parameters')
    
    args = parser.parse_args()
    
    # Load configuration
    if not os.path.exists(args.config):
        print(f"Error: Config file not found: {args.config}")
        sys.exit(1)
    
    try:
        with open(args.config, 'r') as f:
            config = json.load(f)
    except Exception as e:
        print(f"Error loading config: {e}")
        sys.exit(1)
    
    metric_name = config.get('metric_name', 'distance')
    print(f"Generating distance analysis")
    print(f"Metric name: {metric_name}")
    print(f"Atom selection: {config.get('atom', 'N/A')}")
    print(f"Residue selection: {config.get('residue', 'N/A')}")
    print(f"Distance metric: {config.get('distance_metric', 'min')}")
    
    try:
        # Find structure files
        structure_files = []
        for ext in ['.pdb', '.cif', '.mmcif']:
            pattern = os.path.join(args.input_dir, f"*{ext}")
            structure_files.extend(glob.glob(pattern))
        
        if not structure_files:
            print(f"No structure files found in {args.input_dir}")
            sys.exit(1)
        
        print(f"Found {len(structure_files)} structure files")
        
        # Initialize results list for CSV
        csv_results = []
        
        # Process each structure file
        for structure_file in structure_files:
            structure_name = os.path.basename(structure_file)
            structure_id = os.path.splitext(structure_name)[0]
            print(f"Processing: {structure_name}")
            
            # Calculate distance
            distance = calculate_distance(structure_file, config)
            
            if distance is None:
                print(f"  Could not calculate distance - using 999.0")
                distance = 999.0  # Large value indicating no measurement
            else:
                print(f"  {metric_name}={distance:.2f} Å")
            
            # Create CSV row
            csv_row = {
                'id': structure_id,
                'source_structure': structure_name,
                metric_name: distance,
                'atom_selection': config.get('atom', ''),
                'residue_selection': config.get('residue', ''),
                'distance_metric': config.get('distance_metric', 'min')
            }
            
            csv_results.append(csv_row)
        
        # Create output directory if needed
        os.makedirs(os.path.dirname(args.output_csv), exist_ok=True)
        
        # Save results as CSV
        df = pd.DataFrame(csv_results)
        df.to_csv(args.output_csv, index=False)
        
        print(f"\nAnalysis CSV saved: {args.output_csv}")
        print(f"Analyzed {len(csv_results)} structures")
        print(f"Columns: {list(df.columns)}")
        
        # Show summary statistics
        if metric_name in df.columns:
            # Exclude invalid measurements (999.0) from stats
            valid_distances = df[df[metric_name] < 999.0][metric_name]
            if len(valid_distances) > 0:
                distance_stats = valid_distances.describe()
                print(f"\n{metric_name} statistics (valid measurements):")
                print(f"  Count: {len(valid_distances)}/{len(df)}")
                print(f"  Mean: {distance_stats['mean']:.2f} Å")
                print(f"  Min: {distance_stats['min']:.2f} Å")
                print(f"  Max: {distance_stats['max']:.2f} Å")
            else:
                print(f"\nNo valid distance measurements found")
        
    except Exception as e:
        print(f"Error during analysis: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()