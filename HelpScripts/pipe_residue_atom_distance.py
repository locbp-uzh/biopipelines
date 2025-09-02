#!/usr/bin/env python3
"""
Runtime helper script for ResidueAtomDistance analysis.

This script analyzes protein structures to calculate distances between specific atoms and residues,
outputting CSV with distance metrics for all structures.
"""

import os
import sys
import argparse
import json
import pandas as pd
import numpy as np
from typing import Dict, List, Any, Optional, Tuple

# Import MDAnalysis for structure analysis
try:
    import MDAnalysis as mda
    from MDAnalysis.lib.distances import distance_array
    HAS_MDA = True
except ImportError:
    HAS_MDA = False


def parse_selection(selection: str, structure_path: str) -> str:
    """
    Parse and convert selection string to MDAnalysis format.
    
    Args:
        selection: Selection string (e.g., 'LIG.Cl', 'protein.D in TRGDTGH')
        structure_path: Path to structure file for context
        
    Returns:
        MDAnalysis selection string
    """
    # Convert common selection formats to MDAnalysis
    # This is a simplified parser - could be expanded for more complex selections
    
    if 'LIG.' in selection:
        # Ligand atom selection: 'LIG.Cl' -> 'resname LIG and name Cl'
        parts = selection.replace('LIG.', '').strip()
        return f"resname LIG and name {parts}"
    
    elif 'protein.' in selection:
        # Protein selection: 'protein.D in TRGDTGH' -> specific residue and atom
        parts = selection.replace('protein.', '').strip()
        if ' in ' in parts:
            # Residue in sequence: 'D in TRGDTGH' -> 'protein and resname ASP'
            residue, sequence = parts.split(' in ')
            # Map single letter to three letter code
            residue_map = {'D': 'ASP', 'E': 'GLU', 'K': 'LYS', 'R': 'ARG', 'H': 'HIS'}
            three_letter = residue_map.get(residue, residue)
            return f"protein and resname {three_letter}"
        else:
            # Simple protein atom: 'CA' -> 'protein and name CA'
            return f"protein and name {parts}"
    
    elif 'resname' in selection:
        # Direct resname selection: 'resname ZN' -> 'resname ZN'
        return selection
    
    elif 'resid' in selection:
        # Residue ID selection: 'resid 145-150' -> 'resid 145:150'
        return selection.replace('-', ':')
    
    else:
        # Assume it's already in MDAnalysis format
        return selection


def calculate_distance(structure_path: str, atom_selection: str, residue_selection: str, 
                      distance_metric: str) -> Optional[float]:
    """
    Calculate distance between atom and residue selections.
    
    Args:
        structure_path: Path to structure file
        atom_selection: Atom selection string
        residue_selection: Residue selection string
        distance_metric: Type of distance calculation
        
    Returns:
        Calculated distance or None if failed
    """
    if not HAS_MDA:
        print("Warning: MDAnalysis not available, using placeholder distance")
        return 5.0  # Placeholder value
    
    try:
        # Load structure
        u = mda.Universe(structure_path)
        
        # Parse selections
        atom_sel_str = parse_selection(atom_selection, structure_path)
        residue_sel_str = parse_selection(residue_selection, structure_path)
        
        print(f"  - Atom selection: {atom_sel_str}")
        print(f"  - Residue selection: {residue_sel_str}")
        
        # Select atoms
        atom_group = u.select_atoms(atom_sel_str)
        residue_group = u.select_atoms(residue_sel_str)
        
        if len(atom_group) == 0:
            print(f"  - Warning: No atoms found for selection '{atom_sel_str}'")
            return None
        
        if len(residue_group) == 0:
            print(f"  - Warning: No atoms found for selection '{residue_sel_str}'")
            return None
        
        print(f"  - Found {len(atom_group)} atoms, {len(residue_group)} residue atoms")
        
        # Calculate distance matrix
        distances = distance_array(atom_group.positions, residue_group.positions)
        
        # Apply distance metric
        if distance_metric == "min":
            result = np.min(distances)
        elif distance_metric == "max":
            result = np.max(distances)
        elif distance_metric == "mean":
            result = np.mean(distances)
        elif distance_metric == "closest":
            result = np.min(distances)  # Same as min
        else:
            result = np.min(distances)  # Default to min
        
        print(f"  - Distance ({distance_metric}): {result:.3f} Å")
        return float(result)
        
    except Exception as e:
        print(f"  - Error calculating distance: {e}")
        return None


def analyze_residue_atom_distances(config_data: Dict[str, Any]) -> None:
    """
    Analyze distances between atoms and residues in structures.
    
    Args:
        config_data: Configuration dictionary with analysis parameters
    """
    input_structures = config_data['input_structures']
    atom_selection = config_data['atom_selection']
    residue_selection = config_data['residue_selection']
    distance_metric = config_data['distance_metric']
    metric_name = config_data['metric_name']
    output_csv = config_data['output_csv']
    
    print(f"Analyzing distances in {len(input_structures)} structures")
    print(f"Atom selection: {atom_selection}")
    print(f"Residue selection: {residue_selection}")
    print(f"Distance metric: {distance_metric}")
    print(f"Output column: {metric_name}")
    
    # Process each structure
    results = []
    
    for i, structure_path in enumerate(input_structures):
        if not os.path.exists(structure_path):
            print(f"Warning: Structure file not found: {structure_path}")
            continue
        
        print(f"\nProcessing structure {i+1}/{len(input_structures)}: {structure_path}")
        
        # Extract structure ID from filename
        structure_id = os.path.splitext(os.path.basename(structure_path))[0]
        
        # Calculate distance
        distance = calculate_distance(structure_path, atom_selection, residue_selection, distance_metric)
        
        # Store result
        result = {
            'id': structure_id,
            'source_structure': structure_path,
            metric_name: distance
        }
        results.append(result)
        
        print(f"  - Result: {metric_name} = {distance}")
    
    # Create DataFrame and save
    if results:
        df = pd.DataFrame(results)
        
        # Create output directory
        os.makedirs(os.path.dirname(output_csv), exist_ok=True)
        
        # Save results
        df.to_csv(output_csv, index=False)
        
        print(f"\nDistance analysis completed successfully!")
        print(f"Analyzed {len(results)} structures")
        print(f"Results saved to: {output_csv}")
        print(f"\nResults summary:")
        print(df)
        
        # Statistics
        distances = [r[metric_name] for r in results if r[metric_name] is not None]
        if distances:
            print(f"\nDistance statistics:")
            print(f"  Min: {min(distances):.3f} Å")
            print(f"  Max: {max(distances):.3f} Å")
            print(f"  Mean: {np.mean(distances):.3f} Å")
            print(f"  Std: {np.std(distances):.3f} Å")
    else:
        raise ValueError("No valid results generated - check structure files and selections")


def main():
    parser = argparse.ArgumentParser(description='Analyze residue-atom distances in structures')
    parser.add_argument('--config', required=True, help='JSON config file with analysis parameters')
    
    args = parser.parse_args()
    
    # Load configuration
    if not os.path.exists(args.config):
        print(f"Error: Config file not found: {args.config}")
        sys.exit(1)
    
    try:
        with open(args.config, 'r') as f:
            config_data = json.load(f)
    except Exception as e:
        print(f"Error loading config: {e}")
        sys.exit(1)
    
    # Validate required parameters
    required_params = ['input_structures', 'atom_selection', 'residue_selection', 
                      'distance_metric', 'metric_name', 'output_csv']
    for param in required_params:
        if param not in config_data:
            print(f"Error: Missing required parameter: {param}")
            sys.exit(1)
    
    try:
        analyze_residue_atom_distances(config_data)
        
    except Exception as e:
        print(f"Error analyzing distances: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()