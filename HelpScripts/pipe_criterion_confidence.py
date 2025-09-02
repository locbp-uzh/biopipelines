#!/usr/bin/env python3
"""
Runtime helper script for Confidence analysis.

This script analyzes protein structures to extract confidence scores (pLDDT)
and outputs results as CSV for downstream filtering and analysis.
"""

import os
import sys
import argparse
import json
import re
import numpy as np
import pandas as pd
import glob
from typing import Dict, List, Any, Optional, Union

# Try to import structural biology libraries
try:
    from Bio.PDB import PDBParser
    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False

try:
    import MDAnalysis as mda
    MDANALYSIS_AVAILABLE = True
except ImportError:
    MDANALYSIS_AVAILABLE = False


def parse_position_string(position_str: str) -> List[int]:
    """
    Parse position string like "45-50,52,55-57" into list of residue numbers.
    """
    positions = []
    
    for part in str(position_str).split(','):
        part = part.strip()
        if '-' in part and not part.startswith('-'):
            try:
                start, end = map(int, part.split('-'))
                positions.extend(range(start, end + 1))
            except ValueError:
                continue
        else:
            try:
                positions.append(int(part))
            except ValueError:
                continue
    
    return positions


def resolve_datasheet_selection(selection_ref: str, context: Dict[str, Any]) -> Optional[List[int]]:
    """
    Resolve datasheet-based selection to list of residue numbers.
    """
    try:
        if not selection_ref.startswith("input.datasheets."):
            return None
        
        parts = selection_ref.split(".")
        if len(parts) < 4:
            return None
        
        datasheet_name = parts[2]
        column_name = parts[3]
        
        datasheets = context.get('datasheets', {})
        
        if isinstance(datasheets, dict) and datasheet_name in datasheets:
            datasheet_info = datasheets[datasheet_name]
            
            if isinstance(datasheet_info, dict) and 'path' in datasheet_info:
                datasheet_path = datasheet_info['path']
            else:
                datasheet_path = str(datasheet_info)
            
            if os.path.exists(datasheet_path):
                df = pd.read_csv(datasheet_path)
                
                if column_name in df.columns:
                    positions = []
                    for value in df[column_name].dropna():
                        if isinstance(value, str):
                            positions.extend(parse_position_string(value))
                        elif isinstance(value, int):
                            positions.append(value)
                    return positions if positions else None
    
    except Exception as e:
        print(f"Warning: Could not resolve datasheet selection '{selection_ref}': {e}")
    
    return None


def calculate_confidence_score(structure_file: str, config: Dict[str, Any], context: Dict[str, Any]) -> Optional[float]:
    """
    Calculate confidence score from structure file.
    """
    selection = config.get('selection')
    score_metric = config.get('score_metric', 'mean')
    score_source = config.get('score_source', 'bfactor')
    
    try:
        if BIOPYTHON_AVAILABLE:
            return calculate_confidence_biopython(structure_file, selection, score_metric, score_source, context)
        elif MDANALYSIS_AVAILABLE:
            return calculate_confidence_mdanalysis(structure_file, selection, score_metric, score_source, context)
        else:
            print("Warning: Neither BioPython nor MDAnalysis available - using fallback parsing")
            return calculate_confidence_fallback(structure_file, score_source)
    
    except Exception as e:
        print(f"Error calculating confidence for {structure_file}: {e}")
        return None


def calculate_confidence_biopython(structure_file: str, selection: Optional[Union[str, List]], 
                                  score_metric: str, score_source: str, context: Dict[str, Any]) -> Optional[float]:
    """Calculate confidence using BioPython."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("structure", structure_file)
    
    confidence_scores = []
    target_residues = None
    
    # Resolve selection
    if isinstance(selection, str) and selection.startswith("input.datasheets."):
        target_residues = resolve_datasheet_selection(selection, context)
    elif isinstance(selection, list):
        if all(isinstance(x, int) for x in selection):
            target_residues = selection
        elif all(isinstance(x, str) for x in selection):
            # Parse residue identifiers like "A45", "B123"
            target_residues = []
            for res_id in selection:
                if len(res_id) > 1 and res_id[1:].isdigit():
                    target_residues.append(int(res_id[1:]))
    
    # Extract confidence scores
    for model in structure:
        for chain in model:
            for residue in chain:
                res_num = residue.id[1]
                
                # Check if this residue should be included
                if target_residues is not None and res_num not in target_residues:
                    continue
                
                # Get confidence score from atoms
                residue_scores = []
                for atom in residue:
                    if score_source == 'bfactor':
                        score = atom.bfactor
                    elif score_source == 'occupancy':
                        score = atom.occupancy
                    else:
                        score = atom.bfactor
                    
                    if score > 0:  # Valid confidence score
                        residue_scores.append(score)
                
                if residue_scores:
                    if score_metric == 'mean':
                        residue_score = np.mean(residue_scores)
                    elif score_metric == 'min':
                        residue_score = np.min(residue_scores)
                    elif score_metric == 'max':
                        residue_score = np.max(residue_scores)
                    elif score_metric == 'median':
                        residue_score = np.median(residue_scores)
                    else:
                        residue_score = np.mean(residue_scores)
                    
                    confidence_scores.append(residue_score)
    
    if not confidence_scores:
        return None
    
    # Calculate final score
    if score_metric == 'mean':
        return np.mean(confidence_scores)
    elif score_metric == 'min':
        return np.min(confidence_scores)
    elif score_metric == 'max':
        return np.max(confidence_scores)
    elif score_metric == 'median':
        return np.median(confidence_scores)
    else:
        return np.mean(confidence_scores)


def calculate_confidence_mdanalysis(structure_file: str, selection: Optional[Union[str, List]], 
                                   score_metric: str, score_source: str, context: Dict[str, Any]) -> Optional[float]:
    """Calculate confidence using MDAnalysis."""
    u = mda.Universe(structure_file)
    
    # Build selection string
    if selection is None:
        selection_str = "all"
    elif isinstance(selection, str) and selection.startswith("input.datasheets."):
        target_residues = resolve_datasheet_selection(selection, context)
        if target_residues:
            res_str = " ".join(str(r) for r in target_residues)
            selection_str = f"resid {res_str}"
        else:
            selection_str = "all"
    elif isinstance(selection, list):
        if all(isinstance(x, int) for x in selection):
            res_str = " ".join(str(r) for r in selection)
            selection_str = f"resid {res_str}"
        else:
            selection_str = "all"
    else:
        selection_str = "all"
    
    atoms = u.select_atoms(selection_str)
    
    if len(atoms) == 0:
        return None
    
    # Extract confidence scores
    if score_source == 'bfactor':
        scores = atoms.tempfactors
    elif score_source == 'occupancy':
        scores = atoms.occupancies
    else:
        scores = atoms.tempfactors
    
    # Filter valid scores
    valid_scores = scores[scores > 0]
    
    if len(valid_scores) == 0:
        return None
    
    # Calculate final score
    if score_metric == 'mean':
        return np.mean(valid_scores)
    elif score_metric == 'min':
        return np.min(valid_scores)
    elif score_metric == 'max':
        return np.max(valid_scores)
    elif score_metric == 'median':
        return np.median(valid_scores)
    else:
        return np.mean(valid_scores)


def calculate_confidence_fallback(structure_file: str, score_source: str) -> Optional[float]:
    """Fallback confidence calculation using simple text parsing."""
    try:
        confidence_scores = []
        
        with open(structure_file, 'r') as f:
            for line in f:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    try:
                        if score_source == 'bfactor':
                            score = float(line[60:66].strip())
                        elif score_source == 'occupancy':
                            score = float(line[54:60].strip())
                        else:
                            score = float(line[60:66].strip())
                        
                        if score > 0:
                            confidence_scores.append(score)
                    except (ValueError, IndexError):
                        continue
        
        if confidence_scores:
            return np.mean(confidence_scores)
        else:
            return None
    
    except Exception as e:
        print(f"Error in fallback parsing: {e}")
        return None


def main():
    parser = argparse.ArgumentParser(description='Generate confidence analysis CSV')
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
    
    metric_name = config.get('metric_name', 'pLDDT')
    print(f"Generating confidence analysis")
    print(f"Metric name: {metric_name}")
    print(f"Selection: {config.get('selection', 'all residues')}")
    print(f"Score metric: {config.get('score_metric', 'mean')}")
    
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
        
        # Context for datasheet resolution (if needed)
        context = {"datasheets": {}}  # TODO: Get actual context from pipeline
        
        # Process each structure file
        for structure_file in structure_files:
            structure_name = os.path.basename(structure_file)
            structure_id = os.path.splitext(structure_name)[0]
            print(f"Processing: {structure_name}")
            
            # Calculate confidence score
            score = calculate_confidence_score(structure_file, config, context)
            
            if score is None:
                print(f"  Could not calculate score - using 0.0")
                score = 0.0
            else:
                print(f"  {metric_name}={score:.2f}")
            
            # Create CSV row
            csv_row = {
                'id': structure_id,
                'source_structure': structure_name,
                metric_name: score,
                'selection': config.get('selection', 'all_residues'),
                'score_metric': config.get('score_metric', 'mean'),
                'score_source': config.get('score_source', 'bfactor')
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
            score_stats = df[metric_name].describe()
            print(f"\n{metric_name} statistics:")
            print(f"  Mean: {score_stats['mean']:.2f}")
            print(f"  Min: {score_stats['min']:.2f}")
            print(f"  Max: {score_stats['max']:.2f}")
        
    except Exception as e:
        print(f"Error during analysis: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()