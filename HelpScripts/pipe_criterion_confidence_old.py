#!/usr/bin/env python3
"""
Runtime helper script for Confidence criterion evaluation.

This script analyzes protein structures at runtime to extract confidence scores
(pLDDT) from PDB files, following the pipeline pattern.
"""

import os
import sys
import argparse
import json
import re
import numpy as np
import pandas as pd
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
    
    Args:
        position_str: Position specification string
        
    Returns:
        List of residue numbers
    """
    positions = []
    
    for part in str(position_str).split(','):
        part = part.strip()
        if '-' in part and not part.startswith('-'):
            # Range specification
            try:
                start, end = map(int, part.split('-'))
                positions.extend(range(start, end + 1))
            except ValueError:
                continue
        else:
            # Single position
            try:
                positions.append(int(part))
            except ValueError:
                continue
    
    return positions


def resolve_datasheet_selection(selection_ref: str, context: Dict[str, Any]) -> Optional[List[int]]:
    """
    Resolve datasheet-based selection to list of residue numbers.
    
    Args:
        selection_ref: Reference like "input.datasheets.sequences.designed_residues"
        context: Context containing datasheets
        
    Returns:
        List of residue numbers or None if resolution fails
    """
    try:
        # Parse the reference
        if not selection_ref.startswith("input.datasheets."):
            return None
        
        parts = selection_ref.split(".")
        if len(parts) < 4:
            return None
        
        datasheet_name = parts[2]  # e.g., "sequences"
        column_name = parts[3]     # e.g., "designed_residues"
        
        # Look for datasheet in context
        datasheets = context.get('datasheets', {})
        
        if isinstance(datasheets, dict) and datasheet_name in datasheets:
            datasheet_info = datasheets[datasheet_name]
            
            # Get path from datasheet info
            if isinstance(datasheet_info, dict) and 'path' in datasheet_info:
                datasheet_path = datasheet_info['path']
            else:
                datasheet_path = str(datasheet_info)
            
            if os.path.exists(datasheet_path):
                df = pd.read_csv(datasheet_path)
                
                if column_name in df.columns:
                    # Parse position specifications
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


def resolve_selection(selection: Union[str, List, None], context: Dict[str, Any]) -> Optional[List[int]]:
    """
    Resolve selection parameter to list of residue numbers.
    
    Args:
        selection: Selection parameter
        context: Context information
        
    Returns:
        List of residue numbers to analyze, or None for all residues
    """
    if selection is None:
        return None
    
    if isinstance(selection, list):
        if all(isinstance(x, int) for x in selection):
            # List of residue numbers
            return selection
        elif all(isinstance(x, str) for x in selection):
            # List of residue identifiers - extract numbers
            residue_numbers = []
            for res_id in selection:
                # Handle formats like "A45", "B12", or just "45"
                match = re.search(r'(\d+)', res_id)
                if match:
                    residue_numbers.append(int(match.group(1)))
            return residue_numbers if residue_numbers else None
    
    elif isinstance(selection, str):
        # Handle datasheet references
        if selection.startswith("input.datasheets."):
            return resolve_datasheet_selection(selection, context)
        else:
            # Try to parse as residue range or selection
            return parse_selection_string(selection)
    
    return None


def parse_selection_string(selection_str: str) -> Optional[List[int]]:
    """
    Parse selection string into residue numbers.
    
    Args:
        selection_str: Selection string
        
    Returns:
        List of residue numbers or None
    """
    # Handle common formats
    if re.match(r'^\d+-\d+$', selection_str):
        # Range format like "45-50"
        start, end = map(int, selection_str.split('-'))
        return list(range(start, end + 1))
    elif re.match(r'^[\d,\-\s]+$', selection_str):
        # Comma-separated format like "45,46,47" or "45-50,52,55-57"
        return parse_position_string(selection_str)
    
    return None


def extract_scores_mdanalysis(structure_file: str, target_residues: Optional[List[int]] = None,
                             score_source: str = "bfactor") -> List[float]:
    """Extract confidence scores using MDAnalysis."""
    if not MDANALYSIS_AVAILABLE:
        return []
    
    try:
        u = mda.Universe(structure_file)
        
        if target_residues is not None:
            # Select specific residues
            resid_selection = " or ".join([f"resid {res}" for res in target_residues])
            selection = u.select_atoms(f"protein and ({resid_selection})")
        else:
            # Select all protein atoms
            selection = u.select_atoms("protein")
        
        if len(selection) == 0:
            return []
        
        # Group by residue and take mean score per residue
        unique_residues = np.unique(selection.resids)
        residue_scores = []
        
        for resid in unique_residues:
            residue_atoms = selection.select_atoms(f"resid {resid}")
            if score_source == "bfactor":
                residue_score = np.mean(residue_atoms.tempfactors)
            elif score_source == "occupancy":
                residue_score = np.mean(residue_atoms.occupancies)
            else:
                residue_score = np.mean(residue_atoms.tempfactors)
            
            residue_scores.append(residue_score)
        
        return residue_scores
        
    except Exception as e:
        print(f"    Error extracting scores with MDAnalysis: {e}")
        return []


def extract_scores_biopython(structure_file: str, target_residues: Optional[List[int]] = None,
                           score_source: str = "bfactor") -> List[float]:
    """Extract confidence scores using BioPython."""
    if not BIOPYTHON_AVAILABLE:
        return []
    
    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("structure", structure_file)
        
        scores = []
        
        for model in structure:
            for chain in model:
                for residue in chain:
                    # Get residue number
                    res_num = residue.id[1]
                    
                    # Skip if not in target residues
                    if target_residues is not None and res_num not in target_residues:
                        continue
                    
                    # Extract scores from atoms (usually all atoms have same score)
                    residue_scores = []
                    for atom in residue:
                        if score_source == "bfactor":
                            score = atom.bfactor
                        elif score_source == "occupancy":
                            score = atom.occupancy
                        else:
                            score = atom.bfactor  # Default to bfactor
                        
                        residue_scores.append(score)
                    
                    # Use mean score for residue (they should all be the same for pLDDT)
                    if residue_scores:
                        scores.append(np.mean(residue_scores))
        
        return scores
        
    except Exception as e:
        print(f"    Error extracting scores with BioPython: {e}")
        return []


def extract_confidence_scores(structure_file: str, target_residues: Optional[List[int]] = None,
                            score_source: str = "bfactor") -> List[float]:
    """
    Extract confidence scores from structure file.
    
    Args:
        structure_file: Path to structure file
        target_residues: Specific residues to analyze, or None for all
        score_source: Source of scores in PDB ("bfactor", "occupancy")
        
    Returns:
        List of confidence scores
    """
    # Try MDAnalysis first, then BioPython
    if MDANALYSIS_AVAILABLE:
        scores = extract_scores_mdanalysis(structure_file, target_residues, score_source)
        if scores:
            return scores
    
    if BIOPYTHON_AVAILABLE:
        scores = extract_scores_biopython(structure_file, target_residues, score_source)
        if scores:
            return scores
    
    return []


def calculate_aggregate_score(scores: List[float], score_metric: str = "mean") -> Optional[float]:
    """
    Calculate aggregate confidence score based on metric.
    
    Args:
        scores: List of individual residue confidence scores
        score_metric: Aggregation method
        
    Returns:
        Aggregate score or None if no scores
    """
    if not scores:
        return None
    
    if score_metric == "mean":
        return float(np.mean(scores))
    elif score_metric == "min":
        return float(np.min(scores))
    elif score_metric == "max":
        return float(np.max(scores))
    elif score_metric == "median":
        return float(np.median(scores))
    else:
        return float(np.mean(scores))  # Default


def calculate_confidence_score(structure_file: str, config: Dict[str, Any], 
                             context: Dict[str, Any] = None) -> Optional[float]:
    """
    Calculate confidence score for a structure.
    
    Args:
        structure_file: Path to structure file
        config: Configuration with parameters
        context: Additional context (datasheets, etc.)
        
    Returns:
        Confidence score or None if calculation fails
    """
    try:
        selection = config["parameters"].get("selection")
        score_metric = config["parameters"].get("score_metric", "mean")
        score_source = config["parameters"].get("score_source", "bfactor")
        
        # Resolve selection for this structure
        target_residues = resolve_selection(selection, context or {})
        
        if target_residues:
            print(f"    Analyzing {len(target_residues)} selected residues")
        else:
            print(f"    Analyzing all residues")
        
        # Extract confidence scores
        residue_scores = extract_confidence_scores(structure_file, target_residues, score_source)
        
        if not residue_scores:
            print(f"    Could not extract confidence scores")
            return None
        
        # Calculate aggregate score
        aggregate_score = calculate_aggregate_score(residue_scores, score_metric)
        
        if aggregate_score is None:
            print(f"    Could not calculate aggregate score")
            return None
        
        return aggregate_score
        
    except Exception as e:
        print(f"    Error calculating confidence score: {e}")
        return None


def evaluate_expression(score: float, expression: str, variable_name: str) -> bool:
    """
    Evaluate the filtering expression with the calculated score.
    
    Args:
        score: Calculated score for the item
        expression: Boolean expression to evaluate
        variable_name: Variable name to replace in expression
        
    Returns:
        True if the item passes the criterion
    """
    # Replace variable name with actual score
    expr = expression.lower().replace(variable_name.lower(), str(score))
    
    # Replace logical operators for Python evaluation
    expr = expr.replace(' and ', ' & ').replace(' or ', ' | ')
    
    try:
        # Use eval safely with restricted globals
        safe_globals = {
            "__builtins__": {},
            "__name__": "__main__",
            "__doc__": None,
        }
        return bool(eval(expr, safe_globals))
    except Exception as e:
        raise ValueError(f"Error evaluating expression '{expression}' with {variable_name}={score}: {e}")


def main():
    """Main execution function."""
    parser = argparse.ArgumentParser(
        description="Evaluate Confidence criterion on protein structures"
    )
    
    parser.add_argument(
        "--structures_dir",
        required=True,
        help="Directory containing structure files to analyze"
    )
    
    parser.add_argument(
        "--output",
        required=True,
        help="Output JSON file with results"
    )
    
    parser.add_argument(
        "--config",
        required=True,
        help="JSON config file with criterion parameters"
    )
    
    args = parser.parse_args()
    
    try:
        # Check library availability
        if not (BIOPYTHON_AVAILABLE or MDANALYSIS_AVAILABLE):
            print("Error: Neither BioPython nor MDAnalysis is available")
            print("Install with: pip install biopython MDAnalysis")
            sys.exit(1)
        
        print(f"Structure analysis libraries available:")
        print(f"  BioPython: {BIOPYTHON_AVAILABLE}")
        print(f"  MDAnalysis: {MDANALYSIS_AVAILABLE}")
        
        # Load configuration
        with open(args.config, 'r') as f:
            config = json.load(f)
        
        print(f"Evaluating Confidence criterion:")
        print(f"  Selection: {config['parameters'].get('selection', 'all_residues')}")
        print(f"  Expression: {config['expression']}")
        print(f"  Score metric: {config['parameters'].get('score_metric', 'mean')}")
        print(f"  Score source: {config['parameters'].get('score_source', 'bfactor')}")
        
        # Find structure files
        structure_files = []
        for ext in ['.pdb', '.cif', '.mmcif']:
            pattern_files = [f for f in os.listdir(args.structures_dir) if f.lower().endswith(ext)]
            structure_files.extend([os.path.join(args.structures_dir, f) for f in pattern_files])
        
        print(f"Found {len(structure_files)} structure files")
        
        # Process each structure
        results = {
            "criterion_type": config["criterion_type"],
            "criterion_class": config["criterion_class"],
            "expression": config["expression"],
            "variable_name": config["variable_name"],
            "parameters": config["parameters"],
            "kept_items": [],
            "filtered_items": [],
            "item_scores": {},
            "total_input": 0,
            "kept_count": 0,
            "filtered_count": 0,
            "pass_rate": 0.0
        }
        
        # Context for datasheet resolution (if needed)
        context = {}  # Could be extended with actual datasheet context
        
        for structure_file in structure_files:
            if not os.path.exists(structure_file):
                continue
                
            print(f"Processing: {os.path.basename(structure_file)}")
            
            # Calculate confidence score
            score = calculate_confidence_score(structure_file, config, context)
            
            if score is None:
                # Could not calculate score - filter out
                results["filtered_items"].append(structure_file)
                results["item_scores"][structure_file] = 0.0
                print(f"  Could not calculate score")
                continue
            
            results["item_scores"][structure_file] = score
            
            # Evaluate expression
            try:
                passes = evaluate_expression(score, config["expression"], config["variable_name"])
                
                if passes:
                    results["kept_items"].append(structure_file)
                    print(f"  {config['variable_name']}={score:.2f} ✓")
                else:
                    results["filtered_items"].append(structure_file)
                    print(f"  {config['variable_name']}={score:.2f} ✗")
                    
            except Exception as e:
                # Expression evaluation failed - filter out
                results["filtered_items"].append(structure_file)
                print(f"  Expression error: {e}")
        
        # Calculate summary statistics
        results["total_input"] = len(results["kept_items"]) + len(results["filtered_items"])
        results["kept_count"] = len(results["kept_items"])
        results["filtered_count"] = len(results["filtered_items"])
        results["pass_rate"] = results["kept_count"] / results["total_input"] if results["total_input"] > 0 else 0.0
        
        # Save results
        with open(args.output, 'w') as f:
            json.dump(results, f, indent=2)
        
        print(f"Results saved to {args.output}")
        print(f"Summary: {results['kept_count']}/{results['total_input']} structures kept ({results['pass_rate']:.1%})")
        
    except Exception as e:
        print(f"Error during criterion evaluation: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()