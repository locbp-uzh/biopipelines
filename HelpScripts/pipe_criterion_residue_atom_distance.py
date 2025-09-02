#!/usr/bin/env python3
"""
Runtime helper script for ResidueAtomDistance criterion evaluation.

This script analyzes protein structures at runtime to calculate distances between
specific atoms and residues, following the pipeline pattern.
"""

import os
import sys
import argparse
import json
import numpy as np
import pandas as pd
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
    
    Args:
        selection: Selection string
        
    Returns:
        Dictionary with parsed selection criteria
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


def get_coordinates_pymol(structure_file: str, selection_info: Dict[str, Any]) -> List[Tuple[float, float, float]]:
    """
    Get atom coordinates using PyMOL based on selection.
    
    Args:
        structure_file: Path to structure file
        selection_info: Parsed selection information
        
    Returns:
        List of (x, y, z) coordinates
    """
    if not PYMOL_AVAILABLE:
        return []
    
    try:
        # Initialize PyMOL in quiet mode
        pymol.finish_launching(['pymol', '-cq'])
        cmd.reinitialize()
        
        # Load the structure
        obj_name = "struct"
        cmd.load(structure_file, obj_name)
        
        coords = []
        selection_type = selection_info.get("selection_type", "raw")
        
        if selection_type == "entity_atom":
            entity = selection_info["entity"]
            atom_type = selection_info["atom_type"]
            
            if entity == "protein":
                # Select protein atoms
                if atom_type == "CA":
                    selection = f"{obj_name} and polymer and name CA"
                else:
                    selection = f"{obj_name} and polymer and name {atom_type}"
                    
            elif entity == "ligand":
                # Select ligand atoms (non-polymer, non-water)
                if atom_type.upper() == "CL":
                    # Try different chlorine naming conventions
                    selection = f"{obj_name} and not polymer and not solvent and (name CL or name CL1 or name CL01 or name Cl or name Cl1)"
                elif atom_type.lower() == "all":
                    selection = f"{obj_name} and not polymer and not solvent"
                else:
                    selection = f"{obj_name} and not polymer and not solvent and name {atom_type}"
            
            # Get coordinates
            model = cmd.get_model(selection)
            for atom in model.atom:
                coords.append((atom.coord[0], atom.coord[1], atom.coord[2]))
        
        elif selection_type == "sequence_context":
            entity = selection_info["entity"]
            residue_type = selection_info["residue_type"]
            sequence_context = selection_info["sequence_context"]
            
            if entity == "protein":
                # Get sequence and find matches
                chains = cmd.get_chains(obj_name)
                
                for chain in chains:
                    # Get residue sequence for this chain
                    residues = []
                    cmd.iterate(f"{obj_name} and chain {chain} and polymer and name CA",
                              "residues.append((resi, resn))", space={'residues': residues})
                    
                    if not residues:
                        continue
                        
                    # Build sequence string (convert 3-letter to 1-letter)
                    aa_map = {
                        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
                        'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
                        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
                        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
                    }
                    
                    sequence = ""
                    residue_ids = []
                    for resi, resn in residues:
                        if resn in aa_map:
                            sequence += aa_map[resn]
                            residue_ids.append(resi)
                    
                    # Find sequence context matches
                    context_len = len(sequence_context)
                    for i in range(len(sequence) - context_len + 1):
                        if sequence[i:i + context_len] == sequence_context:
                            # Found match, look for the specific residue type
                            target_1letter = residue_type
                            for j in range(context_len):
                                if sequence[i + j] == target_1letter:
                                    # This is the target residue
                                    target_resi = residue_ids[i + j]
                                    # Get all atoms from this residue
                                    selection = f"{obj_name} and chain {chain} and resi {target_resi}"
                                    model = cmd.get_model(selection)
                                    for atom in model.atom:
                                        coords.append((atom.coord[0], atom.coord[1], atom.coord[2]))
        
        # Clean up PyMOL
        cmd.delete("all")
        
        return coords
        
    except Exception as e:
        print(f"Error extracting coordinates with PyMOL from {structure_file}: {e}")
        return []


def get_coordinates_mdanalysis(structure_file: str, selection_info: Dict[str, Any]) -> List[Tuple[float, float, float]]:
    """
    Get atom coordinates using MDAnalysis based on selection.
    
    Args:
        structure_file: Path to structure file
        selection_info: Parsed selection information
        
    Returns:
        List of (x, y, z) coordinates
    """
    if not MDANALYSIS_AVAILABLE:
        return []
    
    try:
        u = mda.Universe(structure_file)
        
        selection_type = selection_info.get("selection_type", "raw")
        
        if selection_type == "entity_atom":
            entity = selection_info["entity"]
            atom_type = selection_info["atom_type"]
            
            if entity == "protein":
                base_sel = "protein"
            elif entity == "ligand":
                base_sel = "not protein and not resname HOH"
            else:
                base_sel = ""
            
            if atom_type.lower() == "all":
                atom_group = u.select_atoms(base_sel)
            else:
                atom_group = u.select_atoms(f"{base_sel} and name {atom_type}")
        
        elif selection_type == "mdanalysis":
            atom_group = u.select_atoms(selection_info["raw_selection"])
        
        else:
            # Try raw selection
            atom_group = u.select_atoms(selection_info.get("raw_selection", "all"))
        
        return [tuple(pos) for pos in atom_group.positions]
        
    except Exception as e:
        print(f"Error extracting coordinates with MDAnalysis from {structure_file}: {e}")
        return []


def get_coordinates(structure_file: str, selection: str) -> List[Tuple[float, float, float]]:
    """
    Get coordinates for a selection using the best available library.
    
    Args:
        structure_file: Path to structure file
        selection: Selection string
        
    Returns:
        List of (x, y, z) coordinates
    """
    selection_info = parse_selection_string(selection)
    
    # Try PyMOL first, then MDAnalysis
    if PYMOL_AVAILABLE:
        coords = get_coordinates_pymol(structure_file, selection_info)
        if coords:
            return coords
    
    if MDANALYSIS_AVAILABLE:
        coords = get_coordinates_mdanalysis(structure_file, selection_info)
        if coords:
            return coords
    
    return []


def calculate_residue_atom_distance(structure_file: str, config: Dict[str, Any]) -> Optional[float]:
    """
    Calculate distance between atoms and residues for a structure.
    
    Args:
        structure_file: Path to structure file
        config: Configuration with atom, residue selections and distance metric
        
    Returns:
        Distance value or None if calculation fails
    """
    try:
        atom_selection = config["parameters"]["atom"]
        residue_selection = config["parameters"]["residue"]
        distance_metric = config["parameters"].get("distance_metric", "min")
        
        # Get atom coordinates
        atom_coords = get_coordinates(structure_file, atom_selection)
        residue_coords = get_coordinates(structure_file, residue_selection)
        
        if not atom_coords or not residue_coords:
            print(f"    Could not get coordinates: atoms={len(atom_coords)}, residues={len(residue_coords)}")
            return None
        
        # Calculate all pairwise distances
        distances = []
        for atom_coord in atom_coords:
            for residue_coord in residue_coords:
                dist = np.linalg.norm(np.array(atom_coord) - np.array(residue_coord))
                distances.append(dist)
        
        if not distances:
            return None
        
        # Apply distance metric
        if distance_metric == "min" or distance_metric == "closest":
            return float(np.min(distances))
        elif distance_metric == "max":
            return float(np.max(distances))
        elif distance_metric == "mean":
            return float(np.mean(distances))
        else:
            return float(np.min(distances))  # Default to minimum
            
    except Exception as e:
        print(f"    Error calculating distance: {e}")
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
        description="Evaluate ResidueAtomDistance criterion on protein structures"
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
        if not (PYMOL_AVAILABLE or MDANALYSIS_AVAILABLE):
            print("Error: Neither PyMOL nor MDAnalysis is available")
            print("Install with: pip install pymol-open-source MDAnalysis")
            sys.exit(1)
        
        print(f"Structure analysis libraries available:")
        print(f"  PyMOL: {PYMOL_AVAILABLE}")
        print(f"  MDAnalysis: {MDANALYSIS_AVAILABLE}")
        
        # Load configuration
        with open(args.config, 'r') as f:
            config = json.load(f)
        
        print(f"Evaluating ResidueAtomDistance criterion:")
        print(f"  Atom: {config['parameters']['atom']}")
        print(f"  Residue: {config['parameters']['residue']}")
        print(f"  Expression: {config['expression']}")
        print(f"  Metric: {config['parameters'].get('distance_metric', 'min')}")
        
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
        
        for structure_file in structure_files:
            if not os.path.exists(structure_file):
                continue
                
            print(f"Processing: {os.path.basename(structure_file)}")
            
            # Calculate distance score
            score = calculate_residue_atom_distance(structure_file, config)
            
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