#!/usr/bin/env python3
"""
PyMOL session creation helper script for BioPipelines.

Creates PyMOL sessions from structure outputs with alignment, coloring,
and multi-tool structure combination support.
"""

import argparse
import json
import os
import sys
from pathlib import Path
from typing import Dict, List, Any, Optional

import pymol
from pymol import cmd


def setup_pymol():
    """Initialize PyMOL in headless mode with standard settings."""
    pymol.pymol_argv = ['pymol', '-c']  # -c for no GUI
    pymol.finish_launching()
    
    # Standard visualization settings
    cmd.do("show cartoon")
    cmd.set("seq_view", 1)
    cmd.set("cartoon_gap_cutoff", 20)
    cmd.set("cartoon_sampling", 20)
    cmd.set("sphere_scale", 0.2)
    cmd.set("ray_trace_mode", 1)
    cmd.set("ray_shadows", 0)
    cmd.set("spec_reflect", 0)
    cmd.set("ray_trace_frames", 1)
    cmd.set("ray_trace_color", "gray20")


def setup_plddt_coloring():
    """Set up standard pLDDT color scheme."""
    # AlphaFold confidence colors
    blue_rgb = [0, 76, 202]
    lightblue_rgb = [73, 196, 238]  
    yellow_rgb = [255, 213, 57]
    orange_rgb = [255, 113, 67]
    
    # Convert to 0-1 range
    blue = [c/255.0 for c in blue_rgb]
    lightblue = [c/255.0 for c in lightblue_rgb]
    yellow = [c/255.0 for c in yellow_rgb]
    orange = [c/255.0 for c in orange_rgb]
    
    cmd.set_color('blue_plddt', blue)
    cmd.set_color('lightblue_plddt', lightblue)
    cmd.set_color('yellow_plddt', yellow)
    cmd.set_color('orange_plddt', orange)


def apply_plddt_coloring(selection="all"):
    """Apply pLDDT-based coloring to selection."""
    # Color ranges based on b-factor (pLDDT scores)
    cmd.color('blue_plddt', f"({selection}) & ! b < 90.0 & ! b > 100.0")
    cmd.color('lightblue_plddt', f"({selection}) & ! b < 70.0 & ! b > 90.0") 
    cmd.color('yellow_plddt', f"({selection}) & ! b < 50.0 & ! b > 70.0")
    cmd.color('orange_plddt', f"({selection}) & b < 50.0")


def calculate_average_plddt(object_name: str) -> float:
    """Calculate average pLDDT for a structure."""
    try:
        atom_iterator = cmd.get_model(f"{object_name} and name CA")
        residues_inspected = []
        plddt_sum = 0
        
        for atom in atom_iterator.atom:
            resi = int(atom.resi)
            if resi in residues_inspected:
                continue
            plddt_sum += atom.b
            residues_inspected.append(resi)
        
        if residues_inspected:
            return plddt_sum / len(residues_inspected)
        return 0.0
    except Exception as e:
        print(f"Error calculating pLDDT for {object_name}: {e}")
        return 0.0


def load_structures_from_source(source_config: Dict[str, Any]) -> List[str]:
    """Load all PDB structures from a source directory."""
    input_folder = source_config['input_folder']
    prefix = source_config['prefix']
    
    loaded_objects = []
    
    # Find all PDB files in the source directory
    pdb_files = []
    if os.path.exists(input_folder):
        for file in os.listdir(input_folder):
            if file.endswith('.pdb'):
                pdb_files.append(file)
    
    if not pdb_files:
        print(f"Warning: No PDB files found in {input_folder}")
        return loaded_objects
    
    # Load each structure with prefixed name
    for pdb_file in pdb_files:
        pdb_path = os.path.join(input_folder, pdb_file)
        structure_name = os.path.splitext(pdb_file)[0]
        object_name = f"{prefix}_{structure_name}"
        
        try:
            cmd.load(pdb_path, object_name)
            loaded_objects.append(object_name)
            print(f"Loaded: {object_name} from {pdb_file}")
        except Exception as e:
            print(f"Error loading {pdb_file}: {e}")
    
    return loaded_objects


def load_datasheet_values(datasheet_path: str, column: str) -> Dict[str, float]:
    """Load values from a datasheet CSV for coloring."""
    if not os.path.exists(datasheet_path):
        print(f"Warning: Datasheet not found: {datasheet_path}")
        return {}
    
    try:
        import pandas as pd
        df = pd.read_csv(datasheet_path)
        
        if 'id' not in df.columns:
            print(f"Warning: 'id' column not found in {datasheet_path}")
            return {}
        
        if column not in df.columns:
            print(f"Warning: Column '{column}' not found in {datasheet_path}")
            return {}
        
        # Convert to dictionary mapping id -> value
        return dict(zip(df['id'], df[column]))
    
    except Exception as e:
        print(f"Error reading datasheet {datasheet_path}: {e}")
        return {}


def apply_custom_coloring(loaded_objects: List[str], color_values: Dict[str, float], prefix: str):
    """Apply custom coloring based on datasheet values."""
    if not color_values:
        print("No color values available, using default coloring")
        return
    
    # Create color gradient based on values
    values = list(color_values.values())
    if not values:
        return
    
    min_val = min(values)
    max_val = max(values)
    val_range = max_val - min_val if max_val != min_val else 1.0
    
    for obj_name in loaded_objects:
        # Extract structure ID from object name (remove prefix)
        if obj_name.startswith(f"{prefix}_"):
            structure_id = obj_name[len(prefix)+1:]
        else:
            structure_id = obj_name
        
        if structure_id in color_values:
            value = color_values[structure_id]
            # Normalize to 0-1 range
            normalized = (value - min_val) / val_range
            
            # Create color based on normalized value (blue to red gradient)
            red = normalized
            blue = 1.0 - normalized
            green = 0.5
            
            color_name = f"custom_{structure_id}"
            cmd.set_color(color_name, [red, green, blue])
            cmd.color(color_name, obj_name)
            
            print(f"Colored {obj_name} with value {value:.2f} (normalized: {normalized:.2f})")


def perform_alignment(loaded_objects: List[str], alignment_method: str, reference_obj: Optional[str] = None):
    """Perform alignment of structures."""
    if alignment_method == "none" or len(loaded_objects) <= 1:
        return
    
    # Determine reference structure
    if reference_obj and reference_obj in loaded_objects:
        ref_structure = reference_obj
    else:
        ref_structure = loaded_objects[0]
    
    print(f"Using {ref_structure} as reference for alignment")
    
    # Align all other structures to reference
    for obj_name in loaded_objects:
        if obj_name == ref_structure:
            continue
        
        try:
            if alignment_method == "align":
                rmsd = cmd.align(f"{obj_name}", f"{ref_structure}")[0]
                print(f"Aligned {obj_name} to {ref_structure}: RMSD = {rmsd:.2f}")
            elif alignment_method == "cealign":
                result = cmd.cealign(f"{ref_structure}", f"{obj_name}")
                rmsd = result.get("RMSD", 0.0) if isinstance(result, dict) else 0.0
                print(f"CE-aligned {obj_name} to {ref_structure}: RMSD = {rmsd:.2f}")
            elif alignment_method == "super":
                rmsd = cmd.super(f"{obj_name}", f"{ref_structure}")[0]
                print(f"Super-aligned {obj_name} to {ref_structure}: RMSD = {rmsd:.2f}")
            elif alignment_method == "alignto":
                cmd.alignto(ref_structure)
                print(f"Aligned all structures to {ref_structure}")
                break  # alignto aligns all at once
                
        except Exception as e:
            print(f"Error aligning {obj_name}: {e}")


def main():
    parser = argparse.ArgumentParser(description='Create PyMOL session from pipeline structures')
    parser.add_argument('--config', required=True, help='JSON configuration file')
    parser.add_argument('--output', required=True, help='Output PyMOL session file (.pse)')
    
    args = parser.parse_args()
    
    # Load configuration
    with open(args.config, 'r') as f:
        config = json.load(f)
    
    # Initialize PyMOL
    setup_pymol()
    setup_plddt_coloring()
    
    all_loaded_objects = []
    color_values = {}
    
    # Load structures from all sources
    for source_config in config['structure_sources']:
        print(f"\nLoading structures from {source_config['prefix']}...")
        loaded_objects = load_structures_from_source(source_config)
        all_loaded_objects.extend(loaded_objects)
    
    if not all_loaded_objects:
        print("Error: No structures were loaded")
        sys.exit(1)
    
    print(f"\nTotal structures loaded: {len(all_loaded_objects)}")
    
    # Apply coloring
    color_by = config.get('color_by')
    if color_by:
        print(f"\nApplying coloring based on: {color_by}")
        # Parse color_by format: tool.output.datasheets.sheet.column
        parts = color_by.split('.')
        if len(parts) >= 5:
            column = parts[-1]
            sheet = parts[-2]
            
            # Try to find corresponding datasheet files
            # This is a simplified approach - in practice, we'd need to resolve
            # the tool reference and find the correct datasheet
            for source_config in config['structure_sources']:
                datasheet_path = os.path.join(source_config['input_folder'], f"{sheet}.csv")
                if os.path.exists(datasheet_path):
                    values = load_datasheet_values(datasheet_path, column)
                    color_values.update(values)
                    break
            
            if color_values:
                apply_custom_coloring(all_loaded_objects, color_values, "")
            else:
                print("No color values found, applying pLDDT coloring")
                for obj in all_loaded_objects:
                    apply_plddt_coloring(obj)
        else:
            print("Invalid color_by format, applying pLDDT coloring")
            for obj in all_loaded_objects:
                apply_plddt_coloring(obj)
    else:
        # Default pLDDT coloring
        print("Applying default pLDDT coloring")
        for obj in all_loaded_objects:
            apply_plddt_coloring(obj)
    
    # Perform alignment
    alignment_method = config.get('alignment', 'align')
    reference_structure = config.get('reference_structure')
    if reference_structure:
        # Find the actual object name with the reference prefix
        reference_obj = None
        for obj in all_loaded_objects:
            if obj.startswith(f"{reference_structure}_"):
                reference_obj = obj
                break
    else:
        reference_obj = None
    
    print(f"\nPerforming alignment using method: {alignment_method}")
    perform_alignment(all_loaded_objects, alignment_method, reference_obj)
    
    # Print summary statistics
    print(f"\nStructure summary:")
    for obj in all_loaded_objects:
        avg_plddt = calculate_average_plddt(obj)
        print(f"  {obj}: Average pLDDT = {avg_plddt:.1f}")
    
    # Save PyMOL session
    print(f"\nSaving PyMOL session: {args.output}")
    cmd.save(args.output)
    
    # Quit PyMOL
    cmd.quit()
    print("PyMOL session creation completed successfully")


if __name__ == "__main__":
    main()