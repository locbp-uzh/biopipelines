#!/usr/bin/env python3
"""
Runtime helper script for concatenating tool outputs.

This script concatenates structures, sequences, compounds and datasheets
from multiple tool outputs into a unified result for iterative optimization cycles.
"""

import os
import sys
import argparse
import json
import pandas as pd
import shutil
from typing import Dict, List, Any, Optional
from glob import glob


def concatenate_structures(pool_configs: List[Dict[str, Any]], output_dir: str) -> List[str]:
    """
    Concatenate structure files from multiple pools.
    
    Args:
        pool_configs: List of pool configurations
        output_dir: Output directory for concatenated structures
        
    Returns:
        List of output structure file paths
    """
    print("Concatenating structure files...")
    os.makedirs(output_dir, exist_ok=True)
    
    output_files = []
    structure_counter = 0
    
    for pool_config in pool_configs:
        output_folder = pool_config.get('output_folder')
        if not output_folder or not os.path.exists(output_folder):
            print(f"Warning: Pool output folder not found: {output_folder}")
            continue
        
        # Search for structure files in predicted directories
        structure_dirs = pool_config.get('structure_dirs', [output_folder])
        found_structures = []
        
        for structure_dir in structure_dirs:
            if os.path.exists(structure_dir):
                # Find PDB files
                pdb_files = glob(os.path.join(structure_dir, "*.pdb"))
                found_structures.extend(pdb_files)
                
                # If structure_dir is the output folder itself, look for PDB files there too
                if structure_dir == output_folder:
                    pdb_files = glob(os.path.join(output_folder, "*.pdb"))
                    found_structures.extend(pdb_files)
                break
        
        print(f"Found {len(found_structures)} structures in {pool_config['prefix']}")
        
        # Copy structures to output directory with new names
        for structure_file in found_structures:
            if os.path.exists(structure_file):
                output_name = f"concatenated_structure_{structure_counter:03d}.pdb"
                output_path = os.path.join(output_dir, output_name)
                shutil.copy2(structure_file, output_path)
                output_files.append(output_path)
                structure_counter += 1
    
    print(f"Concatenated {len(output_files)} structure files")
    return output_files


def concatenate_compounds(pool_configs: List[Dict[str, Any]], output_dir: str) -> List[str]:
    """
    Concatenate compound files from multiple pools.
    
    Args:
        pool_configs: List of pool configurations
        output_dir: Output directory for concatenated compounds
        
    Returns:
        List of output compound file paths
    """
    print("Concatenating compound files...")
    os.makedirs(output_dir, exist_ok=True)
    
    output_files = []
    compound_counter = 0
    
    for pool_config in pool_configs:
        output_folder = pool_config.get('output_folder')
        if not output_folder or not os.path.exists(output_folder):
            continue
        
        # Search for compound files
        compound_files = []
        for ext in ['*.sdf', '*.mol', '*.mol2']:
            compound_files.extend(glob(os.path.join(output_folder, ext)))
            compound_files.extend(glob(os.path.join(output_folder, 'compounds', ext)))
        
        print(f"Found {len(compound_files)} compounds in {pool_config['prefix']}")
        
        # Copy compounds to output directory with new names
        for compound_file in compound_files:
            if os.path.exists(compound_file):
                file_ext = os.path.splitext(compound_file)[1]
                output_name = f"concatenated_compound_{compound_counter:03d}{file_ext}"
                output_path = os.path.join(output_dir, output_name)
                shutil.copy2(compound_file, output_path)
                output_files.append(output_path)
                compound_counter += 1
    
    print(f"Concatenated {len(output_files)} compound files")
    return output_files


def concatenate_datasheets(data_configs: List[Dict[str, Any]], output_csv: str) -> None:
    """
    Concatenate datasheets from multiple data sources.
    
    Args:
        data_configs: List of data source configurations
        output_csv: Output CSV file path
    """
    print("Concatenating datasheets...")
    
    all_dataframes = []
    
    for i, data_config in enumerate(data_configs):
        prefix = data_config['prefix']
        datasheets = data_config.get('datasheets', {})
        
        if not datasheets:
            print(f"Warning: No datasheets found for {prefix}")
            continue
        
        # Find the main datasheet to use
        csv_path = None
        priority_names = ['analysis', 'merged', 'filtered', 'combined', 'datasheet', 'main']
        
        for name in priority_names:
            if name in datasheets:
                csv_path = datasheets[name]
                break
        
        if not csv_path:
            # Use first available datasheet
            csv_path = next(iter(datasheets.values()))
        
        if not os.path.exists(csv_path):
            print(f"Warning: Datasheet file not found: {csv_path}")
            continue
        
        try:
            df = pd.read_csv(csv_path)
            print(f"Loaded {prefix}: {df.shape[0]} rows, {df.shape[1]} columns")
            
            # Add source pool information
            df['source_pool'] = prefix
            df['source_pool_index'] = i
            
            all_dataframes.append(df)
            
        except Exception as e:
            print(f"Error loading datasheet {csv_path}: {e}")
            continue
    
    if not all_dataframes:
        print("Warning: No datasheets could be loaded, creating empty output")
        # Create minimal empty dataframe
        empty_df = pd.DataFrame({
            'id': [],
            'source_pool': [],
            'source_pool_index': []
        })
        empty_df.to_csv(output_csv, index=False)
        return
    
    # Concatenate all dataframes
    print(f"Concatenating {len(all_dataframes)} datasheets...")
    combined_df = pd.concat(all_dataframes, ignore_index=True, sort=False)
    
    # Ensure unique IDs by adding source pool info
    if 'id' in combined_df.columns:
        combined_df['original_id'] = combined_df['id']
        combined_df['id'] = combined_df['source_pool'] + '_' + combined_df['id'].astype(str)
    
    print(f"Final concatenated datasheet: {combined_df.shape[0]} rows, {combined_df.shape[1]} columns")
    
    # Create output directory
    os.makedirs(os.path.dirname(output_csv), exist_ok=True)
    
    # Save concatenated dataframe
    combined_df.to_csv(output_csv, index=False)
    print(f"Concatenated datasheet saved: {output_csv}")


def concatenate_datasheets_only(config_data: Dict[str, Any]) -> None:
    """
    Concatenate datasheets vertically (like SQL UNION).
    
    Args:
        config_data: Configuration dictionary with datasheet inputs and settings
    """
    datasheet_configs = config_data['datasheet_configs']
    fill = config_data.get('fill', None)
    output_csv = config_data['output_csv']
    
    print(f"Concatenating {len(datasheet_configs)} datasheets")
    print(f"Fill strategy: {fill if fill is not None else 'remove non-common columns'}")
    
    # Load all datasheets
    dataframes = []
    for i, config in enumerate(datasheet_configs):
        datasheet_path = config['datasheet_path']
        
        if not os.path.exists(datasheet_path):
            print(f"Warning: Datasheet {i} not found: {datasheet_path}")
            continue
            
        try:
            df = pd.read_csv(datasheet_path)
            df['source_datasheet'] = f"datasheet_{i}"
            dataframes.append(df)
            print(f"Loaded datasheet {i}: {df.shape[0]} rows from {datasheet_path}")
        except Exception as e:
            print(f"Error loading datasheet {i} from {datasheet_path}: {e}")
            continue
    
    if not dataframes:
        print("Error: No datasheets could be loaded")
        return
    
    # Concatenate datasheets
    if fill is None:
        # Remove non-common columns - use intersection of all columns
        common_columns = set(dataframes[0].columns)
        for df in dataframes[1:]:
            common_columns = common_columns.intersection(set(df.columns))
        common_columns = list(common_columns)
        
        print(f"Using common columns: {common_columns}")
        concatenated_df = pd.concat([df[common_columns] for df in dataframes], ignore_index=True)
    else:
        # Fill missing columns with specified value
        concatenated_df = pd.concat(dataframes, ignore_index=True).fillna(fill)
    
    # Save concatenated result
    os.makedirs(os.path.dirname(output_csv), exist_ok=True)
    concatenated_df.to_csv(output_csv, index=False)
    
    print(f"\nConcatenation completed successfully!")
    print(f"Total rows: {len(concatenated_df)}")
    print(f"Total columns: {len(concatenated_df.columns)}")
    print(f"Output: {output_csv}")


def main():
    parser = argparse.ArgumentParser(description='Concatenate datasheets vertically')
    parser.add_argument('--config', required=True, help='JSON config file with concatenation parameters')
    
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
    required_params = ['datasheet_configs', 'output_csv']
    for param in required_params:
        if param not in config_data:
            print(f"Error: Missing required parameter: {param}")
            sys.exit(1)
    
    try:
        concatenate_datasheets_only(config_data)
        
    except Exception as e:
        print(f"Error concatenating datasheets: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()