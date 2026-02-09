#!/usr/bin/env python3
# Copyright (c) 2026 Gianluca Quargnali @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Compound library processing script for the pipeline system.

Handles SMILES library generation, filtering, expansion, and preparation
for structure prediction tools like Boltz2.
"""

import pandas as pd
import argparse
import os
import sys
from typing import List, Dict, Any, Optional

def validate_smiles(smiles: str) -> bool:
    """
    Validate SMILES string.
    
    Args:
        smiles: SMILES string to validate
        
    Returns:
        True if valid, False otherwise
    """
    try:
        # Try to import rdkit for validation
        from rdkit import Chem
        mol = Chem.MolFromSmiles(smiles)
        return mol is not None
    except ImportError:
        # If rdkit not available, do basic validation
        if not smiles or len(smiles.strip()) == 0:
            return False
        # Very basic SMILES validation - check for balanced parentheses
        paren_count = smiles.count('(') - smiles.count(')')
        bracket_count = smiles.count('[') - smiles.count(']')
        return paren_count == 0 and bracket_count == 0


def process_library(input_files: List[str], 
                   output_file: str,
                   library_type: str = "SMILES",
                   primary_key: Optional[str] = None,
                   max_compounds: Optional[int] = None,
                   validate_smiles_flag: bool = True,
                   remove_duplicates: bool = True,
                   ) -> None:
    """
    Process compound library files.
    
    Args:
        input_files: List of input file paths
        output_file: Output file path
        library_type: Type of library (SMILES, SDF, CCD, ChEMBL)
        primary_key: Primary key column name
        max_compounds: Maximum number of compounds
        validate_smiles_flag: Whether to validate SMILES
        remove_duplicates: Whether to remove duplicates
    """
    
    # Read input files
    all_compounds = []
    
    for input_file in input_files:
        print(f"Processing {input_file}...")
        
        if input_file.endswith('.csv'):
            df = pd.read_csv(input_file)
        elif input_file.endswith('.txt'):
            # Assume simple text file with one SMILES per line
            with open(input_file, 'r') as f:
                smiles_list = [line.strip() for line in f if line.strip()]
            df = pd.DataFrame({'SMILES': smiles_list})
        else:
            print(f"Warning: Unsupported file format: {input_file}")
            continue
        
        all_compounds.append(df)
    
    if not all_compounds:
        print("Error: No valid input files found")
        sys.exit(1)
    
    # Combine all dataframes
    combined_df = pd.concat(all_compounds, ignore_index=True)
    
    print(f"Initial compounds: {len(combined_df)}")
    
    # Identify SMILES column
    smiles_col = None
    for col in combined_df.columns:
        if col.lower() in ['smiles', 'smi', 'canonical_smiles']:
            smiles_col = col
            break
    
    if smiles_col is None and 'SMILES' not in combined_df.columns:
        # Add SMILES column if not found
        if len(combined_df.columns) == 1:
            combined_df.rename(columns={combined_df.columns[0]: 'SMILES'}, inplace=True)
            smiles_col = 'SMILES'
        else:
            print("Error: Could not identify SMILES column")
            sys.exit(1)
    else:
        smiles_col = smiles_col or 'SMILES'
    
    # Validate SMILES if requested
    if validate_smiles_flag:
        print("Validating SMILES...")
        valid_mask = combined_df[smiles_col].apply(validate_smiles)
        invalid_count = (~valid_mask).sum()
        if invalid_count > 0:
            print(f"Removed {invalid_count} invalid SMILES")
            combined_df = combined_df[valid_mask]
    
    # Remove duplicates if requested
    if remove_duplicates:
        initial_count = len(combined_df)
        combined_df = combined_df.drop_duplicates(subset=[smiles_col])
        removed_count = initial_count - len(combined_df)
        if removed_count > 0:
            print(f"Removed {removed_count} duplicate compounds")
    
    
    # Apply compound limit
    if max_compounds and len(combined_df) > max_compounds:
        print(f"Limiting to top {max_compounds} compounds")
        combined_df = combined_df.head(max_compounds)
    
    # Add ID column if primary key specified and not present
    if primary_key and primary_key not in combined_df.columns:
        combined_df[primary_key] = range(1, len(combined_df) + 1)
    
    # Ensure standard columns exist
    if 'id' not in combined_df.columns:
        combined_df['id'] = range(1, len(combined_df) + 1)
    
    print(f"Final compounds: {len(combined_df)}")
    
    # Save output
    if output_file.endswith('.csv'):
        combined_df.to_csv(output_file, index=False)
    elif output_file.endswith('.txt'):
        # Save only SMILES column
        combined_df[smiles_col].to_csv(output_file, index=False, header=False)
    else:
        # Default to CSV
        combined_df.to_csv(output_file, index=False)
    
    print(f"Library saved to: {output_file}")

def main():
    parser = argparse.ArgumentParser(description="Process compound libraries")
    parser.add_argument("input_files", nargs="+", help="Input library files")
    parser.add_argument("--output", required=True, help="Output file path")
    parser.add_argument("--type", default="SMILES", choices=["SMILES", "SDF", "CCD", "ChEMBL"], 
                       help="Library type")
    parser.add_argument("--primary-key", help="Primary key column name")
    parser.add_argument("--max-compounds", type=int, help="Maximum number of compounds")
    parser.add_argument("--no-validate", action="store_true", help="Skip SMILES validation")
    parser.add_argument("--keep-duplicates", action="store_true", help="Keep duplicate compounds")
    
    args = parser.parse_args()
    
    process_library(
        input_files=args.input_files,
        output_file=args.output,
        library_type=args.type,
        primary_key=args.primary_key,
        max_compounds=args.max_compounds,
        validate_smiles_flag=not args.no_validate,
        remove_duplicates=not args.keep_duplicates
    )

if __name__ == "__main__":
    main()