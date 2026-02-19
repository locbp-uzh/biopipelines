#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
SMILES properties calculation script for the pipeline system.

Calculates molecular properties for compounds in a library.
"""

import pandas as pd
import argparse
import os
import sys
from typing import Dict, Any, Optional

def calculate_basic_properties(smiles: str) -> Dict[str, Any]:
    """
    Calculate basic molecular properties from SMILES.
    
    Args:
        smiles: SMILES string
        
    Returns:
        Dictionary of molecular properties
    """
    properties = {
        'MW': None,
        'LogP': None,
        'HBD': None,
        'HBA': None,
        'TPSA': None,
        'Rotatable_Bonds': None,
        'Rings': None
    }
    
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors, rdMolDescriptors
        
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return properties
        
        properties['MW'] = Descriptors.MolWt(mol)
        properties['LogP'] = Descriptors.MolLogP(mol)
        properties['HBD'] = Descriptors.NumHDonors(mol)
        properties['HBA'] = Descriptors.NumHAcceptors(mol)
        properties['TPSA'] = Descriptors.TPSA(mol)
        properties['Rotatable_Bonds'] = Descriptors.NumRotatableBonds(mol)
        properties['Rings'] = rdMolDescriptors.CalcNumRings(mol)
        
    except ImportError:
        # If rdkit not available, return None values
        print("Warning: RDKit not available, cannot calculate molecular properties")
    except Exception as e:
        print(f"Warning: Error calculating properties for {smiles}: {e}")
    
    return properties

def calculate_properties_for_library(input_file: str, output_file: str) -> None:
    """
    Calculate properties for all compounds in a library.
    
    Args:
        input_file: Input CSV file with SMILES
        output_file: Output CSV file with properties
    """
    # Read input library
    df = pd.read_csv(input_file)
    
    # Identify SMILES column
    smiles_col = None
    for col in df.columns:
        if col.lower() in ['smiles', 'smi', 'canonical_smiles']:
            smiles_col = col
            break
    
    if smiles_col is None:
        print("Error: Could not identify SMILES column")
        sys.exit(1)
    
    print(f"Calculating properties for {len(df)} compounds...")
    
    # Calculate properties for each compound
    property_data = []
    for idx, row in df.iterrows():
        smiles = row[smiles_col]
        props = calculate_basic_properties(smiles)
        
        # Combine original data with properties
        row_data = row.to_dict()
        row_data.update(props)
        property_data.append(row_data)
        
        if (idx + 1) % 100 == 0:
            print(f"Processed {idx + 1}/{len(df)} compounds")
    
    # Create output dataframe
    output_df = pd.DataFrame(property_data)
    
    # Save to file
    output_df.to_csv(output_file, index=False)
    print(f"Properties saved to: {output_file}")
    
    # Print summary statistics
    numeric_cols = ['MW', 'LogP', 'HBD', 'HBA', 'TPSA', 'Rotatable_Bonds', 'Rings']
    valid_props = [col for col in numeric_cols if col in output_df.columns and output_df[col].notna().any()]
    
    if valid_props:
        print("\nProperty Statistics:")
        print("===================")
        for prop in valid_props:
            values = output_df[prop].dropna()
            if len(values) > 0:
                print(f"{prop}: mean={values.mean():.2f}, min={values.min():.2f}, max={values.max():.2f}")

def main():
    parser = argparse.ArgumentParser(description="Calculate molecular properties for compound library")
    parser.add_argument("input_file", help="Input CSV file with SMILES")
    parser.add_argument("output_file", help="Output CSV file with properties")
    
    args = parser.parse_args()
    
    if not os.path.exists(args.input_file):
        print(f"Error: Input file not found: {args.input_file}")
        sys.exit(1)
    
    calculate_properties_for_library(args.input_file, args.output_file)

if __name__ == "__main__":
    main()