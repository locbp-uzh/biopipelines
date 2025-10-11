#!/usr/bin/env python3
"""
Prepare compounds CSV with 'ligand' column from tool output compounds datasheet.
DynamicBind reads from this file which contains both full compound info and the 'ligand' column.

Usage: pipe_dynamicbind_prepare_ligands.py <input_compounds_csv> <output_compounds_csv>
"""

import os
import sys
import pandas as pd


def main():
    """Main execution function."""
    if len(sys.argv) != 3:
        print("Usage: pipe_dynamicbind_prepare_ligands.py <input_compounds_csv> <output_compounds_csv>")
        sys.exit(1)

    input_compounds_csv = sys.argv[1]
    output_compounds_csv = sys.argv[2]

    # Read input compounds
    try:
        df = pd.read_csv(input_compounds_csv)
    except Exception as e:
        print(f"Error reading compounds datasheet: {e}")
        sys.exit(1)

    # Check for 'smiles' column
    if 'smiles' not in df.columns:
        print(f"Error: Input compounds datasheet must have 'smiles' column. Found columns: {list(df.columns)}")
        sys.exit(1)

    # Add 'ligand' column to compounds datasheet if not already present
    # DynamicBind expects a column named 'ligand' with SMILES strings
    if 'ligand' not in df.columns:
        df['ligand'] = df['smiles']

    # Save compounds datasheet to output folder with 'ligand' column
    df.to_csv(output_compounds_csv, index=False)
    print(f"Created compounds CSV with {len(df)} compounds and 'ligand' column: {output_compounds_csv}")


if __name__ == "__main__":
    main()
