#!/usr/bin/env python3
"""
Prepare ligands CSV with 'ligand' column from tool output compounds datasheet.

Usage: pipe_dynamicbind_prepare_ligands.py <input_compounds_csv> <output_ligands_csv> <output_compounds_csv>
"""

import os
import sys
import pandas as pd


def main():
    """Main execution function."""
    if len(sys.argv) != 4:
        print("Usage: pipe_dynamicbind_prepare_ligands.py <input_compounds_csv> <output_ligands_csv> <output_compounds_csv>")
        sys.exit(1)

    input_compounds_csv = sys.argv[1]
    output_ligands_csv = sys.argv[2]
    output_compounds_csv = sys.argv[3]

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

    # Create ligands CSV with 'ligand' column (rename 'smiles' to 'ligand')
    ligands_df = pd.DataFrame({'ligand': df['smiles']})
    ligands_df.to_csv(output_ligands_csv, index=False)
    print(f"Created ligands CSV with {len(ligands_df)} compounds: {output_ligands_csv}")

    # Copy compounds datasheet to output folder
    df.to_csv(output_compounds_csv, index=False)
    print(f"Copied compounds datasheet to output: {output_compounds_csv}")


if __name__ == "__main__":
    main()
