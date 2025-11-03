#!/usr/bin/env python3
"""
Helper scripts for DynamicBind pipeline processing.

1. pipe_dynamicbind_prepare_ligands.py: Prepares ligands CSV with 'ligand' column from tool output
2. pipe_dynamicbind_table.py: Creates standardized table from DynamicBind outputs
"""

import os
import sys
import glob
import pandas as pd
import re


def prepare_ligands_csv(input_compounds_csv, output_ligands_csv, output_compounds_csv):
    """
    Prepare ligands CSV with 'ligand' column from tool output compounds table.

    Args:
        input_compounds_csv: Path to input compounds table (must have 'smiles' column)
        output_ligands_csv: Path to output ligands CSV for DynamicBind (with 'ligand' column)
        output_compounds_csv: Path to copy of input compounds for output reference
    """
    # Read input compounds
    try:
        df = pd.read_csv(input_compounds_csv)
    except Exception as e:
        print(f"Error reading compounds table: {e}")
        sys.exit(1)

    # Check for 'smiles' column
    if 'smiles' not in df.columns:
        print(f"Error: Input compounds table must have 'smiles' column. Found columns: {list(df.columns)}")
        sys.exit(1)

    # Create ligands CSV with 'ligand' column (rename 'smiles' to 'ligand')
    ligands_df = pd.DataFrame({'ligand': df['smiles']})
    ligands_df.to_csv(output_ligands_csv, index=False)
    print(f"Created ligands CSV with {len(ligands_df)} compounds: {output_ligands_csv}")

    # Copy compounds table to output folder
    df.to_csv(output_compounds_csv, index=False)
    print(f"Copied compounds table to output: {output_compounds_csv}")


def parse_dynamicbind_outputs(output_folder):
    """
    Parse DynamicBind output folder to create standardized table.

    Args:
        output_folder: Path to DynamicBind output folder

    Returns:
        DataFrame with columns: id, ligand_id, structure, affinity, lddt, rank
    """
    # Look for affinity prediction file
    affinity_file = os.path.join(output_folder, "affinity_prediction.csv")

    if not os.path.exists(affinity_file):
        print(f"Warning: affinity_prediction.csv not found in {output_folder}")
        return pd.DataFrame(columns=["id", "ligand_id", "structure", "affinity", "lddt", "rank"])

    # Read affinity predictions
    affinity_df = pd.read_csv(affinity_file)

    # Find all generated SDF files
    # Pattern: rank{N}_ligand_lddt{score}_affinity{score}_relaxed.sdf
    sdf_files = glob.glob(os.path.join(output_folder, "*.sdf"))

    results = []

    # Parse each SDF file to extract metadata
    sdf_pattern = re.compile(r"rank(\d+)_ligand_lddt([\d.]+)_affinity([\d.]+)(?:_relaxed)?\.sdf")

    for sdf_path in sdf_files:
        filename = os.path.basename(sdf_path)
        match = sdf_pattern.match(filename)

        if match:
            rank = int(match.group(1))
            lddt = float(match.group(2))
            affinity = float(match.group(3))

            # Generate ID based on rank
            structure_id = f"rank{rank}"

            results.append({
                "id": structure_id,
                "ligand_id": "ligand",  # Could be extracted from ligand file if needed
                "structure": sdf_path,
                "affinity": affinity,
                "lddt": lddt,
                "rank": rank
            })

    # Create dataframe and sort by rank
    if results:
        df = pd.DataFrame(results)
        df = df.sort_values("rank").reset_index(drop=True)
    else:
        df = pd.DataFrame(columns=["id", "ligand_id", "structure", "affinity", "lddt", "rank"])
        print(f"Warning: No SDF files found matching expected pattern in {output_folder}")

    return df


def main():
    """Main execution function."""
    if len(sys.argv) != 3:
        print("Usage: pipe_dynamicbind_table.py <output_folder> <output_csv>")
        sys.exit(1)

    output_folder = sys.argv[1]
    output_csv = sys.argv[2]

    # Parse DynamicBind outputs
    df = parse_dynamicbind_outputs(output_folder)

    # Save to CSV
    df.to_csv(output_csv, index=False)
    print(f"Created table with {len(df)} structures: {output_csv}")


if __name__ == "__main__":
    main()
