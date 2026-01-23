#!/usr/bin/env python3
"""
Runtime helper script for SASA (Solvent Accessible Surface Area) analysis.

Calculates delta SASA for protein-ligand complexes:
- SASA of ligand alone (fully exposed)
- SASA of ligand in complex (partially buried by protein)
- delta_SASA = SASA_alone - SASA_complex (positive = buried surface)

Uses PyMOL's get_area command for SASA calculations.
"""

import os
import sys
import argparse
import pandas as pd
from typing import List, Tuple

# Import PyMOL
import pymol
from pymol import cmd


def calculate_sasa(structure_path: str, ligand_resn: str, dot_density: int = 4) -> Tuple[float, float, float]:
    """
    Calculate SASA for a ligand in complex and alone.

    Args:
        structure_path: Path to PDB/CIF structure file
        ligand_resn: Ligand residue name (e.g., "X", "LIG", "AMX")
        dot_density: Dot density for SASA calculation (1-4)

    Returns:
        Tuple of (sasa_alone, sasa_complex, delta_sasa)
    """
    # Clear PyMOL state
    cmd.delete("all")

    # Load structure
    cmd.load(structure_path, "complex")

    # Set dot density for SASA calculation
    cmd.set("dot_solvent", 1)
    cmd.set("dot_density", dot_density)

    # Select ligand and protein
    ligand_sel = f"resn {ligand_resn}"
    protein_sel = f"polymer and not {ligand_sel}"

    # Check if ligand exists
    ligand_count = cmd.count_atoms(ligand_sel)
    if ligand_count == 0:
        print(f"Warning: No atoms found for ligand '{ligand_resn}' in {structure_path}")
        return 0.0, 0.0, 0.0

    # Calculate SASA of ligand in complex (with protein present)
    sasa_complex = cmd.get_area(ligand_sel)

    # Remove protein to calculate SASA of ligand alone
    cmd.remove(protein_sel)
    sasa_alone = cmd.get_area(ligand_sel)

    # Calculate delta SASA (positive = surface buried by protein)
    delta_sasa = sasa_alone - sasa_complex

    return sasa_alone, sasa_complex, delta_sasa


def process_structures(structure_paths: List[str], ligand_resn: str,
                       output_csv: str, dot_density: int = 4) -> None:
    """
    Process multiple structures and calculate SASA for each.

    Args:
        structure_paths: List of structure file paths
        ligand_resn: Ligand residue name
        output_csv: Output CSV file path
        dot_density: Dot density for SASA calculation
    """
    results = []

    for struct_path in structure_paths:
        if not os.path.exists(struct_path):
            print(f"Warning: Structure file not found: {struct_path}")
            continue

        struct_id = os.path.splitext(os.path.basename(struct_path))[0]
        print(f"Processing: {struct_id}")

        try:
            sasa_alone, sasa_complex, delta_sasa = calculate_sasa(
                struct_path, ligand_resn, dot_density
            )

            results.append({
                "id": struct_id,
                "structure": struct_path,
                "sasa_ligand_alone": round(sasa_alone, 2),
                "sasa_ligand_complex": round(sasa_complex, 2),
                "delta_sasa": round(delta_sasa, 2)
            })

            print(f"  SASA alone: {sasa_alone:.2f} Å²")
            print(f"  SASA in complex: {sasa_complex:.2f} Å²")
            print(f"  Delta SASA: {delta_sasa:.2f} Å² (buried)")

        except Exception as e:
            print(f"Error processing {struct_path}: {e}")
            results.append({
                "id": struct_id,
                "structure": struct_path,
                "sasa_ligand_alone": 0.0,
                "sasa_ligand_complex": 0.0,
                "delta_sasa": 0.0
            })

    # Save results
    df = pd.DataFrame(results)
    df.to_csv(output_csv, index=False)
    print(f"\nResults saved to: {output_csv}")
    print(f"Processed {len(results)} structures")


def main():
    parser = argparse.ArgumentParser(
        description="Calculate delta SASA for protein-ligand complexes"
    )
    parser.add_argument(
        "--structures",
        required=True,
        help="File containing list of structure file paths (one per line)"
    )
    parser.add_argument(
        "--ligand",
        required=True,
        help="Ligand residue name (e.g., X, LIG, AMX)"
    )
    parser.add_argument(
        "--output_csv",
        required=True,
        help="Output CSV file path"
    )
    parser.add_argument(
        "--dot_density",
        type=int,
        default=4,
        help="Dot density for SASA calculation (1-4, default: 4)"
    )

    args = parser.parse_args()

    # Read structure paths from file (one path per line)
    with open(args.structures, 'r') as f:
        structure_paths = [line.strip() for line in f if line.strip()]

    if not structure_paths:
        print(f"Error: No structures found in: {args.structures}")
        sys.exit(1)

    # Initialize PyMOL in quiet mode
    pymol.finish_launching(['pymol', '-qc'])

    # Process structures
    process_structures(
        structure_paths=structure_paths,
        ligand_resn=args.ligand,
        output_csv=args.output_csv,
        dot_density=args.dot_density
    )

    # Quit PyMOL
    cmd.quit()


if __name__ == "__main__":
    main()
