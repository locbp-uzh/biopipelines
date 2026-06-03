#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

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

# Import unified I/O utilities
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from biopipelines.biopipelines_io import load_datastream, iterate_files
from biopipelines.pdb_parser import relative_accessibility

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

    # PyMOL auto-flags HETATM atoms as "ignore" on load; get_area silently
    # skips ignored atoms and returns 0. Clear the flag on the ligand so its
    # SASA is actually computed.
    cmd.flag("ignore", ligand_sel, "clear")

    # Calculate SASA of ligand in complex (with protein present)
    sasa_complex = cmd.get_area(ligand_sel)

    # Remove protein to calculate SASA of ligand alone
    cmd.remove(protein_sel)
    sasa_alone = cmd.get_area(ligand_sel)

    # Calculate delta SASA (positive = surface buried by protein)
    delta_sasa = sasa_alone - sasa_complex

    return sasa_alone, sasa_complex, delta_sasa


def process_structures(structures_ds, ligand_resn: str,
                       output_csv: str, dot_density: int = 4) -> None:
    """
    Process multiple structures and calculate SASA for each.

    Args:
        structures_ds: DataStream with structure files
        ligand_resn: Ligand residue name
        output_csv: Output CSV file path
        dot_density: Dot density for SASA calculation
    """
    results = []

    # Use iterate_files for proper ID-file matching
    for struct_id, struct_path in iterate_files(structures_ds):
        if not os.path.exists(struct_path):
            print(f"Warning: Structure file not found: {struct_path}")
            continue

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


def calculate_residue_sasa(structure_path: str, dot_density: int = 4) -> List[dict]:
    """Per-residue protein SASA. Returns rows of chain, resi, resn, sasa."""
    cmd.delete("all")
    cmd.load(structure_path, "s")
    cmd.set("dot_solvent", 1)
    cmd.set("dot_density", dot_density)

    # load_b=1 stamps each atom's SASA contribution into its B-factor.
    cmd.get_area("polymer", load_b=1)

    # resv is PyMOL's integer residue number; the iterate namespace has no builtins.
    atom_rows = []
    cmd.iterate("polymer", "atom_rows.append((chain, resv, resn, b))",
                space={"atom_rows": atom_rows})

    per_residue = {}
    order = []
    for chain, resi, resn, b in atom_rows:
        key = (chain, resi, resn)
        if key not in per_residue:
            per_residue[key] = 0.0
            order.append(key)
        per_residue[key] += b

    rows = []
    for chain, resi, resn in order:
        rows.append({
            "chain": chain if chain else "A",
            "resi": resi,
            "resn": resn,
            "sasa": round(per_residue[(chain, resi, resn)], 2),
        })
    return rows


def process_residue_sasa(structures_ds, accessibility_dir: str,
                         accessibility_map_csv: str, dot_density: int = 4) -> None:
    """Per-residue SASA + rsa for each structure, one resi-csv per structure."""
    os.makedirs(accessibility_dir, exist_ok=True)
    cols = ["id", "chain", "resi", "resn", "sasa", "rsa"]
    map_rows = []

    for struct_id, struct_path in iterate_files(structures_ds):
        if not os.path.exists(struct_path):
            print(f"Warning: Structure file not found: {struct_path}")
            continue
        print(f"Processing: {struct_id}")
        residue_rows = calculate_residue_sasa(struct_path, dot_density)
        out_rows = []
        for r in residue_rows:
            rsa = relative_accessibility(r["sasa"], r["resn"])
            out_rows.append({
                "id": struct_id,
                "chain": r["chain"],
                "resi": r["resi"],
                "resn": r["resn"],
                "sasa": r["sasa"],
                "rsa": "" if rsa is None else rsa,
            })
        out_path = os.path.join(accessibility_dir, f"{struct_id}.csv")
        pd.DataFrame(out_rows, columns=cols).to_csv(out_path, index=False)
        map_rows.append({"id": struct_id, "file": out_path})
        print(f"  {len(out_rows)} residues -> {out_path}")

    os.makedirs(os.path.dirname(accessibility_map_csv), exist_ok=True)
    pd.DataFrame(map_rows, columns=["id", "file"]).to_csv(accessibility_map_csv, index=False)
    print(f"\nAccessibility resi-csv map: {accessibility_map_csv} ({len(map_rows)} rows)")
    if not map_rows:
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description="Calculate delta SASA for protein-ligand complexes"
    )
    parser.add_argument(
        "--mode",
        default="ligand",
        choices=["ligand", "residues"],
        help="ligand: delta-SASA per complex; residues: per-residue accessibility"
    )
    parser.add_argument(
        "--structures",
        required=True,
        help="JSON file containing DataStream with structure files"
    )
    parser.add_argument(
        "--ligand",
        default=None,
        help="Ligand residue name (e.g., X, LIG, AMX); required for mode=ligand"
    )
    parser.add_argument(
        "--output_csv",
        default=None,
        help="Output CSV file path (mode=ligand)"
    )
    parser.add_argument(
        "--accessibility-dir",
        default=None,
        help="Directory for per-residue resi-csv files (mode=residues)"
    )
    parser.add_argument(
        "--accessibility-map-csv",
        default=None,
        help="Accessibility stream map CSV (mode=residues)"
    )
    parser.add_argument(
        "--dot_density",
        type=int,
        default=4,
        help="Dot density for SASA calculation (1-4, default: 4)"
    )

    args = parser.parse_args()

    # Load structures DataStream using pipe_biopipelines_io
    structures_ds = load_datastream(args.structures)

    if not structures_ds.ids_expanded:
        print(f"Error: No structures found in: {args.structures}")
        sys.exit(1)

    print(f"Loaded {len(structures_ds.ids_expanded)} structures from DataStream")

    # Initialize PyMOL in quiet mode
    pymol.finish_launching(['pymol', '-qc'])

    if args.mode == "residues":
        process_residue_sasa(
            structures_ds=structures_ds,
            accessibility_dir=args.accessibility_dir,
            accessibility_map_csv=args.accessibility_map_csv,
            dot_density=args.dot_density,
        )
    else:
        if not args.ligand:
            print("Error: --ligand is required for mode=ligand")
            sys.exit(1)
        process_structures(
            structures_ds=structures_ds,
            ligand_resn=args.ligand,
            output_csv=args.output_csv,
            dot_density=args.dot_density
        )

    # Quit PyMOL
    cmd.quit()


if __name__ == "__main__":
    main()
