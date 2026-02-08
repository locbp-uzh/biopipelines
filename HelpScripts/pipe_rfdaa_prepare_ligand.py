#!/usr/bin/env python3
# Copyright (c) 2026 Gianluca Quargnali @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# This software is freely available for use, modification, and redistribution.
# If you use this software or any derivative of it in published work,
# you must cite the original author and this repository.
"""
Helper script for RFDAA_PrepareLigand tool.

Adds a dummy peptide (1PEF) to ligand-only PDB files to satisfy RFdiffusion-AllAtom's
requirement for protein atoms in the input structure. The peptide is positioned away
from the ligand and assigned to chain A with residue numbers starting from 3.
"""

import os
import sys
import argparse
import pandas as pd
from typing import List, Tuple


def download_1pef(pdbs_folder: str) -> str:
    """
    Download 1PEF structure if not already present.

    Args:
        pdbs_folder: Directory to save/check for 1PEF.pdb

    Returns:
        Path to 1PEF.pdb file
    """
    pef_path = os.path.join(pdbs_folder, "1PEF.pdb")

    if os.path.exists(pef_path):
        print(f"Found 1PEF locally: {pef_path}")
        return pef_path

    print(f"Downloading 1PEF from RCSB...")
    try:
        import requests

        url = "https://files.rcsb.org/download/1PEF.pdb"
        headers = {
            'User-Agent': 'BioPipelines-RFDAA/1.0 (https://gitlab.uzh.ch/locbp/public/biopipelines)'
        }

        response = requests.get(url, headers=headers, timeout=30)
        response.raise_for_status()

        # Save to PDBs folder
        os.makedirs(pdbs_folder, exist_ok=True)
        with open(pef_path, 'w') as f:
            f.write(response.text)

        print(f"Downloaded 1PEF to: {pef_path}")
        return pef_path

    except Exception as e:
        print(f"Error downloading 1PEF: {str(e)}")
        sys.exit(1)


def read_pdb_lines(pdb_path: str) -> List[str]:
    """
    Read PDB file and return ATOM/HETATM lines.

    Args:
        pdb_path: Path to PDB file

    Returns:
        List of ATOM/HETATM lines
    """
    lines = []
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                lines.append(line.rstrip('\n'))
    return lines


def renumber_peptide_lines(lines: List[str], target_chain: str = 'A', start_resnum: int = 3) -> List[str]:
    """
    Renumber peptide residues to target chain and starting residue number.

    Args:
        lines: List of ATOM lines from peptide
        target_chain: Target chain ID (default: 'A')
        start_resnum: Starting residue number (default: 3)

    Returns:
        List of renumbered ATOM lines
    """
    renumbered = []

    # Get unique residue numbers from input
    res_nums = []
    for line in lines:
        if line.startswith('ATOM'):
            res_num = int(line[22:26].strip())
            if res_num not in res_nums:
                res_nums.append(res_num)

    # Create mapping from old to new residue numbers
    res_map = {}
    for i, old_res in enumerate(sorted(res_nums)):
        res_map[old_res] = start_resnum + i

    # Renumber lines
    for line in lines:
        if line.startswith('ATOM'):
            old_res = int(line[22:26].strip())
            new_res = res_map[old_res]

            # Reconstruct line with new chain and residue number
            new_line = (
                line[:21] +           # Up to and including atom name
                target_chain +         # New chain ID
                f"{new_res:4d}" +     # New residue number
                line[26:]              # Rest of line
            )
            renumbered.append(new_line)

    return renumbered


def translate_coordinates(lines: List[str], dx: float, dy: float, dz: float) -> List[str]:
    """
    Translate all coordinates by given offsets.

    Args:
        lines: List of ATOM/HETATM lines
        dx, dy, dz: Translation offsets in Angstroms

    Returns:
        List of lines with translated coordinates
    """
    translated = []

    for line in lines:
        if line.startswith(('ATOM', 'HETATM')):
            # Parse coordinates
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())

            # Translate
            x += dx
            y += dy
            z += dz

            # Reconstruct line with new coordinates
            new_line = (
                line[:30] +
                f"{x:8.3f}" +
                f"{y:8.3f}" +
                f"{z:8.3f}" +
                line[54:]
            )
            translated.append(new_line)
        else:
            translated.append(line)

    return translated


def combine_structures(ligand_pdb: str, peptide_pdb: str, output_pdb: str):
    """
    Combine ligand and peptide structures into single PDB file.

    The peptide is renumbered to chain A, residues starting from 3.
    The ligand is kept as-is (typically chain A, residue 1).
    The peptide is translated 20 Angstroms away in the x-direction.

    Args:
        ligand_pdb: Path to ligand PDB file
        peptide_pdb: Path to peptide PDB file (1PEF)
        output_pdb: Path to output combined PDB file
    """
    # Read ligand lines
    print(f"Reading ligand: {ligand_pdb}")
    ligand_lines = read_pdb_lines(ligand_pdb)
    print(f"  Found {len(ligand_lines)} ligand atoms")

    # Read peptide lines
    print(f"Reading peptide: {peptide_pdb}")
    peptide_lines = read_pdb_lines(peptide_pdb)
    print(f"  Found {len(peptide_lines)} peptide atoms")

    # Renumber peptide to chain A, residues starting from 3
    print("Renumbering peptide to chain A, residues 3+")
    peptide_renumbered = renumber_peptide_lines(peptide_lines, target_chain='A', start_resnum=3)

    # Translate peptide away from ligand (20 Angstroms in x-direction)
    print("Translating peptide 20 Angstroms away")
    peptide_translated = translate_coordinates(peptide_renumbered, dx=20.0, dy=0.0, dz=0.0)

    # Combine structures
    print(f"Writing combined structure: {output_pdb}")
    with open(output_pdb, 'w') as f:
        f.write("REMARK   1 Combined ligand + dummy peptide for RFdiffusion-AllAtom\n")
        f.write("REMARK   1 Ligand: chain A, residue 1 (original)\n")
        f.write("REMARK   1 Peptide: chain A, residues 3+ (from 1PEF, translated +20A in x)\n")

        # Write ligand
        for line in ligand_lines:
            f.write(line + '\n')

        # Write peptide
        for line in peptide_translated:
            f.write(line + '\n')

        f.write("END\n")

    print(f"Successfully created combined structure with {len(ligand_lines) + len(peptide_translated)} atoms")


def create_output_table(output_pdb: str, output_csv: str):
    """
    Create structures.csv output table.

    Args:
        output_pdb: Path to prepared PDB file
        output_csv: Path to output CSV file
    """
    df = pd.DataFrame({
        'id': ['prepared_ligand'],
        'file_path': [output_pdb]
    })
    df.to_csv(output_csv, index=False)
    print(f"Created output table: {output_csv}")


def main():
    parser = argparse.ArgumentParser(
        description='Prepare ligand structure for RFdiffusion-AllAtom by adding dummy peptide'
    )
    parser.add_argument('--ligand_pdb', required=True, help='Input ligand PDB file')
    parser.add_argument('--output_pdb', required=True, help='Output combined PDB file')
    parser.add_argument('--output_csv', required=True, help='Output structures CSV file')
    parser.add_argument('--pdbs_folder', required=True, help='PDBs folder for caching 1PEF')

    args = parser.parse_args()

    # Validate input
    if not os.path.exists(args.ligand_pdb):
        print(f"Error: Ligand PDB file not found: {args.ligand_pdb}")
        sys.exit(1)

    try:
        # Download/get 1PEF
        peptide_pdb = download_1pef(args.pdbs_folder)

        # Combine structures
        combine_structures(args.ligand_pdb, peptide_pdb, args.output_pdb)

        # Create output table
        create_output_table(args.output_pdb, args.output_csv)

        print("\nSuccess! Ligand prepared for RFdiffusion-AllAtom")

    except Exception as e:
        print(f"Error: {str(e)}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
