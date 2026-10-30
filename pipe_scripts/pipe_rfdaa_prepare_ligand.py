#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Helper script for RFDAA_PrepareLigand tool.

Adds a dummy peptide (1PEF) to ligand-only PDB files to satisfy RFdiffusion-AllAtom's
requirement for protein atoms in the input structure. The peptide is positioned away
from the ligand and assigned to chain A with residue numbers starting from 3.

Each input ligand is prepared under its own id (<structures_folder>/<id>.pdb), so the output stream is 1:1 with the input stream and keeps its ids.
"""

import os
import sys
import argparse
import pandas as pd
from typing import List, Tuple

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from biopipelines.biopipelines_io import load_datastream, iterate_files  # noqa: E402


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
            'User-Agent': 'BioPipelines-RFDAA/1.0 (https://github.com/locbp-uzh/biopipelines)'
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

        f.write("END".ljust(80) + "\n")

    print(f"Successfully created combined structure with {len(ligand_lines) + len(peptide_translated)} atoms")


def create_output_table(prepared: List[Tuple[str, str]], output_csv: str):
    """
    Create structures.csv output table.

    Args:
        prepared: (id, prepared PDB path) of every ligand
        output_csv: Path to output CSV file
    """
    df = pd.DataFrame({
        'id': [pid for pid, _ in prepared],
        'file_path': [path for _, path in prepared]
    })
    df.to_csv(output_csv, index=False)
    print(f"Created output table: {output_csv}")


def create_structures_map(prepared: List[Tuple[str, str]], structures_map: str):
    """Write the structures stream map_table (id | file)."""
    pd.DataFrame({'id': [pid for pid, _ in prepared], 'file': [path for _, path in prepared]}).to_csv(structures_map, index=False)
    print(f"Created structures map: {structures_map}")


def main():
    parser = argparse.ArgumentParser(
        description='Prepare ligand structure for RFdiffusion-AllAtom by adding dummy peptide'
    )
    parser.add_argument('--ligand_json', required=True, help='Input ligand structures DataStream JSON')
    parser.add_argument('--structures_folder', required=True, help='Folder receiving <id>.pdb for every ligand')
    parser.add_argument('--output_csv', required=True, help='Output structures CSV file')
    parser.add_argument('--structures_map', required=True, help='Output structures map_table CSV')
    parser.add_argument('--pdbs_folder', required=True, help='PDBs folder for caching 1PEF')

    args = parser.parse_args()

    ligands = list(iterate_files(load_datastream(args.ligand_json)))
    if not ligands:
        print("Error: no ligand structures to prepare")
        sys.exit(1)
    absent = [path for _, path in ligands if not os.path.exists(path)]
    if absent:
        print(f"Error: Ligand PDB file not found: {absent[0]}")
        sys.exit(1)

    try:
        # Download/get 1PEF
        peptide_pdb = download_1pef(args.pdbs_folder)

        # Combine structures, one output per input id
        os.makedirs(args.structures_folder, exist_ok=True)
        prepared = []
        for ligand_id, ligand_pdb in ligands:
            output_pdb = os.path.join(args.structures_folder, f"{ligand_id}.pdb")
            combine_structures(ligand_pdb, peptide_pdb, output_pdb)
            prepared.append((ligand_id, output_pdb))

        # Create output table + stream map_table
        create_output_table(prepared, args.output_csv)
        create_structures_map(prepared, args.structures_map)

        print("\nSuccess! Ligand prepared for RFdiffusion-AllAtom")

    except Exception as e:
        print(f"Error: {str(e)}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
