# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.


import argparse
import pandas as pd
import os
import glob
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from biopipelines.pdb_parser import parse_pdb_file

parser = argparse.ArgumentParser(description='Extract pLDDT confidence from ESMFold PDB B-factor columns')
parser.add_argument('output_folder', type=str, help='Path to folder containing ESMFold output PDB files')
parser.add_argument('confidence_csv', type=str, help='Output CSV file path for confidence metrics')

args = parser.parse_args()


def extract_esmfold_confidence(output_folder, confidence_csv):
    """
    Extract pLDDT confidence metrics from ESMFold PDB files.

    ESMFold stores per-residue pLDDT scores in the B-factor column of PDB files.
    This script reads each PDB, computes the mean pLDDT across all CA atoms,
    and writes a summary CSV.

    Args:
        output_folder: Path to folder with ESMFold output PDB files
        confidence_csv: Output CSV file path
    """
    confidence_data = []
    failed = []

    for pdb_file in sorted(glob.glob(os.path.join(output_folder, "*.pdb"))):
        basename = os.path.basename(pdb_file)
        seq_id = os.path.splitext(basename)[0]

        try:
            atoms = parse_pdb_file(pdb_file)

            # Extract B-factors from CA atoms (one per residue)
            ca_bfactors = []
            for atom in atoms:
                if atom.atom_name == "CA":
                    # B-factor is stored at a fixed position in PDB ATOM records
                    # parse_pdb_file returns Atom namedtuples; read B-factor directly from file
                    pass

            # parse_pdb_file doesn't expose B-factors, so read them directly
            ca_bfactors = _extract_bfactors_from_pdb(pdb_file)

            if not ca_bfactors:
                print(f"WARNING: No CA B-factors found in {basename}", file=sys.stderr)
                failed.append(seq_id)
                continue

            avg_plddt = sum(ca_bfactors) / len(ca_bfactors)

            confidence_data.append({
                'id': seq_id,
                'structure': pdb_file,
                'plddt': round(avg_plddt, 2)
            })
        except Exception as e:
            print(f"WARNING: {seq_id} failed: {e}", file=sys.stderr)
            failed.append(seq_id)

    # Always write results (even partial)
    if confidence_data:
        df = pd.DataFrame(confidence_data)
        df.to_csv(confidence_csv, index=False)
        print(f"Extracted confidence metrics for {len(confidence_data)} structures")
    else:
        # Write empty CSV with headers
        pd.DataFrame(columns=['id', 'structure', 'plddt']).to_csv(confidence_csv, index=False)
        print("Warning: No confidence data extracted")

    if failed:
        print(f"Failed {len(failed)}/{len(failed)+len(confidence_data)}: {failed}", file=sys.stderr)

    if not confidence_data:
        sys.exit(1)


def _extract_bfactors_from_pdb(pdb_file):
    """
    Extract B-factor values from CA atoms in a PDB file.

    PDB ATOM record format (columns are 1-indexed):
    - Columns 13-16: Atom name
    - Columns 61-66: B-factor (temperature factor)

    Returns:
        List of B-factor floats for CA atoms.
    """
    bfactors = []
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                atom_name = line[12:16].strip()
                if atom_name == "CA":
                    try:
                        bfactor = float(line[60:66].strip())
                        bfactors.append(bfactor)
                    except (ValueError, IndexError):
                        pass
    return bfactors


if __name__ == "__main__":
    extract_esmfold_confidence(args.output_folder, args.confidence_csv)
