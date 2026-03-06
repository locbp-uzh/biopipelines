# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Generate 2D molecule images (PNG) from a compounds CSV using RDKit.

Usage:
    python pipe_compound_images.py <compounds_csv> <images_folder>
"""

import sys
import csv
import os

def main():
    if len(sys.argv) != 3:
        print("Usage: pipe_compound_images.py <compounds_csv> <images_folder>")
        sys.exit(1)

    compounds_csv = sys.argv[1]
    images_folder = sys.argv[2]

    from rdkit import Chem
    from rdkit.Chem import Draw

    os.makedirs(images_folder, exist_ok=True)

    with open(compounds_csv, newline='') as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    generated = 0
    skipped = 0
    for row in rows:
        cid = row.get('id', '').strip()
        smiles = row.get('smiles', '').strip()
        if not cid or not smiles:
            skipped += 1
            continue
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"  Warning: could not parse SMILES for {cid}: {smiles}")
            skipped += 1
            continue
        img = Draw.MolToImage(mol, size=(300, 300))
        img.save(os.path.join(images_folder, f"{cid}.png"))
        generated += 1

    print(f"Generated {generated} compound images (skipped {skipped})")


if __name__ == "__main__":
    main()
