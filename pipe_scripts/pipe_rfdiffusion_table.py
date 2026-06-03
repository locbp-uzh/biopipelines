# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.


"""
Create RFdiffusion results table from TRB files.

Reads per-design .trb pickle files produced by RFdiffusion / RFdiffusion-AllAtom
to extract fixed/designed region information in a chain-aware way, replacing the
fragile log-file parsing approach.
"""

import argparse
import pickle
import os
import sys
import pandas as pd

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from biopipelines.sele_utils import chain_aware_sele


def pdb_residues(pdb_path):
    """Return ordered list of unique (chain, resnum) tuples from a PDB file."""
    seen = set()
    residues = []
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            chain = line[21]
            resnum = int(line[22:26].strip())
            key = (chain, resnum)
            if key not in seen:
                seen.add(key)
                residues.append(key)
    return residues


def parse_trb(trb_path, pdb_path):
    """Parse a TRB file and return fixed, designed, source_fixed, plddt_mean.

    Args:
        trb_path: path to .trb pickle file
        pdb_path: path to matching .pdb output file (used to derive all residues)

    Returns:
        dict with keys: fixed, designed, source_fixed, plddt_mean
    """
    with open(trb_path, "rb") as f:
        trb = pickle.load(f)

    fixed_residues = trb.get("con_hal_pdb_idx", [])
    source_residues = trb.get("con_ref_pdb_idx", [])

    fixed = chain_aware_sele(fixed_residues)
    source_fixed = chain_aware_sele(source_residues)

    # Designed = all output residues minus fixed
    fixed_set = set(fixed_residues)
    all_residues = pdb_residues(pdb_path)
    designed_residues = [r for r in all_residues if r not in fixed_set]
    designed = chain_aware_sele(designed_residues)

    # plddt: shape (T, N) in RFdiffusion; not present in AllAtom
    plddt_mean = ""
    if "plddt" in trb:
        plddt = trb["plddt"]
        # Use the final timestep row
        final = plddt[-1] if plddt.ndim == 2 else plddt
        plddt_mean = round(float(final.mean()), 4)

    return {
        "fixed": fixed,
        "designed": designed,
        "source_fixed": source_fixed,
        "plddt_mean": plddt_mean,
    }


def main():
    parser = argparse.ArgumentParser(description="Create RFdiffusion results table from TRB files")
    parser.add_argument("output_folder", type=str, help="RFdiffusion output folder")
    parser.add_argument("table_path", type=str, help="Output CSV table path")

    args = parser.parse_args()

    # Scan the structures folder for every produced design. This covers both
    # the single-PDB / unconditional case (one prefix) and the multi-PDB case
    # (one <pdb_id> prefix per input structure) without reconstructing names.
    import glob
    pdb_paths = sorted(glob.glob(os.path.join(args.output_folder, "*.pdb")))

    designs = []
    for pdb_path in pdb_paths:
        pdb_file = os.path.basename(pdb_path)
        design_id = os.path.splitext(pdb_file)[0]
        trb_path = os.path.join(args.output_folder, f"{design_id}.trb")

        if os.path.exists(trb_path):
            trb_data = parse_trb(trb_path, pdb_path)
            status = "ok"
        else:
            trb_data = {"fixed": "", "designed": "", "source_fixed": "", "plddt_mean": ""}
            status = "missing_trb"

        designs.append({
            "id": design_id,
            "pdb": pdb_file,
            "fixed": trb_data["fixed"],
            "designed": trb_data["designed"],
            "source_fixed": trb_data["source_fixed"],
            "plddt_mean": trb_data["plddt_mean"],
            "status": status,
        })

    df = pd.DataFrame(designs, columns=[
        "id", "pdb", "fixed", "designed", "source_fixed", "plddt_mean", "status"])
    os.makedirs(os.path.dirname(args.table_path), exist_ok=True)
    df.to_csv(args.table_path, index=False)


if __name__ == "__main__":
    main()
