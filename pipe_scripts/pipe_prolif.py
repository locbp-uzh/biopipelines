#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""ProLIF runner. For each complex PDB, computes a single-frame interaction
fingerprint between the named ligand residue and the protein, then exports it
as a long-form CSV (one row per residue x interaction-type).

Input PDBs must already carry explicit hydrogens. The ligand stream supplies
the canonical SMILES that's used as a bond-order template; without it,
RDKit's PDB reader produces a chemically broken ligand (spurious charges,
missing aromaticity) and ProLIF detects no interactions."""

import argparse
import os
import re
import sys

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from biopipelines.biopipelines_io import load_datastream, iterate_files, iterate_values, read_upstream_missing, MISSING_COLUMNS, step_id_from_table_path  # noqa: E402
from biopipelines.ligand_utils import resolve_ligand_code, templated_ligand_mol  # noqa: E402

import MDAnalysis as mda  # noqa: E402
import prolif as plf  # noqa: E402
from rdkit import Chem  # noqa: E402


FP_COLS = ["id", "residue", "resn", "chain", "resnum", "resi", "interaction_type", "present"]
RES_ID_RE = re.compile(r"^([A-Za-z0-9]+?)(\d+)(?:\.([A-Za-z0-9]+))?$")


def load_template_smiles(ligand_json: str, ligand_code: str) -> str:
    """Read the SMILES for the ligand identified by `ligand_code` from the
    compounds stream's map_table. Errors if the row is missing or has no SMILES."""
    ds = load_datastream(ligand_json)
    for _cid, values in iterate_values(ds, columns=["code", "smiles"]):
        if str(values.get("code", "")).strip().upper() == ligand_code.upper():
            smi = str(values.get("smiles", "")).strip()
            if not smi:
                raise ValueError(f"ligand stream has no SMILES for code {ligand_code}")
            return smi
    raise ValueError(f"ligand stream has no row with code {ligand_code}")


def parse_res_id(res_id) -> tuple:
    """Decompose a ProLIF ResidueId (e.g. 'ALA138.A') into (resname, resnum, chain)."""
    s = str(res_id)
    m = RES_ID_RE.match(s)
    if not m:
        return s, "", ""
    resname, resnum, chain = m.group(1), m.group(2), m.group(3) or ""
    return resname, resnum, chain


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--structures-json", required=True)
    p.add_argument("--ligand-json", required=True)
    p.add_argument("--extras-dir", required=True)
    p.add_argument("--fingerprints-csv", required=True)
    p.add_argument("--missing-csv", required=True)
    p.add_argument("--upstream-missing", nargs="*", default=None)
    args = p.parse_args()

    os.makedirs(args.extras_dir, exist_ok=True)
    ligand_code = resolve_ligand_code(args.ligand_json)
    template_smiles = load_template_smiles(args.ligand_json, ligand_code)
    if Chem.MolFromSmiles(template_smiles) is None:
        raise RuntimeError(f"RDKit failed to parse template SMILES for {ligand_code}: {template_smiles}")

    ds = load_datastream(args.structures_json)
    rows, failed = [], []
    step_id = step_id_from_table_path(args.missing_csv)
    fp = plf.Fingerprint()

    for sid, pdb_path in iterate_files(ds):
        try:
            u = mda.Universe(pdb_path, guess_bonds=True)
            lig_ag = u.select_atoms(f"resname {ligand_code}")
            prot_ag = u.select_atoms("protein")
            if not lig_ag.atoms or not prot_ag.atoms:
                print(f"WARNING: {sid} missing protein or ligand atoms (code={ligand_code})", file=sys.stderr)
                failed.append({"id": sid, "removed_by": step_id, "kind": "failure",
                               "cause": f"missing protein or ligand atoms (code={ligand_code})"})
                continue

            lig_pdb = os.path.join(args.extras_dir, f"{sid}_lig.pdb")
            lig_ag.write(lig_pdb)
            # Shared bond-order restoration (PDB coords + SMILES template),
            # returning the mol object ProLIF needs (no intermediate SDF).
            lig_with_bonds = templated_ligand_mol(lig_pdb, template_smiles)
            lig_mol = plf.Molecule(lig_with_bonds)
            prot_mol = plf.Molecule.from_mda(prot_ag)

            fp.run_from_iterable([lig_mol], prot_mol, progress=False)
            df = fp.to_dataframe()
            df = df.iloc[0]
            for (_lig_id, res_id, interaction), present in df.items():
                resname, resnum, chain = parse_res_id(res_id)
                rows.append({
                    "id": sid,
                    "residue": resname,
                    "resn": resname,
                    "chain": chain,
                    "resnum": resnum,
                    "resi": resnum,
                    "interaction_type": interaction,
                    "present": int(bool(present)),
                })
            print(f"  {sid}: {len(df)} residue x interaction entries")
        except Exception as e:
            print(f"WARNING: {sid} ProLIF failed: {e}", file=sys.stderr)
            failed.append({"id": sid, "removed_by": step_id, "kind": "failure", "cause": str(e)[:200]})

    all_missing = read_upstream_missing(args.upstream_missing) + failed

    for d in (args.fingerprints_csv, args.missing_csv):
        os.makedirs(os.path.dirname(d), exist_ok=True)
    pd.DataFrame(rows, columns=FP_COLS).to_csv(args.fingerprints_csv, index=False)
    pd.DataFrame(all_missing, columns=MISSING_COLUMNS).to_csv(args.missing_csv, index=False)
    print(f"Fingerprints: {args.fingerprints_csv} ({len(rows)} rows)")
    print(f"Missing: {args.missing_csv} ({len(all_missing)} rows)")

    if failed:
        print(f"Failed: {len(failed)}/{len(failed)+len(set(r['id'] for r in rows))}", file=sys.stderr)
    if not rows:
        sys.exit(1)


if __name__ == "__main__":
    main()
