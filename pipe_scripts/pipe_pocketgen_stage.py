#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
PocketGen input staging.

PocketGen's generate_new.py expects each input pair to live as
    {target}/{name}/{name}.pdb
    {target}/{name}/{name}_ligand.sdf

PocketGen's preprocessing selects the "edit pocket" residues by proximity to
the ligand atoms in 3D space (utils/protein_ligand.PDBProtein.query_residues_ligand).
A free-floating ligand (RDKit-embedded SMILES) returns 0-1 pocket residues
and crashes in collate_mols_block with a 2D/3D tensor mismatch.

To get a docked ligand SDF we mirror ProLIF's pattern (pipe_prolif.py): for
each (scaffold, ligand) pair, extract the HETATM block matching the ligand's
3-letter `code` from the scaffold PDB (keeping the original bound
coordinates), then use the ligand stream's `smiles` column as a bond-order
template via AssignBondOrdersFromTemplate before writing the SDF.

The ligand stream must carry `smiles`; the `code` comes either from
--ligand-code (when the wrapper was given an explicit code) or from the
ligand stream's `code` column (the ligand stream must then resolve to a
single distinct code). Either way, the scaffold PDB must contain HETATM
records for that code. We do NOT fall back to free-floating embedding —
that produces silent garbage out of PocketGen.
"""

import argparse
import os
import shutil
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from biopipelines.biopipelines_io import (  # noqa: E402
    load_datastream,
    iterate_files,
    iterate_values,
)
from biopipelines.ligand_utils import resolve_ligand_code, write_ligand_sdf  # noqa: E402
from biopipelines.pdb_parser import field_res_name, field_chain, field_res_seq  # noqa: E402


def load_template_smiles(ligand_json: str, ligand_code: str) -> str:
    """Read the SMILES for the ligand identified by `ligand_code` from the
    ligand stream's map_table. Errors if the row is missing or has no SMILES.
    Mirrors pipe_prolif.py:load_template_smiles."""
    ds = load_datastream(ligand_json)
    for _cid, values in iterate_values(ds, columns=["code", "smiles"]):
        if str(values.get("code", "") or "").strip().upper() == ligand_code.upper():
            smi = str(values.get("smiles", "") or "").strip()
            if not smi:
                raise ValueError(
                    f"ligand stream has no SMILES for code {ligand_code!r}; "
                    "supply it via Ligand(..., smiles=\"...\")"
                )
            return smi
    raise ValueError(f"ligand stream has no row with code {ligand_code!r}")


def extract_hetatm_block(scaffold_pdb: str, code: str) -> str:
    """Return a mini-PDB string with only HETATM records whose residue name
    (cols 18-20, 1-indexed) matches `code`. If multiple copies of the ligand
    are present, keep the first (same chain, same resSeq) — PocketGen expects
    a single ligand molecule. Raises if no matching HETATM records exist.
    """
    code_u = code.strip().upper()
    selected = []
    first_key = None  # (chain, resnum)
    with open(scaffold_pdb) as f:
        for line in f:
            if not line.startswith("HETATM"):
                continue
            rname = field_res_name(line).upper()
            if rname != code_u:
                continue
            chain = field_chain(line)
            resnum = field_res_seq(line)
            key = (chain, resnum)
            if first_key is None:
                first_key = key
            if key != first_key:
                continue
            selected.append(line.rstrip("\n"))
    if not selected:
        raise ValueError(
            f"No HETATM records with residue name {code_u!r} found in {scaffold_pdb}. "
            "PocketGen requires the ligand to be bound in the scaffold PDB."
        )
    return "\n".join(selected) + "\nEND\n"


def materialise_ligand_sdf(
    pair_id: str,
    scaffold_pdb: str,
    code: str,
    smiles: str,
    dst_sdf: str,
    extras_dir: str,
) -> None:
    """Carve the HETATM block for `code` out of `scaffold_pdb` (keeps docked
    coords), assign bond orders from `smiles`, write the result as an SDF
    PocketGen can consume.

    AssignBondOrdersFromTemplate matches by heavy-atom skeleton. Raw crystal
    PDBs strip hydrogens, so we must NOT call AddHs on the template — a
    protonated template (e.g. 59 atoms vs the PDB's 34 heavy) can't match.
    ProLIF's call uses AddHs because it runs on already-protonated PDBs (the
    Reduce tool runs upstream of it); for PocketGen the scaffold is raw, so
    use the bare SMILES template.

    If the bare-template match fails (e.g. crystal has missing atoms or
    different protonation), fall back to AddHs on both sides — sometimes
    works for ligands where the PDB *does* carry explicit Hs.
    """
    # Carve the bound HETATM (keeps docked coords); the bond-order SDF write is
    # the shared write_ligand_sdf (one implementation across all ligand tools).
    block = extract_hetatm_block(scaffold_pdb, code)
    os.makedirs(extras_dir, exist_ok=True)
    lig_pdb = os.path.join(extras_dir, f"{pair_id}_lig.pdb")
    with open(lig_pdb, "w") as f:
        f.write(block)

    try:
        write_ligand_sdf(lig_pdb, dst_sdf, smiles)
    except Exception as exc:
        raise RuntimeError(
            f"bond-order SDF write failed for {pair_id} (code={code!r}): {exc}. "
            f"Check that the SMILES matches the bound ligand's skeleton; extracted "
            f"PDB left at {lig_pdb} for inspection."
        )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--structures-json", required=True)
    ap.add_argument("--ligand-json", required=True)
    ap.add_argument("--staging-folder", required=True)
    ap.add_argument(
        "--extras-dir",
        required=True,
        help="Tool's _extras/ folder. Per-pair extracted-ligand PDBs land here "
             "so they survive failed staging for inspection.",
    )
    args = ap.parse_args()

    ligand_code = resolve_ligand_code(args.ligand_json)
    template_smiles = load_template_smiles(args.ligand_json, ligand_code)
    print(f"Ligand code: {ligand_code}")
    print(f"Template SMILES: {template_smiles}")

    proteins = list(iterate_files(load_datastream(args.structures_json)))
    if not proteins:
        print("ERROR: no scaffolds to stage", file=sys.stderr)
        sys.exit(1)

    os.makedirs(args.staging_folder, exist_ok=True)

    # Iterate ligand stream rows so per-ligand IDs flow through to pair_id.
    ligand_ds = load_datastream(args.ligand_json)
    ligand_ids = list(ligand_ds.ids_expanded)
    if not ligand_ids:
        print("ERROR: ligand stream is empty", file=sys.stderr)
        sys.exit(1)

    os.makedirs(args.extras_dir, exist_ok=True)

    n_pairs = 0
    failed = []
    for prot_id, prot_path in proteins:
        for lig_id in ligand_ids:
            pair_id = f"{prot_id}+{lig_id}"
            pair_dir = os.path.join(args.staging_folder, pair_id)
            os.makedirs(pair_dir, exist_ok=True)

            dst_pdb = os.path.join(pair_dir, f"{pair_id}.pdb")
            dst_sdf = os.path.join(pair_dir, f"{pair_id}_ligand.sdf")

            shutil.copyfile(prot_path, dst_pdb)
            try:
                materialise_ligand_sdf(
                    pair_id=pair_id,
                    scaffold_pdb=dst_pdb,
                    code=ligand_code,
                    smiles=template_smiles,
                    dst_sdf=dst_sdf,
                    extras_dir=args.extras_dir,
                )
                print(f"  [stage] {pair_id}: extracted HETATM {ligand_code} from scaffold + bond-order template")
                n_pairs += 1
            except Exception as exc:
                print(f"  [stage] {pair_id}: FAILED to materialise docked SDF: {exc}", file=sys.stderr)
                failed.append((pair_id, str(exc)))
                # The extracted ligand PDB (if any) is still under
                # args.extras_dir for inspection. Leave the pair_dir alone —
                # PocketGen's driver picks up only folders that contain a
                # readable {name}_ligand.sdf, so partial staging is harmless
                # downstream.

    print(f"Staged {n_pairs} (scaffold, ligand) pair(s) under {args.staging_folder}")
    if failed:
        print(f"FAILED {len(failed)} pair(s):", file=sys.stderr)
        for pair_id, msg in failed:
            print(f"  - {pair_id}: {msg}", file=sys.stderr)
        if n_pairs == 0:
            sys.exit(1)


if __name__ == "__main__":
    main()
