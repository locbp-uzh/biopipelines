#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""XTB runner. For each complex PDB, splits the structure into a ligand
fragment (resname == ligand_code) and the rest of the structure (protein),
runs three single-point xtb calculations, and writes per-complex interaction
energies E_complex - E_protein - E_ligand."""

import argparse
import os
import re
import shlex
import subprocess
import sys

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from biopipelines.biopipelines_io import load_datastream, iterate_files, read_upstream_missing, step_id_from_table_path  # noqa: E402
from biopipelines.ligand_utils import resolve_ligand_code  # noqa: E402

import MDAnalysis as mda  # noqa: E402


def _shell_join(cmd):
    return " ".join(shlex.quote(c) for c in cmd)


KJ_PER_HARTREE = 2625.4996  # Hartree -> kJ/mol
KCAL_PER_KJ = 1.0 / 4.184

E_COLS = ["id", "e_complex_kj", "e_protein_kj", "e_ligand_kj",
          "e_interaction_kj", "e_interaction_kcal", "charge_complex"]


def split_complex(pdb_path: str, ligand_code: str, out_dir: str) -> tuple:
    """Write protein-only and ligand-only PDB fragments, return their paths."""
    u = mda.Universe(pdb_path)
    lig_ag = u.select_atoms(f"resname {ligand_code}")
    prot_ag = u.select_atoms(f"not resname {ligand_code} and not resname HOH WAT")
    if not lig_ag.atoms:
        raise RuntimeError(f"no atoms found for ligand resname {ligand_code}")
    if not prot_ag.atoms:
        raise RuntimeError("protein fragment is empty after ligand removal")
    prot_pdb = os.path.join(out_dir, "protein.pdb")
    lig_pdb = os.path.join(out_dir, "ligand.pdb")
    prot_ag.write(prot_pdb)
    lig_ag.write(lig_pdb)
    return prot_pdb, lig_pdb


def run_xtb_sp(pdb_path: str, charge: int, method: str, solvent: str,
               opt: bool, scratch_dir: str) -> float:
    """Run an xtb single-point (or optimisation) and return total energy in kJ/mol."""
    method_flag = {"gfn2": "--gfn2", "gfn1": "--gfn1", "gfnff": "--gfnff"}[method]
    cmd = ["xtb", pdb_path, method_flag, "--chrg", str(charge)]
    if solvent:
        cmd += ["--alpb", solvent]
    if opt:
        cmd += ["--opt", "loose"]
    # xtb segfaults (exit -11) on large systems when the OpenMP thread stack
    # is too small; an unlimited heap stack and a generous per-thread stack
    # are the documented fix. Also bump the SCC iteration cap for proteins.
    env = dict(os.environ)
    env.setdefault("OMP_STACKSIZE", "4G")
    env.setdefault("OMP_MAX_ACTIVE_LEVELS", "1")
    cmd += ["--iterations", "500"]
    res = subprocess.run(
        f"ulimit -s unlimited 2>/dev/null; exec {_shell_join(cmd)}",
        shell=True, executable="/bin/bash",
        capture_output=True, text=True, cwd=scratch_dir, env=env,
    )
    if res.returncode != 0:
        raise RuntimeError(f"xtb failed (exit {res.returncode}): {res.stderr.strip()[:200]}")
    # Parse "TOTAL ENERGY  -123.45678 Eh"
    m = re.search(r"TOTAL ENERGY\s+(-?\d+\.\d+)\s+Eh", res.stdout)
    if not m:
        raise RuntimeError("could not parse TOTAL ENERGY from xtb output")
    return float(m.group(1)) * KJ_PER_HARTREE


def infer_ligand_charge(lig_pdb: str, scratch_dir: str) -> int:
    """Use rdkit to read the ligand and compute formal charge. Falls back to 0."""
    try:
        from rdkit import Chem
        mol = Chem.MolFromPDBFile(lig_pdb, sanitize=False, removeHs=False)
        if mol is None:
            return 0
        return int(sum(a.GetFormalCharge() for a in mol.GetAtoms()))
    except Exception:
        return 0


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--structures-json", required=True)
    p.add_argument("--ligand-json", required=True)
    p.add_argument("--method", default="gfn2", choices=["gfn2", "gfn1", "gfnff"])
    p.add_argument("--charge", type=int, default=0)
    p.add_argument("--solvent", default="")
    p.add_argument("--opt", action="store_true")
    p.add_argument("--scratch-dir", required=True)
    p.add_argument("--interaction-csv", required=True)
    p.add_argument("--missing-csv", required=True)
    p.add_argument("--upstream-missing", nargs="*", default=None)
    args = p.parse_args()

    os.makedirs(args.scratch_dir, exist_ok=True)

    ligand_code = resolve_ligand_code(args.ligand_json)

    ds = load_datastream(args.structures_json)
    rows, missing = [], []
    step_id = step_id_from_table_path(args.missing_csv)

    for sid, pdb_path in iterate_files(ds):
        work = os.path.join(args.scratch_dir, sid)
        os.makedirs(work, exist_ok=True)
        try:
            prot_pdb, lig_pdb = split_complex(pdb_path, ligand_code, work)
            lig_q = infer_ligand_charge(lig_pdb, work)
            prot_q = args.charge - lig_q

            e_lig = run_xtb_sp(lig_pdb, lig_q, args.method, args.solvent,
                               args.opt, _mkdir(work, "lig"))
            e_prot = run_xtb_sp(prot_pdb, prot_q, args.method, args.solvent,
                                args.opt, _mkdir(work, "prot"))
            e_cpx = run_xtb_sp(pdb_path, args.charge, args.method, args.solvent,
                               args.opt, _mkdir(work, "cpx"))
            e_int_kj = e_cpx - e_prot - e_lig
            rows.append({
                "id": sid,
                "e_complex_kj": round(e_cpx, 4),
                "e_protein_kj": round(e_prot, 4),
                "e_ligand_kj": round(e_lig, 4),
                "e_interaction_kj": round(e_int_kj, 4),
                "e_interaction_kcal": round(e_int_kj * KCAL_PER_KJ, 4),
                "charge_complex": args.charge,
            })
            print(f"  {sid}: E_int = {e_int_kj:.2f} kJ/mol "
                  f"({e_int_kj * KCAL_PER_KJ:.2f} kcal/mol)")
        except Exception as e:
            print(f"WARNING: {sid} xtb scoring failed: {e}", file=sys.stderr)
            missing.append({"id": sid, "removed_by": step_id, "kind": "failure", "cause": str(e)[:200]})

    all_missing = read_upstream_missing(args.upstream_missing) + missing

    for d in (args.interaction_csv, args.missing_csv):
        os.makedirs(os.path.dirname(d), exist_ok=True)
    pd.DataFrame(rows, columns=E_COLS).to_csv(args.interaction_csv, index=False)
    pd.DataFrame(all_missing, columns=["id", "removed_by", "kind", "cause"]).to_csv(args.missing_csv, index=False)
    print(f"Interactions: {args.interaction_csv} ({len(rows)} rows)")
    print(f"Missing: {args.missing_csv} ({len(all_missing)} rows)")

    if not rows:
        sys.exit(1)


def _mkdir(parent: str, name: str) -> str:
    d = os.path.join(parent, name)
    os.makedirs(d, exist_ok=True)
    return d


if __name__ == "__main__":
    main()
