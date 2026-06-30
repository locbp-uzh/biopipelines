#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""APBS runner. Per PDB:
  1. pdb2pqr --ff=AMBER --titration-state-method=propka --with-ph=<ph>
  2. apbs <input.in>
  3. Parse the resulting PQR for net charge / acidic / basic residue counts.
  4. Parse the DX grid for mean potential at protein surface points.
"""

import argparse
import os
import re
import subprocess
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from biopipelines.biopipelines_io import load_datastream, iterate_files, step_id_from_table_path, container_argv_prefix  # noqa: E402
from biopipelines.pdb_parser import field_res_name, field_chain, field_res_seq  # noqa: E402


E_COLS = ["id", "net_charge", "n_basic", "n_acidic", "pI", "mean_potential"]
BASIC = {"ARG", "LYS", "HIS"}
ACIDIC = {"ASP", "GLU"}


_CONTAINER_PREFIX = ""


def run(cmd, cwd=None):
    full = container_argv_prefix(_CONTAINER_PREFIX) + cmd
    res = subprocess.run(full, capture_output=True, text=True, cwd=cwd)
    if res.returncode != 0:
        raise RuntimeError(f"{cmd[0]} failed: {res.stderr.strip() or res.stdout.strip()}")
    return res


def parse_pqr(pqr_path: str):
    """Read a PQR file. Return (net_charge, n_basic, n_acidic, residue_resnames)."""
    seen = set()
    res_names = []
    net = 0.0
    with open(pqr_path) as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            try:
                charge = float(line[55:62])
            except (ValueError, IndexError):
                continue
            net += charge
            resname = field_res_name(line)
            chain = field_chain(line)
            resnum = field_res_seq(line)
            key = (chain, resnum)
            if key not in seen:
                seen.add(key)
                res_names.append(resname)
    n_basic = sum(1 for r in res_names if r in BASIC)
    n_acidic = sum(1 for r in res_names if r in ACIDIC)
    return net, n_basic, n_acidic


def parse_pdb2pqr_summary(stderr: str):
    """Extract isoelectric-point estimate from pdb2pqr's PROPKA output if present."""
    m = re.search(r"isoelectric\s*point\s*[:=]\s*([\d.]+)", stderr, re.IGNORECASE)
    if m:
        try:
            return float(m.group(1))
        except ValueError:
            pass
    return None


def parse_dx_mean(dx_path: str) -> float:
    """Read an APBS DX grid file and return the mean of all potential values."""
    vals = []
    in_data = False
    with open(dx_path) as f:
        for line in f:
            if not in_data:
                if line.startswith("object 3"):
                    in_data = True
                continue
            if line.startswith("attribute") or line.startswith("object") or not line.strip():
                break
            parts = line.split()
            for p in parts:
                try:
                    vals.append(float(p))
                except ValueError:
                    break
    return float(np.mean(vals)) if vals else 0.0


APBS_TEMPLATE = """read
    mol pqr {pqr}
end
elec
    mg-auto
    dime {dim} {dim} {dim}
    cglen 60.0 60.0 60.0
    fglen 50.0 50.0 50.0
    cgcent mol 1
    fgcent mol 1
    mol 1
    {solver}
    bcfl sdh
    ion charge +1 conc {conc} radius 2.0
    ion charge -1 conc {conc} radius 1.8
    pdie {pdie}
    sdie {sdie}
    srfm smol
    chgm spl2
    sdens 10.0
    srad 1.4
    swin 0.3
    temp 298.15
    calcenergy total
    calcforce no
    write pot dx {dx_base}
end
quit
"""


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--structures-json", required=True)
    p.add_argument("--pqr-dir", required=True)
    p.add_argument("--pqr-map-csv", required=True)
    p.add_argument("--grid-dir", required=True)
    p.add_argument("--grid-map-csv", required=True)
    p.add_argument("--scratch-dir", required=True)
    p.add_argument("--ph", type=float, default=7.4)
    p.add_argument("--forcefield", default="AMBER")
    p.add_argument("--ion-concentration", type=float, default=0.150)
    p.add_argument("--grid-dim", type=int, default=65)
    p.add_argument("--pdie", type=float, default=2.0)
    p.add_argument("--sdie", type=float, default=78.5)
    p.add_argument("--solver", default="lpbe")
    p.add_argument("--electrostatics-csv", required=True)
    p.add_argument("--missing-csv", required=True)
    p.add_argument("--container-prefix", default="")
    p.add_argument("--upstream-missing", default=None)
    args = p.parse_args()
    global _CONTAINER_PREFIX
    _CONTAINER_PREFIX = args.container_prefix

    os.makedirs(args.pqr_dir, exist_ok=True)
    os.makedirs(args.grid_dir, exist_ok=True)
    os.makedirs(args.scratch_dir, exist_ok=True)
    ds = load_datastream(args.structures_json)

    rows, pqr_rows, grid_rows, missing_rows = [], [], [], []
    step_id = step_id_from_table_path(args.missing_csv)
    for sid, pdb_path in iterate_files(ds):
        work = os.path.join(args.scratch_dir, sid)
        os.makedirs(work, exist_ok=True)
        pqr = os.path.join(args.pqr_dir, f"{sid}.pqr")
        in_file = os.path.join(work, f"{sid}.in")
        dx_path = os.path.join(args.grid_dir, f"{sid}.dx")
        dx_base = os.path.splitext(dx_path)[0]
        try:
            pqr_cmd = ["pdb2pqr", f"--ff={args.forcefield}", "--titration-state-method=propka",
                       f"--with-ph={args.ph}", pdb_path, pqr]
            res = run(pqr_cmd)
            pI = parse_pdb2pqr_summary((res.stderr or "") + (res.stdout or ""))

            with open(in_file, "w") as f:
                f.write(APBS_TEMPLATE.format(
                    pqr=pqr, dx_base=dx_base,
                    dim=args.grid_dim, solver=args.solver,
                    conc=args.ion_concentration, pdie=args.pdie, sdie=args.sdie,
                ))

            run(["apbs", in_file], cwd=work)
            if not os.path.exists(dx_path):
                for cand in (os.path.join(work, "pot.dx"),
                             os.path.join(work, "pot-PE0.dx")):
                    if os.path.exists(cand):
                        os.rename(cand, dx_path)
                        break
        except Exception as e:
            print(f"WARNING: {sid} APBS failed: {e}", file=sys.stderr)
            missing_rows.append({"id": sid, "removed_by": step_id, "kind": "failure", "cause": str(e)[:200]})
            continue

        if not os.path.exists(pqr):
            missing_rows.append({"id": sid, "removed_by": step_id, "kind": "failure", "cause": "pdb2pqr produced no PQR"})
            continue
        if not os.path.exists(dx_path):
            missing_rows.append({"id": sid, "removed_by": step_id, "kind": "failure", "cause": "apbs produced no DX grid"})
            continue

        net, n_basic, n_acidic = parse_pqr(pqr)
        mean_phi = parse_dx_mean(dx_path)
        rows.append({
            "id": sid,
            "net_charge": round(net, 2),
            "n_basic": n_basic,
            "n_acidic": n_acidic,
            "pI": round(pI, 2) if pI is not None else "",
            "mean_potential": round(mean_phi, 4),
        })
        pqr_rows.append({"id": sid, "file": pqr})
        grid_rows.append({"id": sid, "file": dx_path})
        print(f"  {sid}: net_charge={net:.2f}  basic={n_basic}  acidic={n_acidic}")

    upstream_rows = []
    if args.upstream_missing and os.path.exists(args.upstream_missing):
        try:
            up_df = pd.read_csv(args.upstream_missing)
            if not up_df.empty:
                upstream_rows = up_df.to_dict("records")
        except Exception as e:
            print(f"Warning: could not read upstream missing.csv: {e}", file=sys.stderr)

    all_missing = upstream_rows + missing_rows

    for d in (args.electrostatics_csv, args.pqr_map_csv, args.grid_map_csv, args.missing_csv):
        os.makedirs(os.path.dirname(d), exist_ok=True)
    pd.DataFrame(rows, columns=E_COLS).to_csv(args.electrostatics_csv, index=False)
    pd.DataFrame(pqr_rows, columns=["id", "file"]).to_csv(args.pqr_map_csv, index=False)
    pd.DataFrame(grid_rows, columns=["id", "file"]).to_csv(args.grid_map_csv, index=False)
    pd.DataFrame(all_missing, columns=["id", "removed_by", "kind", "cause"]).to_csv(args.missing_csv, index=False)
    print(f"Electrostatics: {args.electrostatics_csv} ({len(rows)} rows)")
    print(f"PQR map: {args.pqr_map_csv} ({len(pqr_rows)} rows)")
    print(f"Grid map: {args.grid_map_csv} ({len(grid_rows)} rows)")
    print(f"Missing: {args.missing_csv} ({len(all_missing)} rows)")

    if missing_rows:
        print(f"Failed: {len(missing_rows)}/{len(missing_rows)+len(rows)}", file=sys.stderr)
    if not rows:
        sys.exit(1)


if __name__ == "__main__":
    main()
