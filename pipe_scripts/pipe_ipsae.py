# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Execution half of the IPSAE tool: run upstream ipsae.py per id and parse its score table.

Upstream writes its outputs next to the structure it is given and derives file names by replacing substrings of the whole path (".pdb", ".cif", "pae"), so each id is staged as model.<ext> and pae_model.<ext> inside its own folder and the script runs there. It also exits 0 after printing a message when an input is unusable, so success is judged by the score table it wrote, not by the exit code.
"""

import argparse
import hashlib
import json
import math
import os
import shutil
import subprocess
import sys

import pandas as pd
import yaml

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from biopipelines.biopipelines_io import (  # noqa: E402
    load_datastream,
    iterate_files,
    resolve_file,
    step_id_from_table_path,
)
from biopipelines.ipsae import PAIR_COLUMNS, SCORE_COLUMNS, SUMMARY_COLUMNS  # noqa: E402

NAN = float("nan")

# Header of the upstream score table -> output column. PAE, Dist and Model are constants of the run and are not carried.
UPSTREAM_COLUMNS = {
    "Chn1": "chain1", "Chn2": "chain2", "Type": "type",
    "ipSAE": "ipsae", "ipSAE_d0chn": "ipsae_d0chn", "ipSAE_d0dom": "ipsae_d0dom",
    "ipTM_af": "iptm_af", "ipTM_d0chn": "iptm_d0chn",
    "pDockQ": "pdockq", "pDockQ2": "pdockq2", "LIS": "lis",
    "n0res": "n0res", "n0chn": "n0chn", "n0dom": "n0dom",
    "d0res": "d0res", "d0chn": "d0chn", "d0dom": "d0dom",
    "nres1": "nres1", "nres2": "nres2", "dist1": "dist1", "dist2": "dist2",
}
TEXT_COLUMNS = ("chain1", "chain2", "type")
INT_COLUMNS = ("n0res", "n0chn", "n0dom", "nres1", "nres2", "dist1", "dist2")


def parse_ipsae_output(text):
    """Upstream's <stem>_<pae>_<dist>.txt -> one dict per row, columns renamed per UPSTREAM_COLUMNS."""
    lines = text.splitlines()
    start = next((i for i, line in enumerate(lines) if line.split()[:2] == ["Chn1", "Chn2"]), None)
    if start is None:
        raise ValueError("no score table in the upstream output")
    header = lines[start].split()
    absent = [name for name in UPSTREAM_COLUMNS if name not in header]
    if absent:
        raise ValueError(f"upstream score table has no column(s) {absent}; its header changed")
    rows = []
    for line in lines[start + 1:]:
        if not line.strip():
            continue
        fields = line.split()
        if len(fields) != len(header):
            raise ValueError(f"upstream score row has {len(fields)} fields, header has {len(header)}: {line.strip()[:80]}")
        raw = dict(zip(header, fields))
        row = {}
        for name, column in UPSTREAM_COLUMNS.items():
            if column in TEXT_COLUMNS:
                row[column] = raw[name]
            elif column in INT_COLUMNS:
                row[column] = int(raw[name])
            else:
                row[column] = float(raw[name])
        rows.append(row)
    return rows


class PaeNotFound(Exception):
    """An id with no PAE entry or no PAE file: a lookup that found no partner, recorded as kind="unmatched"."""


def find_pae(pae_ds, sid):
    """The PAE file of an id; PaeNotFound when the stream has no entry or no file for it."""
    if sid not in pae_ds.ids_expanded:
        raise PaeNotFound(f"'{sid}' is not in the pae stream")
    try:
        path = resolve_file(pae_ds, sid)
    except FileNotFoundError as e:
        raise PaeNotFound(str(e)) from e
    if not os.path.exists(path):
        raise PaeNotFound(f"PAE file not found: {path}")
    return path


def structure_extension(path):
    """The extension upstream dispatches on: .pdb or .cif."""
    ext = os.path.splitext(path)[1].lower()
    if ext == ".pdb":
        return ".pdb"
    if ext in (".cif", ".mmcif"):
        return ".cif"
    raise ValueError(f"structure {os.path.basename(path)} is neither PDB nor mmCIF")


def pae_extension(path):
    ext = os.path.splitext(path)[1].lower()
    if ext not in (".npz", ".json"):
        raise ValueError(f"PAE file {os.path.basename(path)} is neither npz nor json")
    return ext


def confidence_available(pae_path, structure_ext):
    """(has_plddt, has_iptm): which upstream inputs the staged PAE file actually carries.

    npz: upstream looks for plddt_*.npz and confidence_*.json beside it, and none is staged.
    json with a PDB: ColabFold layout, plddt and iptm are keys of the file itself.
    json with a CIF: AlphaFold3 layout, plddt is atom_plddts and ipTM would come from a summary file that is not staged.
    """
    if pae_extension(pae_path) == ".npz":
        return False, False
    with open(pae_path) as fh:
        data = json.load(fh)
    if structure_ext == ".pdb":
        return "plddt" in data, "iptm" in data
    return "atom_plddts" in data, False


def mask_unavailable(rows, has_plddt, has_iptm):
    """NaN the columns whose upstream input was absent, which upstream prints as a number computed from zeros."""
    for row in rows:
        if not has_plddt:
            row["pdockq"] = NAN
            row["pdockq2"] = NAN
        if not has_iptm:
            row["iptm_af"] = NAN
    return rows


def summarize(rows, binder_chain, target_chains):
    """The max-type row of binder_chain against the target chain with the highest ipSAE."""
    chains = []
    for row in rows:
        for chain in (row["chain1"], row["chain2"]):
            if chain not in chains:
                chains.append(chain)
    if not chains:
        raise ValueError(f"upstream scored no chain pair, so binder_chain {binder_chain!r} has no partner: the structure holds fewer than two chains")
    if binder_chain not in chains:
        raise ValueError(f"binder_chain {binder_chain!r} is not a chain of the structure (chains: {', '.join(chains)})")
    targets = list(target_chains) if target_chains else [c for c in sorted(chains) if c != binder_chain]
    if not targets:
        raise ValueError(f"the structure has no chain other than binder_chain {binder_chain!r}")
    absent = [c for c in targets if c not in chains]
    if absent:
        raise ValueError(f"target chain(s) {absent} not in the structure (chains: {', '.join(chains)})")

    candidates = []
    for target in targets:
        match = [r for r in rows if r["type"] == "max" and {r["chain1"], r["chain2"]} == {binder_chain, target}]
        if not match:
            raise ValueError(f"upstream printed no max row for chains {binder_chain} and {target}")
        candidates.append((target, match[0]))
    target, best = max(candidates, key=lambda c: -math.inf if math.isnan(c[1]["ipsae"]) else c[1]["ipsae"])
    summary = {"binder_chain": binder_chain, "target_chain": target}
    summary.update({column: best[column] for column in SCORE_COLUMNS})
    return summary


def run_ipsae(script, run_dir, structure_path, pae_path, pae_cutoff, dist_cutoff):
    """Stage one id, run upstream in its folder and return the parsed score rows."""
    structure_ext = structure_extension(structure_path)
    pae_ext = pae_extension(pae_path)
    has_plddt, has_iptm = confidence_available(pae_path, structure_ext)

    shutil.rmtree(run_dir, ignore_errors=True)
    os.makedirs(run_dir)
    staged_structure = f"model{structure_ext}"
    staged_pae = f"pae_model{pae_ext}"
    shutil.copyfile(structure_path, os.path.join(run_dir, staged_structure))
    shutil.copyfile(pae_path, os.path.join(run_dir, staged_pae))

    proc = subprocess.run(
        [sys.executable, os.path.abspath(script), staged_pae, staged_structure, str(pae_cutoff), str(dist_cutoff)],
        cwd=run_dir, capture_output=True, text=True)
    os.remove(os.path.join(run_dir, staged_structure))
    os.remove(os.path.join(run_dir, staged_pae))

    said = (proc.stderr.strip() or proc.stdout.strip()).splitlines()
    tail = said[-1] if said else "no output"
    if proc.returncode != 0:
        raise RuntimeError(f"ipsae.py exited {proc.returncode}: {tail}")
    tables = [f for f in os.listdir(run_dir)
              if f.startswith("model_") and f.endswith(".txt") and not f.endswith("_byres.txt")]
    if len(tables) != 1:
        raise RuntimeError(f"ipsae.py wrote no score table: {tail}")
    with open(os.path.join(run_dir, tables[0])) as fh:
        text = fh.read()
    try:
        rows = parse_ipsae_output(text)
    except ValueError as e:
        raise ValueError(f"{e}; ipsae.py said: {tail}")
    return mask_unavailable(rows, has_plddt, has_iptm)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--config-yaml", required=True)
    p.add_argument("--structures-json", required=True)
    p.add_argument("--pae-json", required=True)
    p.add_argument("--ipsae-script", required=True)
    p.add_argument("--runs-dir", required=True)
    p.add_argument("--ipsae-csv", required=True)
    p.add_argument("--pairs-csv", required=True)
    p.add_argument("--local-missing-csv", required=True)
    args = p.parse_args()

    with open(args.config_yaml) as fh:
        config = yaml.safe_load(fh) or {}
    binder_chain = config["binder_chain"]
    target_chains = config.get("target_chains")
    pae_cutoff = config["pae_cutoff"]
    dist_cutoff = config["dist_cutoff"]

    if not os.path.isfile(args.ipsae_script):
        print(f"ERROR: {args.ipsae_script} not found; run IPSAE.install()", file=sys.stderr)
        sys.exit(1)
    with open(args.ipsae_script, "rb") as fh:
        print(f"ipsae.py sha256 {hashlib.sha256(fh.read()).hexdigest()}")

    structures_ds = load_datastream(args.structures_json)
    pae_ds = load_datastream(args.pae_json)
    step_id = step_id_from_table_path(args.local_missing_csv)

    pair_rows, summary_rows, local_missing_rows = [], [], []
    attempted = failures = 0

    for sid, structure_path in iterate_files(structures_ds):
        attempted += 1
        try:
            pae_path = find_pae(pae_ds, sid)
            rows = run_ipsae(args.ipsae_script, os.path.join(args.runs_dir, sid),
                             structure_path, pae_path, pae_cutoff, dist_cutoff)
            summary = summarize(rows, binder_chain, target_chains)
            pair_rows.extend({"id": sid, "pae_cutoff": float(pae_cutoff), "dist_cutoff": float(dist_cutoff), **r}
                             for r in rows)
            summary_rows.append({"id": sid, **summary})
            print(f"{sid}: ipsae {summary['ipsae']:.4f} ({binder_chain} vs {summary['target_chain']})")
        except Exception as e:
            print(f"WARNING: {sid} failed: {e}", file=sys.stderr)
            failures += 1
            summary_rows.append({"id": sid, "binder_chain": binder_chain})
            kind = "unmatched" if isinstance(e, PaeNotFound) else "failure"
            local_missing_rows.append({"id": sid, "removed_by": step_id,
                                       "kind": kind, "cause": str(e)[:200]})

    for path in (args.ipsae_csv, args.pairs_csv, args.local_missing_csv):
        os.makedirs(os.path.dirname(path), exist_ok=True)
    pd.DataFrame(summary_rows, columns=SUMMARY_COLUMNS).to_csv(args.ipsae_csv, index=False)
    pd.DataFrame(pair_rows, columns=PAIR_COLUMNS).to_csv(args.pairs_csv, index=False)
    pd.DataFrame(local_missing_rows,
                 columns=["id", "removed_by", "kind", "cause"]).to_csv(args.local_missing_csv, index=False)

    print(f"ipsae: {args.ipsae_csv} ({len(summary_rows)} rows, {failures} failed)")
    print(f"ipsae_pairs: {args.pairs_csv} ({len(pair_rows)} rows)")

    if not attempted:
        print("ERROR: no input entities to process (upstream stream is empty or its files are absent)",
              file=sys.stderr)
        sys.exit(1)
    if failures == attempted:
        print(f"ERROR: all {attempted} attempted id(s) failed", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
