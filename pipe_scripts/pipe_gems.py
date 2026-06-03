#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""GEMS runner. Stages a directory of <pair_id>.pdb + <pair_id>.sdf files
from two DataStreams (cross product), then runs the upstream
GEMS_dataprep_workflow.py and inference.py against that directory. Reads
the resulting predictions CSV and writes our standard affinity table."""

import argparse
import os
import shutil
import subprocess
import sys

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from biopipelines.biopipelines_io import load_datastream, iterate_files, iterate_values, read_upstream_missing, step_id_from_table_path  # noqa: E402
# Shared coordinate-ligand -> sanitized SDF (bond-order template) — single
# implementation, also used by OpenBabel and the other ligand-consuming tools.
from biopipelines.ligand_utils import write_ligand_sdf  # noqa: E402


AFF_COLS = ["id", "structures.id", "ligands.id", "pkd_pred"]


def load_smiles_lookup(smiles_json):
    """Build {ligand_id: smiles} from the optional smiles stream."""
    if not smiles_json:
        return {}
    ds = load_datastream(smiles_json)
    out = {}
    for cid, values in iterate_values(ds, columns=["smiles"]):
        smi = str(values.get("smiles", "")).strip()
        if smi:
            out[cid] = smi
    return out


def stage_pairs(prot_pairs, lig_pairs, work_dir, smiles_lookup):
    """Stage each (protein, ligand) pair into work_dir as pair_id.pdb / pair_id.sdf.

    The staged PDB is protein-only (HETATM lines stripped). GEMS's graph builder
    (dataprep/graph_construction.py) tries to bond the SDF ligand to the nearest
    protein residue and silently SkipComplexException's on "Ligand connected to
    unknown protein residue" when the protein PDB still carries a bound HETATM
    (the typical Boltz/co-crystal case). PDBbind convention is `_protein.pdb` =
    protein-only and `_ligand.sdf` = ligand; we honor that here.
    """
    os.makedirs(work_dir, exist_ok=True)
    pair_meta = []
    for prot_id, prot_pdb in prot_pairs:
        for lig_id, lig_file in lig_pairs:
            pair_id = f"{prot_id}__{lig_id}"
            dst_pdb = os.path.join(work_dir, f"{pair_id}.pdb")
            dst_sdf = os.path.join(work_dir, f"{pair_id}.sdf")
            with open(prot_pdb) as src, open(dst_pdb, "w") as dst:
                for line in src:
                    if not line.startswith("HETATM"):
                        dst.write(line)
            write_ligand_sdf(lig_file, dst_sdf, smiles_lookup.get(lig_id))
            pair_meta.append({"pair_id": pair_id, "structures.id": prot_id, "ligands.id": lig_id})
    return pair_meta


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--structures-json", required=True)
    p.add_argument("--compounds-json", required=True)
    p.add_argument("--smiles-json", default="")
    p.add_argument("--gems-repo", required=True)
    p.add_argument("--scratch-dir", required=True)
    p.add_argument("--affinity-csv", required=True)
    p.add_argument("--missing-csv", required=True)
    p.add_argument("--skip-ligand-embedding", action="store_true")
    p.add_argument("--upstream-missing", nargs="*", default=None)
    args = p.parse_args()

    os.makedirs(args.scratch_dir, exist_ok=True)

    prot_ds = load_datastream(args.structures_json)
    lig_ds = load_datastream(args.compounds_json)
    prot_pairs = [(sid, p) for sid, p in iterate_files(prot_ds)]
    lig_pairs = [(cid, l) for cid, l in iterate_files(lig_ds)]
    if not prot_pairs or not lig_pairs:
        raise SystemExit("Empty structures or compounds stream")

    smiles_lookup = load_smiles_lookup(args.smiles_json)

    data_dir = os.path.join(args.scratch_dir, "pairs")
    if os.path.isdir(data_dir):
        shutil.rmtree(data_dir)
    pair_meta = stage_pairs(prot_pairs, lig_pairs, data_dir, smiles_lookup)
    print(f"Staged {len(pair_meta)} pair(s) under {data_dir}")

    # GEMS's dataprep/inference import matplotlib.pyplot. On Colab the inherited
    # MPLBACKEND=module://matplotlib_inline.backend_inline is invalid in this env
    # and crashes the import. Force a headless backend for the subprocesses.
    #
    # PATH: this pipe script already runs under the gems-env interpreter (env
    # activation got it here), but a non-interactive SLURM step does not have
    # that env's bin first on PATH. GEMS_dataprep_workflow.py re-spawns its
    # sub-steps as bare `python -m dataprep.*`, which would otherwise resolve to
    # /usr/bin/python (no ankh). Prepend the gems-env bin so every descendant's
    # bare `python` resolves to the same interpreter.
    env_bin = os.path.dirname(sys.executable)
    sub_env = dict(os.environ,
                   MPLBACKEND="Agg",
                   PATH=env_bin + os.pathsep + os.environ.get("PATH", ""))

    # Step 1: dataprep. Must run with cwd = gems_repo so its relative imports work.
    res = subprocess.run(
        [sys.executable, "GEMS_dataprep_workflow.py", "--data_dir", data_dir],
        cwd=args.gems_repo, capture_output=True, text=True, env=sub_env,
    )
    if res.returncode != 0:
        raise RuntimeError(f"GEMS dataprep failed (exit {res.returncode}): "
                           f"{res.stderr.strip()[:600]}")
    print(res.stdout.strip()[-400:])

    # dataprep writes <basename>_dataset.pt in cwd (gems_repo).
    data_basename = os.path.basename(os.path.normpath(data_dir))
    dataset_pt = os.path.join(args.gems_repo, f"{data_basename}_dataset.pt")
    if not os.path.isfile(dataset_pt):
        raise RuntimeError(f"GEMS dataprep produced no dataset file at {dataset_pt}")

    # Step 2: inference.
    inf_cmd = [sys.executable, "inference.py", "--dataset_path", dataset_pt]
    if args.skip_ligand_embedding:
        inf_cmd += ["--skip_ligand_embedding", "True"]
    res = subprocess.run(inf_cmd, cwd=args.gems_repo, capture_output=True, text=True, env=sub_env)
    if res.returncode != 0:
        raise RuntimeError(f"GEMS inference failed (exit {res.returncode}): "
                           f"{res.stderr.strip()[-1500:]}")
    print(res.stdout.strip()[-400:])

    pred_csv = os.path.join(args.gems_repo, f"{data_basename}_dataset_predictions.csv")
    if not os.path.isfile(pred_csv):
        raise RuntimeError(f"GEMS inference produced no predictions at {pred_csv}")

    preds = pd.read_csv(pred_csv)
    pred_lookup = {str(row["id"]): float(row["y_pred"]) for _, row in preds.iterrows()}

    rows, missing = [], []
    step_id = step_id_from_table_path(args.missing_csv)
    for meta in pair_meta:
        pair_id = meta["pair_id"]
        if pair_id in pred_lookup:
            rows.append({
                "id": pair_id,
                "structures.id": meta["structures.id"],
                "ligands.id": meta["ligands.id"],
                "pkd_pred": round(pred_lookup[pair_id], 4),
            })
            print(f"  {pair_id}: pKd = {pred_lookup[pair_id]:.3f}")
        else:
            missing.append({"id": pair_id, "removed_by": step_id, "kind": "failure", "cause": "GEMS produced no prediction for this pair"})

    all_missing = read_upstream_missing(args.upstream_missing) + missing

    for d in (args.affinity_csv, args.missing_csv):
        os.makedirs(os.path.dirname(d), exist_ok=True)
    pd.DataFrame(rows, columns=AFF_COLS).to_csv(args.affinity_csv, index=False)
    pd.DataFrame(all_missing, columns=["id", "removed_by", "kind", "cause"]).to_csv(args.missing_csv, index=False)
    print(f"Affinity: {args.affinity_csv} ({len(rows)} rows)")
    print(f"Missing: {args.missing_csv} ({len(all_missing)} rows)")

    if not rows:
        sys.exit(1)


if __name__ == "__main__":
    main()
