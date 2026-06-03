#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""RTMScore runner. Iterates (protein, ligand) pairs from two DataStreams,
shells out to the upstream rtmscore.py script with --gen_pocket against the
ligand as the reference, and aggregates the per-pose scores into a single
table."""

import argparse
import os
import subprocess
import sys

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from biopipelines.biopipelines_io import load_datastream, iterate_files, iterate_values, read_upstream_missing, step_id_from_table_path  # noqa: E402
# Shared coordinate-ligand -> sanitized SDF (bond-order template) — single
# implementation, also used by OpenBabel and the other ligand-consuming tools.
from biopipelines.ligand_utils import write_ligand_sdf  # noqa: E402


SCORE_COLS = ["id", "structures.id", "ligands.id", "pose", "rtmscore"]


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


def prepare_protein_pdb(src_pdb, dst_pdb):
    """Write a protein-only PDB for RTMScore's pocket graph. The structures
    stream is a *complex* (the protein source), but upstream prot_to_graph pads
    every residue to RES_MAX_NATOMS atoms — a bound ligand or other HETATM, read
    as one oversized 'residue', overflows that pad ('negative dimensions are not
    allowed'). Keep only ATOM records (drop HETATM/waters/ligands)."""
    with open(src_pdb) as fh, open(dst_pdb, "w") as out:
        for line in fh:
            rec = line[:6]
            if rec in ("ATOM  ", "TER   ", "MODEL ", "ENDMDL") or line.startswith("END"):
                out.write(line)


def _find_babel_libdir():
    """Locate the OpenBabel plugin dir inside the active env.

    Resolve from sys.prefix (the running interpreter's env root) rather than
    $CONDA_PREFIX, which is not reliably exported under `mamba run`."""
    import glob
    candidates = []
    if os.environ.get("CONDA_PREFIX"):
        candidates.append(os.environ["CONDA_PREFIX"])
    candidates.append(sys.prefix)
    for prefix in candidates:
        matches = sorted(glob.glob(os.path.join(prefix, "lib", "openbabel", "*")))
        if matches:
            return matches[-1]
    return None


def run_rtmscore(rtm_script, model_path, prot_pdb, lig_sdf, cutoff, work_dir):
    """Invoke rtmscore.py --gen_pocket; return the path to its output CSV."""
    os.makedirs(work_dir, exist_ok=True)
    out_prefix = os.path.join(work_dir, "scores")
    cmd = [
        "python", rtm_script,
        "-p", prot_pdb,
        "-l", lig_sdf,
        "-rl", lig_sdf,
        "-gen_pocket",
        "-c", str(cutoff),
        "-m", model_path,
        "-o", out_prefix,
    ]
    # RTMScore's pocket extraction loads the OpenBabel C++ format plugins,
    # which need BABEL_LIBDIR pointing at THIS env's plugin dir. We override
    # unconditionally: on Colab the base-Python openbabel-wheel exports a
    # BABEL_LIBDIR for its own (py3.12) build, which would otherwise be
    # inherited here and make the rtmscore env's (py3.10) openbabel try to
    # dlopen ABI-incompatible plugins (undefined-symbol -> no PDB written ->
    # "graph of pocket cannot be generated").
    env = dict(os.environ)
    babel_libdir = _find_babel_libdir()
    if babel_libdir:
        env["BABEL_LIBDIR"] = babel_libdir
    # RTMScore imports MDAnalysis.analysis.dihedrals -> matplotlib.pyplot. On
    # Colab the inherited MPLBACKEND=module://matplotlib_inline.backend_inline is
    # invalid in this env and crashes the import. Force a headless backend.
    env["MPLBACKEND"] = "Agg"
    res = subprocess.run(cmd, capture_output=True, text=True, env=env)
    if res.returncode != 0:
        raise RuntimeError(f"rtmscore.py failed (exit {res.returncode}): "
                           f"{res.stderr.strip()[-1500:]}")
    out_csv = f"{out_prefix}.csv"
    if not os.path.isfile(out_csv):
        raise RuntimeError(f"rtmscore.py produced no output CSV at {out_csv}")
    return out_csv


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--structures-json", required=True)
    p.add_argument("--compounds-json", required=True)
    p.add_argument("--smiles-json", default="")
    p.add_argument("--rtmscore-script", required=True)
    p.add_argument("--model-path", required=True)
    p.add_argument("--cutoff", type=float, default=10.0)
    p.add_argument("--scratch-dir", required=True)
    p.add_argument("--scores-csv", required=True)
    p.add_argument("--missing-csv", required=True)
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

    rows, missing = [], []
    step_id = step_id_from_table_path(args.missing_csv)
    for prot_id, prot_pdb in prot_pairs:
        for lig_id, lig_file in lig_pairs:
            pair_id = f"{prot_id}+{lig_id}"
            work = os.path.join(args.scratch_dir, pair_id)
            try:
                # Stage a bond-order-correct SDF for the reference ligand and a
                # protein-only PDB for the pocket graph (the input is a complex).
                os.makedirs(work, exist_ok=True)
                lig_sdf = os.path.join(work, f"{lig_id}.sdf")
                write_ligand_sdf(lig_file, lig_sdf, smiles_lookup.get(lig_id))
                prot_clean = os.path.join(work, f"{prot_id}_protein.pdb")
                prepare_protein_pdb(prot_pdb, prot_clean)
                out_csv = run_rtmscore(args.rtmscore_script, args.model_path,
                                       prot_clean, lig_sdf, args.cutoff, work)
                df = pd.read_csv(out_csv)
                # rtmscore.py writes columns id,score where id is the
                # SDF molecule name (pose label).
                for _, r in df.iterrows():
                    rows.append({
                        "id": pair_id,
                        "structures.id": prot_id,
                        "ligands.id": lig_id,
                        "pose": str(r.get("id", "")),
                        "rtmscore": float(r["score"]),
                    })
                print(f"  {pair_id}: {len(df)} pose(s), best={df['score'].max():.3f}")
            except Exception as e:
                print(f"WARNING: {pair_id} RTMScore failed: {e}", file=sys.stderr)
                missing.append({"id": pair_id, "removed_by": step_id, "kind": "failure", "cause": str(e)[:200]})

    all_missing = read_upstream_missing(args.upstream_missing) + missing

    for d in (args.scores_csv, args.missing_csv):
        os.makedirs(os.path.dirname(d), exist_ok=True)
    pd.DataFrame(rows, columns=SCORE_COLS).to_csv(args.scores_csv, index=False)
    pd.DataFrame(all_missing, columns=["id", "removed_by", "kind", "cause"]).to_csv(args.missing_csv, index=False)
    print(f"Scores: {args.scores_csv} ({len(rows)} rows)")
    print(f"Missing: {args.missing_csv} ({len(all_missing)} rows)")

    if not rows:
        sys.exit(1)


if __name__ == "__main__":
    main()
