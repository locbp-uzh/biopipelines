#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""AF2BIND runner. For each input structure, builds a ColabDesign `binder`
model on the target chain, runs one AlphaFold2 pass, applies the trained
AF2BIND linear head to the pair representation, and writes per-residue
binding probabilities plus a top-k binding-residue selection.

Mirrors sokrypton/af2bind's af2bind.ipynb: 20 bait residues are appended as a
"binder", the target's pair-representation slices against the baits are
concatenated, normalized with the head's stored mean/std, and projected
through the head's linear weights to a per-residue p_bind.

Weights layout:
    --af2-params-dir/params/params_model_*.npz   AlphaFold2 params (shared cache;
                                                 data_dir passed to colabdesign)
    --head-dir/attempt_7_2k_lam0-03/<model_type>.pickle   AF2BIND linear head
"""

import argparse
import os
import pickle
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from biopipelines.biopipelines_io import load_datastream, iterate_files, step_id_from_table_path  # noqa: E402
from biopipelines.sele_utils import chain_aware_sele  # noqa: E402

# 20-residue bait alphabet (AF2BIND's fixed binder sequence order).
BAIT_SEQ = "ACDEFGHIKLMNPQRSTVWY"
AA_3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}
# colabdesign/AlphaFold integer aatype order (alphabetical 3-letter -> 1-letter).
AATYPE_ORDER = "ARNDCQEGHILKMFPSTWYVX"

BINDING_COLS = ["id", "chain", "resi", "resn", "p_bind"]
SUMMARY_COLS = ["id", "n_residues", "top_resi", "top_p_bind", "binding_residues"]


def sigmoid(x):
    return 1.0 / (1.0 + np.exp(-x))


def load_head(head_dir: str, mask_sidechains: bool, seed: int = 0):
    """Load the AF2BIND linear head params for the chosen masking mode."""
    if mask_sidechains:
        model_type = f"split_nosc_pair_A_split_nosc_pair_B_{seed}"
    else:
        model_type = f"split_pair_A_split_pair_B_{seed}"
    pkl = os.path.join(head_dir, "attempt_7_2k_lam0-03", f"{model_type}.pickle")
    with open(pkl, "rb") as fh:
        params = pickle.load(fh)
    params = dict(**params["~"], **params["linear"])
    return {k: np.asarray(v) for k, v in params.items()}


def af2bind_head(pair, head):
    """Apply the AF2BIND linear head to an AF2 pair representation.

    `pair` is the full [L+20, L+20, C] pair tensor; the last 20 rows/cols are
    the bait residues. Returns p_bind over the L target residues.
    """
    pair_A = pair[:-20, -20:]
    pair_B = pair[-20:, :-20].swapaxes(0, 1)
    pair_A = pair_A.reshape(pair_A.shape[0], -1)
    pair_B = pair_B.reshape(pair_B.shape[0], -1)
    x = np.concatenate([pair_A, pair_B], -1)
    x = (x - head["mean"]) / head["std"]
    x = (x * head["w"][:, 0]) + (head["b"] / x.shape[-1])
    p_bind_aa = x.reshape(x.shape[0], 2, 20, -1).sum((1, 3))
    return sigmoid(p_bind_aa.sum(-1))


def predict_structure(af_model, mk_model_kwargs, pdb_path, chain,
                      mask_sidechains, mask_sequence):
    """Build a fresh binder model for one PDB and return (pair, idx_chain,
    idx_resi, idx_aa). Reusing one model across PDBs is unsafe because
    prep_inputs reshapes the model's static input arrays per target."""
    from colabdesign import mk_afdesign_model  # noqa: F401  (import checked by caller)
    m = mk_afdesign_model(protocol="binder", debug=True, **mk_model_kwargs)
    m.prep_inputs(pdb_filename=pdb_path, chain=chain, binder_len=20,
                  rm_target_sc=mask_sidechains, rm_target_seq=mask_sequence)
    # Spread the 20 bait residue indices far from the target (notebook trick).
    r_idx = m._inputs["residue_index"][-20] + (1 + np.arange(20)) * 50
    m._inputs["residue_index"][-20:] = r_idx.flatten()
    m.set_seq(BAIT_SEQ)
    m.predict(verbose=False)
    pair = np.asarray(m.aux["debug"]["outputs"]["representations"]["pair"])
    chains = m._pdb["idx"]["chain"]
    resis = m._pdb["idx"]["residue"]
    aatype = m._pdb["batch"]["aatype"]
    return pair, chains, resis, aatype


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--structures-json", required=True)
    p.add_argument("--af2-params-dir", required=True,
                   help="dir whose params/ subdir holds params_model_*.npz "
                        "(colabdesign data_dir)")
    p.add_argument("--head-dir", required=True,
                   help="dir holding attempt_7_2k_lam0-03/<model_type>.pickle")
    p.add_argument("--chain", default="A")
    p.add_argument("--top-k", type=int, default=15)
    p.add_argument("--binding-csv", required=True)
    p.add_argument("--binding-dir", required=True)
    p.add_argument("--binding-map-csv", required=True)
    p.add_argument("--summary-csv", required=True)
    p.add_argument("--missing-csv", required=True)
    p.add_argument("--mask-sidechains", action="store_true")
    p.add_argument("--mask-sequence", action="store_true")
    p.add_argument("--upstream-missing", default=None)
    args = p.parse_args()

    # Fail loudly and early if the ML stack is unavailable.
    try:
        import colabdesign  # noqa: F401
    except ImportError as e:
        print(f"ERROR: colabdesign not importable: {e}", file=sys.stderr)
        sys.exit(1)

    head = load_head(args.head_dir, args.mask_sidechains)
    # colabdesign reads AF2 params from <data_dir>/params/.
    mk_model_kwargs = {"data_dir": args.af2_params_dir}

    os.makedirs(args.binding_dir, exist_ok=True)
    ds = load_datastream(args.structures_json)
    binding_rows, summary_rows, missing_rows, binding_map_rows = [], [], [], []
    step_id = step_id_from_table_path(args.missing_csv)

    for sid, pdb_path in iterate_files(ds):
        try:
            pair, chains, resis, aatype = predict_structure(
                None, mk_model_kwargs, pdb_path, args.chain,
                args.mask_sidechains, args.mask_sequence,
            )
            p_bind = af2bind_head(pair, head)
        except Exception as e:
            print(f"WARNING: {sid} AF2BIND failed: {e}", file=sys.stderr)
            missing_rows.append({"id": sid, "removed_by": step_id, "kind": "failure", "cause": str(e)[:200]})
            continue

        target_len = len(p_bind)
        per_res = []
        per_id_rows = []
        for i in range(target_len):
            c = chains[i]
            r = int(resis[i])
            resn = AATYPE_ORDER[int(aatype[i])] if int(aatype[i]) < len(AATYPE_ORDER) else "X"
            pb = float(p_bind[i])
            row = {"id": sid, "chain": c, "resi": r, "resn": resn, "p_bind": pb}
            binding_rows.append(row)
            per_id_rows.append(row)
            per_res.append((c, r, pb))

        if not per_res:
            missing_rows.append({"id": sid, "removed_by": step_id, "kind": "failure", "cause": "no target residues scored"})
            continue

        # One resi-csv file per structure for the `binding` stream.
        bind_path = os.path.join(args.binding_dir, f"{sid}.csv")
        pd.DataFrame(per_id_rows, columns=BINDING_COLS).to_csv(bind_path, index=False)
        binding_map_rows.append({"id": sid, "file": bind_path})

        ranked = sorted(per_res, key=lambda t: t[2], reverse=True)
        top = ranked[: args.top_k]
        top_residues = chain_aware_sele([(c, r) for c, r, _ in top])
        best_c, best_r, best_p = ranked[0]
        summary_rows.append({
            "id": sid,
            "n_residues": target_len,
            "top_resi": f"{best_c}{best_r}",
            "top_p_bind": round(best_p, 4),
            "binding_residues": top_residues,
        })
        print(f"  {sid}: {target_len} residues, top {best_c}{best_r} p_bind={best_p:.3f}")

    upstream_rows = []
    if args.upstream_missing and os.path.exists(args.upstream_missing):
        try:
            up_df = pd.read_csv(args.upstream_missing)
            if not up_df.empty:
                upstream_rows = up_df.to_dict("records")
        except Exception as e:
            print(f"Warning: could not read upstream missing.csv: {e}", file=sys.stderr)

    all_missing = upstream_rows + missing_rows

    for d in (args.binding_csv, args.binding_map_csv, args.summary_csv, args.missing_csv):
        os.makedirs(os.path.dirname(d), exist_ok=True)
    pd.DataFrame(binding_rows, columns=BINDING_COLS).to_csv(args.binding_csv, index=False)
    pd.DataFrame(binding_map_rows, columns=["id", "file"]).to_csv(args.binding_map_csv, index=False)
    pd.DataFrame(summary_rows, columns=SUMMARY_COLS).to_csv(args.summary_csv, index=False)
    pd.DataFrame(all_missing, columns=["id", "removed_by", "kind", "cause"]).to_csv(args.missing_csv, index=False)
    print(f"Binding: {args.binding_csv} ({len(binding_rows)} rows)")
    print(f"Binding resi-csv files: {args.binding_map_csv} ({len(binding_map_rows)} rows)")
    print(f"Summary: {args.summary_csv} ({len(summary_rows)} rows)")
    print(f"Missing: {args.missing_csv} ({len(all_missing)} rows)")

    if missing_rows:
        print(f"Failed: {len(missing_rows)}/{len(missing_rows)+len(summary_rows)}", file=sys.stderr)
    if not summary_rows:
        sys.exit(1)


if __name__ == "__main__":
    main()
