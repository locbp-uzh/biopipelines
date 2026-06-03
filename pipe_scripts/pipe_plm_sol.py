#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""PLM_Sol runtime helper.

Reads the input sequences CSV (columns: id, sequence), then drives PLM_Sol's
two-stage flow inside the cloned repo:

  1. Generate per-residue ProtT5 embeddings directly with transformers
     (T5EncoderModel, Rostlab/prot_t5_xl_uniref50) -> .h5 keyed by id (L x 1024
     per protein) + a matching remapping FASTA. transformers is used directly
     because bio-embeddings pins torch<=1.10, conflicting with PLM_Sol's torch
     2.0.1. PLM_Sol's biLSTM/TextCNN needs the per-residue length dim.
  2. Run the repo's inference.py with a generated config pointing at those,
     which writes 'protTrans_prediction_result.csv' (cols protein_ID, sequence,
     predict_result; predict_result is the solubility probability, binary call
     is predict_result >= 0.5) into its working directory.

Then maps the result back to BioPipelines outputs:

  - solubility.csv : id | sequences.id | solubility | soluble
                     (sorted by solubility descending)
  - missing.csv    : sequences with no prediction (id | removed_by | kind | cause)

The prot_t5_xl_uniref50 weights download on first use (heavy; GPU strongly
recommended). Runs under the `plm_sol` conda env.
"""

import argparse
import os
import re
import subprocess
import sys

import pandas as pd
import yaml

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from biopipelines.biopipelines_io import load_datastream, step_id_from_table_path  # noqa: E402

SOLUBILITY_COLUMNS = ["id", "sequences.id", "solubility", "soluble"]
MISSING_COLUMNS = ["id", "removed_by", "kind", "cause"]

PROT_T5_MODEL = "Rostlab/prot_t5_xl_uniref50"


def read_sequences_csv(csv_path):
    """Read (id, sequence) pairs from the content-bearing sequences CSV."""
    df = pd.read_csv(csv_path)
    if "id" not in df.columns or "sequence" not in df.columns:
        raise ValueError(
            f"sequences CSV must have 'id' and 'sequence' columns. Found: {list(df.columns)}"
        )
    return [(str(i), str(s)) for i, s in zip(df["id"], df["sequence"])]


def remap_id(seq_id):
    """Match PLM_Sol's key_format='fasta_descriptor' id sanitisation: first
    whitespace token, '.'/'/' -> '_'."""
    return str(seq_id).split(" ")[0].replace(".", "_").replace("/", "_")


def run_embeddings(pairs, emb_dir):
    """Produce per-residue ProtT5 embeddings with transformers.

    Writes an .h5 keyed by the remapped id (L x 1024 per-residue array per
    protein) and a remapping FASTA with the same ids. PLM_Sol's biLSTM/TextCNN
    is a sequence model and needs the per-residue length dim (key_format=
    'fasta_descriptor', 'lm' mode). Returns (embeddings_h5, remapped_fasta).
    """
    import h5py
    import numpy as np
    import torch
    from transformers import T5EncoderModel, T5Tokenizer

    os.makedirs(emb_dir, exist_ok=True)
    embeddings_h5 = os.path.join(emb_dir, "reduced_embeddings_file.h5")
    remapped_fasta = os.path.join(emb_dir, "remapped_sequences.fasta")

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    tokenizer = T5Tokenizer.from_pretrained(PROT_T5_MODEL, do_lower_case=False)
    model = T5EncoderModel.from_pretrained(PROT_T5_MODEL).to(device).eval()

    with h5py.File(embeddings_h5, "w") as h5, open(remapped_fasta, "w") as fa:
        for seq_id, seq in pairs:
            rid = remap_id(seq_id)
            # ProtT5 expects space-separated residues; rare AAs -> X.
            spaced = " ".join(re.sub(r"[UZOB]", "X", seq.upper()))
            enc = tokenizer.batch_encode_plus(
                [spaced], add_special_tokens=True, padding="longest", return_tensors="pt"
            )
            input_ids = enc["input_ids"].to(device)
            attention_mask = enc["attention_mask"].to(device)
            with torch.no_grad():
                emb = model(input_ids=input_ids, attention_mask=attention_mask).last_hidden_state
            # Drop the trailing </s>; keep PER-RESIDUE (L x 1024): PLM_Sol's
            # biLSTM/TextCNN is a sequence model and needs the length dim
            # (solver.py keys on shape==2), not a mean-pooled vector.
            seq_len = int(attention_mask[0].sum()) - 1
            per_residue = emb[0, :seq_len].cpu().numpy().astype(np.float32)
            h5.create_dataset(rid, data=per_residue)
            fa.write(f">{rid}\n{seq}\n")

    return embeddings_h5, remapped_fasta


def run_inference(repo_dir, embeddings_h5, remapped_fasta, run_dir):
    """Run PLM_Sol inference.py from run_dir. Returns the prediction DataFrame
    (cols protein_ID, sequence, predict_result)."""
    template = os.path.join(repo_dir, "configs", "inference_Sol_biLSTM_TextCNN.yml")
    with open(template) as f:
        config = yaml.safe_load(f)

    # Repoint the hardcoded paths at our runtime artefacts. Checkpoint stays as
    # the repo's committed model_param; resolve it absolutely so cwd=run_dir
    # doesn't break it.
    config["embeddings"] = embeddings_h5
    config["remapping"] = remapped_fasta
    config["key_format"] = "fasta_descriptor"
    # 'lm' mode keys the h5 by protein id (what we write); 'profiles' keys by
    # sequence. Our reduced per-protein h5 is id-keyed, so force 'lm'.
    config["embedding_mode"] = "lm"
    ckpt = config.get("checkpoints_list") or [os.path.join("model_param", "model_param.t7")]
    config["checkpoints_list"] = [
        c if os.path.isabs(c) else os.path.join(repo_dir, c) for c in ckpt
    ]

    run_config = os.path.join(run_dir, "inference_config.yml")
    with open(run_config, "w") as f:
        yaml.safe_dump(config, f, sort_keys=False)

    # inference.py also reads './model_param/train_arguments.yml' relative to cwd,
    # so expose the repo's model_param/ inside run_dir.
    mp_link = os.path.join(run_dir, "model_param")
    if not os.path.exists(mp_link):
        os.symlink(os.path.join(repo_dir, "model_param"), mp_link)

    # inference.py writes 'protTrans_prediction_result.csv' into its cwd, and
    # imports model code relative to the repo, so put the repo on PYTHONPATH and
    # run from run_dir to capture the output there.
    env = dict(os.environ)
    env["PYTHONPATH"] = repo_dir + os.pathsep + env.get("PYTHONPATH", "")
    subprocess.run(
        [sys.executable, os.path.join(repo_dir, "inference.py"), "--config", run_config],
        check=True, cwd=run_dir, env=env,
    )

    result_csv = os.path.join(run_dir, "protTrans_prediction_result.csv")
    if not os.path.exists(result_csv):
        raise RuntimeError(f"inference.py did not produce {result_csv}")
    return pd.read_csv(result_csv)


def main():
    parser = argparse.ArgumentParser(description="PLM_Sol solubility scorer")
    parser.add_argument("--sequences-json", required=True)
    parser.add_argument("--sequences-csv", required=True)
    parser.add_argument("--repo-dir", required=True)
    parser.add_argument("--work-dir", required=True)
    parser.add_argument("--solubility-csv", required=True)
    parser.add_argument("--missing-csv", required=True)
    args = parser.parse_args()

    # load_datastream validates the input; the content-bearing CSV is the
    # source of (id, sequence) pairs (same contract as VespaG/ESMFold).
    load_datastream(args.sequences_json)
    pairs = read_sequences_csv(args.sequences_csv)
    if not pairs:
        print("Error: no sequences to score", file=sys.stderr)
        sys.exit(1)

    os.makedirs(args.work_dir, exist_ok=True)
    emb_dir = os.path.join(args.work_dir, "embeddings")

    embeddings_h5, remapped_fasta = run_embeddings(pairs, emb_dir)
    pred_df = run_inference(args.repo_dir, embeddings_h5, remapped_fasta, args.work_dir)

    # Join predictions back to original ids by sequence content (the remapped
    # ids are sanitised; the inference CSV preserves the original sequence).
    seq_to_id = {}
    for sid, seq in pairs:
        seq_to_id.setdefault(seq, sid)

    rows = []
    scored_ids = set()
    for _, r in pred_df.iterrows():
        seq = str(r["sequence"])
        prob = float(r["predict_result"])
        sid = seq_to_id.get(seq)
        if sid is None:
            print(f"WARNING: prediction sequence not in inputs (skipped): {seq[:30]}...",
                  file=sys.stderr)
            continue
        scored_ids.add(sid)
        rows.append({
            "id": sid,
            "sequences.id": sid,
            "solubility": prob,
            "soluble": int(prob >= 0.5),
        })

    out_df = pd.DataFrame(rows, columns=SOLUBILITY_COLUMNS)
    if not out_df.empty:
        out_df = out_df.sort_values("solubility", ascending=False)
    out_df.to_csv(args.solubility_csv, index=False)

    step_id = step_id_from_table_path(args.missing_csv)
    missing_rows = [
        {"id": sid, "removed_by": step_id, "kind": "no_prediction",
         "cause": "no solubility prediction returned"}
        for sid, _ in pairs if sid not in scored_ids
    ]
    pd.DataFrame(missing_rows, columns=MISSING_COLUMNS).to_csv(args.missing_csv, index=False)

    print(f"PLM_Sol scored {len(out_df)}/{len(pairs)} sequences")
    if missing_rows:
        print(f"Missing {len(missing_rows)}: {[r['id'] for r in missing_rows]}", file=sys.stderr)
    if out_df.empty:
        print("Error: no sequences produced a solubility prediction", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
