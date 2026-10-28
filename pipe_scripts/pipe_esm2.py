#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Execution half of the ESM2 tool: pseudo-log-likelihood and mean-pooled embedding per sequence.

The numerical core (`masked_log_likelihoods`, `mean_embedding`) works on numpy arrays and takes the model as a callable, so it does not import torch; `Esm2Model` is the thin torch adapter around `EsmForMaskedLM`. Every id that enters gets a row in `scores` (NaN when it failed) and, when it succeeded, a `<id>.npy` and a map row; failures also go to the local missing file.
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd
import yaml

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from biopipelines.biopipelines_io import (  # noqa: E402
    load_datastream,
    iterate_values,
    step_id_from_table_path,
)
from biopipelines.esm2 import RESIDUE_ALPHABET, SCORE_COLUMNS  # noqa: E402

NAN = float("nan")
MISSING_COLUMNS = ["id", "removed_by", "kind", "cause"]


def resolve_device(requested, cuda_available):
    """Resolve auto to cuda when present; an explicit cuda without a GPU is an error, never a silent CPU run."""
    if requested == "auto":
        return "cuda" if cuda_available else "cpu"
    if requested == "cuda" and not cuda_available:
        raise RuntimeError("device 'cuda' was requested but no CUDA device is available")
    if requested not in ("cpu", "cuda"):
        raise ValueError(f"device must be auto, cpu or cuda, got {requested!r}")
    return requested


def check_weights_folder(folder):
    """The pinned folder written by ESM2.install(): config, tokenizer vocabulary and weights."""
    hint = "run ESM2.install() for this model"
    if not os.path.isdir(folder):
        raise FileNotFoundError(f"weights folder {folder} does not exist; {hint}")
    names = set(os.listdir(folder))
    for required in ("config.json", "vocab.txt"):
        if required not in names:
            raise FileNotFoundError(f"{required} not found in {folder}; {hint}")
    if not any(n.endswith((".safetensors", ".bin")) for n in names):
        raise FileNotFoundError(f"no weights file found in {folder}; {hint}")


def validate_sequence(sequence):
    """The sequence unchanged if every letter is in the ESM-2 alphabet; otherwise a ValueError naming the others."""
    if not isinstance(sequence, str) or not sequence:
        raise ValueError(f"empty or non-text sequence: {sequence!r}")
    unknown = sorted(set(sequence) - set(RESIDUE_ALPHABET))
    if unknown:
        raise ValueError(f"letters outside the ESM-2 alphabet {RESIDUE_ALPHABET}: {''.join(unknown)}")
    return sequence


def log_softmax(logits):
    """Row-wise log-softmax in float64."""
    logits = np.asarray(logits, dtype=np.float64)
    peak = logits.max(axis=-1, keepdims=True)
    return logits - peak - np.log(np.exp(logits - peak).sum(axis=-1, keepdims=True))


def masked_log_likelihoods(token_ids, mask_id, logits_fn, batch_size):
    """log p(x_i | all other tokens) for every residue i, each position masked on its own.

    token_ids is [<cls>, x_1 .. x_L, <eos>]. logits_fn maps an int array (B, T) to logits (B, T, V). The L masked copies go through logits_fn batch_size at a time; the log-softmax runs over the whole vocabulary at the masked position and is read at the true token.
    """
    token_ids = np.asarray(token_ids)
    length = len(token_ids) - 2
    if length < 1:
        raise ValueError("token_ids must hold at least one residue between <cls> and <eos>")
    out = np.empty(length, dtype=np.float64)
    for start in range(0, length, batch_size):
        positions = np.arange(start + 1, min(start + batch_size, length) + 1)
        rows = np.arange(len(positions))
        batch = np.tile(token_ids, (len(positions), 1))
        batch[rows, positions] = mask_id
        logits = np.asarray(logits_fn(batch))
        if logits.shape[:2] != batch.shape:
            raise ValueError(f"logits_fn returned shape {logits.shape} for a batch of shape {batch.shape}")
        out[start:start + len(positions)] = log_softmax(logits[rows, positions])[rows, token_ids[positions]]
    return out


def pseudo_log_likelihood(per_residue):
    """(sum, mean) of the per-residue log-likelihoods."""
    return float(np.sum(per_residue)), float(np.mean(per_residue))


def mean_embedding(hidden):
    """Mean over the residue positions of a (T, D) hidden state, <cls> and <eos> excluded, as float32."""
    hidden = np.asarray(hidden)
    if hidden.ndim != 2 or hidden.shape[0] < 3:
        raise ValueError(f"hidden state must be (tokens, dim) with at least one residue, got shape {hidden.shape}")
    return hidden[1:-1].mean(axis=0).astype(np.float32)


class Esm2Model:
    """EsmForMaskedLM and its tokenizer loaded from a pinned local folder, in fp32 and eval mode."""

    def __init__(self, folder, device):
        import torch
        from transformers import AutoTokenizer, EsmForMaskedLM

        self.torch = torch
        self.device = resolve_device(device, torch.cuda.is_available())
        self.tokenizer = AutoTokenizer.from_pretrained(folder, local_files_only=True)
        self.model = EsmForMaskedLM.from_pretrained(folder, local_files_only=True).to(self.device).eval()
        self.mask_id = self.tokenizer.mask_token_id

    def encode(self, sequence):
        """[<cls>, residues, <eos>] as int64; one token per residue is a precondition of the position bookkeeping."""
        ids = np.asarray(self.tokenizer(sequence)["input_ids"], dtype=np.int64)
        if len(ids) != len(sequence) + 2:
            raise ValueError(f"tokenizer produced {len(ids)} tokens for {len(sequence)} residues (expected {len(sequence) + 2})")
        return ids

    def logits(self, batch):
        with self.torch.inference_mode():
            ids = self.torch.as_tensor(batch, dtype=self.torch.long, device=self.device)
            return self.model(input_ids=ids).logits.float().cpu().numpy()

    def hidden(self, token_ids):
        """Last hidden state (T, D) after the encoder's final layer norm."""
        with self.torch.inference_mode():
            ids = self.torch.as_tensor(np.asarray(token_ids)[None], dtype=self.torch.long, device=self.device)
            return self.model.esm(input_ids=ids).last_hidden_state[0].float().cpu().numpy()


def read_revision(folder):
    path = os.path.join(folder, "REVISION")
    return open(path).read().strip() if os.path.isfile(path) else "unrecorded"


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--config-yaml", required=True)
    p.add_argument("--sequences-json", required=True)
    p.add_argument("--weights-dir", required=True)
    p.add_argument("--embeddings-dir", required=True)
    p.add_argument("--embeddings-map-csv", required=True)
    p.add_argument("--scores-csv", required=True)
    p.add_argument("--local-missing-csv", required=True)
    args = p.parse_args()

    with open(args.config_yaml) as fh:
        config = yaml.safe_load(fh) or {}
    compute_pll = config["compute_pll"]
    compute_embeddings = config["compute_embeddings"]
    batch_size = config["batch_size"]

    check_weights_folder(args.weights_dir)
    model = Esm2Model(args.weights_dir, config["device"])
    print(f"ESM-2 weights: {args.weights_dir} ({read_revision(args.weights_dir)}), device {model.device}")

    ds = load_datastream(args.sequences_json)
    step_id = step_id_from_table_path(args.local_missing_csv)
    os.makedirs(args.embeddings_dir, exist_ok=True)

    score_rows, map_rows, local_missing_rows = [], [], []
    attempted = failures = 0

    for sid, values in iterate_values(ds, columns=["sequence"]):
        attempted += 1
        try:
            sequence = validate_sequence(values["sequence"])
            token_ids = model.encode(sequence)
            row = {"id": sid, "length": len(sequence)}
            if compute_pll:
                pll, pll_mean = pseudo_log_likelihood(
                    masked_log_likelihoods(token_ids, model.mask_id, model.logits, batch_size))
                row.update(pll=pll, pll_per_residue=pll_mean)
            if compute_embeddings:
                vector = mean_embedding(model.hidden(token_ids))
                path = os.path.join(args.embeddings_dir, f"{sid}.npy")
                np.save(path, vector)
                map_rows.append({"id": sid, "file": path})
            score_rows.append(row)
            print(f"{sid}: {len(sequence)} residues" + (f", pll {row['pll']:.3f}" if compute_pll else ""))
        except Exception as e:
            print(f"WARNING: {sid} failed: {e}", file=sys.stderr)
            failures += 1
            score_rows.append({"id": sid})
            local_missing_rows.append({"id": sid, "removed_by": step_id, "kind": "failure", "cause": str(e)[:200]})

    for path in (args.scores_csv, args.embeddings_map_csv, args.local_missing_csv):
        os.makedirs(os.path.dirname(path), exist_ok=True)
    if compute_pll:
        scores = pd.DataFrame(score_rows, columns=SCORE_COLUMNS)
        scores["length"] = scores["length"].astype("Int64")
        scores.to_csv(args.scores_csv, index=False)
        print(f"scores: {args.scores_csv} ({len(scores)} rows, {failures} failed)")
    if compute_embeddings:
        pd.DataFrame(map_rows, columns=["id", "file"]).to_csv(args.embeddings_map_csv, index=False)
        print(f"embeddings: {args.embeddings_map_csv} ({len(map_rows)} rows)")
    pd.DataFrame(local_missing_rows, columns=MISSING_COLUMNS).to_csv(args.local_missing_csv, index=False)

    if not attempted:
        print("ERROR: no input entities to process (upstream stream is empty or its files are absent)", file=sys.stderr)
        sys.exit(1)
    if failures == attempted:
        print(f"ERROR: all {attempted} attempted id(s) failed", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
