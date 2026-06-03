#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
ESMFold inference script.

Reads a sequences CSV (columns: id, sequence) produced by the Sequence tool,
runs ESMFold structure prediction, and writes:
  - one <id>.pdb per sequence in --output-dir
  - ESMFold_scores.json with per-sequence pLDDT and pTM scores

This script runs under the 'esmfold' conda environment.

Usage:
    python pipe_esmfold_inference.py \\
        --sequences-csv <path/to/sequences.csv> \\
        --output-dir <output_dir> \\
        [--chunk-size 124] \\
        [--num-recycles 4] \\
        [--max-tokens-per-batch 1024] \\
        [--cpu-only] \\
        [--cpu-offload]
"""

import argparse
import json
import logging
import os
import sys
import typing as T
from pathlib import Path
from timeit import default_timer as timer

import numpy as np
import pandas as pd
import torch
import esm

logger = logging.getLogger()
logger.setLevel(logging.INFO)

formatter = logging.Formatter(
    "%(asctime)s | %(levelname)s | %(name)s | %(message)s",
    datefmt="%y/%m/%d %H:%M:%S",
)
console_handler = logging.StreamHandler(sys.stdout)
console_handler.setLevel(logging.INFO)
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)


def read_sequences_csv(csv_path: str) -> T.List[T.Tuple[str, str]]:
    """Read (id, sequence) pairs from a CSV file produced by the Sequence tool."""
    df = pd.read_csv(csv_path)
    if "id" not in df.columns or "sequence" not in df.columns:
        raise ValueError(
            f"sequences CSV must have 'id' and 'sequence' columns. Found: {list(df.columns)}"
        )
    return list(zip(df["id"].astype(str), df["sequence"].astype(str)))


def perres_plddt_from_pdb(pdb_string: str) -> T.List[float]:
    """Read per-residue pLDDT from the B-factor column of a PDB string."""
    lines = [x.split() for x in pdb_string.splitlines() if x.startswith("ATOM")]
    peratm = [(x[4] + x[5].zfill(8), x[-2]) for x in lines]
    reslist = list(set(r for r, _ in peratm))
    d: T.Dict[str, T.List[float]] = {r: [] for r in reslist}
    for res, plddt in peratm:
        d[res].append(float(plddt))
    return [round(float(np.mean(v)), 2) for _, v in sorted(d.items())]


def enable_cpu_offloading(model):
    import socket
    from torch.distributed.fsdp import CPUOffload, FullyShardedDataParallel
    from torch.distributed.fsdp.wrap import enable_wrap, wrap

    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.bind(("", 0))
        free_port = s.getsockname()[1]

    torch.distributed.init_process_group(
        backend="nccl", init_method=f"tcp://localhost:{free_port}", world_size=1, rank=0
    )
    wrapper_kwargs = dict(cpu_offload=CPUOffload(offload_params=True))
    with enable_wrap(wrapper_cls=FullyShardedDataParallel, **wrapper_kwargs):
        for layer_name, layer in model.layers.named_children():
            setattr(model.layers, layer_name, wrap(layer))
        model = wrap(model)
    return model


def init_model_on_gpu_with_cpu_offloading(model):
    model = model.eval()
    model_esm = enable_cpu_offloading(model.esm)
    del model.esm
    model.cuda()
    model.esm = model_esm
    return model


def create_batched_sequence_dataset(
    sequences: T.List[T.Tuple[str, str]], max_tokens_per_batch: int = 1024
) -> T.Generator[T.Tuple[T.List[str], T.List[str]], None, None]:
    batch_headers, batch_sequences, num_tokens = [], [], 0
    for header, seq in sequences:
        if (len(seq) + num_tokens > max_tokens_per_batch) and num_tokens > 0:
            yield batch_headers, batch_sequences
            batch_headers, batch_sequences, num_tokens = [], [], 0
        batch_headers.append(header)
        batch_sequences.append(seq)
        num_tokens += len(seq)
    yield batch_headers, batch_sequences


def main():
    parser = argparse.ArgumentParser(description="ESMFold inference from sequences CSV")
    parser.add_argument("--sequences-csv", required=True, help="Input sequences CSV (id, sequence columns)")
    parser.add_argument("--output-dir", required=True, help="Output directory for PDB files")
    parser.add_argument("--chunk-size", type=int, default=None,
                        help="Axial attention chunk size to reduce memory (e.g. 128, 64, 32)")
    parser.add_argument("--num-recycles", type=int, default=None,
                        help="Number of recycles (default: 4, as used in training)")
    parser.add_argument("--max-tokens-per-batch", type=int, default=1024,
                        help="Maximum tokens per GPU forward pass for batched prediction")
    parser.add_argument("--cpu-only", action="store_true", help="Run on CPU only")
    parser.add_argument("--cpu-offload", action="store_true", help="Enable CPU offloading")
    args = parser.parse_args()

    if not os.path.exists(args.sequences_csv):
        logger.error(f"Input sequences CSV not found: {args.sequences_csv}")
        sys.exit(1)

    os.makedirs(args.output_dir, exist_ok=True)

    logger.info(f"Reading sequences from {args.sequences_csv}")
    all_sequences = sorted(read_sequences_csv(args.sequences_csv), key=lambda hs: len(hs[1]))
    logger.info(f"Loaded {len(all_sequences)} sequences")

    logger.info("Loading ESMFold model")
    model = esm.pretrained.esmfold_v1()
    model = model.eval()
    model.set_chunk_size(args.chunk_size)

    if args.cpu_only:
        model.cpu()
    elif args.cpu_offload:
        model = init_model_on_gpu_with_cpu_offloading(model)
    else:
        model.cuda()

    logger.info("Starting predictions")
    batches = create_batched_sequence_dataset(all_sequences, args.max_tokens_per_batch)

    num_completed = 0
    num_sequences = len(all_sequences)
    out_dict: T.Dict[str, T.List] = {
        "description": [],
        "plddt": [],
        "ptm": [],
        "perresidue_plddt": [],
    }

    for headers, sequences in batches:
        start = timer()
        try:
            output = model.infer(sequences, num_recycles=args.num_recycles)
        except RuntimeError as e:
            if e.args[0].startswith("CUDA out of memory"):
                if len(sequences) > 1:
                    logger.warning(
                        f"CUDA OOM on batch of {len(sequences)}. "
                        "Try lowering --max-tokens-per-batch."
                    )
                else:
                    logger.warning(f"CUDA OOM on sequence {headers[0]} (length {len(sequences[0])})")
                continue
            raise

        output = {k: v.cpu() for k, v in output.items()}
        pdbs = model.output_to_pdb(output)
        elapsed = timer() - start
        time_str = f"{elapsed / len(headers):.1f}s"
        if len(sequences) > 1:
            time_str += f" (amortized, batch {len(sequences)})"

        for header, seq, pdb_string, mean_plddt, ptm in zip(
            headers, sequences, pdbs, output["mean_plddt"], output["ptm"]
        ):
            out_path = os.path.join(args.output_dir, f"{header}.pdb")
            with open(out_path, "w", encoding="UTF-8") as f:
                f.write(pdb_string)
            num_completed += 1
            logger.info(
                f"Predicted {header} | length {len(seq)} | pLDDT {float(mean_plddt):.1f} | "
                f"pTM {float(ptm):.3f} | {time_str} | {num_completed}/{num_sequences}"
            )
            out_dict["description"].append(header)
            out_dict["plddt"].append(float(mean_plddt))
            out_dict["ptm"].append(float(ptm))
            out_dict["perresidue_plddt"].append(perres_plddt_from_pdb(pdb_string))

    scores_path = os.path.join(args.output_dir, "ESMFold_scores.json")
    with open(scores_path, "w") as f:
        json.dump(out_dict, f)
    logger.info(f"Scores written to {scores_path}")

    if num_completed == 0:
        logger.error("No structures predicted successfully")
        sys.exit(1)

    if num_completed < num_sequences:
        logger.warning(f"{num_sequences - num_completed}/{num_sequences} sequences failed (likely CUDA OOM)")

    logger.info(f"ESMFold complete: {num_completed}/{num_sequences} structures predicted")


if __name__ == "__main__":
    main()
