#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Execution half of the Embedding tool: one vector per input, and the cosine similarity of each requested pair.

The registry in biopipelines/embedding.py says which adapter embeds a model; an adapter turns a weights file into an object with `embed(path) -> (vector, nodes)`. The numerical core (`cosine_similarity`, `pair_rows`) works on numpy arrays and does not import torch; `ProgresGraphEmbedder` is the torch adapter around the `progres` package. Every input gets a `<id>.npy` and a map row when it embeds, and a local missing row when it does not; every pair row gets a similarities row, NaN when a member has no embedding.
"""

import argparse
import hashlib
import os
import sys
from typing import Dict, Set

import numpy as np
import pandas as pd
import yaml

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from biopipelines.biopipelines_io import (  # noqa: E402
    load_datastream,
    iterate_files,
    step_id_from_table_path,
)
from biopipelines.embedding import SIMILARITY_COLUMNS, registry_entry  # noqa: E402
from biopipelines.pair_alignment import pair_id, read_pairs  # noqa: E402

NAN = float("nan")
MISSING_COLUMNS = ["id", "removed_by", "kind", "cause"]
# Progres's tau torsion feature is padded for 3 or more residues; with 2 the node features and coordinates differ in length and torch raises an opaque size error.
MIN_RESIDUES = 3


def resolve_device(requested, cuda_available):
    """Resolve auto to cuda when present; an explicit cuda without a GPU is an error, never a silent CPU run."""
    if requested == "auto":
        return "cuda" if cuda_available else "cpu"
    if requested == "cuda" and not cuda_available:
        raise RuntimeError("device 'cuda' was requested but no CUDA device is available")
    if requested not in ("cpu", "cuda"):
        raise ValueError(f"device must be auto, cpu or cuda, got {requested!r}")
    return requested


def weights_file(entry, environ):
    """The checkpoint path the step exported under the entry's variable; it must name an existing file."""
    path = environ.get(entry.weights_env, "")
    hint = f" Upstream: {entry.upstream_default}." if entry.upstream_default else ""
    if not path:
        raise FileNotFoundError(f"{entry.weights_env} is not set, so there are no weights for model {entry.name!r}.{hint}")
    if not os.path.isfile(path):
        raise FileNotFoundError(f"{entry.weights_env}={path} is not a file.{hint}")
    return path


def weights_digest(path):
    """First 12 hex digits of the SHA-256 of the weights file, printed so a log names the exact weights that ran."""
    sha = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            sha.update(block)
    return sha.hexdigest()[:12]


def cosine_similarity(a, b):
    """Cosine of two vectors in float64, clipped to [-1, 1] against rounding; a zero or non-finite vector is a ValueError."""
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    if a.ndim != 1 or a.shape != b.shape:
        raise ValueError(f"cosine needs two vectors of one length, got shapes {a.shape} and {b.shape}")
    norm_a, norm_b = np.linalg.norm(a), np.linalg.norm(b)
    if not (np.isfinite(norm_a) and np.isfinite(norm_b)) or norm_a == 0 or norm_b == 0:
        raise ValueError("cosine is undefined for a zero or non-finite embedding")
    return float(np.clip(np.dot(a, b) / (norm_a * norm_b), -1.0, 1.0))


def check_embedding(vector, dimension):
    """The vector as float32 if it has the registered length and only finite values; otherwise a ValueError."""
    vector = np.asarray(vector, dtype=np.float32)
    if vector.shape != (dimension,):
        raise ValueError(f"embedding has shape {vector.shape}, the registry declares ({dimension},)")
    if not np.isfinite(vector).all():
        raise ValueError("embedding holds non-finite values")
    return vector


def pair_rows(pairs, vectors: Dict[str, np.ndarray], declared: Set[str], attempted: Set[str]):
    """One similarities row per pair, and (pair id, kind, cause) for every pair that could not be compared.

    `pairs` come from `pair_alignment.read_pairs` (distinct, ids as text, "" for an empty cell). `declared` are the ids of the input stream, `attempted` those whose file was present and embedded or failed; a member without an embedding is named by which of the three it is. An empty member, one outside the stream and one dropped upstream are lookups with no partner (kind "unmatched"); a member that failed to embed, or vectors the cosine rejects, are a "failure".
    """
    rows, failures = [], []
    for first, second in pairs:
        pid = pair_id(first, second)
        kind = "failure"
        try:
            for member in (first, second):
                if not member:
                    kind = "unmatched"
                    raise ValueError("a pair has an empty id")
                if member not in vectors:
                    if member not in declared:
                        kind, why = "unmatched", "is not in the input stream"
                    elif member not in attempted:
                        kind, why = "unmatched", "was dropped upstream (no file)"
                    else:
                        why = "failed to embed"
                    raise ValueError(f"'{member}' {why}")
            similarity = cosine_similarity(vectors[first], vectors[second])
        except ValueError as e:
            rows.append({"id": pid, "first": first, "second": second, "cosine_similarity": NAN})
            failures.append((pid, kind, str(e)))
            continue
        rows.append({"id": pid, "first": first, "second": second, "cosine_similarity": similarity})
    return rows, failures


class ProgresGraphEmbedder:
    """Progres family: the C-alpha graph of the first chain through the EGNN network, one L2-normalized vector per structure.

    Built from the `progres` package's own `Model`, `read_coords` and `coords_to_graph` (what its `read_graph` and `progres embed` chain), with the checkpoint loaded from the path given instead of through `progres.load_trained_model`, which takes no path and first downloads the package's default data and weights. The checkpoint is read with weights_only=True: the Progres and CIRPIN files hold only a state dict and an epoch number. Per-node embeddings are captured with a forward hook on the entry's `node_module`, so the model is run once for both outputs.
    """

    def __init__(self, entry, weights, device, want_nodes):
        import torch
        import progres.progres as progres_module

        self.torch = torch
        self.read_coords = progres_module.read_coords
        self.coords_to_graph = progres_module.coords_to_graph
        self.reader_args = dict(entry.extra_args)
        self.device = resolve_device(device, torch.cuda.is_available())
        checkpoint = torch.load(weights, map_location=self.device, weights_only=True)
        if not isinstance(checkpoint, dict) or "model" not in checkpoint:
            raise ValueError(f"{weights} is not a Progres-family checkpoint: expected a dict with a 'model' state dict")
        self.model = progres_module.Model().to(self.device)
        self.model.load_state_dict(checkpoint["model"])
        self.model.eval()
        self.nodes = None
        if want_nodes:
            self.model.get_submodule(entry.node_module).register_forward_hook(self._capture)

    def _capture(self, module, args, output):
        self.nodes = output.detach()

    def embed(self, path):
        """(vector (D,), node array (residues, D) or None) of the first chain of the structure at `path`."""
        coords = self.read_coords(path, **self.reader_args)
        if len(coords) < MIN_RESIDUES:
            raise ValueError(
                f"{len(coords)} C-alpha atom(s) in the first chain (ATOM records only; at least {MIN_RESIDUES} are needed): "
                "check that the file holds a protein chain and ends in .pdb or .cif")
        graph = self.coords_to_graph(coords)
        self.nodes = None
        with self.torch.no_grad():
            vector = self.model(graph.to(self.device)).squeeze(0).float().cpu().numpy()
        nodes = None
        if self.nodes is not None:
            nodes = self.nodes.squeeze(0).float().cpu().numpy()
        return vector, nodes


ADAPTERS = {"progres_graph": ProgresGraphEmbedder}


def build_embedder(entry, weights, device, want_nodes):
    return ADAPTERS[entry.adapter](entry, weights, device, want_nodes)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--config-yaml", required=True)
    p.add_argument("--inputs-json", required=True)
    p.add_argument("--embeddings-dir", required=True)
    p.add_argument("--embeddings-map-csv", required=True)
    p.add_argument("--node-embeddings-dir", required=True)
    p.add_argument("--node-embeddings-map-csv", required=True)
    p.add_argument("--similarities-csv", required=True)
    p.add_argument("--local-missing-csv", required=True)
    args = p.parse_args()

    with open(args.config_yaml) as fh:
        config = yaml.safe_load(fh) or {}
    entry = registry_entry(config["model"])
    want_nodes = config["node_embeddings"]
    want_pairs = config["first"] is not None

    pairs = read_pairs(config["first"], config["second"]) if want_pairs else []
    weights = weights_file(entry, os.environ)
    embedder = build_embedder(entry, weights, config["device"], want_nodes)
    print(f"Model {entry.name}: weights {weights} (sha256 {weights_digest(weights)}), device {embedder.device}")

    ds = load_datastream(args.inputs_json)
    step_id = step_id_from_table_path(args.local_missing_csv)
    os.makedirs(args.embeddings_dir, exist_ok=True)
    if want_nodes:
        os.makedirs(args.node_embeddings_dir, exist_ok=True)

    vectors: Dict[str, np.ndarray] = {}
    map_rows, node_map_rows, local_missing_rows = [], [], []
    attempted_ids: Set[str] = set()
    attempted = failures = 0

    for sid, path in iterate_files(ds):
        attempted += 1
        attempted_ids.add(sid)
        try:
            vector, nodes = embedder.embed(path)
            vector = check_embedding(vector, entry.dimension)
            if want_nodes:
                if nodes is None or np.ndim(nodes) != 2:
                    raise ValueError("the model returned no per-node embedding of shape (residues, dimension)")
                nodes = np.asarray(nodes, dtype=np.float32)
            out = os.path.join(args.embeddings_dir, f"{sid}.npy")
            np.save(out, vector)
            map_rows.append({"id": sid, "file": out})
            if want_nodes:
                node_out = os.path.join(args.node_embeddings_dir, f"{sid}.npy")
                np.save(node_out, nodes)
                node_map_rows.append({"id": sid, "file": node_out})
            vectors[sid] = vector
            print(f"{sid}: embedding {vector.shape}" + (f", {len(nodes)} residues" if want_nodes else ""))
        except Exception as e:
            print(f"WARNING: {sid} failed: {e}", file=sys.stderr)
            failures += 1
            local_missing_rows.append({"id": sid, "removed_by": step_id, "kind": "failure", "cause": str(e)[:200]})

    written = [args.embeddings_map_csv, args.local_missing_csv]
    written += [args.node_embeddings_map_csv] if want_nodes else []
    written += [args.similarities_csv] if want_pairs else []
    for path in written:
        os.makedirs(os.path.dirname(path), exist_ok=True)
    pd.DataFrame(map_rows, columns=["id", "file"]).to_csv(args.embeddings_map_csv, index=False)
    print(f"embeddings: {args.embeddings_map_csv} ({len(map_rows)} rows)")
    if want_nodes:
        pd.DataFrame(node_map_rows, columns=["id", "file"]).to_csv(args.node_embeddings_map_csv, index=False)
        print(f"node embeddings: {args.node_embeddings_map_csv} ({len(node_map_rows)} rows)")
    if want_pairs:
        rows, pair_failures = pair_rows(pairs, vectors, set(ds.ids_expanded), attempted_ids)
        pd.DataFrame(rows, columns=SIMILARITY_COLUMNS).to_csv(args.similarities_csv, index=False)
        for pid, kind, cause in pair_failures:
            local_missing_rows.append({"id": pid, "removed_by": step_id, "kind": kind, "cause": cause[:200]})
        print(f"similarities: {args.similarities_csv} ({len(rows)} rows, {len(pair_failures)} NaN)")
    pd.DataFrame(local_missing_rows, columns=MISSING_COLUMNS).to_csv(args.local_missing_csv, index=False)
    print(f"Local failures: {failures} input(s)")

    if not attempted:
        print("ERROR: no input entities to process (upstream stream is empty or its files are absent)", file=sys.stderr)
        sys.exit(1)
    if failures == attempted:
        print(f"ERROR: all {attempted} attempted id(s) failed", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
