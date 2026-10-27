#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Post-process RFD4-Proteina outputs into BioPipelines streams and tables.

`rfd4-generate` writes `<run>/inference_samples/<key>_<i>/<key>_<i>_generated.cif` for every spec key and replicate `i` counted from 0. Here each of those becomes the design `<key>_<i+1>`: the annotated CIF is converted to PDB, its polymer sequences are read off, and one row per design lands in the structures, sequences and metrics tables. Upstream's inline metrics (`inference_metrics/per_example_metrics*.csv`) are re-keyed to the design ids.

The ids come from the inputs manifest, not from the files found, so a design that was expected and never written is a failure row rather than a silent gap.

Usage:
    python pipe_rfd4proteina_postprocess.py \\
        --run-dir <rfd4 out_dir> --manifest-json <inputs_manifest.json> --num-designs N \\
        --structures-dir <dir> --structures-map-csv <csv> --sequences-csv <csv> \\
        --structures-csv <csv> --metrics-csv <csv> --local-missing-csv <csv>
"""

import argparse
import glob
import json
import math
import os
import sys
from typing import Dict, List, Optional

import pandas as pd

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from biopipelines.biopipelines_io import step_id_from_table_path  # noqa: E402
# The CIF to PDB conversion and the sequence reader are shared with RFdiffusion3, which solved the same Biopython column and ligand bond problems.
from pipe_rfdiffusion3_postprocess import convert_cif_to_pdb, extract_sequence_from_pdb  # noqa: E402

# Columns of upstream's per-example metrics CSV that are paths or the configuration dump; seed, ckpt_path and ckpt_name stay as the run record.
DROPPED_PREFIXES = ("generation.",)
DROPPED_COLUMNS = ("designed_structure_path", "output_base_dir", "dataset_name", "metrics_tags")


def design_id(source_id: str, replicate: int) -> str:
    """BioPipelines id of replicate `replicate` (counted from 1)."""
    return f"{source_id}_{replicate}"


def upstream_id(source_id: str, replicate: int) -> str:
    """Upstream names replicates from 0: key `binder` gives `binder_0`, `binder_1`, ..."""
    return f"{source_id}_{replicate - 1}"


def locate_generated_cif(run_dir: str, uid: str) -> Optional[str]:
    path = os.path.join(run_dir, "inference_samples", uid, f"{uid}_generated.cif")
    return path if os.path.isfile(path) and os.path.getsize(path) > 0 else None


def read_inference_metrics(run_dir: str) -> Optional[pd.DataFrame]:
    """Upstream's per-example metrics, all job and re-evaluation files together, or None when scoring was off."""
    frames = []
    for path in sorted(glob.glob(os.path.join(run_dir, "inference_metrics", "per_example_metrics*.csv"))):
        try:
            frames.append(pd.read_csv(path, dtype={"example_id": str}))
        except Exception as e:
            print(f"WARNING: could not read {path}: {e}", file=sys.stderr)
    if not frames:
        return None
    merged = pd.concat(frames, ignore_index=True)
    if "example_id" not in merged.columns:
        return None
    merged = merged.drop_duplicates("example_id", keep="last")
    keep = [c for c in merged.columns
            if not c.startswith(DROPPED_PREFIXES) and c not in DROPPED_COLUMNS]
    return merged[keep]


def build_metrics_table(ids: List[str], uid_by_id: Dict[str, str], source_by_id: Dict[str, str],
                        frame: Optional[pd.DataFrame]) -> pd.DataFrame:
    """One row per design id with its source structure; the upstream metric columns follow, NaN where scoring did not reach a design."""
    base = pd.DataFrame({"id": ids, "source_id": [source_by_id[i] for i in ids]})
    if frame is None:
        return base
    id_by_uid = {uid: oid for oid, uid in uid_by_id.items()}
    scored = frame.copy()
    scored["id"] = scored["example_id"].map(id_by_uid)
    scored = scored.dropna(subset=["id"]).drop(columns=["example_id"])
    return base.merge(scored, on="id", how="left")


def main():
    p = argparse.ArgumentParser(description="Post-process RFD4-Proteina outputs")
    p.add_argument("--run-dir", required=True)
    p.add_argument("--manifest-json", required=True)
    p.add_argument("--num-designs", type=int, required=True)
    p.add_argument("--structures-dir", required=True)
    p.add_argument("--structures-map-csv", required=True)
    p.add_argument("--sequences-csv", required=True)
    p.add_argument("--structures-csv", required=True)
    p.add_argument("--metrics-csv", required=True)
    p.add_argument("--local-missing-csv", required=True)
    args = p.parse_args()

    if not os.path.isfile(args.manifest_json):
        print(f"ERROR: inputs manifest not found: {args.manifest_json}", file=sys.stderr)
        sys.exit(1)
    with open(args.manifest_json) as f:
        manifest = json.load(f)
    failed_inputs = manifest.get("failed", {})

    os.makedirs(args.structures_dir, exist_ok=True)
    step_id = step_id_from_table_path(args.local_missing_csv)

    structures_rows, map_rows, sequences_rows, local_missing_rows = [], [], [], []
    ids, uid_by_id, source_by_id = [], {}, {}
    attempted = failures = 0

    for sid in manifest.get("attempted", []):
        for replicate in range(1, args.num_designs + 1):
            oid = design_id(sid, replicate)
            uid = upstream_id(sid, replicate)
            ids.append(oid)
            uid_by_id[oid] = uid
            source_by_id[oid] = sid
            attempted += 1
            row = {"id": oid, "source_id": sid, "design": replicate, "cif": "", "pdb": "",
                   "n_chains": math.nan, "n_residues": math.nan}
            try:
                if sid in failed_inputs:
                    raise RuntimeError(f"input preparation failed: {failed_inputs[sid]}")
                cif = locate_generated_cif(args.run_dir, uid)
                if cif is None:
                    raise FileNotFoundError(f"no generated CIF for {uid} under {args.run_dir}/inference_samples")
                row["cif"] = cif
                pdb_path = os.path.join(args.structures_dir, f"{oid}.pdb")
                if not convert_cif_to_pdb(cif, pdb_path):
                    if os.path.exists(pdb_path):
                        os.remove(pdb_path)
                    raise RuntimeError("CIF to PDB conversion failed (see the step log)")
                chains = extract_sequence_from_pdb(pdb_path)
                if not chains:
                    raise RuntimeError("the designed structure has no polymer chain")
                row["pdb"] = pdb_path
                row["n_chains"] = len(chains)
                row["n_residues"] = sum(len(s) for s in chains.values())
                map_rows.append({"id": oid, "file": pdb_path, "structures.id": sid})
                for chain_id, sequence in chains.items():
                    sequences_rows.append({
                        "id": f"{oid}_{chain_id}" if len(chains) > 1 else oid,
                        "source_id": oid,
                        "source_pdb": pdb_path,
                        "chain": chain_id,
                        "sequence": sequence,
                        "length": len(sequence),
                    })
            except Exception as e:
                print(f"WARNING: {oid} failed: {e}", file=sys.stderr)
                failures += 1
                local_missing_rows.append({
                    "id": oid, "removed_by": step_id, "kind": "failure", "cause": str(e)[:200]})
            structures_rows.append(row)

    for path in (args.structures_map_csv, args.sequences_csv, args.structures_csv,
                 args.metrics_csv, args.local_missing_csv):
        os.makedirs(os.path.dirname(path), exist_ok=True)

    pd.DataFrame(structures_rows, columns=["id", "source_id", "design", "cif", "pdb", "n_chains", "n_residues"]
                 ).to_csv(args.structures_csv, index=False)
    pd.DataFrame(map_rows, columns=["id", "file", "structures.id"]).to_csv(args.structures_map_csv, index=False)
    pd.DataFrame(sequences_rows, columns=["id", "source_id", "source_pdb", "chain", "sequence", "length"]
                 ).to_csv(args.sequences_csv, index=False)
    build_metrics_table(ids, uid_by_id, source_by_id, read_inference_metrics(args.run_dir)).to_csv(args.metrics_csv, index=False)
    pd.DataFrame(local_missing_rows, columns=["id", "removed_by", "kind", "cause"]
                 ).to_csv(args.local_missing_csv, index=False)

    print(f"Structures: {len(map_rows)} of {attempted} design(s) -> {args.structures_map_csv}")
    print(f"Sequences: {len(sequences_rows)} row(s) -> {args.sequences_csv}")
    print(f"Local failures: {failures}")

    # Nothing attempted means no input reached the model; every attempt failing means the step produced nothing.
    if not attempted:
        print("ERROR: no designs expected (no input structure was prepared)", file=sys.stderr)
        sys.exit(1)
    if failures == attempted:
        print(f"ERROR: all {attempted} expected design(s) failed", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
