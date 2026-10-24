#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""BioEmu2 runner.

For each (sequence, ligand set) combination declared by the wrapper, this
script writes a bioemu2 system YAML, invokes ``bioemu2 sample-equilibrium``,
then splits the resulting trajectory into per-conformer PDBs and copies the
compact .xtc + topology.cif into the trajectories stream.

Writes:
  structures_map   per-conformer PDB map
  trajectories_map one xtc + topology per run
  summary_csv      per-run summary (sequence.id, ligand.id, n_samples, n_residues)
  physicality_csv  per-frame physicality checks (bioemu2 physicality.json)
  missing_csv      failures / upstream-missing propagation (via wrapper's
                   generate_missing_propagation; this file writes LOCAL rows
                   into --local-missing-csv per the TEMPLATE contract)
"""

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys

# bioemu2 brings matplotlib via a transitive import. Headless-safe default.
os.environ.setdefault("MPLBACKEND", "Agg")

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from biopipelines.biopipelines_io import (  # noqa: E402
    iterate_values, load_datastream, step_id_from_table_path,
)

import mdtraj as md  # noqa: E402


SUM_COLS = ["id", "sequence.id", "ligand.id", "n_samples", "n_residues"]
PHYS_COLS = ["id", "frame", "invalid_bond", "clash", "chirality", "chain_break"]
MISS_COLS = ["id", "removed_by", "kind", "cause"]


def _bioemu2_cache_key(sequence: str) -> str:
    """Reproduce bioemu2's ColabFold cache key: sha256 of the sorted-keys
    JSON of the request dict. Change only when upstream changes the schema."""
    payload = {"sequences": [sequence], "use_pairing": False,
               "pairing_strategy": "greedy", "use_env": True}
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()


def _a3m_from_file(path: str) -> str:
    """Read an .a3m file and return its contents (bioemu2's cache stores a3m
    text directly)."""
    with open(path) as fh:
        return fh.read()


def _a3m_from_boltz_csv(path: str) -> str:
    """Convert a Boltz2 MSA CSV (``key,sequence``) to an a3m string. First
    row is the query (header ``>101``); subsequent rows are hits."""
    df = pd.read_csv(path)
    lines = []
    for i, row in enumerate(df.itertuples(index=False)):
        header = ">101" if i == 0 else f">hit_{i}"
        lines.append(header)
        lines.append(row.sequence)
    return "\n".join(lines) + "\n"


def _populate_bioemu2_msa_cache(sequence: str, msa_path: str) -> None:
    """Pre-populate ``~/.cache/bioemu2/colabfold_msa/<hash>.json`` for this
    sequence so bioemu2 skips the ColabFold search. Idempotent: a present
    cache entry is left alone."""
    cache_dir = os.path.join(os.path.expanduser("~"), ".cache",
                             "bioemu2", "colabfold_msa")
    os.makedirs(cache_dir, exist_ok=True)
    key = _bioemu2_cache_key(sequence)
    out = os.path.join(cache_dir, f"{key}.json")
    if os.path.exists(out):
        print(f"  cache hit for {key[:12]}... (left in place)")
        return
    if msa_path.endswith(".a3m"):
        a3m = _a3m_from_file(msa_path)
    elif msa_path.endswith(".csv"):
        a3m = _a3m_from_boltz_csv(msa_path)
    else:
        print(f"WARNING: unknown MSA format for {msa_path}; skipping cache populate",
              file=sys.stderr)
        return
    entry = {
        "request": {
            "sequences": [sequence],
            "use_pairing": False,
            "pairing_strategy": "greedy",
            "use_env": True,
        },
        "a3ms": [a3m],
        "sha256": hashlib.sha256(a3m.encode()).hexdigest(),
    }
    with open(out, "w") as fh:
        json.dump(entry, fh)
    hits = a3m.count(">") - 1
    print(f"  wrote cache entry for {key[:12]}... (hits={hits})")


def _write_system_yaml(path: str, sequence: str, ligand_smiles: str,
                       contacts: list, covalent: dict) -> None:
    """Write a bioemu2 system YAML with one protein chain A, optional ligand
    chain B, and optional contacts / covalent linkage."""
    lines = ["entities:",
             "  - entity_type: protein",
             f"    sequence: {sequence}",
             "    asym_ids: [A]"]
    if ligand_smiles:
        lines += ["  - entity_type: ligand",
                  f'    smiles: "{ligand_smiles}"',
                  "    asym_ids: [B]"]
    if contacts:
        lines.append("contacts:")
        for c in contacts:
            lines.append(f"  - {json.dumps(c)}")
    if covalent:
        lines.append("covalent_linkage:")
        for k, v in covalent.items():
            lines.append(f"  {k}: {json.dumps(v)}")
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


def _sample_once(system_yaml: str, out_dir: str, num_samples: int,
                 batch_size: int, sampler: str, physicality_filter: str,
                 seed=None) -> None:
    """Invoke ``bioemu2 sample-equilibrium`` as a subprocess. Running it
    out-of-process keeps bioemu2's lightning context from leaking across
    sequential runs (one call = one checkpoint load, released on exit)."""
    cmd = ["bioemu2", "sample-equilibrium",
           "--system", system_yaml,
           "--out", out_dir,
           "--num-samples", str(num_samples),
           "--batch-size", str(batch_size),
           "--sampler", sampler,
           "--physicality-filter", physicality_filter]
    if seed is not None:
        cmd += ["--seed", str(seed)]
    print(f"  running: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)


def _split_trajectory(topology_cif: str, samples_xtc: str, run_id: str,
                      structures_dir: str) -> list:
    """Write one PDB per sampled frame; returns [(frame_id, pdb_path), ...]."""
    os.makedirs(structures_dir, exist_ok=True)
    traj = md.load(samples_xtc, top=topology_cif)
    out = []
    for i in range(traj.n_frames):
        frame_id = f"{run_id}_{i+1}"
        pdb_path = os.path.join(structures_dir, f"{frame_id}.pdb")
        traj[i].save_pdb(pdb_path)
        out.append((frame_id, pdb_path))
    return out


def _parse_physicality(run_id: str, physicality_json: str) -> list:
    """Flatten bioemu2's physicality.json into per-frame rows."""
    if not os.path.isfile(physicality_json):
        return []
    try:
        data = json.load(open(physicality_json))
    except Exception as e:
        print(f"WARNING: cannot parse {physicality_json}: {e}", file=sys.stderr)
        return []
    rows = []
    frames = data if isinstance(data, list) else data.get("frames", [])
    for i, f in enumerate(frames):
        rows.append({
            "id": run_id,
            "frame": i + 1,
            "invalid_bond": f.get("invalid_bond"),
            "clash": f.get("clash"),
            "chirality": f.get("chirality"),
            "chain_break": f.get("chain_break"),
        })
    return rows


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--sequences-json", required=True)
    p.add_argument("--ligands-json", default=None,
                   help="JSON serialization of the ligands DataStream (omit for apo)")
    p.add_argument("--num-samples", type=int, default=100)
    p.add_argument("--batch-size", type=int, default=16)
    p.add_argument("--sampler", default="edm")
    p.add_argument("--physicality-filter", default="high_quality")
    p.add_argument("--seed", type=int, default=None)
    p.add_argument("--contacts-json", default=None)
    p.add_argument("--covalent-linkage-json", default=None)
    p.add_argument("--precomputed-msa-json", default=None,
                   help="DataStream JSON mapping sequence id -> .a3m / Boltz-MSA .csv path")
    p.add_argument("--structures-dir", required=True)
    p.add_argument("--structures-map", required=True)
    p.add_argument("--trajectories-dir", required=True)
    p.add_argument("--trajectories-map", required=True)
    p.add_argument("--scratch-dir", required=True)
    p.add_argument("--summary-csv", required=True)
    p.add_argument("--physicality-csv", required=True)
    p.add_argument("--missing-csv", required=True)
    p.add_argument("--upstream-missing", default=None)
    args = p.parse_args()

    for d in (args.structures_dir, args.trajectories_dir, args.scratch_dir):
        os.makedirs(d, exist_ok=True)

    seqs = load_datastream(args.sequences_json)
    sequences = {sid: str(v.get("sequence", "")).strip()
                 for sid, v in iterate_values(seqs, columns=["sequence"])}

    if args.ligands_json:
        ligs_ds = load_datastream(args.ligands_json)
        ligands = {lid: str(v.get("smiles", "")).strip()
                   for lid, v in iterate_values(ligs_ds, columns=["smiles"])}
    else:
        ligands = {None: ""}

    contacts = json.loads(args.contacts_json) if args.contacts_json else []
    covalent = json.loads(args.covalent_linkage_json) if args.covalent_linkage_json else None

    if args.precomputed_msa_json:
        msa_ds = load_datastream(args.precomputed_msa_json)
        msa_paths = {mid: str(v.get("file", "")).strip()
                     for mid, v in iterate_values(msa_ds, columns=["file"])}
        for sid, seq in sequences.items():
            path = msa_paths.get(sid)
            if path:
                try:
                    _populate_bioemu2_msa_cache(seq, path)
                except Exception as e:
                    print(f"WARNING: precomputed MSA for {sid} failed to install: {e}",
                          file=sys.stderr)

    step_id = step_id_from_table_path(args.missing_csv)
    struct_rows, traj_rows, sum_rows, phys_rows, missing = [], [], [], [], []

    for sid, seq in sequences.items():
        if not seq:
            missing.append({"id": sid, "removed_by": step_id, "kind": "failure",
                            "cause": "empty sequence"})
            continue
        for lid, lsmi in ligands.items():
            run_id = sid if lid is None else f"{sid}+{lid}"
            work = os.path.join(args.scratch_dir, run_id)
            os.makedirs(work, exist_ok=True)
            system_yaml = os.path.join(work, "system.yaml")
            _write_system_yaml(system_yaml, seq, lsmi, contacts, covalent)
            try:
                _sample_once(system_yaml, work, args.num_samples, args.batch_size,
                             args.sampler, args.physicality_filter, args.seed)
            except subprocess.CalledProcessError as e:
                print(f"WARNING: {run_id} bioemu2 sampling failed: {e}",
                      file=sys.stderr)
                missing.append({"id": run_id, "removed_by": step_id,
                                "kind": "failure", "cause": str(e)[:200]})
                continue

            topology_cif = os.path.join(work, "topology.cif")
            samples_xtc = os.path.join(work, "samples.xtc")
            if not (os.path.isfile(topology_cif) and os.path.isfile(samples_xtc)):
                missing.append({"id": run_id, "removed_by": step_id,
                                "kind": "failure",
                                "cause": "bioemu2 produced no topology/samples"})
                continue

            try:
                frames = _split_trajectory(topology_cif, samples_xtc, run_id,
                                           args.structures_dir)
            except Exception as e:
                print(f"WARNING: {run_id} trajectory split failed: {e}",
                      file=sys.stderr)
                missing.append({"id": run_id, "removed_by": step_id,
                                "kind": "failure",
                                "cause": f"split: {e}"[:200]})
                continue

            for fid, p_pdb in frames:
                struct_rows.append({"id": fid, "file": p_pdb, "run.id": run_id})

            dst_xtc = os.path.join(args.trajectories_dir, f"{run_id}.xtc")
            dst_top = os.path.join(args.trajectories_dir, f"{run_id}_topology.cif")
            shutil.copyfile(samples_xtc, dst_xtc)
            shutil.copyfile(topology_cif, dst_top)
            traj_rows.append({"id": run_id, "file": dst_xtc, "topology": dst_top})

            n_res = int(md.load(topology_cif).n_residues)
            sum_rows.append({"id": run_id, "sequence.id": sid,
                             "ligand.id": lid if lid is not None else "",
                             "n_samples": len(frames), "n_residues": n_res})
            phys_rows += _parse_physicality(run_id,
                                            os.path.join(work, "physicality.json"))
            print(f"  {run_id}: {len(frames)} conformers, {n_res} residues")

    upstream_rows = []
    if args.upstream_missing and os.path.exists(args.upstream_missing):
        try:
            up_df = pd.read_csv(args.upstream_missing)
            if not up_df.empty:
                upstream_rows = up_df.to_dict("records")
        except Exception as e:
            print(f"Warning: could not read upstream missing.csv: {e}", file=sys.stderr)

    for d in (args.structures_map, args.trajectories_map,
              args.summary_csv, args.physicality_csv, args.missing_csv):
        os.makedirs(os.path.dirname(d), exist_ok=True)
    pd.DataFrame(struct_rows, columns=["id", "file", "run.id"]).to_csv(
        args.structures_map, index=False)
    pd.DataFrame(traj_rows, columns=["id", "file", "topology"]).to_csv(
        args.trajectories_map, index=False)
    pd.DataFrame(sum_rows, columns=SUM_COLS).to_csv(args.summary_csv, index=False)
    pd.DataFrame(phys_rows, columns=PHYS_COLS).to_csv(args.physicality_csv, index=False)
    pd.DataFrame(upstream_rows + missing, columns=MISS_COLS).to_csv(
        args.missing_csv, index=False)
    print(f"Structures map: {args.structures_map} ({len(struct_rows)} rows)")
    print(f"Trajectories map: {args.trajectories_map} ({len(traj_rows)} rows)")
    print(f"Summary: {args.summary_csv} ({len(sum_rows)} rows)")
    print(f"Physicality: {args.physicality_csv} ({len(phys_rows)} rows)")
    print(f"Missing: {args.missing_csv} ({len(upstream_rows) + len(missing)} rows)")

    if not struct_rows:
        sys.exit(1)


if __name__ == "__main__":
    main()
