#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""BioEmu runner. For each input sequence, calls bioemu.sample.main to
generate an equilibrium ensemble, optionally reconstructs side chains, then
splits the resulting trajectory into per-conformer PDBs and stages the
compact .xtc + topology. Writes the structures + trajectories maps and a
per-sequence summary."""

import argparse
import os
import shutil
import sys

# bioemu.sample runs in THIS process and transitively imports its vendored
# alphafold -> tensorflow -> keras -> matplotlib.pyplot. On Colab the inherited
# MPLBACKEND=module://matplotlib_inline.backend_inline is invalid outside the
# notebook kernel and crashes that import. Force a headless backend before any
# bioemu import. Harmless elsewhere (Agg is valid on cluster/local too).
os.environ["MPLBACKEND"] = "Agg"

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from biopipelines.biopipelines_io import load_datastream, iterate_values, step_id_from_table_path  # noqa: E402

import MDAnalysis as mda  # noqa: E402


SUM_COLS = ["id", "sequence.id", "n_samples", "n_residues"]


def run_sample(sequence, num_samples, batch_size, out_dir, filter_samples, msa_host_url):
    """Invoke bioemu.sample.main; returns (topology_pdb, samples_xtc).

    The accepted keyword set varies across bioemu releases (e.g. some
    versions take batch_size / msa_host_url, others don't), so pass only
    the kwargs the installed main() actually declares."""
    import inspect
    from bioemu.sample import main as sample_main

    accepted = set(inspect.signature(sample_main).parameters)
    kwargs = {"sequence": sequence, "num_samples": num_samples, "output_dir": out_dir}
    # Newer bioemu renamed batch_size -> batch_size_100 (a length-100 reference
    # batch it scales per sequence length); support both spellings so the
    # user's batch_size isn't silently dropped on current releases.
    if "batch_size" in accepted:
        kwargs["batch_size"] = batch_size
    elif "batch_size_100" in accepted:
        kwargs["batch_size_100"] = batch_size
    if "filter_samples" in accepted:
        kwargs["filter_samples"] = filter_samples
    if "msa_host_url" in accepted and msa_host_url:
        kwargs["msa_host_url"] = msa_host_url
    sample_main(**kwargs)

    topology = os.path.join(out_dir, "topology.pdb")
    xtc = os.path.join(out_dir, "samples.xtc")
    if not (os.path.isfile(topology) and os.path.isfile(xtc)):
        raise RuntimeError(f"bioemu.sample produced no topology/samples in {out_dir}")
    return topology, xtc


def run_sidechain_relax(topology, xtc, out_dir):
    """Run bioemu.sidechain_relax; returns (topology, xtc) of full-atom output."""
    from bioemu.sidechain_relax import main as relax_main

    relax_main(pdb_path=topology, xtc_path=xtc, outpath=out_dir)
    rec_pdb = os.path.join(out_dir, "samples_sidechain_rec.pdb")
    rec_xtc = os.path.join(out_dir, "samples_sidechain_rec.xtc")
    if os.path.isfile(rec_pdb) and os.path.isfile(rec_xtc):
        return rec_pdb, rec_xtc
    raise RuntimeError(f"sidechain_relax produced no reconstructed output in {out_dir}")


def split_trajectory(topology, xtc, seq_id, structures_dir):
    """Write one PDB per frame as <seq_id>_<n>.pdb. Returns list of (model_id, path)."""
    os.makedirs(structures_dir, exist_ok=True)
    u = mda.Universe(topology, xtc)
    out = []
    for i, _ts in enumerate(u.trajectory, start=1):
        model_id = f"{seq_id}_{i}"
        out_pdb = os.path.join(structures_dir, f"{model_id}.pdb")
        u.atoms.write(out_pdb)
        out.append((model_id, out_pdb))
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--sequences-json", required=True)
    p.add_argument("--num-samples", type=int, default=10)
    p.add_argument("--batch-size", type=int, default=10)
    p.add_argument("--structures-dir", required=True)
    p.add_argument("--structures-map", required=True)
    p.add_argument("--trajectories-dir", required=True)
    p.add_argument("--trajectories-map", required=True)
    p.add_argument("--scratch-dir", required=True)
    p.add_argument("--summary-csv", required=True)
    p.add_argument("--missing-csv", required=True)
    p.add_argument("--reconstruct-sidechains", action="store_true")
    p.add_argument("--no-filter-samples", action="store_true")
    p.add_argument("--msa-host-url", default="")
    p.add_argument("--upstream-missing", default=None)
    args = p.parse_args()

    for d in (args.structures_dir, args.trajectories_dir, args.scratch_dir):
        os.makedirs(d, exist_ok=True)

    ds = load_datastream(args.sequences_json)

    struct_rows, traj_rows, sum_rows, missing = [], [], [], []
    step_id = step_id_from_table_path(args.missing_csv)
    for seq_id, values in iterate_values(ds, columns=["sequence"]):
        sequence = str(values.get("sequence", "")).strip()
        if not sequence:
            missing.append({"id": seq_id, "removed_by": step_id, "kind": "failure", "cause": "empty sequence"})
            continue
        work = os.path.join(args.scratch_dir, seq_id)
        os.makedirs(work, exist_ok=True)
        try:
            topology, xtc = run_sample(
                sequence, args.num_samples, args.batch_size, work,
                filter_samples=not args.no_filter_samples,
                msa_host_url=args.msa_host_url,
            )
            if args.reconstruct_sidechains:
                topology, xtc = run_sidechain_relax(topology, xtc, work)

            models = split_trajectory(topology, xtc, seq_id, args.structures_dir)
            for model_id, out_pdb in models:
                struct_rows.append({"id": model_id, "file": out_pdb, "sequence.id": seq_id})

            dst_xtc = os.path.join(args.trajectories_dir, f"{seq_id}.xtc")
            dst_top = os.path.join(args.trajectories_dir, f"{seq_id}_topology.pdb")
            shutil.copyfile(xtc, dst_xtc)
            shutil.copyfile(topology, dst_top)
            traj_rows.append({"id": seq_id, "file": dst_xtc, "topology": dst_top})

            n_res = len(mda.Universe(topology).residues)
            sum_rows.append({"id": seq_id, "sequence.id": seq_id,
                             "n_samples": len(models), "n_residues": n_res})
            print(f"  {seq_id}: {len(models)} conformers, {n_res} residues")
        except Exception as e:
            print(f"WARNING: {seq_id} BioEmu failed: {e}", file=sys.stderr)
            missing.append({"id": seq_id, "removed_by": step_id, "kind": "failure", "cause": str(e)[:200]})

    upstream_rows = []
    if args.upstream_missing and os.path.exists(args.upstream_missing):
        try:
            up_df = pd.read_csv(args.upstream_missing)
            if not up_df.empty:
                upstream_rows = up_df.to_dict("records")
        except Exception as e:
            print(f"Warning: could not read upstream missing.csv: {e}", file=sys.stderr)

    all_missing = upstream_rows + missing

    for d in (args.structures_map, args.trajectories_map, args.summary_csv, args.missing_csv):
        os.makedirs(os.path.dirname(d), exist_ok=True)
    pd.DataFrame(struct_rows, columns=["id", "file", "sequence.id"]).to_csv(args.structures_map, index=False)
    pd.DataFrame(traj_rows, columns=["id", "file", "topology"]).to_csv(args.trajectories_map, index=False)
    pd.DataFrame(sum_rows, columns=SUM_COLS).to_csv(args.summary_csv, index=False)
    pd.DataFrame(all_missing, columns=["id", "removed_by", "kind", "cause"]).to_csv(args.missing_csv, index=False)
    print(f"Structures map: {args.structures_map} ({len(struct_rows)} rows)")
    print(f"Trajectories map: {args.trajectories_map} ({len(traj_rows)} rows)")
    print(f"Summary: {args.summary_csv} ({len(sum_rows)} rows)")
    print(f"Missing: {args.missing_csv} ({len(all_missing)} rows)")

    if not struct_rows:
        sys.exit(1)


if __name__ == "__main__":
    main()
