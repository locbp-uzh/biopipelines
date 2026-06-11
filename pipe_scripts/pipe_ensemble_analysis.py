#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Execution helper for EnsembleAnalysis.

Reads an input structures stream, builds an ensemble per output id, superposes
the conformers (Kabsch least-squares on CA or backbone atoms), and writes:

  - one per-ensemble resi-csv (id, chain, resi, rmsf, rmsd_mean) into rmsf-dir,
    plus the rmsf stream map_table;
  - an `ensemble` summary table (one row per ensemble);
  - a `frames` table (one row per conformer: rmsd_to_ref, radius of gyration).

Two ensemble shapes are auto-detected per input (no user knob):

  - multi-model file: a single input file with >=2 MODEL records -> each MODEL
    is a frame; output id = the input id.
  - conformer set grouped by provenance: single-model files sharing a `<x>.id`
    provenance column in the input map_table -> frames pooled per group; output
    id = the group key. If no provenance column is present, each input id is its
    own one-frame ensemble (RMSF is then zero, but the metrics still emit).

RMSF is the root-mean-square fluctuation of each residue's selected atoms about
their mean position after superposition; rmsd_mean is the RMS deviation of that
residue from the reference. Residues are matched across frames by (chain,
res_num); only residues present in *every* frame are reported.
"""

import os
import sys
import argparse
from collections import defaultdict

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from biopipelines.biopipelines_io import load_datastream, iterate_files
from biopipelines.id_map_utils import get_mapped_ids
from biopipelines.pdb_parser import parse_models_file

# Backbone atom names for the "backbone" selection; "CA" uses just the alpha C.
_BACKBONE = ("N", "CA", "C", "O")


def parse_models(pdb_path):
    """Parse a structure (PDB or mmCIF) into a list of models, each a list of Atom.

    ATOM records only — HETATM/waters excluded, since RMSF is a protein-residue
    metric and HETATM res_num collide across chains.
    """
    return parse_models_file(pdb_path, records=("ATOM",))


def select_coords(model, selection):
    """Return {(chain, res_num): {atom_name: (x,y,z)}} for the selected atoms.

    For "CA": only the alpha carbon per residue. For "backbone": N, CA, C, O.
    """
    wanted = ("CA",) if selection == "CA" else _BACKBONE
    residues = defaultdict(dict)
    for a in model:
        name = a.atom_name.strip()
        if name in wanted:
            residues[(a.chain, a.res_num)][name] = (a.x, a.y, a.z)
    # Keep only residues that carry the full requested atom set, so every frame
    # contributes the same atom count per residue (needed for a clean RMSF).
    complete = {}
    for key, atoms in residues.items():
        if all(w in atoms for w in wanted):
            complete[key] = atoms
    return complete


def kabsch(P, Q):
    """Optimal rotation aligning P onto Q (both centered N×3). Returns R (3×3)."""
    H = P.T @ Q
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1.0, 1.0, d])
    return U @ D @ Vt  # apply as P @ R


def radius_of_gyration(coords):
    """Rg of an N×3 coordinate array about its centroid (unweighted)."""
    c = coords - coords.mean(axis=0)
    return float(np.sqrt((c ** 2).sum() / len(c)))


def analyze_ensemble(frames, selection, reference):
    """Superpose frames and compute per-residue RMSF/RMSD-to-mean + frame metrics.

    frames: list of {(chain, res_num): {atom_name: (x,y,z)}} (one per conformer).

    Returns (resi_rows, frame_rows, summary) or None if no shared residues.
    """
    wanted = ("CA",) if selection == "CA" else _BACKBONE

    # Residues present in EVERY frame, in a stable (chain, res_num) order.
    common = set(frames[0])
    for fr in frames[1:]:
        common &= set(fr)
    if not common:
        return None
    keys = sorted(common)

    # Stack each frame's selected coords into an (n_atoms, 3) array, atom order
    # = residue order × wanted-atom order. fit_idx marks the rows used to fit
    # the superposition (all selected atoms here — same set for CA & backbone).
    def frame_matrix(fr):
        rows = []
        for key in keys:
            for w in wanted:
                rows.append(fr[key][w])
        return np.asarray(rows, dtype=float)

    mats = [frame_matrix(fr) for fr in frames]
    n_atoms = mats[0].shape[0]

    # Reference: first frame, or iterative mean (align-to-first -> mean -> realign).
    def superpose_to(ref):
        ref_c = ref - ref.mean(axis=0)
        out = []
        for m in mats:
            m_c = m - m.mean(axis=0)
            R = kabsch(m_c, ref_c)
            out.append(m_c @ R)
        return out

    aligned = superpose_to(mats[0])
    if reference == "mean":
        mean_coords = np.mean(aligned, axis=0)
        aligned = superpose_to(mean_coords)
        ref_coords = np.mean(aligned, axis=0)
    else:  # "first"
        ref_coords = aligned[0]

    stack = np.stack(aligned, axis=0)              # (n_frames, n_atoms, 3)
    mean_pos = stack.mean(axis=0)                  # (n_atoms, 3)

    # Per-atom fluctuation about the mean, and deviation about the reference.
    fluct_sq = ((stack - mean_pos) ** 2).sum(axis=2).mean(axis=0)   # (n_atoms,)
    dev_sq = ((stack - ref_coords) ** 2).sum(axis=2).mean(axis=0)   # (n_atoms,)

    # Collapse per-atom -> per-residue (mean over the residue's atoms), then RMS.
    n_per_res = len(wanted)
    resi_rows = []
    rmsf_vals = []
    for i, key in enumerate(keys):
        sl = slice(i * n_per_res, (i + 1) * n_per_res)
        rmsf = float(np.sqrt(fluct_sq[sl].mean()))
        rmsd_mean = float(np.sqrt(dev_sq[sl].mean()))
        chain, res_num = key
        resi_rows.append({
            "chain": chain,
            "resi": res_num,
            "rmsf": round(rmsf, 4),
            "rmsd_mean": round(rmsd_mean, 4),
        })
        rmsf_vals.append(rmsf)

    # Per-frame metrics: RMSD to the reference, radius of gyration.
    frame_rows = []
    for fi in range(stack.shape[0]):
        rmsd_to_ref = float(np.sqrt(((stack[fi] - ref_coords) ** 2).sum() / n_atoms))
        rg = radius_of_gyration(stack[fi])
        frame_rows.append({
            "frame": fi + 1,
            "rmsd_to_ref": round(rmsd_to_ref, 4),
            "rg": round(rg, 4),
        })

    rgs = np.array([fr["rg"] for fr in frame_rows], dtype=float)
    summary = {
        "n_frames": len(frames),
        "n_residues": len(keys),
        "mean_rmsf": round(float(np.mean(rmsf_vals)), 4),
        "max_rmsf": round(float(np.max(rmsf_vals)), 4),
        "rg_mean": round(float(rgs.mean()), 4),
        "rg_std": round(float(rgs.std()), 4),
    }
    return resi_rows, frame_rows, summary


def build_ensembles(ds, selection, groups_ds):
    """Yield (group_id, [frame_residue_maps], group_id) per ensemble.

    The groups stream defines the partition (Consensus convention): each input
    id is matched to a group id via framework id matching, and ALL of each
    member's models (one MODEL or many) are pooled into that group's ensemble.
    The output id is always the group id, so it matches the declared ids.
    """
    files = {sid: fp for sid, fp in iterate_files(ds)}
    parsed = {sid: parse_models(fp) for sid, fp in files.items()}

    group_ids = list(groups_ds.ids_expanded)
    mapped = get_mapped_ids(group_ids, list(parsed.keys()), unique=False)
    for gid in group_ids:
        members = mapped.get(gid, [])
        if not members:
            continue
        frames = [select_coords(m, selection)
                  for member in sorted(members)
                  for m in parsed[member]]
        yield gid, frames, gid


def main():
    ap = argparse.ArgumentParser(description="EnsembleAnalysis execution helper")
    ap.add_argument("--structures-json", required=True)
    ap.add_argument("--selection", required=True, choices=["CA", "backbone"])
    ap.add_argument("--reference", required=True, choices=["mean", "first"])
    ap.add_argument("--rmsf-dir", required=True)
    ap.add_argument("--rmsf-map", required=True)
    ap.add_argument("--residues-csv", required=True)
    ap.add_argument("--ensemble-csv", required=True)
    ap.add_argument("--frames-csv", required=True)
    ap.add_argument("--groups-json", required=True,
                    help="Groups DataStream JSON; its ids define the ensemble partition")
    args = ap.parse_args()

    ds = load_datastream(args.structures_json)
    groups_ds = load_datastream(args.groups_json)

    map_rows = []          # rmsf stream map_table
    residue_rows_all = []  # merged per-residue table
    ensemble_rows = []     # summary table
    frame_rows_all = []    # per-conformer table
    failed = []

    for output_id, frames, _prov in build_ensembles(ds, args.selection, groups_ds):
        try:
            if not frames:
                raise ValueError("no frames")
            result = analyze_ensemble(frames, args.selection, args.reference)
            if result is None:
                raise ValueError("no residues shared across all frames")
            resi_rows, frame_rows, summary = result

            # Per-ensemble resi-csv (the resi-csv stream file).
            out_csv = os.path.join(args.rmsf_dir, f"{output_id}.csv")
            resi_df = pd.DataFrame(
                [{"id": output_id, **r} for r in resi_rows],
                columns=["id", "chain", "resi", "rmsf", "rmsd_mean"],
            )
            resi_df.to_csv(out_csv, index=False)
            map_rows.append({"id": output_id, "file": out_csv})
            residue_rows_all.extend(resi_df.to_dict("records"))

            ensemble_rows.append({"id": output_id, **summary})
            for fr in frame_rows:
                frame_rows_all.append({"id": output_id, **fr})
            print(f"  {output_id}: {summary['n_frames']} frames, "
                  f"{summary['n_residues']} residues, "
                  f"mean RMSF {summary['mean_rmsf']} A")
        except Exception as e:
            print(f"WARNING: ensemble '{output_id}' failed: {e}", file=sys.stderr)
            failed.append(output_id)

    # Stream map_table — only ids whose file was actually written.
    pd.DataFrame(map_rows, columns=["id", "file"]).to_csv(args.rmsf_map, index=False)
    pd.DataFrame(
        residue_rows_all, columns=["id", "chain", "resi", "rmsf", "rmsd_mean"]
    ).to_csv(args.residues_csv, index=False)
    pd.DataFrame(
        ensemble_rows,
        columns=["id", "n_frames", "n_residues", "mean_rmsf", "max_rmsf", "rg_mean", "rg_std"],
    ).to_csv(args.ensemble_csv, index=False)
    pd.DataFrame(
        frame_rows_all, columns=["id", "frame", "rmsd_to_ref", "rg"]
    ).to_csv(args.frames_csv, index=False)

    if failed:
        print(f"Failed {len(failed)}/{len(failed) + len(map_rows)}: {failed}", file=sys.stderr)
    if not map_rows:
        print("ERROR: no ensemble produced any output", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
