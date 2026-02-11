#!/usr/bin/env python3
# Copyright (c) 2026 Gianluca Quargnali @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
GNINA docking pipeline helper script.

Orchestrates the full docking workflow:
1. Protein preparation (protonation with OpenBabel)
2. Ligand/conformer preparation (RDKit conformer generation or loading)
3. Docking execution (GNINA binary calls)
4. Pose consistency analysis (RMSD clustering)
5. Results aggregation (docking_results.csv + conformer_ranking.csv)

Usage:
    python pipe_gnina.py <config_json>
"""

import json
import os
import subprocess
import sys
from collections import defaultdict

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign

from biopipelines_io import load_datastream, iterate_files


# ---------------------------------------------------------------------------
# 1. Protein preparation
# ---------------------------------------------------------------------------

def prepare_proteins(structures_ds, output_folder, protonate, pH):
    """
    Prepare protein structures for docking.

    If protonate=True, adds hydrogens using OpenBabel at the specified pH.
    Writes prepared proteins to prepared_proteins/ subfolder.

    Returns:
        dict: protein_id -> path to prepared PDB
    """
    prep_dir = os.path.join(output_folder, "prepared_proteins")
    os.makedirs(prep_dir, exist_ok=True)

    prepared = {}
    for protein_id, protein_file in iterate_files(structures_ds):
        if protonate:
            out_pdb = os.path.join(prep_dir, f"{protein_id}_h.pdb")
            cmd = [
                "obabel", protein_file, "-O", out_pdb,
                "-h", "-p", str(pH)
            ]
            print(f"Protonating {protein_id}: {' '.join(cmd)}")
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode != 0:
                print(f"Warning: obabel failed for {protein_id}: {result.stderr}")
                # Fall back to original file
                prepared[protein_id] = protein_file
            else:
                prepared[protein_id] = out_pdb
        else:
            prepared[protein_id] = protein_file

    print(f"Prepared {len(prepared)} proteins")
    return prepared


# ---------------------------------------------------------------------------
# 2. Ligand / conformer preparation
# ---------------------------------------------------------------------------

def prepare_conformers(compounds_ds, output_folder, config):
    """
    Prepare ligand conformers for docking.

    Three modes:
    - generate_conformers=True: Generate with RDKit ETKDGv3, filter by energy window
    - Pre-computed conformers with optional energy table: Load as-is
    - Single ligand: Use directly as one conformer

    Returns:
        list of dicts: [{ligand_id, conformer_id, sdf_file, energy}, ...]
    """
    conf_dir = os.path.join(output_folder, "conformers")
    os.makedirs(conf_dir, exist_ok=True)

    generate = config.get("generate_conformers", False)
    num_conformers = config.get("num_conformers", 50)
    energy_window = config.get("energy_window", 2.0)
    conformer_energies_ref = config.get("conformer_energies_ref")

    all_conformers = []

    for ligand_id, ligand_file in iterate_files(compounds_ds):
        suppl = Chem.SDMolSupplier(ligand_file, removeHs=False)
        mols = [m for m in suppl if m is not None]

        if not mols:
            print(f"Warning: No valid molecules in {ligand_file}, skipping {ligand_id}")
            continue

        if generate:
            # Generate conformers with RDKit
            mol = mols[0]
            mol = Chem.AddHs(mol)

            params = AllChem.ETKDGv3()
            params.randomSeed = 42
            cids = AllChem.EmbedMultipleConfs(mol, numConfs=num_conformers, params=params)

            if len(cids) == 0:
                print(f"Warning: Could not generate conformers for {ligand_id}")
                continue

            # Optimize and get energies
            results = AllChem.MMFFOptimizeMoleculeConfs(mol)
            energies = []
            valid_cids = []
            for cid_idx, (converged, energy) in enumerate(results):
                if converged == 0:  # 0 = converged successfully
                    energies.append(energy)
                    valid_cids.append(cids[cid_idx])
                else:
                    # Try to get energy even if not fully converged
                    ff = AllChem.MMFFGetMoleculeForceField(mol, AllChem.MMFFGetMoleculeProperties(mol), confId=cids[cid_idx])
                    if ff is not None:
                        energies.append(ff.CalcEnergy())
                        valid_cids.append(cids[cid_idx])

            if not energies:
                print(f"Warning: No valid conformer energies for {ligand_id}")
                continue

            # Filter by energy window relative to minimum
            min_energy = min(energies)
            for i, (cid, energy) in enumerate(zip(valid_cids, energies)):
                rel_energy = energy - min_energy
                if rel_energy <= energy_window:
                    conf_sdf = os.path.join(conf_dir, f"{ligand_id}_conf{i}.sdf")
                    writer = Chem.SDWriter(conf_sdf)
                    writer.write(mol, confId=cid)
                    writer.close()
                    all_conformers.append({
                        "ligand_id": ligand_id,
                        "conformer_id": i,
                        "sdf_file": conf_sdf,
                        "energy": rel_energy,
                    })

            print(f"Generated {sum(1 for c in all_conformers if c['ligand_id'] == ligand_id)} "
                  f"conformers for {ligand_id} (from {len(cids)} initial)")

        elif len(mols) > 1:
            # Multi-molecule SDF — each mol is a conformer
            energies = _load_external_energies(ligand_id, conformer_energies_ref) if conformer_energies_ref else None

            for i, mol in enumerate(mols):
                energy = ""
                if energies is not None and i < len(energies):
                    energy = energies[i]

                # Apply energy window filter if we have energies
                if energies is not None and energy != "":
                    min_e = min(e for e in energies if e != "")
                    if (energy - min_e) > energy_window:
                        continue

                conf_sdf = os.path.join(conf_dir, f"{ligand_id}_conf{i}.sdf")
                writer = Chem.SDWriter(conf_sdf)
                writer.write(mol)
                writer.close()
                all_conformers.append({
                    "ligand_id": ligand_id,
                    "conformer_id": i,
                    "sdf_file": conf_sdf,
                    "energy": energy,
                })
        else:
            # Single ligand — one "conformer"
            all_conformers.append({
                "ligand_id": ligand_id,
                "conformer_id": 0,
                "sdf_file": ligand_file,
                "energy": "",
            })

    # Write conformer energies CSV
    if all_conformers:
        energies_csv = os.path.join(conf_dir, "conformer_energies.csv")
        pd.DataFrame(all_conformers).to_csv(energies_csv, index=False)
        print(f"Total conformers prepared: {len(all_conformers)}")

    return all_conformers


def _load_external_energies(ligand_id, conformer_energies_ref):
    """Load pre-computed conformer energies from a table reference."""
    if not conformer_energies_ref or not conformer_energies_ref.startswith("DATASHEET_REFERENCE:"):
        return None

    from biopipelines_io import load_table
    try:
        table, column = load_table(conformer_energies_ref)
        if column and column in table.columns:
            return table[column].tolist()
    except Exception as e:
        print(f"Warning: Could not load conformer energies: {e}")

    return None


# ---------------------------------------------------------------------------
# 3. Docking execution
# ---------------------------------------------------------------------------

def run_docking(prepared_proteins, conformers, config):
    """
    Run GNINA docking for each (protein, conformer) pair.

    For each pair, runs GNINA num_runs times with incremented seeds.
    Parses output SDF files for Vina and CNN scores.

    Returns:
        list of dicts: all pose records
    """
    output_folder = config["output_folder"]
    gnina_binary = config["gnina_binary"]
    exhaustiveness = config["exhaustiveness"]
    num_modes = config["num_modes"]
    num_runs = config["num_runs"]
    seed = config["seed"]
    cnn_scoring = config["cnn_scoring"]
    cnn_score_threshold = config["cnn_score_threshold"]
    box = config["box"]

    dock_dir = os.path.join(output_folder, "docking")
    os.makedirs(dock_dir, exist_ok=True)

    # Build box arguments
    box_args = _build_box_args(box)

    all_poses = []
    total_pairs = len(prepared_proteins) * len(conformers)
    pair_idx = 0

    for protein_id, protein_file in prepared_proteins.items():
        for conf in conformers:
            pair_idx += 1
            ligand_id = conf["ligand_id"]
            conformer_id = conf["conformer_id"]
            conf_sdf = conf["sdf_file"]

            print(f"Docking pair {pair_idx}/{total_pairs}: "
                  f"{protein_id} + {ligand_id}_conf{conformer_id}")

            for run_idx in range(num_runs):
                out_sdf = os.path.join(
                    dock_dir,
                    f"{protein_id}_{ligand_id}_conf{conformer_id}_run{run_idx}.sdf"
                )

                cmd = [
                    gnina_binary,
                    "-r", protein_file,
                    "-l", conf_sdf,
                ] + box_args + [
                    "--exhaustiveness", str(exhaustiveness),
                    "--num_modes", str(num_modes),
                    "--cnn_scoring", cnn_scoring,
                    "--seed", str(seed + run_idx),
                    "-o", out_sdf,
                ]

                result = subprocess.run(cmd, capture_output=True, text=True)
                if result.returncode != 0:
                    print(f"Warning: GNINA failed for {protein_id}_{ligand_id}_conf{conformer_id}_run{run_idx}")
                    print(f"  stderr: {result.stderr[:500]}")
                    continue

                # Parse output SDF
                poses = _parse_gnina_output(
                    out_sdf, protein_id, ligand_id, conformer_id,
                    run_idx, cnn_score_threshold
                )
                all_poses.extend(poses)

    print(f"Docking complete: {len(all_poses)} accepted poses")
    return all_poses


def _build_box_args(box):
    """Build GNINA command-line arguments for the binding box."""
    args = []
    if "autobox_ligand" in box:
        args.extend(["--autobox_ligand", box["autobox_ligand"]])
        args.extend(["--autobox_add", str(box.get("autobox_add", 4.0))])
    elif "center" in box and "size" in box:
        cx, cy, cz = [s.strip() for s in str(box["center"]).split(",")]
        args.extend(["--center_x", cx, "--center_y", cy, "--center_z", cz])

        size = box["size"]
        if isinstance(size, (int, float)):
            args.extend(["--size_x", str(size), "--size_y", str(size), "--size_z", str(size)])
        else:
            sx, sy, sz = [s.strip() for s in str(size).split(",")]
            args.extend(["--size_x", sx, "--size_y", sy, "--size_z", sz])
    return args


def _parse_gnina_output(sdf_path, protein_id, ligand_id, conformer_id,
                        run_idx, cnn_score_threshold):
    """Parse a GNINA output SDF and return pose records."""
    poses = []
    if not os.path.exists(sdf_path):
        return poses

    suppl = Chem.SDMolSupplier(sdf_path, removeHs=False)
    for pose_idx, mol in enumerate(suppl):
        if mol is None:
            continue

        vina_score = _get_prop(mol, "minimizedAffinity")
        cnn_score = _get_prop(mol, "CNNscore")
        cnn_affinity = _get_prop(mol, "CNNaffinity")

        # Filter by CNN score threshold
        if cnn_score is not None and cnn_score < cnn_score_threshold:
            continue

        pose_id = f"{protein_id}_{ligand_id}_conf{conformer_id}_run{run_idx}_pose{pose_idx}"
        poses.append({
            "id": pose_id,
            "protein_id": protein_id,
            "ligand_id": ligand_id,
            "conformer_id": conformer_id,
            "run": run_idx,
            "pose": pose_idx,
            "vina_score": vina_score if vina_score is not None else "",
            "cnn_score": cnn_score if cnn_score is not None else "",
            "cnn_affinity": cnn_affinity if cnn_affinity is not None else "",
        })

    return poses


def _get_prop(mol, prop_name):
    """Safely get a float property from an RDKit molecule."""
    try:
        if mol.HasProp(prop_name):
            return float(mol.GetProp(prop_name))
    except (ValueError, RuntimeError):
        pass
    return None


# ---------------------------------------------------------------------------
# 4. Pose consistency analysis
# ---------------------------------------------------------------------------

def analyze_pose_consistency(all_poses, config):
    """
    Analyze pose consistency across runs for each (protein, ligand, conformer) group.

    For each group:
    - Collect the best pose (by Vina score) from each run
    - Compute pairwise RMSD between best poses
    - Cluster by RMSD using single-linkage
    - pose_consistency = |largest_cluster| / num_runs

    Returns:
        dict: (protein_id, ligand_id, conformer_id) -> consistency_score
    """
    rmsd_threshold = config["rmsd_threshold"]
    dock_dir = os.path.join(config["output_folder"], "docking")

    # Group poses by (protein, ligand, conformer)
    groups = defaultdict(list)
    for pose in all_poses:
        key = (pose["protein_id"], pose["ligand_id"], pose["conformer_id"])
        groups[key].append(pose)

    consistency = {}
    for key, poses in groups.items():
        protein_id, ligand_id, conformer_id = key

        # Get best pose per run (lowest Vina score)
        best_per_run = {}
        for pose in poses:
            run = pose["run"]
            vina = pose["vina_score"]
            if vina == "":
                continue
            if run not in best_per_run or vina < best_per_run[run]["vina_score"]:
                best_per_run[run] = pose

        if len(best_per_run) < 2:
            consistency[key] = 1.0 if len(best_per_run) == 1 else 0.0
            continue

        # Load best pose molecules for RMSD calculation
        best_mols = {}
        for run, pose in best_per_run.items():
            sdf_file = os.path.join(
                dock_dir,
                f"{protein_id}_{ligand_id}_conf{conformer_id}_run{run}.sdf"
            )
            if not os.path.exists(sdf_file):
                continue
            suppl = Chem.SDMolSupplier(sdf_file, removeHs=False)
            pose_rank = pose["pose"]
            for idx, mol in enumerate(suppl):
                if idx == pose_rank and mol is not None:
                    best_mols[run] = mol
                    break

        runs = sorted(best_mols.keys())
        if len(runs) < 2:
            consistency[key] = 1.0 if len(runs) == 1 else 0.0
            continue

        # Compute pairwise RMSD
        n = len(runs)
        rmsd_matrix = np.zeros((n, n))
        for i in range(n):
            for j in range(i + 1, n):
                try:
                    rmsd = rdMolAlign.GetBestRMS(best_mols[runs[i]], best_mols[runs[j]])
                except Exception:
                    rmsd = float('inf')
                rmsd_matrix[i][j] = rmsd
                rmsd_matrix[j][i] = rmsd

        # Single-linkage clustering
        clusters = _single_linkage_cluster(rmsd_matrix, rmsd_threshold)
        largest_cluster_size = max(len(c) for c in clusters)
        consistency[key] = largest_cluster_size / len(runs)

    return consistency


def _single_linkage_cluster(rmsd_matrix, threshold):
    """Single-linkage clustering on an RMSD distance matrix."""
    n = rmsd_matrix.shape[0]
    visited = [False] * n
    clusters = []

    for i in range(n):
        if visited[i]:
            continue
        cluster = []
        stack = [i]
        while stack:
            node = stack.pop()
            if visited[node]:
                continue
            visited[node] = True
            cluster.append(node)
            for j in range(n):
                if not visited[j] and rmsd_matrix[node][j] < threshold:
                    stack.append(j)
        clusters.append(cluster)

    return clusters


# ---------------------------------------------------------------------------
# 5. Results aggregation
# ---------------------------------------------------------------------------

def aggregate_results(all_poses, conformers, consistency, config):
    """
    Aggregate results into docking_results.csv and conformer_ranking.csv.

    docking_results.csv: All accepted poses with scores.
    conformer_ranking.csv: Per-conformer aggregated statistics.
    """
    output_folder = config["output_folder"]
    docking_results_csv = config["docking_results_csv"]
    conformer_ranking_csv = config["conformer_ranking_csv"]

    # Write docking_results.csv
    if all_poses:
        results_df = pd.DataFrame(all_poses)
        results_df.to_csv(docking_results_csv, index=False)
        print(f"Wrote {len(results_df)} poses to {docking_results_csv}")
    else:
        results_df = pd.DataFrame(columns=[
            "id", "protein_id", "ligand_id", "conformer_id",
            "run", "pose", "vina_score", "cnn_score", "cnn_affinity"
        ])
        results_df.to_csv(docking_results_csv, index=False)
        print("Warning: No poses passed filtering")

    # Build conformer energy lookup
    energy_lookup = {}
    for conf in conformers:
        key = (conf["ligand_id"], conf["conformer_id"])
        energy_lookup[key] = conf.get("energy", "")

    # Build best_poses directory
    best_poses_dir = os.path.join(output_folder, "best_poses")
    os.makedirs(best_poses_dir, exist_ok=True)
    dock_dir = os.path.join(output_folder, "docking")

    # Group poses for conformer ranking
    groups = defaultdict(list)
    for pose in all_poses:
        key = (pose["protein_id"], pose["ligand_id"], pose["conformer_id"])
        groups[key].append(pose)

    ranking_rows = []
    for key, poses in groups.items():
        protein_id, ligand_id, conformer_id = key

        # Get numeric Vina scores
        vina_scores = [p["vina_score"] for p in poses if p["vina_score"] != ""]
        cnn_scores = [p["cnn_score"] for p in poses if p["cnn_score"] != ""]

        # Best pose per run (by Vina)
        best_per_run = {}
        for pose in poses:
            run = pose["run"]
            vina = pose["vina_score"]
            if vina == "":
                continue
            if run not in best_per_run or vina < best_per_run[run]["vina_score"]:
                best_per_run[run] = pose

        run_vinas = [best_per_run[r]["vina_score"] for r in best_per_run]

        best_vina = min(vina_scores) if vina_scores else ""
        mean_vina = float(np.mean(run_vinas)) if run_vinas else ""
        std_vina = float(np.std(run_vinas)) if len(run_vinas) > 1 else ""
        best_cnn = max(cnn_scores) if cnn_scores else ""
        pose_cons = consistency.get(key, "")

        conf_energy = energy_lookup.get((ligand_id, conformer_id), "")

        # Pseudo binding energy
        pseudo_be = ""
        if best_vina != "" and conf_energy != "":
            try:
                pseudo_be = float(best_vina) + float(conf_energy)
            except (ValueError, TypeError):
                pseudo_be = ""

        # Copy best pose SDF
        best_pose_file = ""
        if vina_scores:
            # Find the overall best pose
            best_pose = min(
                [p for p in poses if p["vina_score"] != ""],
                key=lambda p: p["vina_score"]
            )
            src_sdf = os.path.join(
                dock_dir,
                f"{protein_id}_{ligand_id}_conf{conformer_id}_run{best_pose['run']}.sdf"
            )
            dst_sdf = os.path.join(
                best_poses_dir,
                f"{protein_id}_{ligand_id}_conf{conformer_id}_best.sdf"
            )
            if os.path.exists(src_sdf):
                _extract_pose_to_file(src_sdf, best_pose["pose"], dst_sdf)
                best_pose_file = dst_sdf

        ranking_rows.append({
            "id": f"{protein_id}_{ligand_id}_conf{conformer_id}",
            "protein_id": protein_id,
            "ligand_id": ligand_id,
            "conformer_id": conformer_id,
            "best_vina": best_vina,
            "mean_vina": mean_vina,
            "std_vina": std_vina,
            "best_cnn_score": best_cnn,
            "pose_consistency": pose_cons,
            "conformer_energy": conf_energy,
            "pseudo_binding_energy": pseudo_be,
            "best_pose_file": best_pose_file,
        })

    ranking_df = pd.DataFrame(ranking_rows)
    if not ranking_df.empty and "best_vina" in ranking_df.columns:
        # Sort by best_vina (lower is better), handling empty strings
        ranking_df["_sort_key"] = pd.to_numeric(ranking_df["best_vina"], errors="coerce")
        ranking_df = ranking_df.sort_values("_sort_key").drop(columns=["_sort_key"]).reset_index(drop=True)

    ranking_df.to_csv(conformer_ranking_csv, index=False)
    print(f"Wrote {len(ranking_df)} conformer rankings to {conformer_ranking_csv}")


def _extract_pose_to_file(src_sdf, pose_index, dst_sdf):
    """Extract a specific pose from a multi-molecule SDF to a new file."""
    suppl = Chem.SDMolSupplier(src_sdf, removeHs=False)
    for idx, mol in enumerate(suppl):
        if idx == pose_index and mol is not None:
            writer = Chem.SDWriter(dst_sdf)
            writer.write(mol)
            writer.close()
            return
    # Fallback: copy first molecule if exact pose not found
    suppl = Chem.SDMolSupplier(src_sdf, removeHs=False)
    for mol in suppl:
        if mol is not None:
            writer = Chem.SDWriter(dst_sdf)
            writer.write(mol)
            writer.close()
            return


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    """Main execution function."""
    if len(sys.argv) != 2:
        print("Usage: pipe_gnina.py <config_json>")
        sys.exit(1)

    config_path = sys.argv[1]

    with open(config_path, 'r') as f:
        config = json.load(f)

    print("=" * 60)
    print("GNINA Docking Pipeline")
    print("=" * 60)

    # Load DataStreams
    structures_ds = load_datastream(config["structures_json"])
    compounds_ds = load_datastream(config["compounds_json"])

    # 1. Protein preparation
    print("\n--- Step 1: Protein Preparation ---")
    prepared_proteins = prepare_proteins(
        structures_ds, config["output_folder"],
        config["protonate"], config["pH"]
    )

    # 2. Conformer preparation
    print("\n--- Step 2: Conformer Preparation ---")
    conformers = prepare_conformers(
        compounds_ds, config["output_folder"], config
    )

    if not conformers:
        print("Error: No conformers prepared, cannot proceed with docking")
        sys.exit(1)

    # 3. Docking execution
    print("\n--- Step 3: Docking Execution ---")
    all_poses = run_docking(prepared_proteins, conformers, config)

    # 4. Pose consistency analysis
    print("\n--- Step 4: Pose Consistency Analysis ---")
    consistency = analyze_pose_consistency(all_poses, config)

    # 5. Results aggregation
    print("\n--- Step 5: Results Aggregation ---")
    aggregate_results(all_poses, conformers, consistency, config)

    print("\n" + "=" * 60)
    print("GNINA docking pipeline complete")
    print("=" * 60)


if __name__ == "__main__":
    main()
