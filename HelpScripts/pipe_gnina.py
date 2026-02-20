#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
GNINA docking pipeline helper script.

Orchestrates the full docking workflow:
1. Protein preparation (clean ligands/duplicates, protonate with OpenBabel)
2. Ligand/conformer preparation (RDKit conformer generation or loading)
3. Docking execution (GNINA binary calls)
4. Pose consistency analysis (RMSD clustering across runs)
5. Results aggregation (docking_results.csv + conformer_ranking.csv)

Usage:
    python pipe_gnina.py <config_json>
"""

import json
import os
import shutil
import subprocess
import sys
from collections import defaultdict
from contextlib import contextmanager

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdMolAlign
from rdkit.ML.Cluster import Butina

from biopipelines_io import load_datastream, iterate_files, load_table


@contextmanager
def _suppress_rdkit_warnings():
    """Temporarily suppress RDKit warnings (e.g. 2D/3D flag from GNINA output)."""
    logger = RDLogger.logger()
    logger.setLevel(RDLogger.ERROR)
    try:
        yield
    finally:
        logger.setLevel(RDLogger.WARNING)


# ---------------------------------------------------------------------------
# 1. Protein preparation
# ---------------------------------------------------------------------------

def _parse_pdb_coords(line):
    """
    Extract x, y, z from a PDB ATOM/HETATM line.

    Handles column shifts caused by non-standard residue names (e.g. 5-char
    CCD codes that overflow the 3-char PDB residue name field).
    """
    # Standard PDB columns: x=30:38, y=38:46, z=46:54
    try:
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        return x, y, z
    except (ValueError, IndexError):
        pass
    # Shifted by +2 (5-char residue name instead of 3)
    try:
        x = float(line[32:40])
        y = float(line[40:48])
        z = float(line[48:56])
        return x, y, z
    except (ValueError, IndexError):
        pass
    return None


def _sequences_are_similar(seq_a, seq_b, threshold=0.8):
    """
    Check if two residue-name lists represent the same protein chain.

    Slides the shorter sequence along the longer, allowing partial overlap
    at both ends (homodimer chains often differ at partially-resolved termini).
    Returns True if the best ungapped alignment has >= threshold identity
    over the overlapping region.
    """
    if not seq_a or not seq_b:
        return not seq_a and not seq_b

    shorter, longer = (seq_a, seq_b) if len(seq_a) <= len(seq_b) else (seq_b, seq_a)
    min_overlap = int(len(shorter) * 0.5)
    best_identity = 0.0
    for offset in range(len(longer)):
        overlap = min(len(shorter), len(longer) - offset)
        if overlap < min_overlap:
            break
        matches = sum(1 for i in range(overlap) if shorter[i] == longer[offset + i])
        best_identity = max(best_identity, matches / overlap)
    return best_identity >= threshold


def _clean_protein_pdb(input_pdb, output_pdb, crystal_ligand_pdb=None):
    """
    Clean a PDB file for docking.

    Steps:
    1. Strip all HETATM records (ligands, waters, ions).
    2. Detect chains with identical sequences via sliding-window alignment
       and keep only the first (longest) representative of each unique chain.
    3. Write a clean PDB containing only protein ATOM records.
    4. Extract crystal ligand HETATM records closest to the kept chains
       (used for autoboxing when no explicit box is provided).

    Returns the path to the cleaned PDB. Also writes crystal_ligand_pdb
    if non-water HETATM records were found near the kept chains.
    """
    chains = {}
    chain_seqs = {}
    chain_ca_coords = {}
    hetatm_lines = []

    with open(input_pdb, 'r') as f:
        for line in f:
            if line.startswith("ATOM"):
                chain_id = line[21]
                chains.setdefault(chain_id, []).append(line)
                if line[12:16].strip() == "CA":
                    resname = line[17:20].strip()
                    chain_seqs.setdefault(chain_id, []).append(resname)
                    coords = _parse_pdb_coords(line)
                    if coords:
                        chain_ca_coords.setdefault(chain_id, []).append(coords)
            elif line.startswith("HETATM"):
                resname = line[17:20].strip()
                if resname != "HOH":
                    hetatm_lines.append(line)

    if not chains:
        print(f"  Warning: No ATOM records found in {input_pdb}, using file as-is")
        shutil.copy(input_pdb, output_pdb)
        return output_pdb

    # Identify unique chains — homodimers may differ at the termini
    keep_chains = []
    duplicate_chains = []
    for chain_id in sorted(chains.keys()):
        seq = chain_seqs.get(chain_id, [])
        is_dup = False
        for i, kept_id in enumerate(keep_chains):
            kept_seq = chain_seqs.get(kept_id, [])
            if _sequences_are_similar(seq, kept_seq):
                if len(seq) > len(kept_seq):
                    duplicate_chains.append((kept_id, chain_id))
                    keep_chains[i] = chain_id
                else:
                    duplicate_chains.append((chain_id, kept_id))
                is_dup = True
                break
        if not is_dup:
            keep_chains.append(chain_id)

    for dup_id, kept_id in duplicate_chains:
        print(f"  Removing duplicate chain {dup_id} (same sequence as chain {kept_id})")

    print(f"  Keeping chain(s): {', '.join(keep_chains)} "
          f"({sum(len(chains[c]) for c in keep_chains)} atoms)")

    with open(output_pdb, 'w') as f:
        for chain_id in keep_chains:
            for line in chains[chain_id]:
                f.write(line)
            f.write("TER\n")
        f.write("END\n")

    # Extract crystal ligands closest to the kept chains (for autoboxing).
    # HETATM chain IDs can be unreliable when CCD codes overflow the residue
    # name field, so we use spatial proximity to chain centroids instead.
    if crystal_ligand_pdb and hetatm_lines:
        kept_centroids = {}
        for cid in keep_chains:
            coords = chain_ca_coords.get(cid, [])
            if coords:
                kept_centroids[cid] = (
                    sum(c[0] for c in coords) / len(coords),
                    sum(c[1] for c in coords) / len(coords),
                    sum(c[2] for c in coords) / len(coords),
                )

        removed_centroids = {}
        for dup_id, _ in duplicate_chains:
            coords = chain_ca_coords.get(dup_id, [])
            if coords:
                removed_centroids[dup_id] = (
                    sum(c[0] for c in coords) / len(coords),
                    sum(c[1] for c in coords) / len(coords),
                    sum(c[2] for c in coords) / len(coords),
                )

        ligand_lines_kept = []
        for line in hetatm_lines:
            coords = _parse_pdb_coords(line)
            if coords is None:
                continue
            x, y, z = coords
            best_dist = float('inf')
            best_is_kept = False
            for cid, (cx, cy, cz) in kept_centroids.items():
                d = (x - cx)**2 + (y - cy)**2 + (z - cz)**2
                if d < best_dist:
                    best_dist = d
                    best_is_kept = True
            for cid, (cx, cy, cz) in removed_centroids.items():
                d = (x - cx)**2 + (y - cy)**2 + (z - cz)**2
                if d < best_dist:
                    best_dist = d
                    best_is_kept = False
            if best_is_kept:
                ligand_lines_kept.append(line)

        if ligand_lines_kept:
            with open(crystal_ligand_pdb, 'w') as f:
                for line in ligand_lines_kept:
                    f.write(line)
                f.write("END\n")
            print(f"  Extracted crystal ligand ({len(ligand_lines_kept)} atoms) "
                  f"-> {os.path.basename(crystal_ligand_pdb)}")
        else:
            print("  No crystal ligand found near kept chains")
    elif crystal_ligand_pdb and not hetatm_lines:
        print("  No crystal ligand found in PDB")

    return output_pdb


def prepare_proteins(structures_ds, output_folder, protonate, pH):
    """
    Prepare protein structures for docking.

    For each protein:
    1. Clean: strip HETATM records, deduplicate homodimer chains.
    2. Extract crystal ligand HETATM records for the kept chains (autobox).
    3. Protonate with OpenBabel at the specified pH (if requested).

    Returns:
        dict: protein_id -> {"protein": path, "crystal_ligand": path_or_None}
    """
    prep_dir = os.path.join(output_folder, "prepared_proteins")
    os.makedirs(prep_dir, exist_ok=True)

    prepared = {}
    for protein_id, protein_file in iterate_files(structures_ds):
        if not os.path.exists(protein_file):
            print(f"  Warning: protein file not found, skipping: {protein_file}")
            continue
        print(f"Preparing {protein_id}: {protein_file}")

        cleaned_pdb = os.path.join(prep_dir, f"{protein_id}_clean.pdb")
        crystal_lig = os.path.join(prep_dir, f"{protein_id}_crystal_ligand.pdb")
        _clean_protein_pdb(protein_file, cleaned_pdb, crystal_ligand_pdb=crystal_lig)
        crystal_lig = crystal_lig if os.path.exists(crystal_lig) else None

        if protonate:
            out_pdb = os.path.join(prep_dir, f"{protein_id}_h.pdb")
            cmd = ["obabel", cleaned_pdb, "-O", out_pdb, "-h", "-p", str(pH)]
            print(f"  Protonating: {' '.join(cmd)}")
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode != 0:
                print(f"  Warning: obabel failed for {protein_id}: {result.stderr}")
                prepared[protein_id] = {"protein": cleaned_pdb, "crystal_ligand": crystal_lig}
            else:
                prepared[protein_id] = {"protein": out_pdb, "crystal_ligand": crystal_lig}
        else:
            prepared[protein_id] = {"protein": cleaned_pdb, "crystal_ligand": crystal_lig}

    print(f"Prepared {len(prepared)} proteins")
    return prepared


# ---------------------------------------------------------------------------
# 2. Ligand / conformer preparation
# ---------------------------------------------------------------------------

def _load_compound_mols(compounds_ds, conf_dir):
    """
    Load compound molecules from any upstream tool format.

    Handles three DataStream formats:
    - Ligand tool: CSV map_table with file_path column (PDB/SDF files).
    - CompoundLibrary: CSV map_table with smiles column (no files).
    - Direct SDF files: files list points to actual SDF files.

    Returns:
        list of (ligand_id, mol, sdf_path) tuples.
    """
    results = []

    map_df = None
    if compounds_ds.map_table and os.path.exists(compounds_ds.map_table):
        map_df = pd.read_csv(compounds_ds.map_table)

    for ligand_id in compounds_ds.ids:
        mol = None
        sdf_path = None

        # Strategy 1: map_table has file_path column
        if map_df is not None and "file_path" in map_df.columns:
            row = map_df[map_df["id"] == ligand_id]
            if not row.empty:
                file_path = str(row.iloc[0]["file_path"])
                if file_path and os.path.exists(file_path):
                    mol, sdf_path = _load_mol_from_file(file_path, ligand_id, conf_dir)

        # Strategy 2: map_table has smiles column
        if mol is None and map_df is not None and "smiles" in map_df.columns:
            row = map_df[map_df["id"] == ligand_id]
            if not row.empty:
                smiles = str(row.iloc[0]["smiles"])
                if smiles and smiles != "nan":
                    mol, sdf_path = _mol_from_smiles(smiles, ligand_id, conf_dir)

        # Strategy 3: DataStream files list
        if mol is None and compounds_ds.files:
            idx = compounds_ds.ids.index(ligand_id) if ligand_id in compounds_ds.ids else -1
            if idx >= 0:
                if len(compounds_ds.files) == len(compounds_ds.ids):
                    file_path = compounds_ds.files[idx]
                elif len(compounds_ds.files) == 1:
                    file_path = compounds_ds.files[0]
                else:
                    file_path = None
                if file_path and os.path.exists(file_path) and not file_path.endswith(".csv"):
                    mol, sdf_path = _load_mol_from_file(file_path, ligand_id, conf_dir)

        if mol is not None and sdf_path is not None:
            results.append((ligand_id, mol, sdf_path))
        else:
            print(f"Warning: Could not load molecule for {ligand_id}")

    return results


def _load_mol_from_file(file_path, ligand_id, conf_dir):
    """Load an RDKit mol from a PDB or SDF file, converting to SDF if needed."""
    ext = os.path.splitext(file_path)[1].lower()

    if ext in (".sdf", ".mol"):
        suppl = Chem.SDMolSupplier(file_path, removeHs=False)
        for mol in suppl:
            if mol is not None:
                return mol, file_path
        return None, None

    elif ext == ".pdb":
        mol = Chem.MolFromPDBFile(file_path, removeHs=False, sanitize=True)
        if mol is None:
            mol = Chem.MolFromPDBFile(file_path, removeHs=False, sanitize=False)
        if mol is not None:
            if mol.GetNumConformers() > 0:
                mol.GetConformer().Set3D(True)
            sdf_path = os.path.join(conf_dir, f"{ligand_id}.sdf")
            writer = Chem.SDWriter(sdf_path)
            writer.write(mol)
            writer.close()
            return mol, sdf_path
        return None, None

    elif ext == ".mol2":
        mol = Chem.MolFromMol2File(file_path, removeHs=False)
        if mol is not None:
            sdf_path = os.path.join(conf_dir, f"{ligand_id}.sdf")
            writer = Chem.SDWriter(sdf_path)
            writer.write(mol)
            writer.close()
            return mol, sdf_path
        return None, None

    else:
        print(f"Warning: Unsupported file format '{ext}' for {file_path}")
        return None, None


def _mol_from_smiles(smiles, ligand_id, conf_dir):
    """Generate a 3D molecule from SMILES and write to SDF."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"Warning: Could not parse SMILES for {ligand_id}: {smiles}")
        return None, None

    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    status = AllChem.EmbedMolecule(mol, params)
    if status == -1:
        AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    if mol.GetNumConformers() == 0:
        print(f"Warning: Could not embed 3D coordinates for {ligand_id}")
        return None, None

    try:
        AllChem.MMFFOptimizeMolecule(mol)
    except Exception:
        pass

    if mol.GetNumConformers() > 0:
        mol.GetConformer().Set3D(True)

    sdf_path = os.path.join(conf_dir, f"{ligand_id}.sdf")
    writer = Chem.SDWriter(sdf_path)
    writer.write(mol)
    writer.close()
    return mol, sdf_path


def _cluster_conformers(mol, cids, energies, rmsd_cutoff):
    """
    Cluster conformers by heavy-atom RMSD using the Butina algorithm.

    Uses GetBestRMS for symmetry-aware RMSD: tries all automorphism mappings
    so that conformers related by molecular symmetry (e.g. C2 rotation of a
    xanthene core in rhodamine dyes) are correctly identified as equivalent.

    Args:
        mol: RDKit molecule with multiple conformers.
        cids: list of conformer IDs to cluster.
        energies: list of MMFF94 energies corresponding to cids.
        rmsd_cutoff: RMSD threshold in Angstroms.

    Returns:
        list of (cid, energy) tuples for cluster representatives.
    """
    n = len(cids)
    if n <= 1:
        return list(zip(cids, energies))

    mol_noH = Chem.RemoveHs(mol)

    # Condensed pairwise distance matrix expected by Butina:
    # [d(0,1), d(0,2), d(1,2), d(0,3), d(1,3), d(2,3), ...]
    dists = []
    for i in range(n):
        for j in range(i):
            rmsd = AllChem.GetBestRMS(mol_noH, mol_noH, prbId=cids[i], refId=cids[j])
            dists.append(rmsd)

    clusters = Butina.ClusterData(dists, n, rmsd_cutoff, isDistData=True)

    representatives = []
    for cluster in clusters:
        best_idx = min(cluster, key=lambda idx: energies[idx])
        representatives.append((cids[best_idx], energies[best_idx]))

    return representatives


def prepare_conformers(compounds_ds, output_folder, config):
    """
    Prepare ligand conformers for docking.

    With generate_conformers=True: generates RDKit ETKDGv3 conformers, filters
    by energy window, clusters by RMSD (Butina + symmetry-aware GetBestRMS).
    Otherwise: loads existing conformers from the input file directly.

    Returns:
        list of dicts: [{ligand_id, conformer_id, sdf_file, energy}, ...]
    """
    conf_dir = os.path.join(output_folder, "conformers")
    os.makedirs(conf_dir, exist_ok=True)

    generate = config.get("generate_conformers", False)
    num_conformers = config.get("num_conformers", 50)
    energy_window = config.get("energy_window", 2.0)
    conformer_rmsd = config.get("conformer_rmsd", 1.0)
    conformer_energies_ref = config.get("conformer_energies_ref")

    compound_mols = _load_compound_mols(compounds_ds, conf_dir)
    print(f"Loaded {len(compound_mols)} compounds from upstream")

    all_conformers = []

    for ligand_id, mol, sdf_path in compound_mols:

        if generate:
            mol = Chem.AddHs(mol)
            params = AllChem.ETKDGv3()
            params.randomSeed = 42
            cids = AllChem.EmbedMultipleConfs(mol, numConfs=num_conformers, params=params)

            if len(cids) == 0:
                print(f"Warning: Could not generate conformers for {ligand_id}")
                continue

            results = AllChem.MMFFOptimizeMoleculeConfs(mol)
            energies = []
            valid_cids = []
            for cid_idx, (converged, energy) in enumerate(results):
                if converged == 0:
                    energies.append(energy)
                    valid_cids.append(cids[cid_idx])
                else:
                    props = AllChem.MMFFGetMoleculeProperties(mol)
                    if props is not None:
                        ff = AllChem.MMFFGetMoleculeForceField(mol, props, confId=cids[cid_idx])
                        if ff is not None:
                            energies.append(ff.CalcEnergy())
                            valid_cids.append(cids[cid_idx])

            if not energies:
                print(f"Warning: No valid conformer energies for {ligand_id}")
                continue

            min_energy = min(energies)
            filtered_cids, filtered_energies = [], []
            for cid, energy in zip(valid_cids, energies):
                if (energy - min_energy) <= energy_window:
                    filtered_cids.append(cid)
                    filtered_energies.append(energy - min_energy)

            n_after_energy = len(filtered_cids)

            if conformer_rmsd > 0 and len(filtered_cids) > 1:
                representatives = _cluster_conformers(
                    mol, filtered_cids, filtered_energies, conformer_rmsd)
            else:
                representatives = list(zip(filtered_cids, filtered_energies))

            for cid, _ in representatives:
                mol.GetConformer(cid).Set3D(True)

            for i, (cid, energy) in enumerate(representatives):
                conf_sdf = os.path.join(conf_dir, f"{ligand_id}_conf{i}.sdf")
                writer = Chem.SDWriter(conf_sdf)
                writer.write(mol, confId=cid)
                writer.close()
                all_conformers.append({
                    "ligand_id": ligand_id,
                    "conformer_id": i,
                    "sdf_file": conf_sdf,
                    "energy": energy,
                })

            n_final = sum(1 for c in all_conformers if c['ligand_id'] == ligand_id)
            print(f"  {ligand_id}: {len(cids)} generated -> "
                  f"{n_after_energy} after energy filter -> "
                  f"{n_final} after RMSD clustering (cutoff {conformer_rmsd} Å)")

        else:
            if os.path.exists(sdf_path):
                suppl = Chem.SDMolSupplier(sdf_path, removeHs=False)
                mols_in_file = [m for m in suppl if m is not None]
            else:
                mols_in_file = [mol]

            if len(mols_in_file) > 1:
                energies = _load_external_energies(ligand_id, conformer_energies_ref) if conformer_energies_ref else None
                for i, conf_mol in enumerate(mols_in_file):
                    energy = ""
                    if energies is not None and i < len(energies):
                        energy = energies[i]
                    if energies is not None and energy != "":
                        min_e = min(e for e in energies if e != "")
                        if (energy - min_e) > energy_window:
                            continue
                    conf_sdf = os.path.join(conf_dir, f"{ligand_id}_conf{i}.sdf")
                    writer = Chem.SDWriter(conf_sdf)
                    writer.write(conf_mol)
                    writer.close()
                    all_conformers.append({
                        "ligand_id": ligand_id,
                        "conformer_id": i,
                        "sdf_file": conf_sdf,
                        "energy": energy,
                    })
            else:
                all_conformers.append({
                    "ligand_id": ligand_id,
                    "conformer_id": 0,
                    "sdf_file": sdf_path,
                    "energy": "",
                })

    if all_conformers:
        energies_csv = os.path.join(conf_dir, "conformer_energies.csv")
        pd.DataFrame(all_conformers).to_csv(energies_csv, index=False)
        print(f"Total conformers prepared: {len(all_conformers)}")

    return all_conformers


def _load_external_energies(ligand_id, conformer_energies_ref):
    """Load pre-computed conformer energies from a table reference."""
    if not conformer_energies_ref or not conformer_energies_ref.startswith("DATASHEET_REFERENCE:"):
        return None
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
        list of pose record dicts.
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

    explicit_box_args = _build_box_args(box)

    all_poses = []
    total_pairs = len(prepared_proteins) * len(conformers)
    pair_idx = 0

    for protein_id, prot_info in prepared_proteins.items():
        protein_file = prot_info["protein"]
        crystal_ligand = prot_info.get("crystal_ligand")

        if explicit_box_args:
            box_args = explicit_box_args
        elif crystal_ligand:
            print(f"  Using crystal ligand for autobox: {os.path.basename(crystal_ligand)}")
            box_args = ["--autobox_ligand", crystal_ligand,
                        "--autobox_add", str(box.get("autobox_add", 4.0))]
        else:
            print(f"  Warning: No box defined and no crystal ligand for {protein_id}")
            box_args = []

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
                    print(f"Warning: GNINA failed for "
                          f"{protein_id}_{ligand_id}_conf{conformer_id}_run{run_idx}")
                    print(f"  stderr: {result.stderr[:500]}")
                    continue

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
        args.extend(["--autobox_ligand", box["autobox_ligand"],
                     "--autobox_add", str(box.get("autobox_add", 4.0))])
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
    """Parse a GNINA output SDF and return pose records above the CNN threshold."""
    poses = []
    if not os.path.exists(sdf_path):
        return poses

    with _suppress_rdkit_warnings():
        suppl = Chem.SDMolSupplier(sdf_path, removeHs=False)
        mols = list(enumerate(suppl))

    for pose_idx, mol in mols:
        if mol is None:
            continue
        vina_score = _get_prop(mol, "minimizedAffinity")
        cnn_score = _get_prop(mol, "CNNscore")
        cnn_affinity = _get_prop(mol, "CNNaffinity")

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
    - Collect the best pose (lowest Vina score) from each run.
    - Compute pairwise RMSD between best poses.
    - Cluster by RMSD using single-linkage.
    - pose_consistency = largest_cluster_size / num_runs

    Returns:
        dict: (protein_id, ligand_id, conformer_id) -> consistency_score
    """
    rmsd_threshold = config["rmsd_threshold"]
    dock_dir = os.path.join(config["output_folder"], "docking")

    groups = defaultdict(list)
    for pose in all_poses:
        key = (pose["protein_id"], pose["ligand_id"], pose["conformer_id"])
        groups[key].append(pose)

    consistency = {}
    for key, poses in groups.items():
        protein_id, ligand_id, conformer_id = key

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

        best_mols = {}
        with _suppress_rdkit_warnings():
            for run, pose in best_per_run.items():
                sdf_file = os.path.join(
                    dock_dir,
                    f"{protein_id}_{ligand_id}_conf{conformer_id}_run{run}.sdf"
                )
                if not os.path.exists(sdf_file):
                    continue
                suppl = Chem.SDMolSupplier(sdf_file, removeHs=False)
                for idx, mol in enumerate(suppl):
                    if idx == pose["pose"] and mol is not None:
                        best_mols[run] = mol
                        break

        runs = sorted(best_mols.keys())
        if len(runs) < 2:
            consistency[key] = 1.0 if len(runs) == 1 else 0.0
            continue

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

def aggregate_results(all_poses, conformers, consistency, config, prepared_proteins):
    """
    Aggregate results into docking_results.csv and conformer_ranking.csv.

    Also writes best poses as combined protein+ligand PDB files under
    best_poses/ within the output folder.
    """
    output_folder = config["output_folder"]
    docking_results_csv = config["docking_results_csv"]
    conformer_ranking_csv = config["conformer_ranking_csv"]

    if all_poses:
        results_df = pd.DataFrame(all_poses)
        results_df.to_csv(docking_results_csv, index=False)
        print(f"Wrote {len(results_df)} poses to {docking_results_csv}")
    else:
        pd.DataFrame(columns=[
            "id", "protein_id", "ligand_id", "conformer_id",
            "run", "pose", "vina_score", "cnn_score", "cnn_affinity"
        ]).to_csv(docking_results_csv, index=False)
        print("Warning: No poses passed filtering")

    energy_lookup = {
        (c["ligand_id"], c["conformer_id"]): c.get("energy", "")
        for c in conformers
    }

    best_poses_dir = os.path.join(output_folder, "best_poses")
    os.makedirs(best_poses_dir, exist_ok=True)
    dock_dir = os.path.join(output_folder, "docking")

    groups = defaultdict(list)
    for pose in all_poses:
        key = (pose["protein_id"], pose["ligand_id"], pose["conformer_id"])
        groups[key].append(pose)

    ranking_rows = []
    for key, poses in groups.items():
        protein_id, ligand_id, conformer_id = key

        vina_scores = [p["vina_score"] for p in poses if p["vina_score"] != ""]
        cnn_scores = [p["cnn_score"] for p in poses if p["cnn_score"] != ""]

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

        pseudo_be = ""
        if best_vina != "" and conf_energy != "":
            try:
                pseudo_be = float(best_vina) + float(conf_energy)
            except (ValueError, TypeError):
                pass

        best_pose_file = ""
        if vina_scores:
            best_pose = min(
                [p for p in poses if p["vina_score"] != ""],
                key=lambda p: p["vina_score"]
            )
            src_sdf = os.path.join(
                dock_dir,
                f"{protein_id}_{ligand_id}_conf{conformer_id}_run{best_pose['run']}.sdf"
            )
            dst_pdb = os.path.join(
                best_poses_dir,
                f"{protein_id}_{ligand_id}_conf{conformer_id}_best.pdb"
            )
            prot_info = prepared_proteins.get(protein_id, {})
            protein_pdb = prot_info.get("protein") if isinstance(prot_info, dict) else prot_info
            if os.path.exists(src_sdf) and protein_pdb and os.path.exists(protein_pdb):
                _write_complex_pdb(protein_pdb, src_sdf, best_pose["pose"], dst_pdb)
                best_pose_file = dst_pdb

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
        ranking_df["_sort_key"] = pd.to_numeric(ranking_df["best_vina"], errors="coerce")
        ranking_df = (ranking_df.sort_values("_sort_key")
                                .drop(columns=["_sort_key"])
                                .reset_index(drop=True))

    ranking_df.to_csv(conformer_ranking_csv, index=False)
    print(f"Wrote {len(ranking_df)} conformer rankings to {conformer_ranking_csv}")


def _write_complex_pdb(protein_pdb, ligand_sdf, pose_index, dst_pdb):
    """
    Write a combined protein+ligand PDB file.

    Protein ATOM/TER records come first, followed by the docked ligand as
    HETATM records on chain Z with renumbered serial numbers and proper
    CONECT records for bond connectivity.
    """
    with _suppress_rdkit_warnings():
        ligand_mol = None
        suppl = Chem.SDMolSupplier(ligand_sdf, removeHs=False)
        for idx, mol in enumerate(suppl):
            if idx == pose_index and mol is not None:
                ligand_mol = mol
                break
        if ligand_mol is None:
            suppl = Chem.SDMolSupplier(ligand_sdf, removeHs=False)
            for mol in suppl:
                if mol is not None:
                    ligand_mol = mol
                    break

    if ligand_mol is None:
        print(f"  Warning: Could not extract ligand pose from {ligand_sdf}")
        return

    ligand_pdb_block = Chem.MolToPDBBlock(ligand_mol)
    if ligand_pdb_block is None:
        print("  Warning: Could not convert ligand to PDB format")
        return

    max_protein_serial = 0
    with open(protein_pdb, 'r') as prot:
        for line in prot:
            if line.startswith(("ATOM", "HETATM")):
                try:
                    serial = int(line[6:11].strip())
                    max_protein_serial = max(max_protein_serial, serial)
                except ValueError:
                    pass
    serial_offset = max_protein_serial

    ligand_atoms = []
    ligand_conects = []
    for line in ligand_pdb_block.splitlines():
        if line.startswith(("ATOM", "HETATM")):
            ligand_atoms.append(line)
        elif line.startswith("CONECT"):
            ligand_conects.append(line)

    with open(dst_pdb, 'w') as out:
        with open(protein_pdb, 'r') as prot:
            for line in prot:
                if line.startswith(("CONECT", "MASTER", "END")):
                    continue
                out.write(line)

        for line in ligand_atoms:
            old_serial = int(line[6:11].strip())
            new_serial = old_serial + serial_offset
            hetatm = f"HETATM{new_serial:5d}" + line[11:17] + "LIG Z" + line[22:]
            out.write(hetatm + "\n")
        out.write("TER\n")

        for line in ligand_conects:
            parts = line.split()
            conect_line = "CONECT"
            for part in parts[1:]:
                try:
                    conect_line += f"{int(part) + serial_offset:5d}"
                except ValueError:
                    conect_line += part
            out.write(conect_line + "\n")

        out.write("END\n")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    if len(sys.argv) != 2:
        print("Usage: pipe_gnina.py <config_json>")
        sys.exit(1)

    with open(sys.argv[1], 'r') as f:
        config = json.load(f)

    print("=" * 60)
    print("GNINA Docking Pipeline")
    print("=" * 60)

    structures_ds = load_datastream(config["structures_json"])
    compounds_ds = load_datastream(config["compounds_json"])

    missing_csv = config.get("missing_csv")
    if missing_csv and os.path.exists(missing_csv):
        missing_df = pd.read_csv(missing_csv)
        missing_ids = set(missing_df["id"].astype(str))
        if missing_ids:
            pairs = [(id_, f) for id_, f in zip(structures_ds.ids, structures_ds.files)
                     if id_ not in missing_ids]
            skipped = len(structures_ds.ids) - len(pairs)
            if skipped:
                print(f"  Skipping {skipped} IDs from upstream missing table")
            structures_ds.ids = [p[0] for p in pairs]
            structures_ds.files = [p[1] for p in pairs]

    print("\n--- Step 1: Protein Preparation ---")
    prepared_proteins = prepare_proteins(
        structures_ds, config["output_folder"],
        config["protonate"], config["pH"]
    )

    print("\n--- Step 2: Conformer Preparation ---")
    conformers = prepare_conformers(compounds_ds, config["output_folder"], config)
    if not conformers:
        print("Error: No conformers prepared, cannot proceed with docking")
        sys.exit(1)

    print("\n--- Step 3: Docking Execution ---")
    all_poses = run_docking(prepared_proteins, conformers, config)

    print("\n--- Step 4: Pose Consistency Analysis ---")
    consistency = analyze_pose_consistency(all_poses, config)

    print("\n--- Step 5: Results Aggregation ---")
    aggregate_results(all_poses, conformers, consistency, config, prepared_proteins)

    print("\n" + "=" * 60)
    print("GNINA docking pipeline complete")
    print("=" * 60)


if __name__ == "__main__":
    main()
