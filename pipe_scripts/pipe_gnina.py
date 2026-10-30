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
5. Results aggregation (docking_results.csv + docking_summary.csv)

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

# Add repo root to path so biopipelines package is importable
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from biopipelines.biopipelines_io import (
    load_datastream, iterate_files, iterate_values, load_table,
    read_upstream_missing, MISSING_COLUMNS, step_id_from_table_path,
)
from biopipelines.pdb_parser import field_atom_name, field_res_name
from biopipelines.ligand_utils import write_ligand_sdf
from biopipelines.id_map_utils import get_mapped_ids
from biopipelines.combinatorics import generate_multiplied_ids, predict_single_output_id
from biopipelines.gnina import CONFORMER_PREFIX, CONFORMER_RANK_METRICS, refuse_shared_pair_ids

import pipe_vina_backend as vina_backend


def _allocated_cpus():
    """Cores this step may use, from SLURM's own allocation.

    Vina otherwise sizes its thread pool from the visible core count, which on an
    exclusive node is all 288 even when the step was given a slice of them.
    Returns None off SLURM, letting Vina choose.
    """
    value = os.environ.get("SLURM_CPUS_PER_TASK") or os.environ.get("SLURM_CPUS_ON_NODE")
    if not value:
        return None
    try:
        return max(1, int(value))
    except ValueError:
        return None


def _pair_id(protein_id, ligand_id):
    """Id of a (protein, ligand) pair, composed by the shared axis rule the wrapper declares it with."""
    return predict_single_output_id(
        structures=("each", [protein_id], 0, [], False),
        compounds=("each", [ligand_id], 0, [], False),
    )


def _group_id(protein_id, ligand_id, conformer_id):
    """Id of a (protein, ligand, conformer) group, one row of docking_summary."""
    return generate_multiplied_ids([_pair_id(protein_id, ligand_id)], [f"{CONFORMER_PREFIX}{conformer_id}"])[0][0]


def _run_id(group_id, run_idx):
    """Id of one independent docking run of a group, also the stem of its scratch files."""
    return generate_multiplied_ids([group_id], [f"r{run_idx}"])[0][0]


def _pose_id(run_id, pose_idx):
    """Id of one pose of a run, one row of docking_results."""
    return generate_multiplied_ids([run_id], [f"p{pose_idx}"])[0][0]


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
                if field_atom_name(line) == "CA":
                    resname = field_res_name(line)
                    chain_seqs.setdefault(chain_id, []).append(resname)
                    coords = _parse_pdb_coords(line)
                    if coords:
                        chain_ca_coords.setdefault(chain_id, []).append(coords)
            elif line.startswith("HETATM"):
                resname = field_res_name(line)
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
            f.write("TER".ljust(80) + "\n")
        f.write("END".ljust(80) + "\n")

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
                f.write("END".ljust(80) + "\n")
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

    expanded_ids = compounds_ds.ids_expanded
    expanded_files = compounds_ds.files_expanded
    for ligand_id in expanded_ids:
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
        if mol is None and expanded_files:
            idx = expanded_ids.index(ligand_id) if ligand_id in expanded_ids else -1
            if idx >= 0:
                if len(expanded_files) == len(expanded_ids):
                    file_path = expanded_files[idx]
                elif compounds_ds.is_shared_file:
                    file_path = compounds_ds.files
                elif isinstance(compounds_ds.files, list) and len(compounds_ds.files) == 1:
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
    if not conformer_energies_ref or not conformer_energies_ref.startswith("TABLE_REFERENCE:"):
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

def resolve_autobox_references(protein_ids, autobox_by_id, structures_json=None):
    """Map each receptor id to its own autobox reference file.

    Resolved in ONE pass, and through the structures map_table when one is
    available: an upstream rename leaves no string relation between the ids
    (Panda_1 vs design_3), so suffix matching alone silently finds nothing and
    every receptor falls back to its own crystal HETATM.
    """
    if not autobox_by_id:
        return {}
    maps = None
    if structures_json and os.path.exists(structures_json):
        sds = load_datastream(structures_json)
        maps = [sds.map_table] if sds.map_table else None
    mapped = get_mapped_ids(source_ids=list(protein_ids),
                            target_ids=list(autobox_by_id),
                            map_table_paths=maps, unique=True)
    return {pid: (autobox_by_id.get(pid) or autobox_by_id.get(mapped.get(pid)))
            for pid in protein_ids}


def run_docking(prepared_proteins, conformers, config):
    """
    Run GNINA docking for each (protein, conformer) pair.

    For each pair, runs GNINA num_runs times with incremented seeds.
    Parses output SDF files for Vina and CNN scores.

    Returns:
        (list of pose record dicts, {pair id: [why a run of it failed]}); a pair keeps its reasons even when another run of it posed.
    """
    output_folder = config["output_folder"]
    backend = config["backend"]
    gnina_binary = config["gnina_binary"]
    exhaustiveness = config["exhaustiveness"]
    num_modes = config["num_modes"]
    num_runs = config["num_runs"]
    seed = config["seed"]
    cnn_scoring = config["cnn_scoring"]
    cnn_score_threshold = config["cnn_score_threshold"]
    scoring = config["scoring"]
    cpus = _allocated_cpus()
    box = config["box"]
    # Per-run wall-clock cap (seconds). A flexible-ligand search can occasionally
    # never terminate; kill that run and move on rather than hanging the job.
    # None or 0 disables. Set at config time via Gnina(dock_timeout=...).
    dock_timeout = config.get("dock_timeout", 1800) or None

    dock_dir = os.path.join(output_folder, "docking")
    os.makedirs(dock_dir, exist_ok=True)

    explicit_box_args = _build_box_args(box)

    # Per-structure autobox references, when the caller passed a stream with one
    # entry per receptor. Resolving a single index-0 path and reusing it would box
    # every receptor on the first structure's ligand.
    autobox_by_id = {}
    _ab_ds = box.get("autobox_ligand_ds")
    if _ab_ds and os.path.exists(_ab_ds):
        try:
            _ds = load_datastream(_ab_ds)
            autobox_by_id = {i: p for i, p in iterate_files(_ds)}
            print(f"Per-structure autobox references: {len(autobox_by_id)}")
        except Exception as e:
            print(f"Warning: could not load autobox datastream {_ab_ds}: {e}")

    autobox_for = resolve_autobox_references(
        list(prepared_proteins), autobox_by_id, config.get("structures_json"))

    all_poses = []

    # Why a run of a pair failed, so a pair with no pose can say so in missing.csv instead of vanishing from every output table.
    causes = defaultdict(list)
    total_pairs = len(prepared_proteins) * len(conformers)
    pair_idx = 0

    for protein_id, prot_info in prepared_proteins.items():
        protein_file = prot_info["protein"]
        crystal_ligand = prot_info.get("crystal_ligand")
        # An explicit per-structure reference wins over the HETATM extracted from
        # the receptor itself.
        if autobox_by_id:
            _ref = autobox_for.get(protein_id)
            if _ref:
                crystal_ligand = _ref
            else:
                print(f"  Warning: no autobox reference for {protein_id}")

        if backend == "vina":
            # Vina needs a PDBQT receptor and an explicit numeric box; the
            # crystal ligand becomes center+size rather than a reference path.
            receptor_input = os.path.join(dock_dir, f"{protein_id}_receptor.pdbqt")
            try:
                vina_backend.receptor_pdbqt(protein_file, receptor_input)
            except RuntimeError as e:
                print(f"  Warning: receptor conversion failed for {protein_id}: {e}")
                for ligand_id in dict.fromkeys(conf["ligand_id"] for conf in conformers):
                    causes[_pair_id(protein_id, ligand_id)].append(f"receptor PDBQT conversion failed: {e}")
                continue
            vina_box = dict(box)
            if "center" not in vina_box and crystal_ligand:
                print(f"  Using crystal ligand for autobox: {os.path.basename(crystal_ligand)}")
                vina_box["autobox_ligand"] = crystal_ligand
            box_args = None  # resolved per conformer (may fall back to the ligand itself)
        else:
            receptor_input = protein_file
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

            group_id = _group_id(protein_id, ligand_id, conformer_id)
            pair_id = _pair_id(protein_id, ligand_id)
            print(f"Docking pair {pair_idx}/{total_pairs}: {group_id}")

            for run_idx in range(num_runs):
                label = _run_id(group_id, run_idx)
                out_sdf = os.path.join(dock_dir, f"{label}.sdf")

                if backend == "vina":
                    try:
                        lig_pdbqt = vina_backend.ligand_pdbqt(
                            conf_sdf, os.path.join(dock_dir, f"{label}_lig.pdbqt"))
                        run_box = vina_backend.box_args(vina_box, ligand_sdf=conf_sdf,
                                                        receptor_file=protein_file)
                        out_pdbqt = os.path.join(dock_dir, f"{label}.pdbqt")
                        result = vina_backend.run(
                            gnina_binary, receptor_input, lig_pdbqt, out_pdbqt,
                            run_box, scoring=scoring, exhaustiveness=exhaustiveness,
                            num_modes=num_modes, seed=seed + run_idx, cpu=cpus,
                            mode="docking", timeout=dock_timeout,
                            extra_args=config.get("extra_args"),
                        )
                    except subprocess.TimeoutExpired:
                        print(f"Warning: Vina timed out (>{dock_timeout}s), skipping run {label}")
                        causes[pair_id].append(f"Vina timed out (>{dock_timeout}s) in run {label}")
                        continue
                    except RuntimeError as e:
                        print(f"Warning: Vina input preparation failed for {label}: {e}")
                        causes[pair_id].append(f"Vina input preparation failed for {label}: {e}")
                        continue
                    if result.returncode != 0:
                        print(f"Warning: Vina failed for {label}")
                        print(f"  stderr: {result.stderr[:500]}")
                        causes[pair_id].append(f"Vina failed for {label}: {result.stderr[:100].strip()}")
                        continue
                    if vina_backend.pdbqt_to_sdf(out_pdbqt, out_sdf,
                                                 template_sdf=conf_sdf) is None:
                        print(f"Warning: could not convert Vina output for {label}")
                        causes[pair_id].append(f"could not convert the Vina output of {label}")
                        continue
                    poses = _parse_gnina_output(
                        out_sdf, protein_id, ligand_id, conformer_id,
                        run_idx, cnn_score_threshold
                    )
                    all_poses.extend(poses)
                    continue

                cmd = [
                    gnina_binary,
                    "-r", receptor_input,
                    "-l", conf_sdf,
                ] + box_args + [
                    "--exhaustiveness", str(exhaustiveness),
                    "--num_modes", str(num_modes),
                    "--cnn_scoring", cnn_scoring,
                    "--seed", str(seed + run_idx),
                    "-o", out_sdf,
                # Forwarded kwargs as argv tokens; no shell is involved, so nothing re-parses them.
                ] + list(config.get("extra_args") or [])

                try:
                    result = subprocess.run(cmd, capture_output=True, text=True,
                                            timeout=dock_timeout)
                except subprocess.TimeoutExpired:
                    print(f"Warning: GNINA timed out (>{dock_timeout}s), skipping run {label}")
                    causes[pair_id].append(f"GNINA timed out (>{dock_timeout}s) in run {label}")
                    continue
                if result.returncode != 0:
                    print(f"Warning: GNINA failed for {label}")
                    print(f"  stderr: {result.stderr[:500]}")
                    causes[pair_id].append(f"GNINA failed for {label}: {result.stderr[:100].strip()}")
                    continue

                poses = _parse_gnina_output(
                    out_sdf, protein_id, ligand_id, conformer_id,
                    run_idx, cnn_score_threshold
                )
                all_poses.extend(poses)

    print(f"Docking complete: {len(all_poses)} accepted poses")
    return all_poses, dict(causes)


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
    # An empty file is a run that posed nothing; RDKit raises on it instead of yielding no molecules.
    if not os.path.exists(sdf_path) or os.path.getsize(sdf_path) == 0:
        return poses

    group_id = _group_id(protein_id, ligand_id, conformer_id)
    run_id = _run_id(group_id, run_idx)

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

        poses.append({
            "id": _pose_id(run_id, pose_idx),
            "group_id": group_id,
            "sdf_file": sdf_path,
            "structures.id": protein_id,
            "compounds.id": ligand_id,
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
        dict: group id (the docking_summary id) -> consistency_score
    """
    rmsd_threshold = config["rmsd_threshold"]

    groups = defaultdict(list)
    for pose in all_poses:
        groups[pose["group_id"]].append(pose)

    consistency = {}
    for key, poses in groups.items():
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
                sdf_file = pose["sdf_file"]
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

def _pair_winners(summary_rows, best_poses, rank_by="auto"):
    """The group id whose best pose stands for each pair, and the pairs a named metric could not decide.

    Returns (winners, unrankable): `winners` maps a pair id to its winning group id, `unrankable` maps a pair id to the group ids that lack the metric. "auto" ranks by pseudo_binding_energy (best Vina plus conformer strain) when every posed group of the pair has one, else by best Vina score. A named metric is used as given, in the direction CONFORMER_RANK_METRICS records, and never replaced: a pair with a posed group that lacks it has no winner. The first group wins a tie.
    """
    if rank_by != "auto" and rank_by not in CONFORMER_RANK_METRICS:
        raise ValueError(f"rank_conformers_by must be 'auto' or one of {list(CONFORMER_RANK_METRICS)}, got {rank_by!r}")
    posed = defaultdict(list)
    for row in summary_rows:
        if row["id"] in best_poses:
            posed[_pair_id(row["structures.id"], row["compounds.id"])].append(row)
    winners, unrankable = {}, {}
    for pair_id, rows in posed.items():
        if rank_by == "auto":
            column = ("pseudo_binding_energy"
                      if all(r["pseudo_binding_energy"] != "" for r in rows) else "best_vina")
        else:
            column = rank_by
            lacking = [r["id"] for r in rows if r.get(column, "") == ""]
            if lacking:
                unrankable[pair_id] = lacking
                continue
        pick = min if CONFORMER_RANK_METRICS.get(column, "lower") == "lower" else max
        winners[pair_id] = pick(rows, key=lambda r: r[column])["id"]
    return winners, unrankable


def aggregate_results(all_poses, conformers, consistency, config, prepared_proteins):
    """
    Aggregate results into docking_results.csv and docking_summary.csv.

    Also writes each group's best pose as a combined protein+ligand PDB file. With top_only the pair's winning conformer goes into the structures folder as <pair id>_best.pdb and every other conformer under conformer_poses/; without it every conformer goes there as <summary id>_best.pdb. A pair that a named rank_conformers_by metric cannot decide gets no stream entry and its conformers go under conformer_poses/.

    Returns (the pair ids that have an entry in the structures stream, {pair id: [why it has none]} for the pairs a named metric could not decide).
    """
    output_folder = config["output_folder"]
    docking_results_csv = config["docking_results_csv"]
    docking_summary_csv = config["docking_summary_csv"]
    is_vina = config["backend"] == "vina"
    top_only = config["top_only"]
    rank_by = config["rank_conformers_by"]
    cnn_pose_cols = [] if is_vina else ["cnn_score", "cnn_affinity"]

    results_cols = ["id", "structures.id", "compounds.id", "conformer_id",
                    "run", "pose", "vina_score"] + cnn_pose_cols
    if all_poses:
        results_df = pd.DataFrame(all_poses)[results_cols]
        results_df.to_csv(docking_results_csv, index=False)
        print(f"Wrote {len(results_df)} poses to {docking_results_csv}")
    else:
        pd.DataFrame(columns=results_cols).to_csv(docking_results_csv, index=False)
        print("Warning: No poses passed filtering")

    energy_lookup = {
        (c["ligand_id"], c["conformer_id"]): c.get("energy", "")
        for c in conformers
    }

    best_poses_dir = config.get("best_poses_dir") or os.path.join(output_folder, "best_poses")
    os.makedirs(best_poses_dir, exist_ok=True)
    conformer_poses_dir = os.path.join(output_folder, "conformer_poses")

    groups = defaultdict(list)
    for pose in all_poses:
        groups[pose["group_id"]].append(pose)

    summary_rows = []
    best_poses = {}
    for group_id, poses in groups.items():
        protein_id = poses[0]["structures.id"]
        ligand_id = poses[0]["compounds.id"]
        conformer_id = poses[0]["conformer_id"]

        vina_scores = [p["vina_score"] for p in poses if p["vina_score"] != ""]
        cnn_scores = [p["cnn_score"] for p in poses if p["cnn_score"] != ""]
        cnn_affinities = [p["cnn_affinity"] for p in poses if p["cnn_affinity"] != ""]

        # Best pose per run (by Vina score)
        best_per_run = {}
        for pose in poses:
            run = pose["run"]
            vina = pose["vina_score"]
            if vina == "":
                continue
            if run not in best_per_run or vina < best_per_run[run]["vina_score"]:
                best_per_run[run] = pose

        run_vinas = [best_per_run[r]["vina_score"] for r in best_per_run]
        run_cnn_affinities = [best_per_run[r]["cnn_affinity"] for r in best_per_run
                              if best_per_run[r]["cnn_affinity"] != ""]

        best_vina = min(vina_scores) if vina_scores else ""
        mean_vina = float(np.mean(run_vinas)) if run_vinas else ""
        std_vina = float(np.std(run_vinas)) if len(run_vinas) > 1 else ""
        best_cnn = max(cnn_scores) if cnn_scores else ""
        mean_cnn_affinity = float(np.mean(run_cnn_affinities)) if run_cnn_affinities else ""
        std_cnn_affinity = float(np.std(run_cnn_affinities)) if len(run_cnn_affinities) > 1 else ""
        pose_cons = consistency.get(group_id, "")
        conf_energy = energy_lookup.get((ligand_id, conformer_id), "")

        pseudo_be = ""
        if best_vina != "" and conf_energy != "":
            try:
                pseudo_be = float(best_vina) + float(conf_energy)
            except (ValueError, TypeError):
                pass

        if vina_scores:
            best_poses[group_id] = min(
                [p for p in poses if p["vina_score"] != ""],
                key=lambda p: p["vina_score"]
            )

        row = {
            "id": group_id,
            "structures.id": protein_id,
            "compounds.id": ligand_id,
            "conformer_id": conformer_id,
            "best_vina": best_vina,
            "mean_vina": mean_vina,
            "std_vina": std_vina,
        }
        if not is_vina:
            row.update({
                "best_cnn_score": best_cnn,
                "mean_cnn_affinity": mean_cnn_affinity,
                "std_cnn_affinity": std_cnn_affinity,
            })
        row.update({
            "pose_consistency": pose_cons,
            "conformer_energy": conf_energy,
            "pseudo_binding_energy": pseudo_be,
            "best_pose_file": "",
        })
        summary_rows.append(row)

    winners, unrankable = _pair_winners(summary_rows, best_poses, rank_by) if top_only else ({}, {})
    map_rows = []
    for row in summary_rows:
        best_pose = best_poses.get(row["id"])
        if best_pose is None:
            continue
        pair_id = _pair_id(row["structures.id"], row["compounds.id"])
        is_winner = winners.get(pair_id) == row["id"] if top_only else True
        stream_id = pair_id if top_only else row["id"]
        dst_pdb = (os.path.join(best_poses_dir, f"{stream_id}_best.pdb") if is_winner
                   else os.path.join(conformer_poses_dir, f"{row['id']}_best.pdb"))
        prot_info = prepared_proteins.get(row["structures.id"], {})
        protein_pdb = prot_info.get("protein") if isinstance(prot_info, dict) else prot_info
        if os.path.exists(best_pose["sdf_file"]) and protein_pdb and os.path.exists(protein_pdb):
            os.makedirs(os.path.dirname(dst_pdb), exist_ok=True)
            _write_complex_pdb(protein_pdb, best_pose["sdf_file"], best_pose["pose"], dst_pdb)
            row["best_pose_file"] = dst_pdb
            if is_winner:
                map_rows.append({"id": stream_id, "file": dst_pdb,
                                 "structures.id": row["structures.id"],
                                 "compounds.id": row["compounds.id"]})

    summary_df = pd.DataFrame(summary_rows)
    if not summary_df.empty and "mean_vina" in summary_df.columns:
        summary_df["_sort_key"] = pd.to_numeric(summary_df["mean_vina"], errors="coerce")
        summary_df = (summary_df.sort_values("_sort_key")
                                .drop(columns=["_sort_key"])
                                .reset_index(drop=True))

    summary_df.to_csv(docking_summary_csv, index=False)
    print(f"Wrote {len(summary_df)} docking summaries to {docking_summary_csv}")

    structures_map = config.get("structures_map")
    if structures_map:
        os.makedirs(os.path.dirname(structures_map), exist_ok=True)
        pd.DataFrame(map_rows, columns=["id", "file", "structures.id", "compounds.id"]).to_csv(
            structures_map, index=False
        )
        print(f"Wrote structures map: {structures_map} ({len(map_rows)} rows)")

    written_pairs = {_pair_id(row["structures.id"], row["compounds.id"]) for row in map_rows}
    rank_causes = {pair_id: [f"rank_conformers_by={rank_by!r} has no value for the posed conformer(s) "
                             f"{', '.join(ids)}"] for pair_id, ids in unrankable.items()}
    return written_pairs, rank_causes


def unprepared_causes(protein_ids, ligand_ids, prepared_proteins, conformers):
    """{pair id: [why]} for the pairs that never reached the docking: a protein that was not prepared, or a ligand with no conformer."""
    with_conformers = {conf["ligand_id"] for conf in conformers}
    causes = defaultdict(list)
    for protein_id in protein_ids:
        for ligand_id in ligand_ids:
            if protein_id not in prepared_proteins:
                causes[_pair_id(protein_id, ligand_id)].append("the protein could not be prepared")
            if ligand_id not in with_conformers:
                causes[_pair_id(protein_id, ligand_id)].append("no conformer could be prepared for the ligand")
    return dict(causes)


def local_failures(protein_ids, ligand_ids, written_pairs, causes, step_id):
    """One `failure` row, under the pair id, for every (protein, ligand) pair that has no entry in the structures stream.

    The cause is what the step recorded for the pair, else that no run of it gave an accepted pose.
    """
    rows = []
    for protein_id in protein_ids:
        for ligand_id in ligand_ids:
            pair_id = _pair_id(protein_id, ligand_id)
            if pair_id in written_pairs:
                continue
            cause = "; ".join(causes.get(pair_id, [])) or "no accepted pose"
            rows.append({"id": pair_id, "removed_by": step_id, "kind": "failure", "cause": cause[:200]})
    return rows


def write_local_missing(path, rows):
    """This step's own failures, which the propagation merges with the upstream manifests into tables/missing.csv."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    pd.DataFrame(rows, columns=MISSING_COLUMNS).to_csv(path, index=False)
    print(f"Wrote local missing manifest: {path} ({len(rows)} rows)")


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
        wrote_protein = False
        with open(protein_pdb, 'r') as prot:
            for line in prot:
                if line.startswith(("CONECT", "MASTER", "END")):
                    continue
                if line.startswith("TER"):
                    continue  # re-emit a single TER after the protein block below
                out.write(line)
                if line.startswith(("ATOM", "HETATM")):
                    wrote_protein = True
        # Terminate the protein chain before the ligand. OpenBabel's -h re-emit
        # drops the TER record, and without it viewers/parsers bond the protein
        # C-terminus straight into the chain-Z ligand (the protein renders broken).
        if wrote_protein:
            out.write("TER".ljust(80) + "\n")

        for line in ligand_atoms:
            old_serial = int(line[6:11].strip())
            new_serial = old_serial + serial_offset
            hetatm = f"HETATM{new_serial:5d}" + line[11:17] + "LIG Z" + line[22:]
            out.write(hetatm + "\n")
        out.write("TER".ljust(80) + "\n")

        for line in ligand_conects:
            parts = line.split()
            conect_line = "CONECT"
            for part in parts[1:]:
                try:
                    conect_line += f"{int(part) + serial_offset:5d}"
                except ValueError:
                    conect_line += part
            out.write(conect_line + "\n")

        out.write("END".ljust(80) + "\n")


# ---------------------------------------------------------------------------
# Score / minimize mode (no docking search)
# ---------------------------------------------------------------------------

def _load_score_chemistry(compounds_ds):
    """Map each compound id to its (code, smiles).

    The code locates the in-pocket HETATM residue to score; the SMILES restores
    that pose's bond orders before scoring.
    """
    chem = {}
    for cid, values in iterate_values(compounds_ds, columns=["code", "smiles"]):
        code = str(values.get("code", "")).strip()
        if not code or code.lower() == "nan":
            raise ValueError(f"compound {cid} has no residue code to locate its in-pocket pose")
        smi = str(values.get("smiles", "")).strip()
        if not smi or smi.lower() == "nan":
            raise ValueError(f"compound {cid} (code {code}) has no SMILES for bond-order templating")
        chem[cid] = (code, smi)
    return chem


def _extract_ligand_hetatm(complex_pdb, ligand_code, out_pdb):
    """Write the HETATM records whose residue name matches `ligand_code` to out_pdb.

    Returns out_pdb if any were found, else None.
    """
    lines, serials = [], set()
    conect = []
    with open(complex_pdb, 'r') as f:
        for line in f:
            if line.startswith("HETATM") and field_res_name(line).upper() == ligand_code.upper():
                lines.append(line)
                serials.add(line[6:11].strip())
            elif line.startswith("CONECT"):
                conect.append(line)
    if not lines:
        return None
    # Carry the declared bonds across. Without them the ligand is a bag of atoms and
    # perception falls back to distance, which misses long bonds like Si-C (1.87 A);
    # the molecule then fragments and OpenBabel writes one PDBQT ROOT per piece.
    kept = [c for c in conect
            if c[6:11].strip() in serials
            and any(c[i:i + 5].strip() in serials for i in range(11, min(len(c.rstrip()), 31), 5))]
    with open(out_pdb, 'w') as f:
        f.writelines(lines)
        f.writelines(kept)
        f.write("END".ljust(80) + "\n")
    return out_pdb


# SDF property tags emitted by gnina 1.3.2 --score_only/--minimize (verified
# against real output). CNN_VS is the CNN virtual-screening score; the affinity
# variance is CNNaffinity_variance (NOT CNN_VS). No intramolecular-energy term
# is written in these modes, so none is collected.
_SCORE_TERMS = {
    "vina_affinity": ("minimizedAffinity", "Affinity"),
    "cnn_score": ("CNNscore",),
    "cnn_affinity": ("CNNaffinity",),
    "cnn_vs": ("CNN_VS",),
    "cnn_affinity_variance": ("CNNaffinity_variance",),
}


def _parse_score_output(sdf_path):
    """Parse the first molecule of a gnina --score_only/--minimize SDF.

    Returns a dict of the available score terms (missing ones omitted) and the
    RDKit mol (for emitting the scored pose), or (None, None) on failure.
    """
    if not os.path.exists(sdf_path):
        return None, None
    with _suppress_rdkit_warnings():
        suppl = Chem.SDMolSupplier(sdf_path, removeHs=False)
        mol = next((m for m in suppl if m is not None), None)
    if mol is None:
        return None, None
    scores = {}
    for col, prop_names in _SCORE_TERMS.items():
        for prop in prop_names:
            val = _get_prop(mol, prop)
            if val is not None:
                scores[col] = val
                break
    return scores, mol


def _run_vina_score(vina_binary, clean_protein, lig_sdf, scratch, structure_id,
                    mode, scoring, cpus, timeout, out_sdf, autobox_add=4.0, extra_args=None):
    """Score or locally minimize an in-pocket pose with Vina.

    ``--score_only`` prints its affinity and writes no file, so the input pose is
    re-emitted as the scored pose; ``--local_only`` writes the refined pose.
    Returns (scores_dict, pose_sdf).
    """
    receptor = vina_backend.receptor_pdbqt(
        clean_protein, os.path.join(scratch, f"{structure_id}_receptor.pdbqt"))
    lig_pdbqt = vina_backend.ligand_pdbqt(
        lig_sdf, os.path.join(scratch, f"{structure_id}_lig.pdbqt"))
    out_pdbqt = os.path.join(scratch, f"{structure_id}_{mode}.pdbqt")

    # Vina needs an explicit box even for --score_only, unlike gnina which autoboxes
    # on the input ligand. Box the pose being scored.
    box_args = vina_backend.box_args({"autobox_add": autobox_add}, ligand_sdf=lig_sdf,
                                     receptor_file=clean_protein)

    try:
        result = vina_backend.run(
            vina_binary, receptor, lig_pdbqt, out_pdbqt, box_args, scoring=scoring,
            cpu=cpus, mode=mode, timeout=timeout,
            extra_args=extra_args,
        )
    except subprocess.TimeoutExpired:
        raise RuntimeError(f"vina {mode} timed out (>{timeout}s)")
    if result.returncode != 0:
        raise RuntimeError(f"vina {mode} failed: {result.stderr[:300]}")

    if mode == "score":
        affinity = vina_backend.parse_score_stdout(result.stdout)
        if affinity is None:
            raise RuntimeError("could not parse vina --score_only output")
        return {"vina_affinity": affinity}, lig_sdf

    if vina_backend.pdbqt_to_sdf(out_pdbqt, out_sdf, template_sdf=lig_sdf) is None:
        raise RuntimeError(f"could not convert vina output {out_pdbqt}")
    # gnina tags the SDF with its scores; vina --local_only prints the affinity to
    # stdout exactly as --score_only does and writes no such tags.
    scores, _ = _parse_score_output(out_sdf)
    if not scores:
        affinity = vina_backend.parse_score_stdout(result.stdout)
        if affinity is None:
            raise RuntimeError(f"could not parse vina output {out_sdf}")
        scores = {"vina_affinity": affinity}
    return scores, out_sdf


def run_score_mode(structures_ds, compounds_ds, config):
    """Score (or locally minimize then score) the in-pocket ligand of each complex.

    For each complex: strip it to a clean protein, extract the ligand HETATM
    (identified by the compounds residue code), restore bond orders from the
    SMILES template, then run gnina --score_only/--minimize autoboxed on that
    ligand. Writes scores.csv, the scored-complex PDBs, and the structures map.
    """
    mode = config["mode"]
    backend = config["backend"]
    gnina_binary = config["gnina_binary"]
    cnn_scoring = config["cnn_scoring"]
    scoring = config["scoring"]
    cpus = _allocated_cpus()
    autobox_add = config.get("box", {}).get("autobox_add", 4.0)
    dock_timeout = config.get("dock_timeout", 1800) or None
    output_folder = config["output_folder"]
    best_poses_dir = config["best_poses_dir"]
    scores_csv = config["scores_csv"]
    structures_map = config["structures_map"]

    scratch = os.path.join(output_folder, "scoring")
    os.makedirs(scratch, exist_ok=True)
    os.makedirs(best_poses_dir, exist_ok=True)

    chem = _load_score_chemistry(compounds_ds)
    step_id = step_id_from_table_path(scores_csv)
    # Vina emits a single affinity; the CNN terms have no counterpart.
    score_terms = ["vina_affinity"] if backend == "vina" else list(_SCORE_TERMS)

    # Pair each complex with the ONE compound it carries. A complex already holds
    # a single posed ligand; scoring another compound's SMILES against those
    # coordinates is meaningless. One compound broadcasts to every complex;
    # otherwise the structure→compound link comes from the structures map's
    # provenance (compounds.id), resolved by get_mapped_ids.
    structure_ids = list(structures_ds.ids_expanded)
    compound_ids = list(chem.keys())
    if len(compound_ids) == 1:
        structure_to_compound = {sid: compound_ids[0] for sid in structure_ids}
    else:
        structure_to_compound = get_mapped_ids(
            source_ids=structure_ids,
            target_ids=compound_ids,
            map_table_paths=[structures_ds.map_table] if structures_ds.map_table else None,
            unique=True,
        )

    gnina_flag = "--score_only" if mode == "score" else "--minimize"

    rows, map_rows, failed = [], [], []
    for structure_id, complex_pdb in iterate_files(structures_ds):
        try:
            compound_id = structure_to_compound.get(structure_id)
            if compound_id is None:
                raise RuntimeError(
                    f"no compound paired with complex {structure_id}; the structures map "
                    f"carries no compounds.id provenance and the compounds stream has "
                    f"{len(compound_ids)} compounds")
            ligand_code, template_smiles = chem[compound_id]

            clean_protein = os.path.join(scratch, f"{structure_id}_protein.pdb")
            _clean_protein_pdb(complex_pdb, clean_protein)

            lig_pdb = os.path.join(scratch, f"{structure_id}_lig.pdb")
            if _extract_ligand_hetatm(complex_pdb, ligand_code, lig_pdb) is None:
                raise RuntimeError(f"no HETATM records with code {ligand_code} in {complex_pdb}")

            lig_sdf = os.path.join(scratch, f"{structure_id}_lig.sdf")
            write_ligand_sdf(lig_pdb, lig_sdf, smiles=template_smiles)

            out_sdf = os.path.join(scratch, f"{structure_id}_{mode}.sdf")

            if backend == "vina":
                scores, out_sdf = _run_vina_score(
                    gnina_binary, clean_protein, lig_sdf, scratch, structure_id,
                    mode, scoring, cpus, dock_timeout, out_sdf,
                    autobox_add=autobox_add, extra_args=config.get("extra_args"),
                )
            else:
                cmd = [
                    gnina_binary,
                    "-r", clean_protein,
                    "-l", lig_sdf,
                    gnina_flag,
                    "--autobox_ligand", lig_sdf,
                    "--autobox_add", str(autobox_add),
                    "--cnn_scoring", cnn_scoring,
                    "-o", out_sdf,
                # Forwarded kwargs as argv tokens; no shell is involved, so nothing re-parses them.
                ] + list(config.get("extra_args") or [])
                try:
                    result = subprocess.run(cmd, capture_output=True, text=True, timeout=dock_timeout)
                except subprocess.TimeoutExpired:
                    raise RuntimeError(f"gnina {mode} timed out (>{dock_timeout}s)")
                if result.returncode != 0:
                    raise RuntimeError(f"gnina {mode} failed: {result.stderr[:300]}")

                scores, scored_mol = _parse_score_output(out_sdf)
            if scores is None:
                raise RuntimeError(f"could not parse {backend} output {out_sdf}")

            # "score" leaves the pose untouched; re-emit the input ligand SDF so the
            # complex reflects exactly what was scored. "minimize" emits the refined
            # pose gnina wrote.
            pose_sdf = lig_sdf if mode == "score" else out_sdf
            dst_pdb = os.path.join(best_poses_dir, f"{structure_id}_{mode}.pdb")
            _write_complex_pdb(clean_protein, pose_sdf, 0, dst_pdb)

            row = {"id": structure_id, "structures.id": structure_id, "compounds.id": compound_id}
            for col in score_terms:
                row[col] = scores.get(col, "")
            rows.append(row)
            map_rows.append({"id": structure_id, "file": dst_pdb,
                             "structures.id": structure_id, "compounds.id": compound_id})
            print(f"  {structure_id} ({compound_id}): " + ", ".join(f"{k}={scores[k]}" for k in scores))
        except Exception as e:
            print(f"WARNING: {structure_id} {mode} failed: {e}", file=sys.stderr)
            failed.append({"id": structure_id, "removed_by": step_id, "kind": "failure",
                           "cause": str(e)[:200]})

    score_cols = ["id", "structures.id", "compounds.id"] + list(score_terms)
    os.makedirs(os.path.dirname(scores_csv), exist_ok=True)
    pd.DataFrame(rows, columns=score_cols).to_csv(scores_csv, index=False)
    print(f"Wrote {len(rows)} scores to {scores_csv}")

    if structures_map:
        os.makedirs(os.path.dirname(structures_map), exist_ok=True)
        pd.DataFrame(map_rows, columns=["id", "file", "structures.id", "compounds.id"]).to_csv(
            structures_map, index=False)
        print(f"Wrote structures map: {structures_map} ({len(map_rows)} rows)")

    return failed


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    if len(sys.argv) != 2:
        print("Usage: pipe_gnina.py <config_json>")
        sys.exit(1)

    with open(sys.argv[1], 'r') as f:
        config = json.load(f)

    mode = config.get("mode", "docking")
    engine = "Vina" if config["backend"] == "vina" else "GNINA"
    print("=" * 60)
    print(f"{engine} {mode} pipeline")
    print("=" * 60)

    structures_ds = load_datastream(config["structures_json"])
    compounds_ds = load_datastream(config["compounds_json"])

    upstream_missing_rows = read_upstream_missing(config.get("upstream_missing"))
    missing_ids = {str(r["id"]) for r in upstream_missing_rows}
    if missing_ids:
        expanded_ids = structures_ds.ids_expanded
        expanded_files = structures_ds.files_expanded
        pairs = [(id_, f) for id_, f in zip(expanded_ids, expanded_files)
                 if id_ not in missing_ids]
        skipped = len(expanded_ids) - len(pairs)
        if skipped:
            print(f"  Skipping {skipped} structure IDs from upstream missing table")
        structures_ds.ids = [p[0] for p in pairs]
        structures_ds.files = [p[1] for p in pairs]
        structures_ds._ids_expanded = None
        structures_ds._files_expanded = None

        comp_all = compounds_ds.ids_expanded
        comp_kept = [cid for cid in comp_all if cid not in missing_ids]
        if len(comp_all) - len(comp_kept):
            print(f"  Skipping {len(comp_all) - len(comp_kept)} compound IDs from upstream missing table")
        compounds_ds.ids = comp_kept
        compounds_ds._ids_expanded = comp_kept

    if mode in ("score", "minimize"):
        print(f"\n--- GNINA {mode} (no docking search) ---")
        failed = run_score_mode(structures_ds, compounds_ds, config)
        write_local_missing(config["local_missing_csv"], failed)
        print("\n" + "=" * 60)
        print(f"{engine} {mode} pipeline complete")
        print("=" * 60)
        return

    print("\n--- Step 1: Protein Preparation ---")
    protein_ids = [protein_id for protein_id, _file in iterate_files(structures_ds)]
    ligand_ids = list(compounds_ds.ids_expanded)
    refuse_shared_pair_ids(protein_ids, ligand_ids)
    prepared_proteins = prepare_proteins(
        structures_ds, config["output_folder"],
        config["protonate"], config["pH"]
    )

    print("\n--- Step 2: Conformer Preparation ---")
    conformers = prepare_conformers(compounds_ds, config["output_folder"], config)
    causes = defaultdict(list)
    for pair_id, why in unprepared_causes(protein_ids, ligand_ids, prepared_proteins, conformers).items():
        causes[pair_id] += why

    written_pairs = set()
    if conformers:
        print("\n--- Step 3: Docking Execution ---")
        all_poses, run_causes = run_docking(prepared_proteins, conformers, config)

        print("\n--- Step 4: Pose Consistency Analysis ---")
        consistency = analyze_pose_consistency(all_poses, config)

        print("\n--- Step 5: Results Aggregation ---")
        written_pairs, rank_causes = aggregate_results(
            all_poses, conformers, consistency, config, prepared_proteins)
        for found in (run_causes, rank_causes):
            for pair_id, why in found.items():
                causes[pair_id] += why

    step_id = step_id_from_table_path(config["local_missing_csv"])
    write_local_missing(config["local_missing_csv"],
                        local_failures(protein_ids, ligand_ids, written_pairs, causes, step_id))

    if not written_pairs:
        print(f"ERROR: no (protein, ligand) pair gave a structure; see {config['local_missing_csv']}",
              file=sys.stderr)
        sys.exit(1)

    print("\n" + "=" * 60)
    print(f"{engine} docking pipeline complete")
    print("=" * 60)


if __name__ == "__main__":
    main()
