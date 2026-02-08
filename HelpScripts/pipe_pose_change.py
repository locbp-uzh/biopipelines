#!/usr/bin/env python3
# Copyright (c) 2026 Gianluca Quargnali @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# This software is freely available for use, modification, and redistribution.
# If you use this software or any derivative of it in published work,
# you must cite the original author and this repository.
"""
Helper script for PoseChange - executed during SLURM runtime.

Calculates ligand pose RMSD and distance metrics between reference
and target structures using PyMOL.
"""

import argparse
import json
import os
import sys
import pandas as pd

# Import unified I/O utilities
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from biopipelines_io import load_datastream, iterate_files


def calculate_pose_change(cmd, reference_pdb, reference_ligand, target_pdb, target_id,
                            ligand, alignment_selection, calculate_centroid,
                            calculate_orientation):
    """
    Calculate pose distance metrics for a single target structure.

    Args:
        cmd: PyMOL command object (already initialized)
        reference_pdb: Path to reference structure
        reference_ligand: Ligand residue name in reference
        target_pdb: Path to target structure
        target_id: ID for target structure
        ligand: Ligand residue name in target
        alignment_selection: PyMOL selection for protein alignment
        calculate_centroid: Whether to calculate centroid distance
        calculate_orientation: Whether to calculate orientation angle

    Returns:
        Dictionary with pose distance metrics
    """
    try:
        # Load structures
        cmd.load(reference_pdb, 'reference')
        cmd.load(target_pdb, 'target')

        # Check ligands exist
        ref_lig_count = cmd.count_atoms(f'reference and resn {reference_ligand}')
        if ref_lig_count == 0:
            raise ValueError(
                f"Ligand '{reference_ligand}' not found in reference structure {reference_pdb}"
            )

        tgt_lig_count = cmd.count_atoms(f'target and resn {ligand}')
        if tgt_lig_count == 0:
            raise ValueError(
                f"Ligand '{ligand}' not found in target structure {target_pdb}"
            )

        # Align proteins (not ligands) first
        try:
            # Align entire structures by their object names
            # This is more robust than using selection keywords like "protein"
            alignment_result = cmd.align('target', 'reference')
            alignment_rmsd = alignment_result[0]  # RMSD after alignment

        except Exception as e:
            raise ValueError(f"Protein alignment failed: {e}")

        # Calculate ligand RMSD after protein alignment
        try:
            # Select ligands
            ref_ligand_sel = f'reference and resn {reference_ligand}'
            tgt_ligand_sel = f'target and resn {ligand}'

            # Calculate RMSD between ligands (after protein alignment)
            ligand_rmsd = cmd.rms_cur(tgt_ligand_sel, ref_ligand_sel, matchmaker=-1)

            num_ligand_atoms = cmd.count_atoms(tgt_ligand_sel)

        except Exception as e:
            raise ValueError(f"Ligand RMSD calculation failed: {e}")

        # Initialize result dictionary
        result = {
            'id': target_id,
            'target_structure': os.path.basename(target_pdb),
            'reference_structure': os.path.basename(reference_pdb),
            'ligand_rmsd': round(ligand_rmsd, 3),
            'alignment_rmsd': round(alignment_rmsd, 3),
            'num_ligand_atoms': num_ligand_atoms,
            'alignment_method': 'whole_structure'
        }

        # Calculate centroid distance if requested
        if calculate_centroid:
            try:
                # Get centroid coordinates
                ref_center = cmd.centerofmass(ref_ligand_sel)
                tgt_center = cmd.centerofmass(tgt_ligand_sel)

                # Calculate Euclidean distance
                import math
                centroid_dist = math.sqrt(
                    sum((r - t)**2 for r, t in zip(ref_center, tgt_center))
                )
                result['centroid_distance'] = round(centroid_dist, 3)

            except Exception as e:
                print(f"Warning: Centroid calculation failed for {target_id}: {e}")
                result['centroid_distance'] = None

        # Calculate orientation angle if requested
        if calculate_orientation:
            try:
                # Calculate principal axes for both ligands
                ref_inertia = cmd.intra_fit(ref_ligand_sel, 0)
                tgt_inertia = cmd.intra_fit(tgt_ligand_sel, 0)

                # Note: This is a simplified orientation metric
                # For more sophisticated analysis, use moment of inertia tensors
                result['orientation_angle'] = None  # Placeholder
                result['orientation_axis'] = None   # Placeholder
                print(f"Warning: Orientation calculation not fully implemented yet")

            except Exception as e:
                print(f"Warning: Orientation calculation failed for {target_id}: {e}")
                result['orientation_angle'] = None
                result['orientation_axis'] = None

        return result

    finally:
        # Clean up loaded structures for next iteration
        cmd.delete('all')


def main():
    """Main entry point for pose distance analysis."""
    parser = argparse.ArgumentParser(description='Calculate ligand pose distances')
    parser.add_argument('--config', required=True, help='Path to JSON config file')
    args = parser.parse_args()

    # Load configuration
    with open(args.config, 'r') as f:
        config = json.load(f)

    reference_pdb = config['reference_pdb']
    reference_ligand = config['reference_ligand']
    ligand = config['ligand']
    alignment_selection = config['alignment_selection']
    calculate_centroid = config['calculate_centroid']
    calculate_orientation = config['calculate_orientation']
    output_csv = config['output_csv']

    # Load sample structures DataStream using pipe_biopipelines_io
    samples_ds = load_datastream(config['samples_json'])

    # Initialize PyMOL once for all structures
    try:
        import pymol
        from pymol import cmd
    except ImportError:
        print("Error: PyMOL is required for PoseChange")
        print("Install with: conda install -c conda-forge pymol-open-source")
        sys.exit(1)

    print("Initializing PyMOL...")
    pymol.finish_launching(['pymol', '-cq'])
    print("PyMOL initialized")

    print(f"Analyzing pose distances for {len(samples_ds.ids)} structures")
    print(f"Reference: {os.path.basename(reference_pdb)}")
    print(f"Reference ligand: {reference_ligand}")
    print(f"Target ligand: {ligand}")
    print(f"Alignment: {alignment_selection}")

    # Process each target structure using iterate_files for proper ID-file matching
    results = []
    failed = []
    target_items = list(iterate_files(samples_ds))
    total = len(target_items)

    for idx, (target_id, target_pdb) in enumerate(target_items, 1):
        print(f"[{idx}/{total}] Processing {target_id}...", end=" ")
        try:
            result = calculate_pose_change(
                cmd=cmd,
                reference_pdb=reference_pdb,
                reference_ligand=reference_ligand,
                target_pdb=target_pdb,
                target_id=target_id,
                ligand=ligand,
                alignment_selection=alignment_selection,
                calculate_centroid=calculate_centroid,
                calculate_orientation=calculate_orientation
            )
            results.append(result)
            print(f"OK (RMSD: {result['ligand_rmsd']:.2f} A)")

        except Exception as e:
            print(f"FAILED: {e}")
            failed.append({'id': target_id, 'error': str(e)})

    # Save results to CSV
    if results:
        df = pd.DataFrame(results)
        df.to_csv(output_csv, index=False)
        print(f"Successfully analyzed {len(results)} structures")
        print(f"Results saved to: {output_csv}")

        # Print summary statistics
        print(f"Ligand RMSD statistics:")
        print(f"  Mean:   {df['ligand_rmsd'].mean():.3f} A")
        print(f"  Median: {df['ligand_rmsd'].median():.3f} A")
        print(f"  Min:    {df['ligand_rmsd'].min():.3f} A")
        print(f"  Max:    {df['ligand_rmsd'].max():.3f} A")

        if calculate_centroid and 'centroid_distance' in df.columns:
            print(f"Centroid distance statistics:")
            print(f"  Mean:   {df['centroid_distance'].mean():.3f} A")
            print(f"  Median: {df['centroid_distance'].median():.3f} A")

    else:
        print("All structures failed analysis")
        sys.exit(1)

    if failed:
        print(f"{len(failed)} structures failed:")
        for item in failed:
            print(f"  {item['id']}: {item['error']}")

    # Clean up PyMOL at the very end
    print("Shutting down PyMOL...")
    cmd.quit()


if __name__ == "__main__":
    main()
