#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Runtime helper script for PLIP (Protein-Ligand Interaction Profiler) analysis.

Loads structures via DataStream, runs PLIP on each PDB file to detect
non-covalent protein-ligand interactions, and writes two output CSVs:
  - interactions.csv: one row per interaction (detailed)
  - summary.csv: one row per ligand per structure (counts by type)
"""

import os
import sys
import argparse
import json
import pandas as pd
from typing import Dict, List, Any

# Import unified I/O utilities
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from biopipelines.biopipelines_io import load_datastream, iterate_files


def analyze_structure(structure_id: str, structure_path: str) -> tuple:
    """
    Run PLIP on a single PDB structure.

    Returns:
        Tuple of (interactions_rows, summary_rows) where each is a list of dicts.
        Returns ([], []) on failure.
    """
    from plip.structure.preparation import PDBComplex

    interactions_rows = []
    summary_rows = []

    try:
        mol = PDBComplex()
        mol.load_pdb(structure_path)
        mol.analyze()

        for bsid, interaction_set in mol.interaction_sets.items():
            # Parse binding site ID: typically "LIG:A:1" (resname:chain:resnum)
            parts = str(bsid).split(":")
            if len(parts) >= 2:
                ligand_name = parts[0]
                ligand_chain = parts[1]
            else:
                ligand_name = str(bsid)
                ligand_chain = ""

            counts = {
                "hydrophobic": 0,
                "hbond": 0,
                "water_bridge": 0,
                "salt_bridge": 0,
                "pi_stacking": 0,
                "pi_cation": 0,
                "halogen_bond": 0,
                "metal_complex": 0,
            }

            # --- Hydrophobic contacts ---
            for contact in interaction_set.hydrophobic_contacts:
                counts["hydrophobic"] += 1
                interactions_rows.append({
                    "id": structure_id,
                    "ligand_name": ligand_name,
                    "ligand_chain": ligand_chain,
                    "interaction_type": "hydrophobic",
                    "protein_residue": contact.restype,
                    "protein_chain": contact.reschain,
                    "protein_residue_number": contact.resnr,
                    "distance": round(contact.distance, 2),
                })

            # --- Hydrogen bonds (protein is donor) ---
            for hb in interaction_set.hbonds_pdon:
                counts["hbond"] += 1
                interactions_rows.append({
                    "id": structure_id,
                    "ligand_name": ligand_name,
                    "ligand_chain": ligand_chain,
                    "interaction_type": "hbond",
                    "protein_residue": hb.restype,
                    "protein_chain": hb.reschain,
                    "protein_residue_number": hb.resnr,
                    "distance": round(hb.distance_ah, 2),
                })

            # --- Hydrogen bonds (ligand is donor) ---
            for hb in interaction_set.hbonds_ldon:
                counts["hbond"] += 1
                interactions_rows.append({
                    "id": structure_id,
                    "ligand_name": ligand_name,
                    "ligand_chain": ligand_chain,
                    "interaction_type": "hbond",
                    "protein_residue": hb.restype,
                    "protein_chain": hb.reschain,
                    "protein_residue_number": hb.resnr,
                    "distance": round(hb.distance_ah, 2),
                })

            # --- Water bridges ---
            for wb in interaction_set.water_bridges:
                counts["water_bridge"] += 1
                interactions_rows.append({
                    "id": structure_id,
                    "ligand_name": ligand_name,
                    "ligand_chain": ligand_chain,
                    "interaction_type": "water_bridge",
                    "protein_residue": wb.restype,
                    "protein_chain": wb.reschain,
                    "protein_residue_number": wb.resnr,
                    "distance": round(wb.distance_aw, 2),
                })

            # --- Pi-stacking ---
            for ps in interaction_set.pistacking:
                counts["pi_stacking"] += 1
                interactions_rows.append({
                    "id": structure_id,
                    "ligand_name": ligand_name,
                    "ligand_chain": ligand_chain,
                    "interaction_type": "pi_stacking",
                    "protein_residue": ps.restype,
                    "protein_chain": ps.reschain,
                    "protein_residue_number": ps.resnr,
                    "distance": round(ps.distance, 2),
                })

            # --- Pi-cation (ligand aromatic) ---
            for pc in interaction_set.pication_laro:
                counts["pi_cation"] += 1
                interactions_rows.append({
                    "id": structure_id,
                    "ligand_name": ligand_name,
                    "ligand_chain": ligand_chain,
                    "interaction_type": "pi_cation",
                    "protein_residue": pc.restype,
                    "protein_chain": pc.reschain,
                    "protein_residue_number": pc.resnr,
                    "distance": round(pc.distance, 2),
                })

            # --- Pi-cation (protein aromatic) ---
            for pc in interaction_set.pication_paro:
                counts["pi_cation"] += 1
                interactions_rows.append({
                    "id": structure_id,
                    "ligand_name": ligand_name,
                    "ligand_chain": ligand_chain,
                    "interaction_type": "pi_cation",
                    "protein_residue": pc.restype,
                    "protein_chain": pc.reschain,
                    "protein_residue_number": pc.resnr,
                    "distance": round(pc.distance, 2),
                })

            # --- Salt bridges (ligand negative) ---
            for sb in interaction_set.saltbridge_lneg:
                counts["salt_bridge"] += 1
                interactions_rows.append({
                    "id": structure_id,
                    "ligand_name": ligand_name,
                    "ligand_chain": ligand_chain,
                    "interaction_type": "salt_bridge",
                    "protein_residue": sb.restype,
                    "protein_chain": sb.reschain,
                    "protein_residue_number": sb.resnr,
                    "distance": round(sb.distance, 2),
                })

            # --- Salt bridges (protein negative) ---
            for sb in interaction_set.saltbridge_pneg:
                counts["salt_bridge"] += 1
                interactions_rows.append({
                    "id": structure_id,
                    "ligand_name": ligand_name,
                    "ligand_chain": ligand_chain,
                    "interaction_type": "salt_bridge",
                    "protein_residue": sb.restype,
                    "protein_chain": sb.reschain,
                    "protein_residue_number": sb.resnr,
                    "distance": round(sb.distance, 2),
                })

            # --- Halogen bonds ---
            for hx in interaction_set.halogen_bonds:
                counts["halogen_bond"] += 1
                interactions_rows.append({
                    "id": structure_id,
                    "ligand_name": ligand_name,
                    "ligand_chain": ligand_chain,
                    "interaction_type": "halogen_bond",
                    "protein_residue": hx.restype,
                    "protein_chain": hx.reschain,
                    "protein_residue_number": hx.resnr,
                    "distance": round(hx.distance, 2),
                })

            # --- Metal complexes ---
            for mc in interaction_set.metal_complexes:
                counts["metal_complex"] += 1
                interactions_rows.append({
                    "id": structure_id,
                    "ligand_name": ligand_name,
                    "ligand_chain": ligand_chain,
                    "interaction_type": "metal_complex",
                    "protein_residue": mc.restype,
                    "protein_chain": mc.reschain,
                    "protein_residue_number": mc.resnr,
                    "distance": round(mc.distance, 2),
                })

            total = sum(counts.values())
            summary_rows.append({
                "id": structure_id,
                "ligand_name": ligand_name,
                "ligand_chain": ligand_chain,
                **counts,
                "total_interactions": total,
            })

            print(f"  - Binding site {bsid}: {total} interactions")

    except Exception as e:
        print(f"  - ERROR: {e}")
        import traceback
        traceback.print_exc()

    return interactions_rows, summary_rows


def run_plip(config_data: Dict[str, Any]) -> None:
    """
    Run PLIP analysis on all structures.

    Args:
        config_data: Configuration dictionary with paths
    """
    structures_ds = load_datastream(config_data['structures_json'])
    interactions_csv = config_data['interactions_csv']
    summary_csv = config_data['summary_csv']

    print(f"Running PLIP interaction profiling")
    print(f"Structures: {len(structures_ds.ids_expanded)}")

    all_interactions = []
    all_summaries = []
    structure_items = list(iterate_files(structures_ds))
    total = len(structure_items)

    for i, (structure_id, structure_path) in enumerate(structure_items):
        if not os.path.exists(structure_path):
            print(f"Warning: Structure file not found: {structure_path}")
            continue

        # Skip non-PDB files (CIF not supported by PLIP)
        if not structure_path.lower().endswith(".pdb"):
            print(f"Skipping non-PDB file: {structure_path}")
            continue

        print(f"\nProcessing structure {i+1}/{total}: {structure_path}")
        print(f"  - ID: {structure_id}")

        interactions_rows, summary_rows = analyze_structure(structure_id, structure_path)
        all_interactions.extend(interactions_rows)
        all_summaries.extend(summary_rows)

    # Write results
    output_dir = os.path.dirname(interactions_csv)
    os.makedirs(output_dir, exist_ok=True)

    # Interactions table
    if all_interactions:
        df_interactions = pd.DataFrame(all_interactions)
        df_interactions.to_csv(interactions_csv, index=False)
        print(f"\nInteractions written to: {interactions_csv}")
        print(f"Total interactions: {len(df_interactions)}")
        print(df_interactions)
    else:
        # Write empty CSV with correct headers
        columns = ["id", "ligand_name", "ligand_chain", "interaction_type",
                    "protein_residue", "protein_chain", "protein_residue_number",
                    "distance"]
        pd.DataFrame(columns=columns).to_csv(interactions_csv, index=False)
        print(f"\nNo interactions found. Empty table written to: {interactions_csv}")

    # Summary table
    if all_summaries:
        df_summary = pd.DataFrame(all_summaries)
        df_summary.to_csv(summary_csv, index=False)
        print(f"\nSummary written to: {summary_csv}")
        print(df_summary)
    else:
        columns = ["id", "ligand_name", "ligand_chain",
                    "hydrophobic", "hbond", "water_bridge", "salt_bridge",
                    "pi_stacking", "pi_cation", "halogen_bond", "metal_complex",
                    "total_interactions"]
        pd.DataFrame(columns=columns).to_csv(summary_csv, index=False)
        print(f"\nNo binding sites found. Empty summary written to: {summary_csv}")

    print(f"\nPLIP analysis completed successfully!")


def main():
    parser = argparse.ArgumentParser(description='PLIP protein-ligand interaction profiling')
    parser.add_argument('--config', required=True, help='JSON config file with analysis parameters')

    args = parser.parse_args()

    if not os.path.exists(args.config):
        print(f"Error: Config file not found: {args.config}")
        sys.exit(1)

    try:
        with open(args.config, 'r') as f:
            config_data = json.load(f)
    except Exception as e:
        print(f"Error loading config: {e}")
        sys.exit(1)

    required_params = ['structures_json', 'interactions_csv', 'summary_csv']
    for param in required_params:
        if param not in config_data:
            print(f"Error: Missing required parameter: {param}")
            sys.exit(1)

    try:
        run_plip(config_data)
    except Exception as e:
        print(f"Error running PLIP analysis: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
