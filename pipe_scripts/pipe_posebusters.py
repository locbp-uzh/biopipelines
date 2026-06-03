#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Runtime helper script for PoseBusters validation.

Extracts ligands and proteins from predicted structures, runs PoseBusters
quality checks, and outputs a CSV with boolean pass/fail columns.
"""

import os
import sys
import argparse
import json
import tempfile
import pandas as pd
from typing import Dict, List, Any, Optional, Tuple

# Import unified I/O utilities
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from biopipelines.biopipelines_io import load_datastream, iterate_files
from biopipelines.ligand_utils import resolve_ligand_code, resolve_ligand_smiles as _resolve_ligand_smiles


def extract_ligand_and_protein(structure_path, ligand_name, output_dir, structure_id,
                               ligand_smiles=None):
    """
    Extract ligand to SDF and protein to PDB from a PDB/CIF structure.

    Uses gemmi for structure parsing (works for both PDB and CIF) and RDKit
    for bond perception and SDF writing.

    Args:
        structure_path: Path to input structure (PDB or CIF)
        ligand_name: 3-letter residue code for the ligand
        output_dir: Directory to write extracted files
        structure_id: Identifier for naming output files
        ligand_smiles: Optional SMILES template for bond-order assignment. Large
            conjugated/charged ligands (e.g. Si-rhodamine dyes) cannot be bond-
            perceived from coordinates alone (rdDetermineBonds fails on charge),
            so a template is required for PoseBusters to validate them.

    Returns:
        List of (ligand_sdf_path, protein_pdb_path, suffix) tuples.
        Multiple tuples if multiple copies of the ligand exist.
    """
    import gemmi
    from rdkit import Chem
    from rdkit.Chem import rdDetermineBonds

    # Read structure with gemmi (handles PDB and CIF)
    structure = gemmi.read_structure(structure_path)
    if len(structure) == 0:
        print(f"  [WARNING] Empty structure: {structure_path}")
        return []

    model = structure[0]

    # Collect ligand residues and protein atoms
    ligand_residues = []
    protein_lines = []

    for chain in model:
        for residue in chain:
            if residue.name == ligand_name:
                ligand_residues.append((chain.name, residue))
            else:
                # Write all non-ligand residues to protein PDB
                for atom in residue:
                    record = "HETATM" if not residue.entity_type == gemmi.EntityType.Polymer else "ATOM  "
                    protein_lines.append(
                        f"{record}{atom.serial:5d} {atom.name:<4s} {residue.name:>3s} "
                        f"{chain.name:1s}{residue.seqid.num:4d}    "
                        f"{atom.pos.x:8.3f}{atom.pos.y:8.3f}{atom.pos.z:8.3f}"
                        f"{atom.occ:6.2f}{atom.b_iso:6.2f}          "
                        f"{atom.element.name:>2s}\n"
                    )

    if not ligand_residues:
        print(f"  [WARNING] No ligand '{ligand_name}' found in {structure_path}")
        return []

    # Write protein PDB
    protein_pdb_path = os.path.join(output_dir, f"{structure_id}_protein.pdb")
    with open(protein_pdb_path, 'w') as f:
        f.writelines(protein_lines)
        f.write("END\n")

    # Extract each ligand copy
    results = []
    for lig_idx, (chain_name, residue) in enumerate(ligand_residues):
        suffix = f"_lig{lig_idx + 1}" if len(ligand_residues) > 1 else ""

        # Build PDB block for this ligand residue
        lig_lines = []
        for atom in residue:
            lig_lines.append(
                f"HETATM{atom.serial:5d} {atom.name:<4s} {residue.name:>3s} "
                f"{chain_name:1s}{residue.seqid.num:4d}    "
                f"{atom.pos.x:8.3f}{atom.pos.y:8.3f}{atom.pos.z:8.3f}"
                f"{atom.occ:6.2f}{atom.b_iso:6.2f}          "
                f"{atom.element.name:>2s}\n"
            )
        lig_pdb_block = "".join(lig_lines) + "END\n"

        ligand_sdf_path = os.path.join(output_dir, f"{structure_id}{suffix}_ligand.sdf")
        # Prefer the shared SMILES-templated conversion (bond orders from the
        # template, coordinates kept) — coordinate-only perception fails on
        # charged/conjugated dyes. Write the extracted block to a temp PDB and
        # route through write_ligand_sdf; fall back to rdDetermineBonds only if
        # no SMILES is available.
        lig_tmp_pdb = os.path.join(output_dir, f"{structure_id}{suffix}_ligand_raw.pdb")
        with open(lig_tmp_pdb, "w") as _f:
            _f.write(lig_pdb_block)
        try:
            from biopipelines.ligand_utils import write_ligand_sdf
            write_ligand_sdf(lig_tmp_pdb, ligand_sdf_path, smiles=ligand_smiles)
        except Exception as e:
            print(f"  [WARNING] templated SDF write failed for {ligand_name} copy {lig_idx + 1}: {e}; "
                  f"falling back to coordinate bond perception")
            mol = Chem.MolFromPDBBlock(lig_pdb_block, removeHs=False, sanitize=False)
            if mol is None:
                print(f"  [WARNING] RDKit could not parse ligand {ligand_name} (copy {lig_idx + 1})")
                continue
            try:
                rdDetermineBonds.DetermineBonds(mol)
            except Exception as e2:
                print(f"  [WARNING] Bond perception failed for {ligand_name} copy {lig_idx + 1}: {e2}")
            writer = Chem.SDWriter(ligand_sdf_path)
            writer.write(mol)
            writer.close()

        results.append((ligand_sdf_path, protein_pdb_path, suffix))

    return results


def extract_reference_ligand(reference_path, ligand_code, output_dir):
    """
    Extract the reference ligand from a reference structure for redock mode.

    Args:
        reference_path: Path to reference structure (PDB or CIF)
        ligand_code: Residue code of the ligand in the reference
        output_dir: Directory to write extracted file

    Returns:
        Path to reference ligand SDF, or None if extraction failed.
    """
    import gemmi
    from rdkit import Chem
    from rdkit.Chem import rdDetermineBonds

    # Check if it's already an SDF file
    ext = os.path.splitext(reference_path)[1].lower()
    if ext in ('.sdf', '.mol'):
        return reference_path

    # Parse structure and extract ligand
    structure = gemmi.read_structure(reference_path)
    if len(structure) == 0:
        print(f"  [WARNING] Empty reference structure: {reference_path}")
        return None

    model = structure[0]

    # Find the first matching residue
    for chain in model:
        for residue in chain:
            if residue.name == ligand_code:
                lig_lines = []
                for atom in residue:
                    lig_lines.append(
                        f"HETATM{atom.serial:5d} {atom.name:<4s} {residue.name:>3s} "
                        f"{chain.name:1s}{residue.seqid.num:4d}    "
                        f"{atom.pos.x:8.3f}{atom.pos.y:8.3f}{atom.pos.z:8.3f}"
                        f"{atom.occ:6.2f}{atom.b_iso:6.2f}          "
                        f"{atom.element.name:>2s}\n"
                    )
                lig_pdb_block = "".join(lig_lines) + "END\n"

                mol = Chem.MolFromPDBBlock(lig_pdb_block, removeHs=False, sanitize=False)
                if mol is None:
                    print(f"  [WARNING] RDKit could not parse reference ligand {ligand_code}")
                    return None

                try:
                    rdDetermineBonds.DetermineBonds(mol)
                except Exception as e:
                    print(f"  [WARNING] Bond perception failed for reference ligand: {e}")

                ref_sdf_path = os.path.join(output_dir, "reference_ligand.sdf")
                writer = Chem.SDWriter(ref_sdf_path)
                writer.write(mol)
                writer.close()
                return ref_sdf_path

    print(f"  [WARNING] Reference ligand '{ligand_code}' not found in {reference_path}")
    return None


def run_posebusters(config_data: Dict[str, Any]) -> None:
    """
    Run PoseBusters validation on structures.

    Args:
        config_data: Configuration dictionary with validation parameters
    """
    from posebusters import PoseBusters

    # Load parameters
    structures_ds = load_datastream(config_data['structures_json'])
    ligand_name = resolve_ligand_code(config_data['ligand_json'])
    # SMILES template (if the compounds stream carries one) — needed for ligands
    # whose bond order/charge can't be perceived from coordinates alone.
    ligand_smiles = _resolve_ligand_smiles(config_data['ligand_json'])
    mode = config_data['mode']
    reference_pdb_json = config_data.get('reference_pdb')
    # The reference-structure residue code defaults to the ligand code.
    reference_ligand_code = ligand_name
    output_csv = config_data['output_csv']
    execution_dir = config_data['execution_dir']
    check_columns = list(config_data.get('check_columns') or [])

    print(f"PoseBusters validation")
    print(f"Structures: {len(structures_ds.ids_expanded)}")
    print(f"Ligand: {ligand_name}")
    print(f"Mode: {mode}")

    # Load reference structure for redock mode
    reference_sdf = None
    if mode == "redock" and reference_pdb_json:
        ref_ds = load_datastream(reference_pdb_json)
        ref_files = list(iterate_files(ref_ds))
        if ref_files:
            ref_id, ref_path = ref_files[0]
            print(f"Reference structure: {ref_path}")
            reference_sdf = extract_reference_ligand(ref_path, reference_ligand_code, execution_dir)
            if reference_sdf is None:
                print(f"[ERROR] Could not extract reference ligand '{reference_ligand_code}' from {ref_path}")
                sys.exit(1)
            print(f"Reference ligand SDF: {reference_sdf}")

    # Process structures
    results = []
    structure_items = list(iterate_files(structures_ds))
    total = len(structure_items)
    for i, (structure_id, structure_path) in enumerate(structure_items):
        if not os.path.exists(structure_path):
            print(f"Warning: Structure file not found: {structure_path}")
            continue

        print(f"\nProcessing structure {i + 1}/{total}: {structure_path}")
        print(f"  ID: {structure_id}")

        # Extract ligand(s) and protein
        struct_work_dir = os.path.join(execution_dir, structure_id)
        os.makedirs(struct_work_dir)

        extractions = extract_ligand_and_protein(structure_path, ligand_name, struct_work_dir, structure_id,
                                                  ligand_smiles=ligand_smiles)

        if not extractions:
            print(f"  [WARNING] No ligand extracted, skipping structure")
            results.append({
                'id': structure_id,
                'source_structure': structure_path,
                'all_pass': False
            })
            continue

        for ligand_sdf, protein_pdb, suffix in extractions:
            entry_id = f"{structure_id}{suffix}"
            print(f"  Running PoseBusters ({mode}) for {entry_id}")
            print(f"    Ligand SDF: {ligand_sdf}")
            print(f"    Protein PDB: {protein_pdb}")

            try:
                buster = PoseBusters(config=mode)

                bust_kwargs = {
                    "mol_pred": ligand_sdf,
                    "mol_cond": protein_pdb,
                }
                if mode == "redock" and reference_sdf:
                    bust_kwargs["mol_true"] = reference_sdf

                df_result = buster.bust(**bust_kwargs)

                # Convert results to a flat dict
                result_row = {
                    'id': entry_id,
                    'source_structure': structure_path,
                }

                # PoseBusters returns a DataFrame with boolean columns
                if not df_result.empty:
                    for col in df_result.columns:
                        # Flatten MultiIndex columns if present
                        col_name = col if isinstance(col, str) else "_".join(str(c) for c in col if c)
                        val = df_result.iloc[0][col]
                        result_row[col_name] = val

                # Compute all_pass over the user-selected check subset only.
                # Selected columns absent from this row (e.g. an identity check that
                # PoseBusters did not emit) are ignored rather than treated as failures.
                # Values come from a pandas DataFrame, so they are numpy.bool_, NOT
                # Python bool — accept both. A check is "present" if it is a real
                # boolean (True/False); pandas NA / None / non-bool are skipped.
                import numpy as _np
                def _is_boolish(v):
                    return isinstance(v, (bool, _np.bool_))
                selected = [c for c in check_columns
                            if c in result_row and _is_boolish(result_row[c])]
                if selected:
                    result_row['all_pass'] = all(bool(result_row[c]) for c in selected)
                else:
                    result_row['all_pass'] = False

                results.append(result_row)
                print(f"    all_pass: {result_row['all_pass']} (over {len(selected)} selected checks)")

            except Exception as e:
                print(f"    [ERROR] PoseBusters failed for {entry_id}: {e}")
                import traceback
                traceback.print_exc()
                results.append({
                    'id': entry_id,
                    'source_structure': structure_path,
                    'all_pass': False
                })

    # Assemble and save results
    if results:
        df = pd.DataFrame(results)

        # Reorder: id, source_structure, included checks, all_pass, excluded checks, anything else.
        leading = [c for c in ("id", "source_structure") if c in df.columns]
        included = [c for c in check_columns if c in df.columns and c not in leading]
        all_pass_col = ["all_pass"] if "all_pass" in df.columns else []
        check_set = set(check_columns)
        # Excluded = any other boolean check column PoseBusters emitted that wasn't selected.
        excluded = [c for c in df.columns
                    if c not in leading and c not in included and c not in all_pass_col
                    and c != "file" and c not in check_set]
        ordered = leading + included + all_pass_col + excluded
        # Append anything we somehow missed (defensive — pd preserves any leftover cols).
        leftover = [c for c in df.columns if c not in ordered]
        df = df[ordered + leftover]

        df.to_csv(output_csv, index=False)

        print(f"\nPoseBusters validation completed!")
        print(f"Validated {len(results)} ligand poses from {total} structures")
        print(f"Results saved to: {output_csv}")
        print(f"\nResults summary:")
        print(df.to_string())

        # Summary statistics
        if 'all_pass' in df.columns:
            n_pass = df['all_pass'].sum()
            n_total = len(df)
            print(f"\nOverall: {n_pass}/{n_total} poses passed all checks ({100 * n_pass / n_total:.1f}%)")
    else:
        raise ValueError("No valid results generated - check structure files and ligand name")


def main():
    parser = argparse.ArgumentParser(description='Run PoseBusters validation on structures')
    parser.add_argument('--config', required=True, help='JSON config file with validation parameters')

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

    # Validate required parameters
    required_params = ['structures_json', 'ligand_json', 'mode', 'output_csv']
    for param in required_params:
        if param not in config_data:
            print(f"Error: Missing required parameter: {param}")
            sys.exit(1)

    try:
        run_posebusters(config_data)
    except Exception as e:
        print(f"Error running PoseBusters validation: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
