#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Runtime helper script for ConformationalChange analysis.

This script analyzes protein structures to quantify conformational changes by computing
RMSD using PyMOL's align, super, or cealign methods.
"""

import os
import sys
import argparse
import json
import pandas as pd
from typing import Dict, Any, Optional

# Import PyMOL for structure analysis
import pymol
from pymol import cmd

# Import unified I/O utilities
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from biopipelines.biopipelines_io import load_datastream, iterate_files
from biopipelines.id_map_utils import get_mapped_ids
from biopipelines.sele_utils import sele_group_by_chain


def _sele_to_pymol(obj_name: str, selection_str: str) -> str:
    """Build a PyMOL selection string from a chain-aware selection.

    Converts e.g. ``A1-117+A172-225`` into
    ``obj_name and (chain A and resi 1-117+172-225)``.
    Chainless input like ``1-50+60-80`` becomes
    ``obj_name and resi 1-50+60-80``.
    """
    groups = sele_group_by_chain(selection_str)
    if not groups:
        return f"{obj_name} and resi {selection_str}"

    chain_parts = []
    for chain, resi_str in groups:
        if chain:
            chain_parts.append(f"(chain {chain} and resi {resi_str})")
        else:
            chain_parts.append(f"resi {resi_str}")

    if len(chain_parts) == 1:
        return f"{obj_name} and {chain_parts[0]}"
    return f"{obj_name} and ({' or '.join(chain_parts)})"


def resolve_atoms(atoms: str) -> Optional[str]:
    """
    Resolve an atoms specification to a PyMOL atom name selection string.

    Args:
        atoms: Atom specification. "all" returns None (no filter).
            "backbone" expands to "CA+C+N+O". Any '+'-separated names
            are joined with '+' for PyMOL's ``name`` selector.

    Returns:
        PyMOL name selection string (e.g. "CA+C+N+O"), or None for all atoms
    """
    if atoms == "all":
        return None
    if atoms == "backbone":
        return "CA+C+N+O"
    return atoms


def align_and_compute_rmsd(ref_obj: str, target_obj: str, selection: str,
                           alignment_method: str,
                           atoms: str = "all",
                           pairing: str = "sequence",
                           cycles: int = 5,
                           cutoff: float = 2.0,
                           frame: str = None) -> Dict[str, Any]:
    """
    Superpose target on reference and return RMSD over `selection`.

    Three things here are not PyMOL's defaults, because PyMOL's defaults answer a
    different question than "did this structure fold the way it was designed":

    * `pairing` decides how atoms are put into correspondence. "sequence" uses
      align/super/cealign, which pair by sequence similarity and silently drop
      residues they cannot match — right for homologues, wrong for a design and the
      refold of an inverse-folded sequence, where the sequences differ BY DESIGN.
      "ordered" (cmd.fit matchmaker=-1) pairs the Nth atom with the Nth atom.
    * `cycles` is PyMOL's outlier rejection. At its default of 5 the reported RMSD
      describes only the atoms that survived, so a badly-folded region can be trimmed
      away until what remains fits well. Measured on real designs: a segment 11.4 A
      from its design reported 1.89 A after refinement discarded 95 of 200 atoms.
    * `frame` superposes on one selection and measures another in that frame. A
      segment allowed its own superposition can fit itself well while sitting in
      completely the wrong place: one such segment scored 11.4 A on its own best fit
      and 18.9 A once the cores were aligned.

    Returns RMSD both after and before refinement, so a large gap between them is
    visible in the output rather than silent.
    """
    def _sel(obj, sele):
        base = obj if (sele is None or sele == "all") else _sele_to_pymol(obj, sele)
        names = resolve_atoms(atoms)
        return f"({base}) and name {names}" if names else base

    ref_sel = _sel(ref_obj, selection)
    target_sel = _sel(target_obj, selection)

    # Superpose on `frame` when given, otherwise on the measured selection itself.
    fit_ref = _sel(ref_obj, frame) if frame else ref_sel
    fit_tgt = _sel(target_obj, frame) if frame else target_sel

    n_before = cmd.count_atoms(fit_tgt)
    residues = cmd.count_atoms(fit_tgt + " and name CA")

    if pairing == "sequence":
        if alignment_method == "align":
            r = cmd.align(fit_tgt, fit_ref, cycles=cycles, cutoff=cutoff)
            rmsd, n_after, rmsd_before, n_before = r[0], r[1], r[3], r[4]
            residues = r[6]
        elif alignment_method == "super":
            r = cmd.super(fit_tgt, fit_ref, cycles=cycles, cutoff=cutoff)
            rmsd, n_after, rmsd_before, n_before = r[0], r[1], r[3], r[4]
            residues = r[6]
        elif alignment_method == "cealign":
            r = cmd.cealign(fit_ref, fit_tgt)
            rmsd = rmsd_before = r["RMSD"]
            n_after = n_before = r["alignment_length"]
        else:
            raise ValueError(f"Unknown alignment method: {alignment_method}")
    elif pairing in ("ordered", "identifier"):
        mm = -1 if pairing == "ordered" else 0
        if cmd.count_atoms(fit_tgt) != cmd.count_atoms(fit_ref) and pairing == "ordered":
            raise ValueError(
                f"pairing='ordered' needs the same atom count on both sides, got "
                f"{cmd.count_atoms(fit_tgt)} vs {cmd.count_atoms(fit_ref)}. The structures "
                f"do not share an atom order — use pairing='sequence' or fix the inputs.")
        rmsd = cmd.fit(fit_tgt, fit_ref, matchmaker=mm, cycles=cycles, cutoff=cutoff)
        # cmd.fit reports no post-rejection count and does not narrow the selection,
        # so how many atoms its cycles actually rejected is not observable here --
        # unlike cmd.align, which returns it. n_after therefore equals n_before by
        # construction and atoms_dropped_pct is 0 in this mode; do not read it as
        # evidence that nothing was rejected.
        rmsd_before = cmd.fit(fit_tgt, fit_ref, matchmaker=mm, cycles=0) if cycles else rmsd
        n_after = cmd.count_atoms(fit_tgt)
    else:
        raise ValueError(f"Unknown pairing: {pairing}")

    # With a frame, the fit above moved the target; now measure the requested
    # selection where it landed, without refitting it.
    if frame:
        mm = -1 if pairing in ("sequence", "ordered") else 0
        if mm == -1:
            # matchmaker=-1 pairs the Nth atom with the Nth atom. The count guard
            # in the 'ordered' branch above covers the FRAME selections, not these
            # measured ones, so without this a mismatched selection is paired
            # positionally across two structures that differ by design -- and
            # rms_cur answers 0.000 on unequal counts rather than raising, which is
            # indistinguishable from a perfect superposition.
            n_tgt = cmd.count_atoms(target_sel)
            n_ref = cmd.count_atoms(ref_sel)
            if n_tgt != n_ref:
                raise ValueError(
                    f"selection has {n_tgt} atoms in the target and {n_ref} in the "
                    f"reference, and pairing={pairing!r} pairs them by position. Use "
                    f"pairing='identifier' to pair by chain/residue/atom name, or "
                    f"narrow `selection` to a region both structures share.")
        rmsd = cmd.rms_cur(target_sel, ref_sel, matchmaker=mm)
        rmsd_before = rmsd
        n_after = n_before = cmd.count_atoms(target_sel)
        residues = cmd.count_atoms(target_sel + " and name CA")

    dropped = 100.0 * (1 - n_after / n_before) if n_before else 0.0
    if dropped >= 10.0:
        print(f"  ! refinement dropped {dropped:.0f}% of atoms "
              f"({n_before} -> {n_after}); RMSD {rmsd_before:.2f} -> {rmsd:.2f} A. "
              f"The reported RMSD describes only the atoms that survived.")

    print(f"  - pairing={pairing}, cycles={cycles}, RMSD={rmsd:.3f}, atoms={n_after}")

    return {
        'RMSD': rmsd,
        'num_aligned_atoms': n_after,
        'RMSD_before': rmsd_before,
        'num_atoms_before': n_before,
        'num_residues_aligned': residues,
        'atoms_dropped_pct': round(dropped, 1),
    }


def load_selection_from_table(table_path: str, column_name: str) -> Dict[str, str]:
    """
    Load selection specifications from table CSV file.

    Args:
        table_path: Path to CSV file
        column_name: Column containing selection specifications

    Returns:
        Dictionary mapping structure IDs to selection strings
    """
    if not os.path.exists(table_path):
        raise FileNotFoundError(f"Table file not found: {table_path}")

    df = pd.read_csv(table_path)
    if column_name not in df.columns:
        raise ValueError(f"Column '{column_name}' not found in table. Available columns: {list(df.columns)}")

    # Assuming the first column contains IDs
    id_column = df.columns[0]
    selection_map = {}

    for _, row in df.iterrows():
        structure_id = row[id_column]
        selection_value = row[column_name]
        selection_map[str(structure_id)] = str(selection_value)

    print(f"Loaded selections for {len(selection_map)} structures from {table_path}")

    return selection_map



def analyze_conformational_change(ref_path: str, target_path: str, selection: str,
                                  alignment_method: str,
                                  atoms: str = "all",
                                  pairing: str = "sequence",
                                  cycles: int = 5,
                                  cutoff: float = 2.0,
                                  frame: str = None) -> Optional[Dict[str, Any]]:
    """
    Analyze conformational change between reference and target structures.

    Args:
        ref_path: Path to reference structure file
        target_path: Path to target structure file
        selection: Selection specification (e.g., '10-20+30-40')
        alignment_method: Alignment method ("align", "super", or "cealign")
        atoms: Atom specification ("all", "CA", "backbone", or "CA+CB" etc.)

    Returns:
        Dictionary with RMSD and num_aligned_atoms, or None if failed
    """
    try:
        # Extract structure IDs from filenames for PyMOL object names
        ref_id = os.path.splitext(os.path.basename(ref_path))[0]
        target_id = os.path.splitext(os.path.basename(target_path))[0]

        ref_obj = f"ref_{ref_id}"
        target_obj = f"target_{target_id}"

        # Load structures into PyMOL
        cmd.load(ref_path, ref_obj)
        cmd.load(target_path, target_obj)

        print(f"  - Loaded reference: {ref_obj}")
        print(f"  - Loaded target: {target_obj}")
        print(f"  - Selection: {selection}")
        print(f"  - Atoms: {atoms}")

        # Align and get RMSD from PyMOL
        metrics = align_and_compute_rmsd(ref_obj, target_obj, selection, alignment_method,
                                         atoms, pairing, cycles, cutoff, frame)

        # Clean up PyMOL objects
        cmd.delete(ref_obj)
        cmd.delete(target_obj)

        return metrics

    except Exception as e:
        print(f"  - Error analyzing conformational change: {e}")
        import traceback
        traceback.print_exc()
        return None


def analyze_all_conformational_changes(config_data: Dict[str, Any]) -> None:
    """
    Analyze conformational changes for all structure pairs.

    Args:
        config_data: Configuration dictionary with analysis parameters
    """
    # Load DataStreams using pipe_biopipelines_io
    reference_ds = load_datastream(config_data['reference_structures_json'])
    target_ds = load_datastream(config_data['target_structures_json'])

    selection_config = config_data['selection']
    alignment_method = config_data['alignment_method']
    pairing = config_data.get('pairing', 'sequence')
    cycles = config_data.get('cycles', 5)
    cutoff = config_data.get('cutoff', 2.0)
    frame_config = config_data.get('frame')
    atoms = config_data.get('atoms', 'all')
    output_csv = config_data['output_csv']

    print(f"Analyzing conformational changes")
    print(f"Reference structures: {len(reference_ds.ids_expanded)}")
    print(f"Target structures: {len(target_ds.ids_expanded)}")
    print(f"Selection: {selection_config}")
    print(f"Alignment method: {alignment_method}")
    print(f"Pairing: {pairing} | cycles: {cycles} | cutoff: {cutoff}")
    if frame_config:
        print(f"Frame: {frame_config.get('value', frame_config.get('column_name'))}")
    print(f"Atoms: {atoms}")

    # Initialize PyMOL in headless mode
    pymol.pymol_argv = ['pymol', '-c']
    pymol.finish_launching()
    cmd.set("cartoon_gap_cutoff", 0)

    # Build reference lookup by ID for efficient matching
    reference_files_by_id = {}
    for ref_id, ref_file in iterate_files(reference_ds):
        reference_files_by_id[ref_id] = ref_file

    # Handle selection
    selection_map = {}
    if selection_config['type'] == 'all':
        # Use all CA atoms (whole structure RMSD)
        print("Using all atoms (whole structure RMSD)")
        for target_id in target_ds.ids_expanded:
            selection_map[target_id] = "all"
    elif selection_config['type'] == 'fixed':
        # Fixed selection for all structures
        fixed_selection = selection_config['value']
        print(f"Using fixed selection: {fixed_selection}")
        for target_id in target_ds.ids_expanded:
            selection_map[target_id] = fixed_selection
    else:
        # Load from table
        table_path = selection_config['table_path']
        column_name = selection_config['column_name']
        selection_map = load_selection_from_table(table_path, column_name)

    # Pre-compute ID mappings using get_mapped_ids (handles +components, suffixes, siblings)
    target_ids = list(target_ds.ids_expanded)
    target_to_ref_id = get_mapped_ids(
        source_ids=target_ids,
        target_ids=list(reference_files_by_id.keys()),
        unique=True
    )
    if selection_config['type'] not in ('all', 'fixed'):
        if len(selection_map) == 1:
            single_sele_id, single_sele_value = next(iter(selection_map.items()))
            print(f"Using single selection for all structures: {single_sele_value}")
            target_to_sele_id = None
            for target_id in target_ids:
                selection_map[target_id] = single_sele_value
        else:
            target_to_sele_id = get_mapped_ids(
                source_ids=target_ids,
                target_ids=list(selection_map.keys()),
                unique=True
            )
    else:
        target_to_sele_id = None

    # Per-structure frame selections, mapped to target ids the same way: a frame
    # column typically lives on a design-level table while the targets are the
    # per-sequence folds derived from it.
    frame_selections = {}
    target_to_frame_id = None
    if frame_config and frame_config['type'] == 'table_column':
        frame_selections = load_selection_from_table(frame_config['table_path'],
                                                     frame_config['column_name'])
        if len(frame_selections) == 1:
            single_frame_value = next(iter(frame_selections.values()))
            for target_id in target_ids:
                frame_selections[target_id] = single_frame_value
        else:
            target_to_frame_id = get_mapped_ids(
                source_ids=target_ids,
                target_ids=list(frame_selections.keys()),
                unique=True
            )

    # Determine if reference is single or multiple
    use_single_reference = len(reference_ds.ids_expanded) == 1
    if use_single_reference:
        single_ref_id, single_ref_path = next(iterate_files(reference_ds))
        print(f"Using single reference structure: {single_ref_path}")
    else:
        print(f"Using paired reference structures")

    # Ensure we have compatible number of structures
    if not use_single_reference and len(reference_ds.ids_expanded) != len(target_ds.ids_expanded):
        print(f"Warning: Reference structures ({len(reference_ds.ids_expanded)}) and target structures ({len(target_ds.ids_expanded)}) count mismatch")

    # Process structure pairs using iterate_files for proper ID-file matching
    results = []
    target_items = list(iterate_files(target_ds))

    for i, (target_id, target_path) in enumerate(target_items):
        if not os.path.exists(target_path):
            print(f"Warning: Target structure file not found: {target_path}")
            continue

        # Get reference structure (single or paired by ID)
        if use_single_reference:
            ref_path = single_ref_path
        else:
            # Match by ID using pre-computed mapping (handles +components, suffixes)
            matched_ref_id = target_to_ref_id.get(target_id)
            if matched_ref_id is not None:
                ref_path = reference_files_by_id[matched_ref_id]
                if matched_ref_id != target_id:
                    print(f"  - Matched target '{target_id}' to reference '{matched_ref_id}'")
            else:
                print(f"Warning: No matching reference for target ID: {target_id}")
                continue

        if not os.path.exists(ref_path):
            print(f"Warning: Reference structure file not found: {ref_path}")
            continue

        print(f"\nProcessing structure pair {i+1}/{len(target_items)}")
        print(f"Reference: {ref_path}")
        print(f"Target: {target_path}")
        print(f"ID: {target_id}")

        # Get selection for this structure
        if target_to_sele_id is not None:
            matched_sele_id = target_to_sele_id.get(target_id)
            if matched_sele_id is not None:
                selection = selection_map[matched_sele_id]
                if matched_sele_id != target_id:
                    print(f"  - Matched selection ID: {target_id} -> {matched_sele_id}")
            else:
                selection = None
        else:
            selection = selection_map.get(target_id)
        if not selection:
            print(f"  - Warning: No selection found for structure ID: {target_id}")
            print(f"    Available IDs in selection map: {list(selection_map.keys())[:5]}...")
            continue

        # Analyze conformational change
        frame_sel = None
        if frame_config:
            if frame_config['type'] == 'fixed':
                frame_sel = frame_config['value']
            else:
                frame_key = (target_to_frame_id.get(target_id, target_id)
                             if target_to_frame_id else target_id)
                frame_sel = frame_selections.get(frame_key)
                if frame_sel is None:
                    print(f"Warning: no frame selection for {target_id}; skipping")
                    continue
        metrics = analyze_conformational_change(ref_path, target_path, selection,
                                                alignment_method, atoms,
                                                pairing, cycles, cutoff, frame_sel)

        if metrics is None:
            continue

        # Store result using the proper ID from DataStream
        result = {
            'id': target_id,
            'reference_structure': ref_path,
            'target_structure': target_path,
            'selection': selection,
            'num_aligned_atoms': metrics['num_aligned_atoms'],
            'RMSD': metrics['RMSD']
        }
        results.append(result)

        print(f"  - Result stored")

    # Create DataFrame and save
    if results:
        df = pd.DataFrame(results)

        # Create output directory
        output_dir = os.path.dirname(output_csv)
        print(f"Creating output directory: {output_dir}")
        os.makedirs(output_dir, exist_ok=True)

        # Save results
        print(f"Writing results to: {output_csv}")
        df.to_csv(output_csv, index=False)

        print(f"\nConformational change analysis completed successfully!")
        print(f"Analyzed {len(results)} structure pairs")
        print(f"Results saved to: {output_csv}")
        print(f"\nResults summary:")
        print(df)

        # Statistics
        values = df['RMSD'].dropna()
        if len(values) > 0:
            print(f"\nRMSD statistics:")
            print(f"  Min: {values.min():.3f}")
            print(f"  Max: {values.max():.3f}")
            print(f"  Mean: {values.mean():.3f}")
            print(f"  Std: {values.std():.3f}")
    else:
        raise ValueError("No valid results generated - check structure files and selections")

    # Clean up PyMOL
    cmd.quit()


def main():
    parser = argparse.ArgumentParser(description='Analyze conformational changes between reference and target structures')
    parser.add_argument('--config', required=True, help='JSON config file with analysis parameters')

    args = parser.parse_args()

    # Load configuration
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
    required_params = ['reference_structures_json', 'target_structures_json', 'selection', 'alignment_method', 'output_csv']
    for param in required_params:
        if param not in config_data:
            print(f"Error: Missing required parameter: {param}")
            sys.exit(1)

    try:
        analyze_all_conformational_changes(config_data)

    except Exception as e:
        print(f"Error analyzing conformational changes: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
