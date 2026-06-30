# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.


import json
import os
import sys

# Import PDB parser and I/O utilities
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from biopipelines.pdb_parser import parse_pdb_file
from biopipelines.biopipelines_io import load_datastream, iterate_files, load_table, lookup_table_value

# Establish residues for which to generate a sequence with ProteinMPNN.
# <config_json> keys: structures_json, FIXED, DESIGNED, FIXED_CHAIN,
# fixed_jsonl_file, sele_csv_file. FIXED/DESIGNED are a selection or '-'.
if len(sys.argv) != 2:
    print("Usage: python pipe_pmpnn_fixed_positions.py <config_json>", file=sys.stderr)
    sys.exit(1)

with open(sys.argv[1]) as _cfg_f:
    cfg = json.load(_cfg_f)
structures_json = cfg["structures_json"]
FIXED = cfg["FIXED"]
DESIGNED = cfg["DESIGNED"]
FIXED_CHAIN = cfg["FIXED_CHAIN"]
fixed_jsonl_file = cfg["fixed_jsonl_file"]
sele_csv_file = cfg["sele_csv_file"]

from biopipelines.sele_utils import sele_to_list as _sele_to_list_chain_aware, list_to_sele

def sele_to_list(s):
    """Convert selection string to flat list of residue numbers (strips chain info)."""
    return [r for _, r in _sele_to_list_chain_aware(s)]

def sele_to_dict(s):
    """Convert selection string to dict of chain -> sorted residue list.

    For chainless input, all residues go under key ''.
    """
    pairs = _sele_to_list_chain_aware(s)
    d = {}
    for chain, resnum in pairs:
        d.setdefault(chain, []).append(resnum)
    return d

def get_protein_chains_from_pdb(pdb_path):
    """
    Get all protein chain identifiers from a PDB file.

    Args:
        pdb_path: Path to PDB file

    Returns:
        Sorted list of unique chain identifiers that contain protein residues
    """
    standard_residues = {
        'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
        'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'
    }

    atoms = parse_pdb_file(pdb_path)
    chains = set()

    for atom in atoms:
        if atom.res_name in standard_residues:
            chains.add(atom.chain)

    return sorted(list(chains))


def get_all_chains_from_pdb(pdb_path):
    """
    Get every chain identifier present in a PDB file, regardless of residue type.

    ProteinMPNN's parse_multiple_chains.py records every chain in the structure
    (protein, DNA, RNA, hetero), and tied_featurize requires fixed_positions_dict
    to have an entry for each. Restricting to protein chains here would raise
    KeyError at inference time when a structure carries e.g. a bound DNA strand.

    Args:
        pdb_path: Path to PDB file

    Returns:
        Sorted list of unique chain identifiers
    """
    atoms = parse_pdb_file(pdb_path)
    return sorted({atom.chain for atom in atoms})


def get_residues_in_chain(pdb_path, chain):
    """
    Get every residue number present in the given chain, regardless of residue type.

    Used to fully fix non-protein chains so ProteinMPNN does not redesign them.

    Args:
        pdb_path: Path to PDB file
        chain: Chain identifier

    Returns:
        Sorted list of residue numbers in chain
    """
    atoms = parse_pdb_file(pdb_path)
    residues = {atom.res_num for atom in atoms if atom.chain == chain}
    return sorted(residues)


def get_protein_residues_from_pdb(pdb_path, chain):
    """
    Get all protein residue numbers for a specific chain from PDB file.

    Args:
        pdb_path: Path to PDB file
        chain: Chain identifier (e.g., 'A')

    Returns:
        Sorted list of residue numbers for the specified chain
    """
    # Standard amino acid names
    standard_residues = {
        'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
        'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'
    }

    atoms = parse_pdb_file(pdb_path)
    residues = set()

    for atom in atoms:
        if atom.chain == chain and atom.res_name in standard_residues:
            residues.add(atom.res_num)

    return sorted(list(residues))

def compute_complement(all_residues, redesigned_residues):
    """
    Compute complement: all_residues - redesigned_residues.

    Args:
        all_residues: List of all protein residue numbers
        redesigned_residues: List of residues to redesign

    Returns:
        Sorted list of residues NOT in redesigned_residues
    """
    redesigned_set = set(redesigned_residues)
    complement = [res for res in all_residues if res not in redesigned_set]
    return sorted(complement)

def resolve_table_reference(reference, design_ids):
    """
    Resolve table reference to per-design selections.

    Args:
        reference: Either a table reference like "TABLE_REFERENCE:path:column" or direct PyMOL selection
        design_ids: List of design IDs from DataStream

    Returns:
        Dictionary mapping design IDs to position lists
    """
    if not reference.startswith("TABLE_REFERENCE:"):
        # Direct PyMOL selection - same for all designs
        return {design_id: sele_to_list(reference) for design_id in design_ids}

    # Use pipe_biopipelines_io to load table and column
    table, column_name = load_table(reference)
    positions_per_design = {}

    for design_id in design_ids:
        try:
            selection_value = lookup_table_value(table, design_id, column_name)
            if pd.notna(selection_value) and selection_value != '':
                positions_per_design[design_id] = sele_to_list(str(selection_value))
            else:
                positions_per_design[design_id] = []
        except KeyError:
            print(f"Warning: No table entry found for {design_id} in column {column_name}")
            positions_per_design[design_id] = []

    return positions_per_design

import os
import pandas as pd

# Load DataStream and get (id, file) pairs
structures_ds = load_datastream(structures_json)
design_entries = list(iterate_files(structures_ds))  # List of (design_id, pdb_file) tuples
design_ids = [entry[0] for entry in design_entries]
design_files = {entry[0]: entry[1] for entry in design_entries}  # Map id -> file path

# Auto-detect chain if requested
if FIXED_CHAIN == "auto":
    # Use the first structure to detect chains
    first_pdb = design_files[design_ids[0]]
    all_chains = get_protein_chains_from_pdb(first_pdb)
    if len(all_chains) == 1:
        FIXED_CHAIN = all_chains[0]
        print(f"Auto-detected single protein chain: {FIXED_CHAIN}")
    elif len(all_chains) > 1:
        FIXED_CHAIN = all_chains[0]
        print(f"Multiple protein chains found: {', '.join(all_chains)}. Using first chain: {FIXED_CHAIN}")
    else:
        FIXED_CHAIN = "A"
        print(f"Warning: No protein chains detected, defaulting to chain A")

# Sanitize '-' placeholders (used when no positions are specified)
FIXED = '' if FIXED == '-' else FIXED
DESIGNED = '' if DESIGNED == '-' else DESIGNED

fixed_dict = dict()
mobile_dict = dict()

# Resolve table references if present
fixed_per_design = resolve_table_reference(FIXED, design_ids) if FIXED else {design_id: [] for design_id in design_ids}
designed_per_design = resolve_table_reference(DESIGNED, design_ids) if DESIGNED else {design_id: [] for design_id in design_ids}

for design_id in design_ids:
    fixed_dict[design_id] = dict()
    mobile_dict[design_id] = dict()

    # Store original mobile/designed positions for documentation
    mobile_dict[design_id][FIXED_CHAIN] = designed_per_design[design_id]

    # Compute what ProteinMPNN should keep fixed:
    # Union of explicit fixed + complement of redesigned
    pdb_path = design_files[design_id]

    # Start with explicit fixed positions
    final_fixed = list(fixed_per_design[design_id])

    # Get all protein residues from PDB
    all_residues = get_protein_residues_from_pdb(pdb_path, FIXED_CHAIN)

    if designed_per_design[design_id]:
        # If redesigned positions are specified, add their complement to fixed
        complement = compute_complement(all_residues, designed_per_design[design_id])
        final_fixed = sorted(list(set(final_fixed + complement)))

        print(f"Design: {design_id}, Explicit Fixed: {list_to_sele(fixed_per_design[design_id]) if fixed_per_design[design_id] else ''}, Redesigned: {list_to_sele(designed_per_design[design_id])}, Final Fixed (to ProteinMPNN): {list_to_sele(final_fixed)}")
    elif final_fixed:
        # No redesigned specified but fixed specified, use fixed as-is
        print(f"Design: {design_id}, Fixed: {list_to_sele(final_fixed)}, Redesigned: all")
    else:
        # Neither redesigned nor fixed specified - redesign everything (no fixed positions)
        print(f"Design: {design_id}, No fixed/redesigned specified, redesigning all residues")

    # Store final fixed positions (what ProteinMPNN will use)
    fixed_dict[design_id][FIXED_CHAIN] = final_fixed

# Ensure every chain present in the structure (including DNA/RNA/hetero)
# has an entry in fixed_dict — ProteinMPNN's tied_featurize requires it.
#
# Semantics for unmentioned chains depend on which selector the user passed:
#   * fixed selection given      -> unmentioned chains = [] (fully redesigned)
#   * redesigned selection given -> unmentioned chains = all their residues fixed
#                                   (so ProteinMPNN leaves them untouched)
#   * neither given              -> all chains = [] (everything redesigned)
designed_was_passed = bool(DESIGNED)
for design_id in design_ids:
    pdb_path = design_files[design_id]
    all_chains = get_all_chains_from_pdb(pdb_path)
    for ch in all_chains:
        if ch in fixed_dict[design_id]:
            continue
        if designed_was_passed:
            fixed_dict[design_id][ch] = get_residues_in_chain(pdb_path, ch)
        else:
            fixed_dict[design_id][ch] = []

with open(fixed_jsonl_file,"w") as jsonl_file:
    #Python converts dictionaries to string having keys inside '', json only recognises ""
    jsonl_file.write(str(fixed_dict).replace("\'","\""))

with open(sele_csv_file,"w") as csv_file:
    csv_file.write("id,fixed,mobile")
    for id in fixed_dict.keys():
        fixed = list_to_sele(fixed_dict[id][FIXED_CHAIN]) if FIXED_CHAIN in fixed_dict[id].keys() else ""
        mobile = list_to_sele(mobile_dict[id][FIXED_CHAIN]) if FIXED_CHAIN in fixed_dict[id].keys() else ""
        csv_file.write("\n")
        csv_file.write(f"{id},{fixed},{mobile}")
        
print(f"Fixed positions written to: {fixed_jsonl_file}")
print(f"Selections summary written to: {sele_csv_file}")