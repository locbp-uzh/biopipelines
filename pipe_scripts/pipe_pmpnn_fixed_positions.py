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

from biopipelines.sele_utils import sele_to_list as _sele_to_list_chain_aware, list_to_sele, chain_aware_sele

sele_to_list = _sele_to_list_chain_aware

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

def resolve_table_reference(reference, design_ids, map_table_paths=None):
    """
    Resolve table reference to per-design selections.

    Args:
        reference: Either a table reference like "TABLE_REFERENCE:path:column" or direct PyMOL selection
        design_ids: List of design IDs from DataStream
        map_table_paths: Input stream's map_table, so ids renamed upstream still
            resolve against the table's original id space

    Returns:
        Dictionary mapping design IDs to lists of (chain, resnum) tuples. Chain is '' for
        a selection written without one.
    """
    if not reference.startswith("TABLE_REFERENCE:"):
        # Direct PyMOL selection - same for all designs
        return {design_id: sele_to_list(reference) for design_id in design_ids}

    # Use pipe_biopipelines_io to load table and column
    table, column_name = load_table(reference)
    positions_per_design = {}

    for design_id in design_ids:
        try:
            selection_value = lookup_table_value(table, design_id, column_name,
                                                 map_table_paths=map_table_paths)
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

# Sanitize '-' placeholders (used when no positions are specified)
FIXED = '' if FIXED == '-' else FIXED
DESIGNED = '' if DESIGNED == '-' else DESIGNED


def as_per_chain(selection):
    """A selection as {key chain or '': selection}, whichever form the user wrote.

    A dict already says which chain each part belongs to; anything else is one selection
    under '', whose own chain prefixes (if any) decide where its residues land.
    """
    if isinstance(selection, dict):
        return {c: v for c, v in selection.items() if v not in ('', '-', None)}
    return {'': selection} if selection else {}


structures_maps = [structures_ds.map_table] if structures_ds.map_table else None


def bucket_by_chain(selection, design_ids, maps):
    """{chain: {design_id: [resnum]}} plus whether any residue named no chain.

    A residue's own prefix wins; failing that the dict key it was written under; failing
    that the chain resolved from the structure. Because the parse is chain-aware, a
    qualified string ("A10-20+B5") needs no dict.
    """
    buckets, chainless = {}, False
    for key_chain, sel in as_per_chain(selection).items():
        for design_id, pairs in resolve_table_reference(sel, design_ids, maps).items():
            for chain, resnum in pairs:
                target = chain or key_chain
                if not target:
                    chainless = True
                buckets.setdefault(target, {}).setdefault(design_id, []).append(resnum)
    return buckets, chainless


fixed_raw, fixed_chainless = bucket_by_chain(FIXED, design_ids, structures_maps)
designed_raw, designed_chainless = bucket_by_chain(DESIGNED, design_ids, structures_maps)

# The chain an unqualified residue attaches to. "auto" means read it off the structure;
# several protein chains then make it genuinely ambiguous, which is an error rather than a
# reason to pick the first. A qualified selection never reaches this.
if FIXED_CHAIN == "auto":
    protein_chains = get_protein_chains_from_pdb(design_files[design_ids[0]])
    named = sorted((set(fixed_raw) | set(designed_raw)) - {''})
    if len(protein_chains) == 1:
        FIXED_CHAIN = protein_chains[0]
    elif len(protein_chains) > 1:
        if fixed_chainless or designed_chainless:
            raise ValueError(
                f"{len(protein_chains)} protein chains ({'+'.join(protein_chains)}) and a "
                f"position selection that names no chain: which chain do the residues "
                f"belong to? Qualify them (\"{protein_chains[0]}10-20\"), name them per "
                f"chain (redesigned={{'{protein_chains[0]}': ...}}), or restrict the step "
                f"with chains=\"{protein_chains[0]}\".")
        FIXED_CHAIN = named[0] if named else protein_chains[0]
    else:
        FIXED_CHAIN = "A"
        print(f"Warning: No protein chains detected, defaulting to chain A")


def resolve_chainless(buckets):
    """Fold the unqualified bucket into the chain resolved for it."""
    unqualified = buckets.pop('', None)
    if unqualified:
        target = buckets.setdefault(FIXED_CHAIN, {})
        for design_id, resnums in unqualified.items():
            target.setdefault(design_id, []).extend(resnums)
    return {c: {d: sorted(set(per.get(d, []))) for d in design_ids}
            for c, per in buckets.items()}


fixed_by_chain = resolve_chainless(fixed_raw)
designed_by_chain = resolve_chainless(designed_raw)
selected_chains = sorted(set(fixed_by_chain) | set(designed_by_chain)) or [FIXED_CHAIN]
empty = {design_id: [] for design_id in design_ids}

fixed_dict = dict()
mobile_dict = dict()

for design_id in design_ids:
    fixed_dict[design_id] = dict()
    mobile_dict[design_id] = dict()
    pdb_path = design_files[design_id]

    for chain in selected_chains:
        fixed_here = fixed_by_chain.get(chain, empty)[design_id]
        designed_here = designed_by_chain.get(chain, empty)[design_id]

        # Store original mobile/designed positions for documentation
        mobile_dict[design_id][chain] = designed_here

        # Compute what ProteinMPNN should keep fixed:
        # Union of explicit fixed + complement of redesigned
        final_fixed = list(fixed_here)
        all_residues = get_protein_residues_from_pdb(pdb_path, chain)

        if designed_here:
            complement = compute_complement(all_residues, designed_here)
            final_fixed = sorted(list(set(final_fixed + complement)))
            print(f"Design: {design_id} chain {chain}, Explicit Fixed: {list_to_sele(fixed_here) if fixed_here else ''}, Redesigned: {list_to_sele(designed_here)}, Final Fixed (to ProteinMPNN): {list_to_sele(final_fixed)}")
        elif final_fixed:
            print(f"Design: {design_id} chain {chain}, Fixed: {list_to_sele(final_fixed)}, Redesigned: all")
        else:
            print(f"Design: {design_id} chain {chain}, No fixed/redesigned specified, redesigning all residues")

        # Store final fixed positions (what ProteinMPNN will use)
        fixed_dict[design_id][chain] = final_fixed

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

# ProteinMPNN indexes a mask sized to the chain's residue COUNT:
#     fixed_position_mask = np.ones(chain_length)
#     fixed_position_mask[np.array(fixed_pos_list)-1] = 0.0
# so a position is the residue's 1-based rank within its chain, not its PDB number. The two
# coincide only for a chain numbered 1..N with no gaps, which is what Boltz2/RFD3 emit — hence
# this never bit until a trimmed chain arrived numbered 35-174. Translate here, from the parsed
# residue order, so offsets AND gaps are both handled.
# The summary CSV reports residue numbers, so keep them before converting.
fixed_resnums = {d: dict(chains) for d, chains in fixed_dict.items()}

for design_id in design_ids:
    pdb_path = design_files[design_id]
    for ch, positions in fixed_dict[design_id].items():
        if not positions:
            continue
        rank = {resnum: i + 1 for i, resnum in enumerate(get_residues_in_chain(pdb_path, ch))}
        converted = [rank[p] for p in positions if p in rank]
        dropped = [p for p in positions if p not in rank]
        if dropped:
            print(f"Warning: {design_id} chain {ch}: {len(dropped)} fixed position(s) "
                  f"absent from the structure, ignored: {dropped[:10]}")
        fixed_dict[design_id][ch] = converted

with open(fixed_jsonl_file,"w") as jsonl_file:
    #Python converts dictionaries to string having keys inside '', json only recognises ""
    jsonl_file.write(str(fixed_dict).replace("\'","\""))

def summarize(per_chain):
    """One selection string over the chains actually selected.

    Chain-aware ("A10-20+B5") as soon as more than one chain is in play, which is what an
    inter-tool selection column has to be; a single-chain step keeps the bare form it has
    always written, so existing consumers are unaffected.
    """
    if len(selected_chains) == 1:
        return list_to_sele(per_chain.get(selected_chains[0], []))
    return chain_aware_sele([(c, r) for c in selected_chains for r in per_chain.get(c, [])])


with open(sele_csv_file,"w") as csv_file:
    csv_file.write("id,fixed,mobile")
    for id in fixed_dict.keys():
        fixed = summarize(fixed_resnums[id])
        mobile = summarize(mobile_dict[id])
        csv_file.write("\n")
        csv_file.write(f"{id},{fixed},{mobile}")
        
print(f"Fixed positions written to: {fixed_jsonl_file}")
print(f"Selections summary written to: {sele_csv_file}")