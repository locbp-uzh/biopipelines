#!/usr/bin/env python3
"""
BoltzGenImport helper script.

Converts external structures into BoltzGen filesystem format:
1. Converts PDB files to CIF format
2. Reassigns chains (protein to A, ligand to B)
3. Optionally assigns sequences and zeros sidechain coordinates
4. Generates design_spec.yaml from ligand info and binder spec
5. Generates NPZ metadata files required by BoltzGen dataloader

This enables use of downstream BoltzGen steps (folding, analysis, filtering)
on structures from external tools like RFdiffusion.
"""

import os
import sys
import re
import argparse
import json
import numpy as np
import pandas as pd
from pathlib import Path
import gemmi

# Backbone atom names (coordinates preserved)
BACKBONE_ATOMS = {'N', 'CA', 'C', 'O', 'CB'}

# Standard amino acid 3-letter to 1-letter mapping
AA_3TO1 = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
    'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
    'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
    'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
}

AA_1TO3 = {v: k for k, v in AA_3TO1.items()}


def extract_structure_id(filepath: str) -> str:
    """Extract structure ID from filepath."""
    name = os.path.basename(filepath)
    for ext in ['.pdb', '.cif', '.ent']:
        if name.endswith(ext):
            name = name[:-len(ext)]
            break
    return name


def is_amino_acid(resname: str) -> bool:
    """Check if residue is a standard amino acid."""
    return resname in AA_3TO1

def write_cif(structure, output_path):
    st = gemmi.Structure()
    st.name = "imported"
    model = gemmi.Model("1")
    st.add_model(model)
    for chain in structure[0]:
        gchain = gemmi.Chain(chain.id)
        for i, res in enumerate(chain, start=1):
            gres = gemmi.Residue()
            gres.name = res.resname
            gres.seqid = gemmi.SeqId(str(i)) #no chain break, start at 1
            for atom in res:
                a = gemmi.Atom()
                a.name = atom.name
                a.pos = gemmi.Position(*atom.coord)
                a.element = gemmi.Element(atom.element)
                gres.add_atom(a)
            gchain.add_residue(gres)
        model.add_chain(gchain)

    st.setup_entities()  # ← THIS creates _entity_poly_seq
    doc = st.make_mmcif_document()          # ← this creates a complete mmCIF document
    doc.write_file(output_path)

def convert_and_reassign_chains(
    input_path: str,
    output_path: str,
    protein_chain: str = "A",
    ligand_chain: str = "B",
    new_sequence: str = None,
    ligand_name: str = "LIG"
) -> tuple:

    print(f"→ Processing: {input_path}")

    st_in = gemmi.read_structure(input_path)
    if not st_in:
        print("  ERROR: gemmi.read_structure() returned None/empty")
        return (False, 0, 0)

    if not st_in[0]:
        print("  ERROR: No model 1 in structure")
        return (False, 0, 0)

    chains = list(st_in[0])
    print(f"  Found {len(chains)} chains in input")

    if not chains:
        print("  → No chains at all → empty output")
        return (False, 0, 0)

    def get_chain_id(ch):
        return ch.name if hasattr(ch, 'name') else ch.label_asym_id if hasattr(ch, 'label_asym_id') else '?'
    # ── very important ──
    for i, ch in enumerate(chains, 1):
        res_count = len(ch)
        first_res = ch[0].name if ch else "—"
        chain_id = get_chain_id(ch)
        print(f"    chain {i}: id={chain_id:>2}, {res_count:3d} residues, first res={first_res}")
    
    # Sort by number of residues (descending)
    chains.sort(key=lambda ch: -len(ch))

    
    longest_chain_id = get_chain_id(chains[0])
        
    print(f"  → Longest chain selected as protein: chain {longest_chain_id}, "
          f"{len(chains[0])} residues")

    g_protein = gemmi.Chain(protein_chain)
    g_ligand  = gemmi.Chain(ligand_chain)

    n_protein_res = 0
    n_ligand_atoms = 0

    # =================================================================
    protein_input_chain = chains[0]
    ligand_residues = []

    # Heuristic: last residues that are NOT standard amino acids → ligand
    protein_residues = []
    for res in protein_input_chain:
        if is_amino_acid(res.name.strip()):          # ← use your existing is_amino_acid()
            protein_residues.append(res)
        else:
            ligand_residues.append(res)

    print(f"  Detected {len(protein_residues)} standard AA residues → protein")
    print(f"  Detected {len(ligand_residues)} non-standard residues → ligand")

    # Now build protein from protein_residues only
    n_protein_res = 0
    for i, res in enumerate(protein_residues, 1):
        gres = gemmi.Residue()
        gres.name = res.name
        gres.seqid = gemmi.SeqId(str(i))

        if new_sequence and i-1 < len(new_sequence):
            new_aa = new_sequence[i-1].upper()
            new_3letter = AA_1TO3.get(new_aa, '???')
            gres.name = new_3letter
            # only backbone
            for atom in res:
                a = atom.clone()               # copy name, element, occupancy, etc.
                if not atom.name in BACKBONE_ATOMS:
                    a.pos = gemmi.Position(0.0, 0.0, 0.0)
                gres.add_atom(atom.clone())
        else:
            for atom in res:
                gres.add_atom(atom.clone())

        g_protein.add_residue(gres)
        n_protein_res += 1

    # Ligand part
    n_ligand_atoms = 0
    if ligand_residues:
        gres_lig = gemmi.Residue()           # many tools put ligand as SINGLE residue
        gres_lig.name = "LIG"                # "LIG"
        gres_lig.seqid = gemmi.SeqId("1")

        for orig_res in ligand_residues:
            for atom in orig_res:
                gres_lig.add_atom(atom.clone())
                n_ligand_atoms += 1

        g_ligand.add_residue(gres_lig)

    # Final check before writing
    total_res = n_protein_res + len(g_ligand)
    total_atoms = sum(len(res) for res in g_protein) + n_ligand_atoms

    print(f"  Final stats → protein: {n_protein_res} res, ligand: {n_ligand_atoms} atoms")

    if n_protein_res == 0 and n_ligand_atoms == 0:
        print("  !!! NOTHING WAS COPIED → WILL WRITE EMPTY CIF !!!")

    # ── At this point you should already have ──
    #    g_protein  (Chain "A" with 170 residues)
    #    g_ligand   (Chain "B" with 1 residue, 25 atoms)

        # ── Debug: show what we have ─────────────────────────────────────────────
    print(f"  Before writing:")
    print(f"    Protein chain A: {len(g_protein)} residues, "
          f"{sum(1 for res in g_protein for _ in res)} atoms")
    print(f"    Ligand chain B:  {len(g_ligand)} residues, "
          f"{sum(1 for res in g_ligand for _ in res)} atoms")

    if len(g_protein) == 0 and len(g_ligand) == 0:
        print("  !!! No residues to write → will produce empty file")
        return False, 0, 0

    # After building g_protein and g_ligand as before...

    # ── New writing block ────────────────────────────────────────────────────────

    st_out = gemmi.Structure()
    st_out.name = "imported"
    st_out.spacegroup_hm = "P 1"           # minimal, avoids some complaints

    model = gemmi.Model("1")
    
    model.add_chain(g_protein)
    model.add_chain(g_ligand)

    st_out.add_model(model)
    print(f"    Entities: {len(st_out.entities)}")

    # Now setup entities – crucial order
    st_out.setup_entities()
    # Optional: ensure minimal required categories
    st_out.make_mmcif_headers()   # adds _entry, _cell etc if missing

    # ── Insert this ───────────────────────────────────────
    import tempfile
    import os

    with tempfile.NamedTemporaryFile(suffix='.pdb', delete=False) as tmp:
        tmp_path = tmp.name
        gemmi.write_pdb(st_out, tmp_path)          # or st_out.make_pdb_string() → write yourself

    # Read back → Gemmi now "knows" the flat atom list
    st_round = gemmi.read_structure(tmp_path)
    os.unlink(tmp_path)  # clean up

    doc = st_round.make_mmcif_document()

    block = doc.sole_block()
    loop = block.find_loop("_atom_site.")
    num_atoms = len(loop)/len(loop.tags) if loop else 0
    print(f"→ now has {num_atoms} atoms")
    doc.write_file(output_path)


    return num_atoms>0, len(g_protein), sum(len(r) for r in g_ligand)


def load_sequences(sequences_csv: str) -> dict:
    """
    Load sequences from CSV file.

    Args:
        sequences_csv: Path to CSV with 'id' and 'sequence' columns

    Returns:
        Dictionary mapping structure ID to sequence
    """
    df = pd.read_csv(sequences_csv)

    if 'id' not in df.columns:
        raise ValueError(f"sequences CSV must have 'id' column, found: {list(df.columns)}")

    if 'sequence' not in df.columns:
        raise ValueError(f"sequences CSV must have 'sequence' column, found: {list(df.columns)}")

    sequences = {}
    for _, row in df.iterrows():
        seq_id = str(row['id'])
        sequence = str(row['sequence'])
        sequences[seq_id] = sequence

    return sequences


def apply_id_map(structure_id: str, id_map: dict) -> list:
    """
    Apply id_map pattern to generate possible sequence ID patterns.

    The id_map pattern {"*": "*_<N>"} means structure ID can have
    recursive numeric suffixes added to match sequence IDs.

    Args:
        structure_id: ID of the structure
        id_map: Mapping pattern, e.g., {"*": "*_<N>"}

    Returns:
        List of matching patterns (regex patterns or exact strings)
    """
    patterns = [structure_id]  # Always include exact match

    # Handle the default recursive pattern {"*": "*_<N>"}
    if "*" in id_map:
        pattern = id_map["*"]
        if pattern == "*_<N>":
            # Recursive numeric suffix: matches structure_id_1, structure_id_1_1, etc.
            # Build regex pattern: structure_id followed by one or more _<number> suffixes
            escaped_id = re.escape(structure_id)
            patterns.append(f"^{escaped_id}(_\\d+)+$")
        elif pattern == "*":
            # Exact match only (no transformation)
            pass
        else:
            # Custom pattern - replace * with structure_id
            custom_pattern = pattern.replace("*", structure_id)
            patterns.append(custom_pattern)

    return patterns


def find_sequence_for_structure(structure_id: str, sequences: dict, id_map: dict = None) -> str:
    """
    Find the sequence for a given structure ID using id_map pattern.

    Handles various naming conventions:
    - Exact match: structure_id in sequences
    - id_map pattern match: using recursive suffix matching
    - Fallback prefix match: sequence_id starts with structure_id

    Args:
        structure_id: ID of the structure
        sequences: Dictionary mapping IDs to sequences
        id_map: ID mapping pattern, e.g., {"*": "*_<N>"} for recursive suffixes.
               Default: {"*": "*_<N>"}

    Returns:
        Sequence string or None if not found
    """
    if id_map is None:
        id_map = {"*": "*_<N>"}

    # Get matching patterns from id_map
    patterns = apply_id_map(structure_id, id_map)

    # First: exact match (always first priority)
    if structure_id in sequences:
        return sequences[structure_id]

    # Second: try regex patterns from id_map
    for pattern in patterns:
        if pattern == structure_id:
            continue  # Skip exact match, already checked
        try:
            regex = re.compile(pattern)
            for seq_id, sequence in sequences.items():
                if regex.match(seq_id):
                    return sequence
        except re.error:
            # Not a regex, try exact match
            if pattern in sequences:
                return sequences[pattern]

    # Third: fallback prefix match (e.g., "design_1" matches "design_1_seq_0")
    for seq_id, sequence in sequences.items():
        if seq_id.startswith(structure_id + "_"):
            return sequence

    # Fourth: suffix match (structure might have extra prefix)
    for seq_id, sequence in sequences.items():
        if structure_id.endswith("_" + seq_id) or structure_id == seq_id:
            return sequence

    return None


def load_ligand_info(ligand_csv: str) -> dict:
    """
    Load ligand information from CSV.

    Args:
        ligand_csv: Path to compounds CSV with id, smiles/ccd columns

    Returns:
        Dictionary with ligand info (smiles or ccd)
    """
    df = pd.read_csv(ligand_csv)

    ligand_info = {}

    if 'smiles' in df.columns and pd.notna(df['smiles'].iloc[0]):
        ligand_info['smiles'] = str(df['smiles'].iloc[0])
    elif 'ccd' in df.columns and pd.notna(df['ccd'].iloc[0]):
        ligand_info['ccd'] = str(df['ccd'].iloc[0])

    return ligand_info


def generate_design_spec(
    output_dir: str,
    binder_min: int,
    binder_max: int,
    protein_chain: str,
    ligand_chain: str,
    ligand_info: dict = None
) -> str:
    """
    Generate design_spec.yaml file.

    Args:
        output_dir: Output directory
        binder_min: Minimum binder length
        binder_max: Maximum binder length
        protein_chain: Protein chain ID
        ligand_chain: Ligand chain ID
        ligand_info: Ligand information (smiles or ccd)

    Returns:
        Path to generated design_spec.yaml
    """
    design_spec_path = os.path.join(output_dir, "design_spec.yaml")

    # Build entities
    lines = ["entities:"]

    # Protein entity
    if binder_min == binder_max:
        seq_spec = str(binder_min)
    else:
        seq_spec = f"{binder_min}..{binder_max}"

    lines.append(f"- protein:")
    lines.append(f"    id: {protein_chain}")
    lines.append(f"    sequence: {seq_spec}")

    # Ligand entity (if provided)
    if ligand_info:
        lines.append(f"- ligand:")
        lines.append(f"    id: {ligand_chain}")
        if 'smiles' in ligand_info:
            lines.append(f"    smiles: {ligand_info['smiles']}")
        elif 'ccd' in ligand_info:
            lines.append(f"    ccd: {ligand_info['ccd']}")

    content = "\n".join(lines) + "\n"

    with open(design_spec_path, 'w') as f:
        f.write(content)

    return design_spec_path


def generate_npz_metadata(
    output_path: str,
    n_protein_residues: int,
    n_ligand_atoms: int
) -> bool:
    """
    Generate NPZ metadata file required by BoltzGen dataloader.

    BoltzGen expects each CIF structure file to have a corresponding NPZ file
    with per-token metadata arrays.

    Args:
        output_path: Output NPZ file path
        n_protein_residues: Number of protein residues (each is one token)
        n_ligand_atoms: Number of ligand atoms (each is one token)

    Returns:
        True if successful, False otherwise
    """
    try:
        n_tokens = n_protein_residues + n_ligand_atoms

        if n_tokens == 0:
            print(f"Warning: No tokens to write for {output_path}")
            return False

        # design_mask: 1.0 for protein (designed), 0.0 for ligand (fixed)
        design_mask = np.concatenate([
            np.ones(n_protein_residues, dtype=np.float32),
            np.zeros(n_ligand_atoms, dtype=np.float32)
        ])

        # inverse_fold_design_mask: same as design_mask
        inverse_fold_design_mask = design_mask.copy()

        # ss_type: all UNSPECIFIED (0)
        ss_type = np.zeros(n_tokens, dtype=np.int64)

        # binding_type: all UNSPECIFIED (0)
        binding_type = np.zeros(n_tokens, dtype=np.int64)

        # token_resolved_mask: all resolved (1.0)
        token_resolved_mask = np.ones(n_tokens, dtype=np.float32)

        np.savez(
            output_path,
            design_mask=design_mask,
            inverse_fold_design_mask=inverse_fold_design_mask,
            ss_type=ss_type,
            binding_type=binding_type,
            token_resolved_mask=token_resolved_mask
        )
        return True

    except Exception as e:
        print(f"Error generating NPZ metadata {output_path}: {e}")
        import traceback
        traceback.print_exc()
        return False


def import_structures(
    structures_file: str,
    output_dir: str,
    mode: str,
    protein_chain: str,
    ligand_chain: str,
    binder_min: int,
    binder_max: int,
    sequences_csv: str = None,
    ligand_csv: str = None,
    id_map: dict = None,
    ligand_name: str = "LIG"
) -> dict:
    """
    Import structures into BoltzGen format.

    Args:
        structures_file: File containing list of structure paths
        output_dir: Output directory
        mode: "design" or "inverse_folding"
        protein_chain: Target chain ID for protein
        ligand_chain: Target chain ID for ligand
        binder_min: Minimum binder length
        binder_max: Maximum binder length
        sequences_csv: Path to sequences CSV (for inverse_folding mode)
        ligand_csv: Path to ligand compounds CSV
        id_map: ID mapping pattern for matching structure IDs to sequence IDs.
               Default {"*": "*_<N>"} handles recursive numeric suffixes.
        ligand_name: Residue name for ligand in output (default: "LIG")

    Returns:
        Dictionary with import statistics
    """
    if id_map is None:
        id_map = {"*": "*_<N>"}
    stats = {
        'processed': 0,
        'success': 0,
        'failed': 0,
        'missing_sequences': 0
    }

    # Read structure paths
    with open(structures_file, 'r') as f:
        structures = [line.strip() for line in f if line.strip()]

    print(f"Found {len(structures)} structures to import")

    # Load sequences if provided
    sequences = {}
    if sequences_csv:
        sequences = load_sequences(sequences_csv)
        print(f"Loaded {len(sequences)} sequences from CSV")

    # Load ligand info if provided
    ligand_info = {}
    if ligand_csv and os.path.exists(ligand_csv):
        ligand_info = load_ligand_info(ligand_csv)
        print(f"Loaded ligand info: {ligand_info}")

    # Generate design_spec.yaml
    design_spec_path = generate_design_spec(
        output_dir, binder_min, binder_max,
        protein_chain, ligand_chain, ligand_info
    )
    print(f"Generated design_spec.yaml: {design_spec_path}")

    # Determine output directory
    if mode == "design":
        out_subdir = os.path.join(output_dir, "intermediate_designs")
    else:  # inverse_folding
        out_subdir = os.path.join(output_dir, "intermediate_designs_inverse_folded")

    os.makedirs(out_subdir, exist_ok=True)

    # Process each structure
    for struct_path in structures:
        stats['processed'] += 1
        struct_id = extract_structure_id(struct_path)

        # BoltzGen naming: design_spec_<N>.cif
        # Extract number from struct_id if possible, or use index
        output_name = f"design_spec_{stats['processed'] - 1}.cif"
        output_cif = os.path.join(out_subdir, output_name)

        # Find sequence for this structure (if in inverse_folding mode)
        sequence = None
        if mode == "inverse_folding" and sequences:
            sequence = find_sequence_for_structure(struct_id, sequences, id_map)
            if sequence is None:
                stats['missing_sequences'] += 1
                print(f"  Warning: No sequence found for {struct_id}")

        # Convert and reassign chains
        success, n_protein, n_ligand = convert_and_reassign_chains(
            struct_path, output_cif,
            protein_chain=protein_chain,
            ligand_chain=ligand_chain,
            new_sequence=sequence,
            ligand_name=ligand_name
        )

        if success:
            # Generate NPZ metadata file alongside CIF
            output_npz = output_cif.replace('.cif', '.npz')
            npz_success = generate_npz_metadata(output_npz, n_protein, n_ligand)
            if npz_success:
                stats['success'] += 1
            else:
                stats['failed'] += 1
                print(f"  Failed NPZ generation: {struct_id}")
        else:
            stats['failed'] += 1
            print(f"  Failed: {struct_id}")

        # Progress update every 100 structures
        if stats['processed'] % 100 == 0:
            print(f"  Processed {stats['processed']}/{len(structures)} structures...")

    return stats

def main():
    parser = argparse.ArgumentParser(
        description='Import structures into BoltzGen filesystem format'
    )
    parser.add_argument(
        '--structures',
        required=True,
        help='File containing list of structure paths (one per line)'
    )
    parser.add_argument(
        '--output',
        required=True,
        help='Output directory for BoltzGen format'
    )
    parser.add_argument(
        '--mode',
        choices=['design', 'inverse_folding'],
        default='design',
        help='Import mode: "design" for backbone only, '
             '"inverse_folding" for with sequences'
    )
    parser.add_argument(
        '--protein-chain',
        default='A',
        help='Chain ID for protein in output (default: A)'
    )
    parser.add_argument(
        '--ligand-chain',
        default='B',
        help='Chain ID for ligand in output (default: B)'
    )
    parser.add_argument(
        '--binder-min',
        type=int,
        default=100,
        help='Minimum binder length for design_spec'
    )
    parser.add_argument(
        '--binder-max',
        type=int,
        default=200,
        help='Maximum binder length for design_spec'
    )
    parser.add_argument(
        '--sequences',
        help='CSV file with id and sequence columns (for inverse_folding mode)'
    )
    parser.add_argument(
        '--ligand-csv',
        help='CSV file with ligand info (smiles or ccd columns)'
    )
    parser.add_argument(
        '--id-map',
        default='{"*": "*_<N>"}',
        help='JSON string for ID mapping pattern (default: {"*": "*_<N>"}). '
             'Maps structure IDs to sequence IDs with recursive numeric suffixes.'
    )
    parser.add_argument(
        '--ligand-name',
        default='LIG',
        help='Residue name for ligand in output (default: LIG). '
             'BoltzGen expects all ligands to be named "LIG".'
    )

    args = parser.parse_args()

    # Parse id_map from JSON string
    try:
        id_map = json.loads(args.id_map)
    except json.JSONDecodeError as e:
        print(f"Error parsing --id-map JSON: {e}")
        sys.exit(1)

    # Validate arguments
    if args.mode == 'inverse_folding' and not args.sequences:
        print("Warning: --sequences not provided for inverse_folding mode")

    print("=" * 60)
    print("BoltzGen Import")
    print("=" * 60)
    print(f"Mode: {args.mode}")
    print(f"Output: {args.output}")
    print(f"Chains: protein={args.protein_chain}, ligand={args.ligand_chain}")
    print(f"Binder spec: {args.binder_min}-{args.binder_max}")
    print(f"Ligand name: {args.ligand_name}")
    print(f"ID map: {id_map}")
    if args.sequences:
        print(f"Sequences: {args.sequences}")
    if args.ligand_csv:
        print(f"Ligand CSV: {args.ligand_csv}")

    stats = import_structures(
        structures_file=args.structures,
        output_dir=args.output,
        mode=args.mode,
        protein_chain=args.protein_chain,
        ligand_chain=args.ligand_chain,
        binder_min=args.binder_min,
        binder_max=args.binder_max,
        sequences_csv=args.sequences,
        ligand_csv=args.ligand_csv,
        id_map=id_map,
        ligand_name=args.ligand_name
    )

    print("\n" + "=" * 60)
    print("Import Summary")
    print("=" * 60)
    print(f"Structures processed: {stats['processed']}")
    print(f"Successfully imported: {stats['success']}")
    print(f"Failed: {stats['failed']}")
    if args.mode == 'inverse_folding':
        print(f"Missing sequences: {stats['missing_sequences']}")
    print("=" * 60)

    if stats['failed'] > 0:
        print("\nWarning: Some structures failed to import")
        sys.exit(1)

    return 0


if __name__ == "__main__":
    sys.exit(main())
