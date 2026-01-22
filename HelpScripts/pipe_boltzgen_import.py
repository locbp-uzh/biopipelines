#!/usr/bin/env python3
"""
BoltzGenImport helper script.

Converts external structures into BoltzGen filesystem format:
1. Converts PDB files to CIF format
2. Reassigns chains (protein to A, ligand to B)
3. Optionally assigns sequences and zeros sidechain coordinates
4. Generates design_spec.yaml from ligand info and binder spec

This enables use of downstream BoltzGen steps (folding, analysis, filtering)
on structures from external tools like RFdiffusion.
"""

import os
import sys
import re
import argparse
import json
import pandas as pd
from pathlib import Path


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


def convert_and_reassign_chains(
    input_path: str,
    output_path: str,
    protein_chain: str = "A",
    ligand_chain: str = "B",
    new_sequence: str = None,
    ligand_name: str = "LIG"
) -> bool:
    """
    Convert structure to CIF and reassign chains.

    Input structures (e.g., from RFdiffusion) may have protein+ligand in chain A.
    This function separates them: protein -> protein_chain, ligand -> ligand_chain.

    Args:
        input_path: Input structure file (PDB or CIF)
        output_path: Output CIF file path
        protein_chain: Target chain ID for protein (default: "A")
        ligand_chain: Target chain ID for ligand (default: "B")
        new_sequence: Optional new sequence to assign (zeros sidechains)
        ligand_name: Residue name for ligand in output (default: "LIG").
                    BoltzGen expects ligand to be named "LIG".

    Returns:
        True if successful, False otherwise
    """
    try:
        from Bio.PDB import PDBParser, MMCIFParser, MMCIFIO, Structure, Model, Chain
        from Bio.PDB.Polypeptide import is_aa
    except ImportError:
        print("Error: BioPython is required. Install with: pip install biopython")
        return False

    try:
        # Parse input structure
        if input_path.endswith('.cif'):
            parser = MMCIFParser(QUIET=True)
        else:
            parser = PDBParser(QUIET=True)

        structure = parser.get_structure('structure', input_path)

        # Create new structure with reassigned chains
        new_structure = Structure.Structure('imported')
        new_model = Model.Model(0)
        new_structure.add(new_model)

        protein_chain_obj = Chain.Chain(protein_chain)
        ligand_chain_obj = Chain.Chain(ligand_chain)

        protein_residues = []
        ligand_residues = []

        # Separate protein and ligand residues from all chains
        for model in structure:
            for chain in model:
                for residue in chain:
                    resname = residue.get_resname()
                    # Check if it's a standard amino acid
                    if is_aa(residue, standard=True):
                        protein_residues.append(residue.copy())
                    elif resname not in ['HOH', 'WAT']:  # Skip water
                        # It's a ligand/hetero atom
                        ligand_residues.append(residue.copy())
            break  # Only process first model

        # Add protein residues to protein chain with sequential numbering
        for i, residue in enumerate(protein_residues):
            # Create new residue ID with sequential numbering
            new_id = (' ', i + 1, ' ')
            residue.id = new_id

            # If new sequence provided, mutate and zero sidechains
            if new_sequence and i < len(new_sequence):
                new_aa = new_sequence[i]
                new_resname = AA_1TO3.get(new_aa, 'ALA')
                residue.resname = new_resname

                # Zero sidechain coordinates
                for atom in residue:
                    if atom.name not in BACKBONE_ATOMS:
                        atom.set_coord((0.0, 0.0, 0.0))

            protein_chain_obj.add(residue)

        # Add ligand residues to ligand chain
        for i, residue in enumerate(ligand_residues):
            # Keep original residue ID or renumber
            new_id = residue.id
            if new_id[1] in [r.id[1] for r in ligand_chain_obj]:
                # Renumber to avoid conflicts
                new_id = (' ', i + 1, ' ')
            residue.id = new_id
            # Rename ligand residue to the specified name (default: LIG)
            # BoltzGen expects all ligands to be named "LIG"
            residue.resname = ligand_name
            ligand_chain_obj.add(residue)

        # Add chains to model
        if len(protein_chain_obj) > 0:
            new_model.add(protein_chain_obj)
        if len(ligand_chain_obj) > 0:
            new_model.add(ligand_chain_obj)

        # Write output CIF
        io = MMCIFIO()
        io.set_structure(new_structure)
        io.save(output_path)
        return True

    except Exception as e:
        print(f"Error processing {input_path}: {e}")
        import traceback
        traceback.print_exc()
        return False


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
        success = convert_and_reassign_chains(
            struct_path, output_cif,
            protein_chain=protein_chain,
            ligand_chain=ligand_chain,
            new_sequence=sequence,
            ligand_name=ligand_name
        )

        if success:
            stats['success'] += 1
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
