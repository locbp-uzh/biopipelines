#!/usr/bin/env python3
"""
BoltzGenImport helper script.

Converts external structures into BoltzGen filesystem format:
1. Reads PDB files manually
2. Reassigns chains (protein to A, ligand to B)
3. Optionally assigns sequences and zeros sidechain coordinates
4. Generates design_spec.yaml from ligand info and binder spec
5. Generates NPZ metadata files required by BoltzGen dataloader
6. Writes output as PDB files (not CIF)

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
from collections import defaultdict

# Backbone atom names (coordinates preserved)
BACKBONE_ATOMS = {'N', 'CA', 'C', 'O'}

# Standard amino acid 3-letter to 1-letter mapping
AA_3TO1 = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
    'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
    'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
    'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
}

AA_1TO3 = {v: k for k, v in AA_3TO1.items()}

# Standard amino acid atom templates (backbone + sidechain atoms in order)
# Backbone atoms will get real coordinates, sidechain atoms get (0,0,0)
AA_ATOMS = {
    'ALA': ['N', 'CA', 'C', 'O', 'CB'],
    'ARG': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'NE', 'CZ', 'NH1', 'NH2'],
    'ASN': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'OD1', 'ND2'],
    'ASP': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'OD1', 'OD2'],
    'CYS': ['N', 'CA', 'C', 'O', 'CB', 'SG'],
    'GLN': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'OE1', 'NE2'],
    'GLU': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'OE1', 'OE2'],
    'GLY': ['N', 'CA', 'C', 'O'],
    'HIS': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'ND1', 'CD2', 'CE1', 'NE2'],
    'ILE': ['N', 'CA', 'C', 'O', 'CB', 'CG1', 'CG2', 'CD1'],
    'LEU': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2'],
    'LYS': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'CE', 'NZ'],
    'MET': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'SD', 'CE'],
    'PHE': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
    'PRO': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD'],
    'SER': ['N', 'CA', 'C', 'O', 'CB', 'OG'],
    'THR': ['N', 'CA', 'C', 'O', 'CB', 'OG1', 'CG2'],
    'TRP': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2', 'NE1', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'],
    'TYR': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'OH'],
    'VAL': ['N', 'CA', 'C', 'O', 'CB', 'CG1', 'CG2'],
}

# Element mapping for common atom names
ATOM_ELEMENTS = {
    'N': 'N', 'CA': 'C', 'C': 'C', 'O': 'O',
    'CB': 'C', 'CG': 'C', 'CG1': 'C', 'CG2': 'C',
    'CD': 'C', 'CD1': 'C', 'CD2': 'C',
    'CE': 'C', 'CE1': 'C', 'CE2': 'C', 'CE3': 'C',
    'CZ': 'C', 'CZ2': 'C', 'CZ3': 'C', 'CH2': 'C',
    'ND1': 'N', 'ND2': 'N', 'NE': 'N', 'NE1': 'N', 'NE2': 'N', 'NZ': 'N', 'NH1': 'N', 'NH2': 'N',
    'OD1': 'O', 'OD2': 'O', 'OE1': 'O', 'OE2': 'O', 'OG': 'O', 'OG1': 'O', 'OH': 'O',
    'SD': 'S', 'SG': 'S',
}


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


class PDBAtom:
    """Represents a PDB ATOM/HETATM record."""
    def __init__(self, line):
        self.record = line[0:6].strip()
        self.serial = int(line[6:11].strip())
        self.name = line[12:16].strip()
        self.altloc = line[16:17].strip()
        self.resname = line[17:20].strip()
        self.chain = line[21:22].strip()
        self.resseq = int(line[22:26].strip())
        self.icode = line[26:27].strip()
        self.x = float(line[30:38].strip())
        self.y = float(line[38:46].strip())
        self.z = float(line[46:54].strip())
        self.occupancy = float(line[54:60].strip()) if line[54:60].strip() else 1.0
        self.tempfactor = float(line[60:66].strip()) if line[60:66].strip() else 0.0
        self.element = line[76:78].strip() if len(line) > 76 else self.name[0]

    def to_pdb_line(self, serial, chain, resname, resseq):
        """Convert atom back to PDB format line."""
        return (
            f"{self.record:<6}{serial:>5} "
            f"{self.name:<4}{self.altloc:1}{resname:>3} "
            f"{chain:1}{resseq:>4}{self.icode:1}   "
            f"{self.x:>8.3f}{self.y:>8.3f}{self.z:>8.3f}"
            f"{self.occupancy:>6.2f}{self.tempfactor:>6.2f}          "
            f"{self.element:>2}\n"
        )


def read_pdb(pdb_path: str):
    """
    Read PDB file manually and organize by chain and residue.

    Returns:
        dict: {chain_id: {resseq: [atoms]}}
    """
    structure = defaultdict(lambda: defaultdict(list))

    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith(('ATOM  ', 'HETATM')):
                atom = PDBAtom(line)
                structure[atom.chain][atom.resseq].append(atom)

    return structure


def write_pdb(output_path: str, protein_atoms, ligand_atoms, protein_chain='A', ligand_chain='B'):
    """
    Write PDB file with protein and ligand chains.

    Args:
        output_path: Output PDB file path
        protein_atoms: List of (resname, resseq, [atoms]) for protein
        ligand_atoms: List of atoms for ligand
        protein_chain: Chain ID for protein (default: A)
        ligand_chain: Chain ID for ligand (default: B)
    """
    with open(output_path, 'w') as f:
        serial = 1

        # Write protein chain
        for resname, resseq, atoms in protein_atoms:
            for atom in atoms:
                f.write(atom.to_pdb_line(serial, protein_chain, resname, resseq))
                serial += 1

        # Write TER card after protein
        if protein_atoms:
            f.write(f"TER   {serial:>5}      {resname:>3} {protein_chain}{resseq:>4}\n")
            serial += 1

        # Write ligand chain
        if ligand_atoms:
            for atom in ligand_atoms:
                f.write(atom.to_pdb_line(serial, ligand_chain, 'LIG', 1))
                serial += 1

            # Write TER card after ligand
            f.write(f"TER   {serial:>5}      LIG {ligand_chain}   1\n")

        # Write END card
        f.write("END\n")


def convert_and_reassign_chains(
    input_path: str,
    output_path: str,
    protein_chain: str = "A",
    ligand_chain: str = "B",
    new_sequence: str = None,
    ligand_name: str = "LIG"
) -> tuple:
    """
    Convert PDB structure: reassign chains, optionally apply sequence.

    Args:
        input_path: Input PDB file path
        output_path: Output PDB file path
        protein_chain: Target chain ID for protein
        ligand_chain: Target chain ID for ligand
        new_sequence: Optional new sequence to apply (inverse folding mode)
        ligand_name: Residue name for ligand (default: LIG)

    Returns:
        (success, n_protein_residues, n_ligand_atoms)
    """
    print(f"→ Processing: {input_path}")

    # Read structure
    structure = read_pdb(input_path)

    if not structure:
        print("  ERROR: No atoms found in structure")
        return (False, 0, 0)

    print(f"  Found {len(structure)} chains in input")

    # Organize residues by chain
    all_chains = []
    for chain_id, residues in structure.items():
        chain_residues = []
        for resseq in sorted(residues.keys()):
            atoms = residues[resseq]
            resname = atoms[0].resname
            chain_residues.append((resname, resseq, atoms))
        all_chains.append((chain_id, chain_residues))
        print(f"    chain {len(all_chains)}: id={chain_id:>2}, {len(chain_residues):3d} residues, first res={chain_residues[0][0] if chain_residues else '—'}")

    # Sort by number of residues (descending) - longest is protein
    all_chains.sort(key=lambda x: -len(x[1]))

    longest_chain_id, longest_chain_residues = all_chains[0]
    print(f"  → Longest chain selected as protein: chain {longest_chain_id}, {len(longest_chain_residues)} residues")

    # Separate protein and ligand residues
    protein_residues = []
    ligand_residues = []

    for resname, resseq, atoms in longest_chain_residues:
        if is_amino_acid(resname.strip()):
            protein_residues.append((resname, resseq, atoms))
        else:
            ligand_residues.extend(atoms)

    print(f"  Detected {len(protein_residues)} standard AA residues → protein")
    print(f"  Detected {len(longest_chain_residues) - len(protein_residues)} non-standard residues → ligand")

    # Build output protein chain
    output_protein = []
    n_protein_res = 0

    for i, (resname, resseq, atoms) in enumerate(protein_residues, 1):
        # Determine output residue name
        if new_sequence and i - 1 < len(new_sequence):
            new_aa = new_sequence[i - 1].upper()
            out_resname = AA_1TO3.get(new_aa, resname)
        else:
            out_resname = resname

        # Build atom map from input atoms
        atom_map = {atom.name: atom for atom in atoms}

        # Process atoms based on the NEW residue type template
        out_atoms = []

        if new_sequence and i - 1 < len(new_sequence):
            # Use template for the NEW amino acid type
            atom_template = AA_ATOMS.get(out_resname, ['N', 'CA', 'C', 'O'])

            for atom_name in atom_template:
                # Create new atom
                out_atom = PDBAtom.__new__(PDBAtom)

                # Set basic properties
                out_atom.record = 'ATOM'
                out_atom.serial = 1  # Will be set during writing
                out_atom.name = atom_name
                out_atom.altloc = ''
                out_atom.resname = out_resname
                out_atom.chain = ''  # Will be set during writing
                out_atom.resseq = i
                out_atom.icode = ''
                out_atom.occupancy = 1.0
                out_atom.tempfactor = 0.0
                out_atom.element = ATOM_ELEMENTS.get(atom_name, atom_name[0])

                # Set coordinates: backbone from input, sidechain at (0,0,0)
                if atom_name in BACKBONE_ATOMS and atom_name in atom_map:
                    # Copy backbone coordinates from input
                    input_atom = atom_map[atom_name]
                    out_atom.x = input_atom.x
                    out_atom.y = input_atom.y
                    out_atom.z = input_atom.z
                else:
                    # Sidechain or missing backbone: set to (0,0,0)
                    out_atom.x = 0.0
                    out_atom.y = 0.0
                    out_atom.z = 0.0

                out_atoms.append(out_atom)
        else:
            # No mutation: copy all atoms as-is
            for atom in atoms:
                out_atom = PDBAtom.__new__(PDBAtom)
                out_atom.__dict__.update(atom.__dict__)
                out_atoms.append(out_atom)

        output_protein.append((out_resname, i, out_atoms))
        n_protein_res += 1

    # Build output ligand chain
    output_ligand = []
    n_ligand_atoms = 0
    n_ligand_atoms_with_h = 0
    ligand_atom_elements = []

    for atom in ligand_residues:
        # Clone atom
        out_atom = PDBAtom.__new__(PDBAtom)
        out_atom.__dict__.update(atom.__dict__)

        ligand_atom_elements.append(out_atom.element)
        n_ligand_atoms_with_h += 1

        # BoltzGen tokenizer excludes hydrogen atoms
        if out_atom.element != 'H':
            output_ligand.append(out_atom)
            n_ligand_atoms += 1

    print(f"  Ligand atom elements: {ligand_atom_elements}")
    from collections import Counter
    element_counts = Counter(ligand_atom_elements)
    print(f"  Ligand element counts: {dict(element_counts)}")
    print(f"  Ligand atoms (excluding H): {n_ligand_atoms} (with H: {n_ligand_atoms_with_h})")

    print(f"  Final stats → protein: {n_protein_res} res, ligand: {n_ligand_atoms} atoms")

    if n_protein_res == 0 and n_ligand_atoms == 0:
        print("  !!! NOTHING TO WRITE !!!")
        return False, 0, 0

    # Write output PDB
    write_pdb(output_path, output_protein, output_ligand, protein_chain, ligand_chain)

    # Verify output
    print(f"  Before writing:")
    print(f"    Protein chain {protein_chain}: {n_protein_res} residues, {sum(len(atoms) for _, _, atoms in output_protein)} atoms")
    print(f"    Ligand chain {ligand_chain}:  1 residues, {n_ligand_atoms} atoms")

    # Read back and verify
    verify_structure = read_pdb(output_path)
    total_residues = sum(len(residues) for residues in verify_structure.values())
    total_atoms = sum(len(atoms) for residues in verify_structure.values() for atoms in residues.values())

    print(f"  [PDB VERIFY] Written PDB file: {os.path.basename(output_path)}")
    print(f"    Total chains: {len(verify_structure)}")
    print(f"    Total residues: {total_residues}")
    print(f"    Total atoms: {total_atoms}")

    for chain_id, residues in sorted(verify_structure.items()):
        n_res = len(residues)
        n_atoms = sum(len(atoms) for atoms in residues.values())
        first_res = list(residues.values())[0][0].resname if residues else "—"
        print(f"      Chain {chain_id}: {n_res} residues, {n_atoms} atoms, first={first_res}")

        # For ligand chain, show detailed atom info
        if chain_id == ligand_chain and n_res > 0:
            lig_atoms = list(residues.values())[0]
            lig_elements = [atom.element for atom in lig_atoms]
            lig_element_counts = Counter(lig_elements)
            print(f"        Ligand atoms in PDB: {lig_elements}")
            print(f"        Ligand element counts in PDB: {dict(lig_element_counts)}")
            non_h_count = sum(1 for atom in lig_atoms if atom.element != 'H')
            print(f"        Non-hydrogen atoms: {non_h_count}")

    return True, n_protein_res, n_ligand_atoms


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

    BoltzGen expects each PDB structure file to have a corresponding NPZ file
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

        print(f"  [NPZ DEBUG] Creating NPZ metadata:")
        print(f"    n_protein_residues: {n_protein_residues}")
        print(f"    n_ligand_atoms: {n_ligand_atoms}")
        print(f"    n_tokens (total): {n_tokens}")
        print(f"    Output: {os.path.basename(output_path)}")

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

        print(f"    design_mask shape: {design_mask.shape}")
        print(f"    design_mask nonzero (protein tokens): {np.sum(design_mask)}")

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

        # BoltzGen naming: design_spec_<N>.pdb
        output_name = f"design_spec_{stats['processed'] - 1}.pdb"
        output_pdb = os.path.join(out_subdir, output_name)

        # Find sequence for this structure (if in inverse_folding mode)
        sequence = None
        if mode == "inverse_folding" and sequences:
            sequence = find_sequence_for_structure(struct_id, sequences, id_map)
            if sequence is None:
                stats['missing_sequences'] += 1
                print(f"  Warning: No sequence found for {struct_id}")

        # Convert and reassign chains
        success, n_protein, n_ligand = convert_and_reassign_chains(
            struct_path, output_pdb,
            protein_chain=protein_chain,
            ligand_chain=ligand_chain,
            new_sequence=sequence,
            ligand_name=ligand_name
        )

        if success:
            # Generate NPZ metadata file alongside PDB
            # NOTE: Exclude ligand from NPZ - BoltzGen doesn't tokenize ligand chain
            output_npz = output_pdb.replace('.pdb', '.npz')
            npz_success = generate_npz_metadata(output_npz, n_protein, 0)  # ligand atoms = 0
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

def generate_design_spec_only(args):
    """Generate only the design_spec.yaml file (for molecule pickle generation)."""
    print("=" * 60)
    print("Generating design_spec.yaml")
    print("=" * 60)

    # Load ligand info
    ligand_info = {}
    if args.ligand_csv and os.path.exists(args.ligand_csv):
        ligand_info = load_ligand_info(args.ligand_csv)
        print(f"Loaded ligand info: {ligand_info}")

    # Generate the design_spec.yaml
    output_path = generate_design_spec(
        output_dir=os.path.dirname(args.output),
        binder_min=args.binder_min,
        binder_max=args.binder_max,
        protein_chain=args.protein_chain,
        ligand_chain=args.ligand_chain,
        ligand_info=ligand_info
    )

    # If output path differs from default, rename
    if output_path != args.output:
        import shutil
        shutil.move(output_path, args.output)

    print(f"Generated: {args.output}")
    return 0


def main():
    parser = argparse.ArgumentParser(
        description='Import structures into BoltzGen filesystem format'
    )

    # Mode flag for design_spec generation only
    parser.add_argument(
        '--generate-design-spec',
        action='store_true',
        help='Only generate design_spec.yaml (for molecule pickle generation)'
    )

    parser.add_argument(
        '--structures',
        help='File containing list of structure paths (one per line)'
    )
    parser.add_argument(
        '--output',
        required=True,
        help='Output directory for BoltzGen format (or output file for --generate-design-spec)'
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

    # Handle design_spec generation mode
    if args.generate_design_spec:
        return generate_design_spec_only(args)

    # Regular import mode requires --structures
    if not args.structures:
        parser.error("--structures is required when not using --generate-design-spec")

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
