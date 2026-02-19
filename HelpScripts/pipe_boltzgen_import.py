#!/usr/bin/env python3
# Copyright (c) 2026 Gianluca Quargnali @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
BoltzGenImport helper script.

Converts external structures into BoltzGen filesystem format:
1. Reads PDB files manually
2. Reassigns chains (protein to A, ligand to B)
3. Optionally assigns sequences and zeros sidechain coordinates
4. Generates design_spec.yaml from ligand info and binder spec
5. Generates NPZ metadata files required by BoltzGen dataloader
6. Writes output as CIF files (not PDB) because:
   - CIF allows longer residue names (LIG1) unlike PDB's 3-char limit
   - BoltzGen's parse_mmcif accepts the mols dictionary for ligand resolution,
     while parse_pdb does not

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


def write_cif(output_path: str, protein_atoms, ligand_atoms, protein_chain='A', ligand_chain='B', ligand_name='LIG1'):
    """
    Write mmCIF file with protein and ligand chains in BoltzGen-compatible format.

    CIF format allows longer residue names (like LIG1) unlike PDB's 3-char limit.
    Also, BoltzGen's parse_mmcif accepts the mols dictionary for ligand resolution,
    while parse_pdb does not.

    Args:
        output_path: Output CIF file path
        protein_atoms: List of (resname, resseq, [atoms]) for protein
        ligand_atoms: List of atoms for ligand
        protein_chain: Chain ID for protein (default: A)
        ligand_chain: Chain ID for ligand (default: B)
        ligand_name: Residue name for ligand (default: LIG1)
    """
    # Build one-letter sequence from protein residues
    sequence = ""
    for resname, resseq, atoms in protein_atoms:
        one_letter = AA_3TO1.get(resname, 'X')
        sequence += one_letter

    # Collect unique atom type symbols
    atom_types = set()
    for resname, resseq, atoms in protein_atoms:
        for atom in atoms:
            element = atom.element if atom.element else atom.name[0]
            atom_types.add(element)
    if ligand_atoms:
        for atom in ligand_atoms:
            element = atom.element if atom.element else atom.name[0]
            atom_types.add(element)
    atom_types = sorted(atom_types)

    with open(output_path, 'w') as f:
        # Write header
        f.write("data_model\n")
        f.write("_entry.id model\n")
        f.write("\n")

        # Write cell parameters
        f.write("_cell.entry_id model\n")
        f.write("_cell.length_a 1\n")
        f.write("_cell.length_b 1\n")
        f.write("_cell.length_c 1\n")
        f.write("_cell.angle_alpha 90\n")
        f.write("_cell.angle_beta 90\n")
        f.write("_cell.angle_gamma 90\n")
        f.write("\n")

        # Write symmetry
        f.write("_symmetry.entry_id model\n")
        f.write("_symmetry.space_group_name_H-M ''\n")
        f.write("\n")

        # Write entity loop
        f.write("loop_\n")
        f.write("_entity.id\n")
        f.write("_entity.type\n")
        f.write("1 polymer\n")
        if ligand_atoms:
            f.write("2 non-polymer\n")
        f.write("\n")

        # Write entity_poly loop
        f.write("loop_\n")
        f.write("_entity_poly.entity_id\n")
        f.write("_entity_poly.type\n")
        f.write("_entity_poly.pdbx_strand_id\n")
        f.write("_entity_poly.pdbx_seq_one_letter_code\n")
        f.write(f"1 polypeptide(L) ? {sequence}\n")
        f.write("\n")
        f.write("\n")
        f.write("\n")

        # Write struct_asym loop
        f.write("loop_\n")
        f.write("_struct_asym.id\n")
        f.write("_struct_asym.entity_id\n")
        f.write(f"{protein_chain} 1\n")
        if ligand_atoms:
            f.write(f"{ligand_chain} 2\n")
        f.write("\n")
        f.write("\n")
        f.write("\n")
        f.write("\n")

        # Write atom_type loop
        f.write("loop_\n")
        f.write("_atom_type.symbol\n")
        for atype in atom_types:
            f.write(f"{atype}\n")
        f.write("\n")

        # Write entity_poly_seq loop
        f.write("loop_\n")
        f.write("_entity_poly_seq.entity_id\n")
        f.write("_entity_poly_seq.num\n")
        f.write("_entity_poly_seq.mon_id\n")
        for resname, resseq, atoms in protein_atoms:
            f.write(f"1 {resseq} {resname}\n")
        f.write("\n")

        # Write atom_site loop
        f.write("loop_\n")
        f.write("_atom_site.group_PDB\n")
        f.write("_atom_site.id\n")
        f.write("_atom_site.type_symbol\n")
        f.write("_atom_site.label_atom_id\n")
        f.write("_atom_site.label_alt_id\n")
        f.write("_atom_site.label_comp_id\n")
        f.write("_atom_site.label_asym_id\n")
        f.write("_atom_site.label_entity_id\n")
        f.write("_atom_site.label_seq_id\n")
        f.write("_atom_site.pdbx_PDB_ins_code\n")
        f.write("_atom_site.Cartn_x\n")
        f.write("_atom_site.Cartn_y\n")
        f.write("_atom_site.Cartn_z\n")
        f.write("_atom_site.occupancy\n")
        f.write("_atom_site.B_iso_or_equiv\n")
        f.write("_atom_site.pdbx_formal_charge\n")
        f.write("_atom_site.auth_seq_id\n")
        f.write("_atom_site.auth_asym_id\n")
        f.write("_atom_site.pdbx_PDB_model_num\n")

        serial = 1

        # Write protein atoms (entity 1)
        for resname, resseq, atoms in protein_atoms:
            for atom in atoms:
                group = "ATOM"
                element = atom.element if atom.element else atom.name[0]
                alt_id = "."
                ins_code = "?"
                # B_iso_or_equiv: 1 for backbone (real coords), 100 for sidechain (zeroed)
                b_factor = 1 if atom.name in BACKBONE_ATOMS else 100
                f.write(f"{group} {serial} {element} {atom.name} {alt_id} {resname} "
                       f"{protein_chain} 1 {resseq} {ins_code} "
                       f"{atom.x} {atom.y} {atom.z} "
                       f"1 {b_factor} {ins_code} {resseq} {protein_chain} 1\n")
                serial += 1

        # Write ligand atoms (entity 2)
        if ligand_atoms:
            for atom in ligand_atoms:
                group = "ATOM"
                element = atom.element if atom.element else atom.name[0]
                alt_id = "."
                ins_code = "?"
                # Ligand atoms: label_seq_id is 1, B_iso is 0
                f.write(f"{group} {serial} {element} {atom.name} {alt_id} {ligand_name} "
                       f"{ligand_chain} 2 1 {ins_code} "
                       f"{atom.x} {atom.y} {atom.z} "
                       f"1 0 {ins_code} 1 {ligand_chain} 1\n")
                serial += 1

        f.write("\n")

        # Write pdbx_poly_seq_scheme loop
        f.write("loop_\n")
        f.write("_pdbx_poly_seq_scheme.asym_id\n")
        f.write("_pdbx_poly_seq_scheme.entity_id\n")
        f.write("_pdbx_poly_seq_scheme.seq_id\n")
        f.write("_pdbx_poly_seq_scheme.mon_id\n")
        f.write("_pdbx_poly_seq_scheme.pdb_seq_num\n")
        f.write("_pdbx_poly_seq_scheme.auth_seq_num\n")
        f.write("_pdbx_poly_seq_scheme.pdb_mon_id\n")
        f.write("_pdbx_poly_seq_scheme.auth_mon_id\n")
        f.write("_pdbx_poly_seq_scheme.pdb_strand_id\n")
        f.write("_pdbx_poly_seq_scheme.pdb_ins_code\n")
        f.write("_pdbx_poly_seq_scheme.hetero\n")
        for resname, resseq, atoms in protein_atoms:
            f.write(f"{protein_chain} 1 {resseq} {resname} {resseq} {resseq} {resname} {resname} {protein_chain} . n\n")
        f.write("\n")

        # Write ma_qa_metric loop (pLDDT metadata)
        f.write("loop_\n")
        f.write("_ma_qa_metric.id\n")
        f.write("_ma_qa_metric.name\n")
        f.write("_ma_qa_metric.description\n")
        f.write("_ma_qa_metric.type\n")
        f.write("_ma_qa_metric.mode\n")
        f.write("_ma_qa_metric.type_other_details\n")
        f.write("_ma_qa_metric.software_group_id\n")
        f.write("1 pLDDT 'Predicted lddt' pLDDT local . .\n")
        f.write("\n")

        # Write ma_qa_metric_local loop (per-residue pLDDT)
        f.write("loop_\n")
        f.write("_ma_qa_metric_local.ordinal_id\n")
        f.write("_ma_qa_metric_local.model_id\n")
        f.write("_ma_qa_metric_local.label_asym_id\n")
        f.write("_ma_qa_metric_local.label_seq_id\n")
        f.write("_ma_qa_metric_local.label_comp_id\n")
        f.write("_ma_qa_metric_local.metric_id\n")
        f.write("_ma_qa_metric_local.metric_value\n")
        ordinal = 1
        # Protein residues get pLDDT of 60.000
        for resname, resseq, atoms in protein_atoms:
            f.write(f"{ordinal} 1 {protein_chain} {resseq} {resname} 1 60.000\n")
            ordinal += 1
        # Ligand gets pLDDT of 80.000
        if ligand_atoms:
            f.write(f"{ordinal} 1 {ligand_chain} 1 {ligand_name} 1 80.000\n")
        f.write("\n")


def convert_and_reassign_chains(
    input_path: str,
    output_path: str,
    protein_chain: str = "A",
    ligand_chain: str = "B",
    new_sequence: str = None,
    ligand_name: str = "LIG1"
) -> tuple:
    """
    Convert PDB structure: reassign chains, optionally apply sequence.

    Args:
        input_path: Input PDB file path
        output_path: Output PDB file path
        protein_chain: Target chain ID for protein
        ligand_chain: Target chain ID for ligand
        new_sequence: Optional new sequence to apply (inverse folding mode)
        ligand_name: Residue name for ligand (default: LIG1, must match BoltzGen pattern ^LIG\\d+)

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

    # Write output CIF (not PDB - CIF allows longer residue names like LIG1,
    # and BoltzGen's parse_mmcif accepts the mols dictionary for ligand resolution)
    write_cif(output_path, output_protein, output_ligand, protein_chain, ligand_chain, ligand_name)

    # Log what was written
    print(f"  Written CIF file: {os.path.basename(output_path)}")
    print(f"    Protein chain {protein_chain}: {n_protein_res} residues, {sum(len(atoms) for _, _, atoms in output_protein)} atoms")
    print(f"    Ligand chain {ligand_chain}: 1 residue ({ligand_name}), {n_ligand_atoms} atoms")

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

    The id_map pattern {"*": "*_<S>"} means structure ID can have
    recursive suffixes added to match sequence IDs.

    Args:
        structure_id: ID of the structure
        id_map: Mapping pattern, e.g., {"*": "*_<S>"}

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
        id_map: ID mapping pattern, e.g., {"*": "*_<S>"} for recursive suffixes.
               Default: {"*": "*_<S>"}

    Returns:
        Sequence string or None if not found
    """
    if id_map is None:
        id_map = {"*": "*_<S>"}

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
               Default {"*": "*_<S>"} handles recursive suffix stripping.
        ligand_name: Residue name for ligand in output (default: "LIG")

    Returns:
        Dictionary with import statistics
    """
    if id_map is None:
        id_map = {"*": "*_<S>"}
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
        # Using CIF format because:
        # 1. CIF allows longer residue names (LIG1) unlike PDB's 3-char limit
        # 2. BoltzGen's parse_mmcif accepts the mols dictionary for ligand resolution
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
            # Include ligand atoms in NPZ - BoltzGen tokenizes ligand atoms for co-folding
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
        default='{"*": "*_<S>"}',
        help='JSON string for ID mapping pattern (default: {"*": "*_<S>"}). '
             'Maps structure IDs to sequence IDs with recursive suffix stripping.'
    )
    parser.add_argument(
        '--ligand-name',
        default='LIG1',
        help='Residue name for ligand in output (default: LIG1). '
             'BoltzGen expects ligands to match pattern ^LIG\\d+ (e.g., LIG1, LIG2).'
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
