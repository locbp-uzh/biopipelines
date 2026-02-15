#!/usr/bin/env python3
# Copyright (c) 2026 Gianluca Quargnali @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
PDB file parser for BioPipelines helper scripts.

Lightweight PDB parsing functionality without external dependencies.
Provides atom selection and distance calculation utilities.
"""

import math
from typing import List, Dict, Any, Optional, Tuple, NamedTuple, Set


# Standard amino acid residue names — single source of truth
STANDARD_RESIDUES = {
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
    'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'
}

class Atom(NamedTuple):
    """Represents an atom from PDB file."""
    x: float
    y: float
    z: float
    atom_name: str
    res_name: str
    res_num: int
    chain: str
    element: str = ""

def parse_pdb_file(pdb_path: str) -> List[Atom]:
    """
    Parse PDB file and extract atom information.
    
    Args:
        pdb_path: Path to PDB file
        
    Returns:
        List of Atom objects
    """
    atoms = []
    
    with open(pdb_path, 'r') as f:
        for line in f:
            # Only parse ATOM and HETATM records
            if line.startswith(('ATOM', 'HETATM')):
                try:
                    # Parse PDB format according to specification
                    atom_name = line[12:16].strip()
                    res_name = line[17:20].strip()
                    chain = line[21:22].strip()
                    res_num = int(line[22:26].strip())
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    element = line[76:78].strip() if len(line) > 76 else atom_name[0]
                    
                    atom = Atom(
                        x=x, y=y, z=z,
                        atom_name=atom_name,
                        res_name=res_name,
                        res_num=res_num,
                        chain=chain,
                        element=element
                    )
                    atoms.append(atom)
                    
                except (ValueError, IndexError) as e:
                    # Skip malformed lines
                    continue
                    
    return atoms

def get_protein_sequence(atoms: List[Atom]) -> Dict[str, str]:
    """
    Extract protein sequence from atoms.
    
    Args:
        atoms: List of Atom objects
        
    Returns:
        Dictionary mapping chain -> sequence
    """
    # Standard amino acid mapping (three-letter to one-letter)
    aa_map = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
    }

    sequences = {}

    # Group by chain and get unique residues
    chain_residues = {}
    for atom in atoms:
        if atom.res_name in STANDARD_RESIDUES:
            if atom.chain not in chain_residues:
                chain_residues[atom.chain] = {}
            chain_residues[atom.chain][atom.res_num] = atom.res_name
    
    # Convert to sequences
    for chain, residues in chain_residues.items():
        sorted_residues = sorted(residues.items())
        sequence = ''.join([aa_map[res_name] for _, res_name in sorted_residues])
        sequences[chain] = sequence
    
    return sequences

def select_atoms_by_ligand(atoms: List[Atom], ligand_name: str, atom_name: str = None) -> List[Atom]:
    """
    Select atoms from a specific ligand.
    
    Args:
        atoms: List of all atoms
        ligand_name: Ligand residue name (e.g., 'LIG', 'HAL')
        atom_name: Specific atom name (e.g., 'Cl', 'Br'), None for all atoms
        
    Returns:
        List of selected atoms
    """
    selected = []
    for atom in atoms:
        if atom.res_name == ligand_name:
            if atom_name is None:
                selected.append(atom)
            else:
                # Check for exact match or element-based match (e.g., 'Cl' matches 'CL59')
                atom_clean = atom.atom_name.strip()
                if (atom_clean == atom_name or 
                    atom_clean.startswith(atom_name.upper()) or
                    atom.element.strip().upper() == atom_name.upper()):
                    selected.append(atom)
    return selected

def debug_ligand_atoms(atoms: List[Atom]) -> None:
    """Debug function to show all ligand residues and atoms."""
    ligand_residues = {}
    for atom in atoms:
        if atom.res_name not in STANDARD_RESIDUES:
            if atom.res_name not in ligand_residues:
                ligand_residues[atom.res_name] = set()
            ligand_residues[atom.res_name].add(atom.atom_name.strip())
    
    if ligand_residues:
        print("  - Found ligand residues:")
        for res_name, atom_names in ligand_residues.items():
            print(f"    {res_name}: {sorted(atom_names)}")
    else:
        print("  - No ligand residues found")

def select_atoms_by_residue_number(atoms: List[Atom], residue_numbers: List[int], chain: str = None) -> List[Atom]:
    """
    Select atoms from specific residue numbers.
    Supports negative indexing: -1 for last residue, -2 for second-to-last, etc.

    Args:
        atoms: List of all atoms
        residue_numbers: List of residue numbers to select (supports negative indexing)
        chain: Optional chain filter for negative indexing

    Returns:
        List of selected atoms
    """
    # Check if any negative indices are used
    has_negative = any(num < 0 for num in residue_numbers)

    if has_negative:
        # Build a list of unique protein residue numbers (sorted)
        protein_residues = {}
        for atom in atoms:
            if atom.res_name in STANDARD_RESIDUES:
                key = (atom.chain, atom.res_num)
                if key not in protein_residues:
                    protein_residues[key] = atom.res_num

        # Sort by chain and residue number
        sorted_residues = sorted(protein_residues.items(), key=lambda x: (x[0][0], x[0][1]))

        # If chain is specified, filter to that chain
        if chain:
            sorted_residues = [(k, v) for k, v in sorted_residues if k[0] == chain]

        # Convert negative indices to actual residue numbers
        actual_residue_numbers = []
        for num in residue_numbers:
            if num < 0:
                # Negative indexing: -1 is last, -2 is second-to-last, etc.
                idx = len(sorted_residues) + num
                if 0 <= idx < len(sorted_residues):
                    actual_res_num = sorted_residues[idx][1]
                    actual_residue_numbers.append(actual_res_num)
                    print(f"  - Negative index {num} -> residue {actual_res_num}")
                else:
                    print(f"  - Warning: Negative index {num} out of range (total residues: {len(sorted_residues)})")
            else:
                actual_residue_numbers.append(num)
    else:
        actual_residue_numbers = residue_numbers

    # Select atoms with the resolved residue numbers
    selected = []
    for atom in atoms:
        if atom.res_num in actual_residue_numbers:
            selected.append(atom)
    return selected

def select_atoms_by_sequence_context(atoms: List[Atom], target_residue: str, sequence_context: str) -> List[Atom]:
    """
    Select atoms by finding residue in sequence context.
    
    Args:
        atoms: List of all atoms
        target_residue: Single letter amino acid code (e.g., 'D')
        sequence_context: Sequence context (e.g., 'IGDWG')
        
    Returns:
        List of selected atoms
    """
    # Map single letter to three letter code
    residue_map = {
        'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS',
        'Q': 'GLN', 'E': 'GLU', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
        'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO',
        'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL'
    }
    
    target_three_letter = residue_map.get(target_residue, target_residue)
    
    # Get sequences and find matching positions
    sequences = get_protein_sequence(atoms)
    target_positions = []
    
    print(f"  - Looking for '{target_residue}' in context '{sequence_context}'")
    
    for chain, sequence in sequences.items():
        print(f"  - Chain {chain}: {len(sequence)} residues")
        print(f"    Sequence preview: {sequence[:50]}{'...' if len(sequence) > 50 else ''}")
        
        # Find all occurrences of sequence context
        start = 0
        while True:
            pos = sequence.find(sequence_context, start)
            if pos == -1:
                break
            
            print(f"    Found context '{sequence_context}' at position {pos}")
            
            # Find position of target residue within context
            target_pos_in_context = sequence_context.find(target_residue)
            if target_pos_in_context != -1:
                # Get ordered residue numbers for this chain
                chain_residues = [(a.res_num, a.res_name) for a in atoms 
                                if a.chain == chain and a.res_name in residue_map.values()]
                # Remove duplicates and sort by residue number
                unique_residues = list(set(chain_residues))
                unique_residues.sort(key=lambda x: x[0])
                
                if pos + target_pos_in_context < len(unique_residues):
                    target_res_num = unique_residues[pos + target_pos_in_context][0]
                    target_positions.append(target_res_num)
                    print(f"    Target '{target_residue}' at sequence pos {pos + target_pos_in_context} = residue {target_res_num}")
            
            start = pos + 1
    
    print(f"  - Found target positions: {target_positions}")
    
    # If no sequence context match, raise error - no fallbacks
    if not target_positions:
        raise ValueError(f"Sequence context '{sequence_context}' not found in protein sequence. "
                        f"Target residue '{target_residue}' could not be located in the specified context.")
    
    return select_atoms_by_residue_number(atoms, target_positions)

def resolve_selection(selection: str, atoms: List[Atom]) -> List[Atom]:
    """
    Parse selection string and return matching atoms.

    Supports:
    - Residue.atom: ``10.CA``, ``-1.C`` (numeric before dot)
    - Ligand.atom:  ``LIG.Cl`` (non-numeric before dot)
    - Sequence context: ``D in IGDWG``
    - Residue numbers: ``145``, ``-1``
    - Ranges: ``145-150``
    - Multiple: ``145+147+150``, ``1+-1``
    - Atom name fallback

    Args:
        selection: Selection string
        atoms: List of all atoms

    Returns:
        List of selected atoms
    """
    if '.' in selection:
        parts = selection.split('.', 1)
        prefix = parts[0]
        atom_name = parts[1] if len(parts) > 1 else None

        if prefix.lstrip('-').isdigit():
            # Residue.atom: '10.CA', '-1.C'
            residue_atoms = select_atoms_by_residue_number(atoms, [int(prefix)])
            if atom_name is None:
                return residue_atoms
            selected = []
            for atom in residue_atoms:
                atom_clean = atom.atom_name.strip()
                if atom_clean == atom_name or atom_clean.upper() == atom_name.upper():
                    selected.append(atom)
            return selected
        else:
            # Ligand.atom: 'LIG.Cl'
            return select_atoms_by_ligand(atoms, prefix, atom_name)

    elif ' in ' in selection:
        # Residue in sequence context: 'D in IGDWG'
        parts = selection.split(' in ')
        target_residue = parts[0].strip()
        sequence_context = parts[1].strip()
        return select_atoms_by_sequence_context(atoms, target_residue, sequence_context)

    elif selection.lstrip('-').isdigit():
        # Simple residue number: '145' or '-1' (negative indexing)
        res_num = int(selection)
        return select_atoms_by_residue_number(atoms, [res_num])

    elif '-' in selection and selection.count('-') == 1:
        # Could be residue range '145-150' or negative number '-1'
        parts = selection.split('-')
        # Check if it's a range (two positive numbers)
        if parts[0] and parts[1] and parts[0].isdigit() and parts[1].isdigit():
            # Residue range: '145-150'
            start, end = int(parts[0]), int(parts[1])
            res_nums = list(range(start, end + 1))
            return select_atoms_by_residue_number(atoms, res_nums)
        else:
            # Single negative number: '-1'
            res_num = int(selection)
            return select_atoms_by_residue_number(atoms, [res_num])

    elif '+' in selection:
        # Multiple residues: '145+147+150' or '1+-1' (supports negative)
        parts = selection.split('+')
        res_nums = []
        for part in parts:
            part = part.strip()
            if part.lstrip('-').isdigit():
                res_nums.append(int(part))
        return select_atoms_by_residue_number(atoms, res_nums)

    else:
        # Try as atom name
        selected = []
        for atom in atoms:
            if atom.atom_name == selection:
                selected.append(atom)
        return selected

def calculate_distance(atom1: Atom, atom2: Atom) -> float:
    """
    Calculate Euclidean distance between two atoms.
    
    Args:
        atom1: First atom
        atom2: Second atom
        
    Returns:
        Distance in Angstroms
    """
    dx = atom1.x - atom2.x
    dy = atom1.y - atom2.y
    dz = atom1.z - atom2.z
    return math.sqrt(dx*dx + dy*dy + dz*dz)

def calculate_distances(atoms1: List[Atom], atoms2: List[Atom], metric: str = "min") -> Optional[float]:
    """
    Calculate distance metric between two sets of atoms.
    
    Args:
        atoms1: First set of atoms
        atoms2: Second set of atoms  
        metric: Distance metric ("min", "max", "mean", "closest")
        
    Returns:
        Calculated distance or None if no atoms
    """
    if not atoms1 or not atoms2:
        return None
    
    distances = []
    for a1 in atoms1:
        for a2 in atoms2:
            dist = calculate_distance(a1, a2)
            distances.append(dist)
    
    if not distances:
        return None
    
    if metric == "min" or metric == "closest":
        return min(distances)
    elif metric == "max":
        return max(distances)
    elif metric == "mean":
        return sum(distances) / len(distances)
    else:
        return min(distances)  # Default to min


# ---------------------------------------------------------------------------
# PyMOL range string helpers (pure string <-> tuple conversion)
# ---------------------------------------------------------------------------

def parse_pymol_ranges(selection: str) -> List[Tuple[int, int]]:
    """
    Parse a PyMOL-style range string into ``(start, end)`` tuples.

    Pure string parsing — no atom/structure dependency.

    Examples:
        ``"3-45+58-60"``  -> ``[(3, 45), (58, 60)]``
        ``"10"``          -> ``[(10, 10)]``

    Args:
        selection: PyMOL selection string

    Returns:
        List of (start, end) tuples
    """
    if not selection or selection.strip() == "":
        return []

    ranges = []
    parts = selection.split('+')

    for part in parts:
        part = part.strip()
        if not part:
            continue

        if '-' in part:
            range_parts = part.split('-')
            if len(range_parts) != 2:
                raise ValueError(f"Invalid range format: {part}")
            start = int(range_parts[0])
            end = int(range_parts[1])
            if start > end:
                raise ValueError(f"Invalid range {part}: start > end")
            ranges.append((start, end))
        else:
            res_num = int(part)
            ranges.append((res_num, res_num))

    return ranges


def format_pymol_ranges(residue_numbers: List[int]) -> str:
    """
    Format a list of residue numbers as a compact ``start-end+...`` string.

    Inverse of :func:`parse_pymol_ranges` (but takes a flat list of ints
    rather than tuples).

    Args:
        residue_numbers: List of residue numbers

    Returns:
        Selection string with merged ranges (e.g., ``"10-11+15"``)
    """
    if not residue_numbers:
        return ""

    res_nums = sorted(set(residue_numbers))
    ranges: List[str] = []
    start = res_nums[0]
    end = res_nums[0]

    for i in range(1, len(res_nums)):
        if res_nums[i] == end + 1:
            end = res_nums[i]
        else:
            if start == end:
                ranges.append(f"{start}")
            else:
                ranges.append(f"{start}-{end}")
            start = end = res_nums[i]

    # Final range
    if start == end:
        ranges.append(f"{start}")
    else:
        ranges.append(f"{start}-{end}")

    return "+".join(ranges)