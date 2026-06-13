#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
PDB file parser for BioPipelines helper scripts.

Lightweight PDB parsing functionality without external dependencies.
Provides atom selection and distance calculation utilities.
"""

import math
import re
from typing import List, Dict, Any, Optional, Tuple, NamedTuple, Set


# Standard amino acid residue names — single source of truth
STANDARD_RESIDUES = {
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
    'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'
}

THREE_TO_ONE = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
}

# Theoretical max solvent-accessible surface area per residue (one-letter -> Å²).
MAX_ACC_TIEN = {
    'A': 129.0, 'R': 274.0, 'N': 195.0, 'D': 193.0, 'C': 167.0,
    'E': 223.0, 'Q': 225.0, 'G': 104.0, 'H': 224.0, 'I': 197.0,
    'L': 201.0, 'K': 236.0, 'M': 224.0, 'F': 240.0, 'P': 159.0,
    'S': 155.0, 'T': 172.0, 'W': 285.0, 'Y': 263.0, 'V': 174.0,
}


def relative_accessibility(acc: float, restype: str) -> Optional[float]:
    """rsa = acc / max-ACC; None when restype is not a standard residue."""
    if len(restype) == 3:
        restype = THREE_TO_ONE.get(restype.upper(), restype)
    max_acc = MAX_ACC_TIEN.get(restype)
    if not max_acc:
        return None
    return round(acc / max_acc, 4)

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


# ---------------------------------------------------------------------------
# Fixed-column field accessors — single source of truth for the PDB layout.
# Every pipe script reads ATOM/HETATM fields through these instead of inlining
# offsets, so a future format change is a one-file edit, not a repo-wide sweep.
# ---------------------------------------------------------------------------

def field_atom_name(line: str) -> str:
    """Atom name (cols 13-16)."""
    return line[12:16].strip()


def field_res_name(line: str) -> str:
    """Residue / CCD code (cols 18-20 in fixed-column PDB)."""
    return line[17:20].strip()


def field_chain(line: str) -> str:
    """Chain identifier (col 22)."""
    return line[21:22].strip()


def field_res_seq(line: str) -> str:
    """Residue sequence number as a string (cols 23-26), icode excluded."""
    return line[22:26].strip()


def field_coords(line: str) -> Tuple[float, float, float]:
    """(x, y, z) coordinates (cols 31-54)."""
    return (float(line[30:38]), float(line[38:46]), float(line[46:54]))

def _tokenize_cif_line(line: str) -> List[str]:
    """Split one mmCIF data line into tokens by CIF rules, not shell rules.

    A quote opens a quoted value only at a token start (after whitespace); it
    closes only when the matching quote is followed by whitespace or EOL. So an
    unquoted ``O5'`` is a single literal token (the prime is not a closing quote),
    which ``shlex`` would instead reject as an unterminated quote.
    """
    tokens: List[str] = []
    i = 0
    n = len(line)
    while i < n:
        if line[i].isspace():
            i += 1
            continue
        ch = line[i]
        if ch in ("'", '"'):
            j = i + 1
            while j < n:
                if line[j] == ch and (j + 1 >= n or line[j + 1].isspace()):
                    break
                j += 1
            tokens.append(line[i + 1:j])
            i = j + 1
        else:
            j = i
            while j < n and not line[j].isspace():
                j += 1
            tokens.append(line[i:j])
            i = j
    return tokens


def _cif_value(rec: Dict[str, str], *keys: str) -> Optional[str]:
    """First present, non-null value across ``keys`` for the auth->label fallback.

    mmCIF nulls (``.`` inapplicable, ``?`` unknown) count as absent so a missing
    ``auth_*`` value falls through to its ``label_*`` counterpart instead of
    reaching ``int(".")`` and dropping the atom.
    """
    for k in keys:
        v = rec.get(k)
        if v is not None and v not in (".", "?"):
            return v
    return None


def read_cif_loop(text: str, prefix: str) -> List[Dict[str, str]]:
    """Parse a single mmCIF ``loop_`` into a list of column->value dicts.

    ``prefix`` is the category tag (e.g. ``_atom_site.``). Tokens are split by
    CIF rules (``_tokenize_cif_line``) so an unquoted prime-bearing atom id like
    ``O5'`` stays one token. A single loop row may span several physical lines;
    tokens accumulate until one full record (``len(columns)`` values) is complete.
    """
    columns = []
    records = []
    reading_data = False
    pending: List[str] = []
    for line in text.splitlines():
        s = line.strip()
        if s.startswith(prefix):
            columns.append(s[len(prefix):])
            reading_data = False
            continue
        if columns and not reading_data:
            if not s or s.startswith("#") or s == "loop_" or s.startswith("data_"):
                continue
            reading_data = True
        if reading_data:
            if not s or s.startswith("_") or s.startswith("#") or s == "loop_" or s.startswith("data_"):
                break
            pending.extend(_tokenize_cif_line(s))
            while len(pending) >= len(columns):
                row = pending[:len(columns)]
                del pending[:len(columns)]
                records.append(dict(zip(columns, row)))
    return records


def parse_cif_file(cif_path: str) -> List[Atom]:
    """Parse an mmCIF file's ``_atom_site`` loop into Atom objects.

    Uses auth_* numbering (auth_asym_id / auth_seq_id) so deposited
    crystallographic chain ids and residue numbers are preserved; falls back
    to label_* only when the auth columns are absent.
    """
    with open(cif_path, 'r') as f:
        text = f.read()

    records = read_cif_loop(text, "_atom_site.")
    atoms = []
    for rec in records:
        try:
            atom_name = _cif_value(rec, "auth_atom_id", "label_atom_id")
            res_name = _cif_value(rec, "label_comp_id")
            chain = _cif_value(rec, "auth_asym_id", "label_asym_id") or ""
            res_num = int(_cif_value(rec, "auth_seq_id", "label_seq_id"))
            x = float(rec["Cartn_x"])
            y = float(rec["Cartn_y"])
            z = float(rec["Cartn_z"])
            element = (_cif_value(rec, "type_symbol") or atom_name[0])
            if atom_name is None or res_name is None:
                continue
        except (KeyError, ValueError, IndexError, TypeError):
            continue
        atoms.append(Atom(
            x=x, y=y, z=z,
            atom_name=atom_name,
            res_name=res_name,
            res_num=res_num,
            chain=chain,
            element=element,
        ))
    return atoms


def _is_mmcif(path: str, head: str) -> bool:
    """mmCIF detection: .cif extension or an ``_atom_site.`` loop in the head."""
    return path.lower().endswith((".cif", ".mmcif")) or "_atom_site." in head


def parse_models_file(path: str, records: Tuple[str, ...] = ("ATOM", "HETATM")) -> List[List[Atom]]:
    """Parse a structure file into a list of models, each a ``List[Atom]``.

    Handles multi-model ensembles for both formats: PDB MODEL/ENDMDL records and
    mmCIF's ``_atom_site.pdbx_PDB_model_num`` column. A single-model file returns
    a one-element list. ``records`` filters which record types to keep ("ATOM",
    "HETATM"); pass ``("ATOM",)`` for protein-only ensembles.
    """
    keep_atom = "ATOM" in records
    keep_hetatm = "HETATM" in records

    with open(path, 'r') as f:
        head = f.read(4096)
    if _is_mmcif(path, head):
        with open(path, 'r') as f:
            text = f.read()
        recs = read_cif_loop(text, "_atom_site.")
        models: Dict[str, List[Atom]] = {}
        order: List[str] = []
        for rec in recs:
            group = rec.get("group_PDB", "ATOM")
            if group == "ATOM" and not keep_atom:
                continue
            if group == "HETATM" and not keep_hetatm:
                continue
            try:
                atom_name = _cif_value(rec, "auth_atom_id", "label_atom_id")
                res_name = _cif_value(rec, "label_comp_id")
                chain = _cif_value(rec, "auth_asym_id", "label_asym_id") or ""
                res_num = int(_cif_value(rec, "auth_seq_id", "label_seq_id"))
                x = float(rec["Cartn_x"]); y = float(rec["Cartn_y"]); z = float(rec["Cartn_z"])
                element = (_cif_value(rec, "type_symbol") or atom_name[0])
                if atom_name is None or res_name is None:
                    continue
            except (KeyError, ValueError, IndexError, TypeError):
                continue
            model_key = rec.get("pdbx_PDB_model_num", "1")
            if model_key not in models:
                models[model_key] = []
                order.append(model_key)
            models[model_key].append(Atom(
                x=x, y=y, z=z, atom_name=atom_name, res_name=res_name,
                res_num=res_num, chain=chain, element=element,
            ))
        return [models[k] for k in order]

    starts = tuple(records)
    pdb_models: List[List[Atom]] = []
    current: List[Atom] = []
    seen_model_record = False
    with open(path, 'r') as f:
        for line in f:
            if line.startswith("MODEL"):
                seen_model_record = True
                current = []
            elif line.startswith("ENDMDL"):
                if current:
                    pdb_models.append(current)
                current = []
            elif line.startswith(starts):
                try:
                    atom_name = field_atom_name(line)
                    x, y, z = field_coords(line)
                    current.append(Atom(
                        x=x, y=y, z=z,
                        atom_name=atom_name,
                        res_name=field_res_name(line),
                        res_num=int(field_res_seq(line)),
                        chain=field_chain(line),
                        element=line[76:78].strip() if len(line) > 76 else atom_name[0],
                    ))
                except (ValueError, IndexError):
                    continue
    if current:
        pdb_models.append(current)
    if not seen_model_record and pdb_models:
        return pdb_models[:1]
    return pdb_models


def parse_pdb_file(pdb_path: str) -> List[Atom]:
    """
    Parse a structure file (PDB or mmCIF) into Atom objects.

    mmCIF input (``.cif`` extension or a file carrying an ``_atom_site.`` loop)
    is parsed via :func:`parse_cif_file`; everything else is read as
    fixed-column PDB. Both return the same ``List[Atom]``.

    Args:
        pdb_path: Path to a PDB or mmCIF file

    Returns:
        List of Atom objects
    """
    with open(pdb_path, 'r') as f:
        head = f.read(4096)
    if _is_mmcif(pdb_path, head):
        return parse_cif_file(pdb_path)

    atoms = []
    with open(pdb_path, 'r') as f:
        for line in f:
            # Only parse ATOM and HETATM records
            if line.startswith(('ATOM', 'HETATM')):
                try:
                    # Parse PDB format according to specification
                    atom_name = field_atom_name(line)
                    res_name = field_res_name(line)
                    chain = field_chain(line)
                    res_num = int(field_res_seq(line))
                    x, y, z = field_coords(line)
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
    
    # Convert to sequences, padding gaps with 'X' so the string index stays
    # aligned with residue number (a missing residue is one 'X', not a silent
    # frame shift that breaks position-based lookups and downstream tools).
    for chain, residues in chain_residues.items():
        nums = sorted(residues)
        sequence = ''.join(
            aa_map[residues[n]] if n in residues else 'X'
            for n in range(nums[0], nums[-1] + 1)
        )
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

    # Select atoms with the resolved residue numbers (optionally chain-restricted).
    selected = []
    for atom in atoms:
        if atom.res_num in actual_residue_numbers:
            if chain and atom.chain != chain:
                continue
            selected.append(atom)
    return selected

def select_atoms_by_sequence_context(atoms: List[Atom], target_residue: str, sequence_context: str) -> List[Atom]:
    """
    Select atoms by finding residue in sequence context.

    Disambiguation: the context is matched against every occurrence in the
    chain, and within each matched occurrence the FIRST instance of
    target_residue is taken. So pick a context that occurs once in the sequence
    (callers that need exactly one residue, e.g. covalent_linkage, enforce that),
    and place the intended residue ahead of any other copy of the same letter in
    the window (e.g. for the catalytic S in 'FAMCSTSKV', the first S is selected).

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

def _resolve_keyword_residue(keyword: str, atoms: List[Atom], chain: str = None) -> int:
    """
    Resolve 'first' or 'last' keyword to actual PDB residue number.

    Args:
        keyword: 'first' or 'last'
        atoms: List of all atoms
        chain: Optional chain filter

    Returns:
        Actual PDB residue number
    """
    protein_residues = sorted(
        {(a.chain, a.res_num) for a in atoms if a.res_name in STANDARD_RESIDUES}
    )
    if chain:
        protein_residues = [(c, r) for c, r in protein_residues if c == chain]
    if not protein_residues:
        raise ValueError(f"No protein residues found to resolve '{keyword}'")
    if keyword == "first":
        res_num = protein_residues[0][1]
    else:
        res_num = protein_residues[-1][1]
    print(f"  - Keyword '{keyword}' -> residue {res_num}")
    return res_num


def resolve_selection(selection: str, atoms: List[Atom]) -> List[Atom]:
    """
    Parse selection string and return matching atoms.

    Supports:
    - Residue.atom: ``10.CA``, ``-1.C`` (numeric before dot)
    - Chain-qualified residue.atom: ``A141.CB``, ``B1.SI81`` (single chain letter
      then a residue number — restricts to that chain; also disambiguates a ligand
      copy on a specific chain, e.g. ``B1.SI81`` for the ligand at chain B residue 1)
    - Keyword.atom: ``first.CA``, ``last.C``
    - Ligand.atom:  ``LIG.Cl`` (non-numeric before dot)
    - Sequence context: ``D in IGDWG``
    - Keywords: ``first``, ``last``
    - Residue numbers: ``145``, ``-1``
    - Ranges: ``145-150``
    - Multiple: ``145+147+150``, ``1+-1``, ``first+last``
    - Atom name fallback

    Args:
        selection: Selection string
        atoms: List of all atoms

    Returns:
        List of selected atoms
    """
    # Multi-term union for atom-set references: 'LIG.B41+LIG.B42', 'LIG.O3+LIG.Cl'.
    # Only when terms carry a '.' (atom/ligand refs); pure-number '+' (e.g. '1+-1',
    # 'first+last') stays with the dedicated residue-number branch below.
    if '+' in selection and '.' in selection:
        seen = set()
        union: List[Atom] = []
        for term in selection.split('+'):
            term = term.strip()
            if not term:
                continue
            for atom in resolve_selection(term, atoms):
                key = id(atom)
                if key not in seen:
                    seen.add(key)
                    union.append(atom)
        return union

    if '.' in selection:
        parts = selection.split('.', 1)
        prefix = parts[0]
        atom_name = parts[1] if len(parts) > 1 else None

        # Detect a chain-qualified residue prefix: one chain letter then a
        # (possibly negative) residue number, e.g. 'A141', 'B1', 'A-1'. A
        # ligand resname of this shape (B12, K21, T3) takes precedence, so the
        # chain reading only applies when the prefix isn't a present resname.
        chain_match = re.match(r'^([A-Za-z])(-?\d+)$', prefix)
        if chain_match and any(a.res_name == prefix for a in atoms):
            chain_match = None

        def _filter_by_atom_name(residue_atoms):
            if atom_name is None:
                return residue_atoms
            selected = []
            for atom in residue_atoms:
                atom_clean = atom.atom_name.strip()
                if atom_clean == atom_name or atom_clean.upper() == atom_name.upper():
                    selected.append(atom)
            return selected

        if prefix.lower() in ('first', 'last'):
            # Keyword.atom: 'first.CA', 'last.C'
            res_num = _resolve_keyword_residue(prefix.lower(), atoms)
            return _filter_by_atom_name(select_atoms_by_residue_number(atoms, [res_num]))
        elif chain_match:
            # Chain-qualified residue.atom: 'A141.CB', 'B1.SI81'. Restrict to that
            # chain so a residue number / ligand copy is disambiguated per chain.
            chain = chain_match.group(1)
            res_num = int(chain_match.group(2))
            return _filter_by_atom_name(
                select_atoms_by_residue_number(atoms, [res_num], chain=chain)
            )
        elif prefix.lstrip('-').isdigit():
            # Residue.atom: '10.CA', '-1.C' (any chain)
            return _filter_by_atom_name(select_atoms_by_residue_number(atoms, [int(prefix)]))
        else:
            # Ligand.atom: 'LIG.Cl'
            return select_atoms_by_ligand(atoms, prefix, atom_name)

    elif ' in ' in selection:
        # Residue in sequence context: 'D in IGDWG'
        parts = selection.split(' in ')
        target_residue = parts[0].strip()
        sequence_context = parts[1].strip()
        return select_atoms_by_sequence_context(atoms, target_residue, sequence_context)

    elif selection.lower() in ('first', 'last'):
        # Keyword: 'first' or 'last'
        res_num = _resolve_keyword_residue(selection.lower(), atoms)
        return select_atoms_by_residue_number(atoms, [res_num])

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
        # Multiple residues: '145+147+150', '1+-1', 'first+last' (supports negative & keywords)
        parts = selection.split('+')
        res_nums = []
        for part in parts:
            part = part.strip()
            if part.lower() in ('first', 'last'):
                res_nums.append(_resolve_keyword_residue(part.lower(), atoms))
            elif part.lstrip('-').isdigit():
                res_nums.append(int(part))
        return select_atoms_by_residue_number(atoms, res_nums)

    else:
        # Try as atom name
        selected = []
        for atom in atoms:
            if atom.atom_name == selection:
                selected.append(atom)
        return selected


_ONE_TO_THREE = {
    'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS',
    'Q': 'GLN', 'E': 'GLU', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
    'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO',
    'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL',
}


def resolve_selection_in_sequence(selection: str, sequence: str, chain: str = "A") -> List[int]:
    """Resolve a selection string against a raw sequence, returning residue numbers.

    Builds one pseudo-atom per residue (numbered 1..N) so the common
    ``resolve_selection`` grammar (``"C in ILIPCH"``, ``"145"``, ``"10-20"``,
    ``"first"``/``"last"`` …) works without a structure — for callers that only
    have a sequence (e.g. building a Boltz constraint from a designed sequence)."""
    pseudo = [
        Atom(0.0, 0.0, 0.0, "CA", _ONE_TO_THREE.get(aa, "UNK"), i + 1, chain, "C")
        for i, aa in enumerate(sequence)
        if aa in _ONE_TO_THREE
    ]
    return sorted({a.res_num for a in resolve_selection(selection, pseudo)})

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