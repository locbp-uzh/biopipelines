#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Runtime helper script for PDB tool.

Fetches protein structures with priority-based lookup: local_folder -> pdbs/ -> RCSB download.
Downloads are saved to both pdbs/ folder (for reuse) and tool output folder.
"""

import os
import sys
import argparse
import json
import re
import pandas as pd
from datetime import datetime
from typing import Dict, List, Any, Tuple, Optional
from pathlib import Path

# Add repo root to path so biopipelines package is importable
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from biopipelines.pdb_parser import (
    get_protein_sequence, parse_pdb_file,
    field_atom_name, field_res_name, field_chain, field_res_seq,
)
from biopipelines.id_patterns import expand_ids, contains_pattern, expand_file_pattern


def _is_rcsb_pdb_code(s: str) -> bool:
    """Return True if s looks like a 4-character alphanumeric RCSB PDB code."""
    return isinstance(s, str) and len(s) == 4 and s.isalnum()


def _chain_is_multi(chain) -> bool:
    """True when the chain selector covers more than one chain.

    chain="all"  -> all chains in the structure
    List[str]    -> the listed chain letters
    """
    return isinstance(chain, list) or chain == "all"


def _chain_filters_structure(chain) -> bool:
    """True when chain restricts the structure file to a single chain on disk.

    Only an explicit single chain letter triggers filter_chain_from_content.
    "auto" / "longest" keep the file as-is, "all" / List[str] keep all chains
    (per-chain splitting, when requested, is handled separately via
    split_chains).
    """
    if isinstance(chain, list):
        return False
    return chain not in ("longest", "auto", "all")


def fetch_rcsb_chain_sequences(pdb_code: str) -> List[Tuple[str, str]]:
    """
    Fetch per-chain sequences from the RCSB FASTA endpoint.

    Mirrors biopipelines.sequence._load_from_pdb_code: one entity per FASTA
    record, possibly listing multiple chain letters in 'Chains A, B, ...'.
    Returns one (chain_letter, sequence) pair per chain letter — the same
    sequence appears once per chain letter the entity covers, so downstream
    code can index by chain letter without resolving entity numbering.

    Args:
        pdb_code: 4-character RCSB code (case-insensitive)

    Returns:
        List of (chain_letter, sequence) tuples.

    Raises:
        ValueError if RCSB returns no protein chains for this code.
    """
    try:
        import requests
    except ImportError:
        raise ImportError("'requests' module is required to fetch sequences from RCSB")

    pdb_id = pdb_code.upper()
    url = f"https://www.rcsb.org/fasta/entry/{pdb_id}/display"
    response = requests.get(url, timeout=15,
                            headers={"User-Agent": "BioPipelines-PDB/1.0"})
    response.raise_for_status()

    pairs: List[Tuple[str, str]] = []
    current_chains: List[str] = []
    current_seq_parts: List[str] = []

    def flush():
        if current_chains and current_seq_parts:
            seq = "".join(current_seq_parts)
            for ch in current_chains:
                pairs.append((ch, seq))

    for line in response.text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            flush()
            current_chains = []
            current_seq_parts = []
            parts = line[1:].split("|")
            if len(parts) > 1:
                chain_field = parts[1].strip()
                chain_str = chain_field.replace("Chains", "").replace("Chain", "").strip()
                current_chains = [c.strip() for c in chain_str.split(",") if c.strip()]
        else:
            current_seq_parts.append(line)
    flush()

    if not pairs:
        raise ValueError(f"No chains found in RCSB FASTA for {pdb_id}")
    return pairs


def convert_cif_to_pdb(cif_content: str) -> str:
    """
    Convert CIF format to PDB format using BioPython.

    Args:
        cif_content: CIF file content

    Returns:
        PDB format content

    Raises:
        Exception: If conversion fails
    """
    import tempfile
    from Bio.PDB import MMCIFParser, PDBIO

    # Write CIF content to temporary file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.cif', delete=False) as tmp_cif:
        tmp_cif.write(cif_content)
        tmp_cif.flush()
        cif_path = tmp_cif.name

    try:
        # Parse CIF file
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure("structure", cif_path)

        # Write to PDB format
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as tmp_pdb:
            pdb_path = tmp_pdb.name

        io = PDBIO()
        io.set_structure(structure)
        io.save(pdb_path)

        # Read PDB content
        with open(pdb_path, 'r') as f:
            pdb_content = f.read()

        # Clean up temporary files
        os.unlink(cif_path)
        os.unlink(pdb_path)

        return pdb_content

    except Exception as e:
        # Clean up on error
        if os.path.exists(cif_path):
            os.unlink(cif_path)
        raise Exception(f"CIF to PDB conversion failed: {str(e)}")


def convert_pdb_to_cif(pdb_content: str) -> str:
    """
    Convert PDB format to mmCIF format using BioPython.

    Mirror of ``convert_cif_to_pdb`` for the opposite direction so PDB(convert="cif")
    works for local files and for RCSB downloads that arrived as .pdb.

    Args:
        pdb_content: PDB file content

    Returns:
        mmCIF format content

    Raises:
        Exception: If conversion fails
    """
    import tempfile
    from Bio.PDB import PDBParser, MMCIFIO

    with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as tmp_pdb:
        tmp_pdb.write(pdb_content)
        tmp_pdb.flush()
        pdb_path = tmp_pdb.name

    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("structure", pdb_path)

        with tempfile.NamedTemporaryFile(mode='w', suffix='.cif', delete=False) as tmp_cif:
            cif_path = tmp_cif.name

        io = MMCIFIO()
        io.set_structure(structure)
        io.save(cif_path)

        with open(cif_path, 'r') as f:
            cif_content = f.read()

        os.unlink(pdb_path)
        os.unlink(cif_path)

        return cif_content

    except Exception as e:
        if os.path.exists(pdb_path):
            os.unlink(pdb_path)
        raise Exception(f"PDB to CIF conversion failed: {str(e)}")


def build_rename_mapping(operations: List[Dict[str, Any]]) -> Dict[str, str]:
    """
    Build a mapping from old residue names to new names based on rename operations.

    Args:
        operations: List of operation dictionaries

    Returns:
        Dict mapping old_name -> new_name for rename operations
    """
    mapping = {}
    for op in operations:
        if op.get("op") == "rename":
            mapping[op["old"]] = op["new"]
    return mapping


def apply_operations(content: str, format: str, operations: List[Dict[str, Any]],
                     structure_id: Optional[str] = None) -> str:
    """
    Apply a sequence of operations to structure content.

    Args:
        content: Structure file content
        format: File format ("pdb" or "cif")
        operations: List of operation dictionaries with 'op' key and parameters
        structure_id: ID of the structure being processed (used to resolve a
                      per-structure TableReference selection for the remove op).

    Returns:
        Modified structure content
    """
    for op in operations:
        op_type = op.get("op")
        if op_type == "rename":
            content = apply_rename_operation(content, format, op["old"], op["new"])
        elif op_type == "remove":
            # Back-compat: accept the old key name if an older config carries it.
            remove_het = op.get("remove_hetatm", op.get("include_hetatm", True))
            content = apply_remove_operation(content, format, op["selection"], structure_id,
                                             remove_hetatm=remove_het)
        elif op_type == "break_bond":
            content = apply_break_bond_operation(content, format, op["atom1"], op["atom2"])
        elif op_type == "rotate_bond":
            content = apply_rotate_bond_operation(content, format, op["atom1"], op["atom2"], op["angle"])
        else:
            print(f"Warning: Unknown operation type '{op_type}', skipping")

    return content


def apply_rename_operation(content: str, format: str, old_name: str, new_name: str) -> str:
    """
    Rename a residue/ligand in the structure.

    For PDB format, modifies columns 18-20 (residue name) in ATOM/HETATM records.
    The new name is padded/truncated to fit the 3-character field, unless it's
    a special format like ":L:" which uses the full field width.

    Args:
        content: Structure file content
        format: File format ("pdb" or "cif")
        old_name: Current residue name (e.g., "LIG")
        new_name: New residue name (e.g., ":L:")

    Returns:
        Structure content with renamed residues
    """
    if format != "pdb":
        print(f"Warning: Rename operation not yet implemented for CIF format")
        return content

    lines = content.split('\n')
    modified_lines = []
    rename_count = 0

    for line in lines:
        if line.startswith(('ATOM', 'HETATM')):
            current_res_name = field_res_name(line)

            if current_res_name == old_name:
                # Format new name to fit in 3-character field
                # Right-justify, but for special names like ":L:" use as-is
                if len(new_name) <= 3:
                    formatted_name = new_name.rjust(3)
                else:
                    # For longer names, truncate (though this shouldn't normally happen)
                    formatted_name = new_name[:3]
                    print(f"Warning: Truncating residue name '{new_name}' to '{formatted_name}'")

                # Reconstruct the line with new residue name
                modified_line = line[:17] + formatted_name + line[20:]
                modified_lines.append(modified_line)
                rename_count += 1
            else:
                modified_lines.append(line)
        else:
            modified_lines.append(line)

    if rename_count > 0:
        print(f"  Renamed {rename_count} atoms: {old_name} → {new_name}")
    else:
        print(f"  Warning: No atoms found with residue name '{old_name}'")

    return '\n'.join(modified_lines)


def _parse_residue_selection(selection: str, present_names=None):
    """Parse a residue selection into (ranges, names).

    A token shaped like an (optionally chain-qualified) residue number or range
    is a **range**; anything else is a residue **name** to remove. So 'A1-83',
    '84', '1-83' are ranges while 'NAP', 'EDO', and digit-bearing CCD names like
    '9DP' or 'A7ZK' are names. Names and ranges may mix ('NAP+A1-83'). Returns
    ``(ranges, names)`` where ranges is a list of ``(chain, lo, hi)`` (chain None
    when unqualified) and names is a set of uppercased residue names.

    Forms: "1-83", "1-83+90-95", "A1-83", "84", "A84+A90-95", "NAP", "9DP+EDO".

    A token shaped like one chain letter + digits ('B12', 'K21') is ambiguous —
    chain B residue 12, or the CCD ligand named B12. `present_names` (the residue
    names actually in the structure) breaks the tie the same way
    ``resolve_selection`` does: a token equal to a present residue name is a name,
    else a range. With no `present_names` such a token defaults to a range.
    """
    # Optional leading chain letter(s), then a number or lo-hi range.
    range_re = re.compile(r'^([A-Za-z]*)(\d+)(?:-(\d+))?$')
    # One chain letter + digits, no range: the only chain/name-ambiguous shape.
    ambiguous_re = re.compile(r'^[A-Za-z]\d+$')
    present = {n.upper() for n in (present_names or ())}
    ranges = []
    names = set()
    for tok in str(selection).replace(" ", "").split("+"):
        if not tok:
            continue
        # A ligand name present in the structure wins over the chain reading
        # (matches resolve_selection's 'B12 takes precedence' rule).
        if ambiguous_re.match(tok) and tok.upper() in present:
            names.add(tok.upper())
            continue
        m = range_re.match(tok)
        if not m:
            # A '-' marks an intended range; if it didn't parse, it's malformed
            # (CCD names never contain '-'). Otherwise it's a residue name
            # (including digit-bearing CCD codes like '9DP', 'A7ZK').
            if "-" in tok:
                raise ValueError(
                    f"PDB.remove selection {selection!r} has an unparseable token {tok!r}. "
                    f"Use residue-number ranges ('1-83', 'A1-83', '1-83+90-95') or "
                    f"residue names ('NAP', 'EDO', 'NAP+EDO').")
            names.add(tok.upper())
            continue
        chain = m.group(1) or None
        lo = int(m.group(2))
        hi = int(m.group(3)) if m.group(3) is not None else lo
        ranges.append((chain, lo, hi))
    return ranges, names


def _resolve_remove_selection(selection: str, structure_id: Optional[str]) -> Optional[str]:
    """Resolve a remove-op selection to a literal selection string.

    Plain selection strings pass through. A TABLE_REFERENCE:path:column resolves
    to the cell for `structure_id`, matched via the framework's hierarchical ID
    mapping (exact/provenance/parent/child/sibling). A single-row table therefore
    broadcasts to every structure. Returns None if no row matches (caller skips).
    """
    if not isinstance(selection, str) or not selection.startswith("TABLE_REFERENCE:"):
        return selection
    from biopipelines.biopipelines_io import load_table
    from biopipelines.id_map_utils import get_mapped_ids
    table, column = load_table(selection)
    id_col = "id" if "id" in table.columns else table.columns[0]
    target_ids = [str(x) for x in table[id_col].tolist()]
    if structure_id is None:
        matched = target_ids[0] if len(target_ids) == 1 else None
    else:
        matched = get_mapped_ids([structure_id], target_ids).get(structure_id)
    if matched is None:
        print(f"  Warning: remove — no table row matched structure '{structure_id}'; "
              f"leaving structure unchanged")
        return None
    row = table[table[id_col].astype(str) == matched]
    cell = row.iloc[0][column]
    # An empty/NaN cell means "nothing to remove for this structure" — skip.
    if pd.isna(cell):
        return None
    text = str(cell).strip()
    if not text or text.lower() == "nan":
        return None
    return text


def apply_remove_operation(content: str, format: str, selection, structure_id: Optional[str] = None,
                           remove_hetatm: bool = True) -> str:
    """Remove residues matching `selection` from the structure.

    A range-shaped token (an optionally chain-qualified residue number or range)
    is a (chain,resseq) **range**; any other token is a residue **name** (a HETATM
    group like 'NAP', 'EDO', or a digit-bearing CCD code like '9DP'); the two may
    mix ('NAP+A1-83').

    With `remove_hetatm=True` (default), any residue in the selected ranges is
    dropped regardless of record type — a HETATM (ligand/ion) numbered inside the
    range is removed along with the protein residues. With `remove_hetatm=False`,
    only ATOM records are dropped by range and all HETATM records are preserved.
    A **named** residue is always removed wherever it appears (it is itself the
    HETATM target), so `remove_hetatm` does not gate name removal. Dependent
    records that reference a removed atom (ANISOU, CONECT, LINK) are dropped or
    rewritten so the result stays self-consistent. `selection` may be a literal
    string or a resolved TABLE_REFERENCE (per-structure by ID match)."""
    sel = _resolve_remove_selection(selection, structure_id)
    if sel is None:
        return content
    if format != "pdb":
        print("Warning: remove operation not yet implemented for CIF format")
        return content

    # Residue names present in the structure disambiguate 'B12'-shaped tokens
    # (chain B res 12 vs the ligand named B12) the way resolve_selection does.
    present_names = {field_res_name(l).upper() for l in content.split('\n')
                     if l.startswith(("ATOM", "HETATM"))}
    ranges, names = _parse_residue_selection(sel, present_names)
    if not ranges and not names:
        return content

    def _in_range(chain_id: str, resseq: int) -> bool:
        for ch, lo, hi in ranges:
            if ch is not None and ch != chain_id:
                continue
            if lo <= resseq <= hi:
                return True
        return False

    range_prefixes = ("ATOM", "HETATM") if remove_hetatm else ("ATOM",)

    def _line_residue(line):
        """(chain, resseq) of an ATOM/HETATM/ANISOU line, resseq None if unparseable."""
        chain_id = line[21] if len(line) > 21 else " "
        try:
            return chain_id, int(field_res_seq(line))
        except ValueError:
            return chain_id, None

    def _is_removed_atom(line):
        """True if this atom line matches the name/range selection and is dropped."""
        if names and field_res_name(line).upper() in names:
            return True
        if line.startswith(range_prefixes):
            _, resseq = _line_residue(line)
            if resseq is not None and _in_range(line[21] if len(line) > 21 else " ", resseq):
                return True
        return False

    # Two passes: first collect the serials of removed atoms so dependent
    # records (CONECT, ANISOU) that reference them can be dropped too.
    removed_serials = set()
    for line in content.split('\n'):
        if line.startswith(("ATOM", "HETATM")) and _is_removed_atom(line):
            serial = line[6:11].strip()
            if serial:
                removed_serials.add(serial)

    kept, removed = [], 0
    for line in content.split('\n'):
        if line.startswith(("ATOM", "HETATM")):
            if _is_removed_atom(line):
                removed += 1
                continue
        elif line.startswith("ANISOU"):
            # ANISOU shares its atom's serial; drop it when the atom is gone.
            if line[6:11].strip() in removed_serials:
                continue
        elif line.startswith("CONECT"):
            # CONECT lists bonded atom serials in 5-col fields from col 7 on.
            refs = [line[i:i + 5].strip() for i in range(6, len(line.rstrip()), 5)]
            refs = [r for r in refs if r]
            if any(r in removed_serials for r in refs):
                kept_refs = [r for r in refs if r not in removed_serials]
                # Drop the record if its base atom is gone, or if fewer than two
                # serials remain (a CONECT needs the base plus >=1 bonded partner;
                # a lone base records no bond). Otherwise rewrite without dead refs.
                if refs[0] in removed_serials or len(kept_refs) < 2:
                    continue
                line = "CONECT" + "".join(f"{r:>5}" for r in kept_refs)
        elif line.startswith("LINK"):
            # LINK names two residues by (resname, chain, resseq) at cols 18-27
            # and 48-57. Drop the record if either endpoint was removed.
            def _link_removed(rname, chain, resnum):
                if rname.upper() in names:
                    return True
                try:
                    return _in_range(chain, int(resnum))
                except ValueError:
                    return False
            r1 = _link_removed(line[17:20].strip(), line[21:22].strip(), line[22:26].strip())
            r2 = _link_removed(line[47:50].strip(), line[51:52].strip(), line[52:56].strip())
            if r1 or r2:
                continue
        kept.append(line)

    het_note = "ATOM+HETATM" if remove_hetatm else "ATOM only"
    print(f"  Removed {removed} atoms ({het_note}) matching selection '{sel}'")
    return '\n'.join(kept)


def _parse_atom_selection(sel: str):
    """Parse a BioPipelines '<residue>.<atom>' atom selection.

    Returns (chain or None, resnum or None, resname or None, atom_name). The
    residue part may be a chain-prefixed number ("A145"), a bare number ("145"),
    or a residue name ("LIG"). Examples: "A145.SG", "145.SG", "LIG.C12".
    """
    if "." not in sel:
        raise ValueError(f"atom selection must be '<residue>.<atom>', got {sel!r}")
    res_part, atom_name = sel.rsplit(".", 1)
    chain = resnum = resname = None
    # Optional leading chain letter(s).
    i = 0
    while i < len(res_part) and not (res_part[i].isdigit() or res_part[i] == "-"):
        i += 1
    if i > 0 and i < len(res_part) and (res_part[i].isdigit() or res_part[i] == "-"):
        chain = res_part[:i]
        res_part = res_part[i:]
    if res_part.lstrip("-").isdigit():
        resnum = int(res_part)
    else:
        resname = res_part  # e.g. "LIG"
    return chain, resnum, resname, atom_name.strip()


def _atom_matches(line: str, sel) -> bool:
    """True if a PDB ATOM/HETATM line matches a parsed atom selection tuple."""
    chain, resnum, resname, atom_name = sel
    if field_atom_name(line) != atom_name:
        return False
    if chain is not None and (len(line) <= 21 or line[21] != chain):
        return False
    if resname is not None and field_res_name(line) != resname:
        return False
    if resnum is not None:
        try:
            if int(field_res_seq(line)) != resnum:
                return False
        except ValueError:
            return False
    return True


def apply_break_bond_operation(content: str, format: str, atom1: str, atom2: str) -> str:
    """Remove CONECT/LINK records joining the two named atoms (PyMOL `unbond`).

    Coordinates are untouched; only connectivity records are stripped so a
    covalently-tethered ligand becomes a separate non-covalent entity for
    downstream design. Atom selections use '<residue>.<atom>' syntax."""
    if format != "pdb":
        print("Warning: break_bond not yet implemented for CIF format")
        return content

    sel1, sel2 = _parse_atom_selection(atom1), _parse_atom_selection(atom2)
    # Resolve each selection to the matching atom serial numbers.
    serials1, serials2 = set(), set()
    for line in content.split('\n'):
        if line.startswith(("ATOM", "HETATM")):
            try:
                serial = int(line[6:11])
            except ValueError:
                continue
            if _atom_matches(line, sel1):
                serials1.add(serial)
            if _atom_matches(line, sel2):
                serials2.add(serial)

    if not serials1 or not serials2:
        print(f"  Warning: break_bond — no atoms matched "
              f"({atom1!r}->{len(serials1)}, {atom2!r}->{len(serials2)}); nothing broken")
        return content
    # Unlike rotate_bond, break_bond legitimately strips every matching
    # connectivity record, but a too-broad selection silently breaks more bonds
    # than intended — warn so the user can qualify it with chain/resnum.
    if len(serials1) > 1 or len(serials2) > 1:
        print(f"  Warning: break_bond — ambiguous selection "
              f"({atom1!r}->{len(serials1)} atoms, {atom2!r}->{len(serials2)} atoms); "
              f"breaking every matching bond. Qualify with chain/resnum to narrow.")

    kept, removed = [], 0
    for line in content.split('\n'):
        if line.startswith("CONECT"):
            nums = [int(x) for x in line[6:].split()]
            if not nums:
                kept.append(line); continue
            base, bonded = nums[0], nums[1:]
            # Drop references that join an atom in set1 to one in set2 (either order).
            drop = (base in serials1 and any(b in serials2 for b in bonded)) or \
                   (base in serials2 and any(b in serials1 for b in bonded))
            if drop:
                remaining = [b for b in bonded
                             if not ((base in serials1 and b in serials2) or
                                     (base in serials2 and b in serials1))]
                removed += len(bonded) - len(remaining)
                if remaining:
                    kept.append("CONECT" + "".join(f"{base:5d}" + "".join(f"{r:5d}" for r in remaining)))
                # else: drop the now-empty CONECT entirely
            else:
                kept.append(line)
        elif line.startswith("LINK"):
            # LINK names atoms by name/resname/chain/resseq in fixed columns.
            a1 = line[12:16].strip(); r1 = line[17:20].strip(); c1 = line[21:22]; n1 = line[22:26].strip()
            a2 = line[42:46].strip(); r2 = line[47:50].strip(); c2 = line[51:52]; n2 = line[52:56].strip()
            def _link_end_matches(an, rn, ch, ns, sel):
                cs, rns, rname, atn = sel
                if atn != an: return False
                if cs is not None and cs != ch: return False
                if rname is not None and rname != rn: return False
                if rns is not None and str(rns) != ns: return False
                return True
            j = (_link_end_matches(a1,r1,c1,n1,sel1) and _link_end_matches(a2,r2,c2,n2,sel2)) or \
                (_link_end_matches(a1,r1,c1,n1,sel2) and _link_end_matches(a2,r2,c2,n2,sel1))
            if j:
                removed += 1
                continue
            kept.append(line)
        else:
            kept.append(line)

    print(f"  break_bond: removed {removed} connectivity record(s) between {atom1} and {atom2}")
    return '\n'.join(kept)


def _atom_coords(line: str):
    """Parse (x, y, z) from an ATOM/HETATM line."""
    return (float(line[30:38]), float(line[38:46]), float(line[46:54]))


def apply_rotate_bond_operation(content: str, format: str, atom1: str, atom2: str,
                                angle: float) -> str:
    """Rotate the fragment on atom2's side of the atom1-atom2 bond by `angle` deg.

    atom1's side stays fixed; the moving fragment is the set of atoms reachable
    from atom2 (within atom2's residue) without crossing the atom1-atom2 bond,
    found by a connectivity BFS over inferred bonds (distance <= 1.9 A). Atoms
    are rotated about the atom1->atom2 axis through atom2 (Rodrigues' formula).
    Atom selections use '<residue>.<atom>' syntax."""
    import math

    if format != "pdb":
        print("Warning: rotate_bond not yet implemented for CIF format")
        return content

    sel1, sel2 = _parse_atom_selection(atom1), _parse_atom_selection(atom2)
    lines = content.split('\n')

    # Index atom lines: collect serial, coords, residue key, and match flags.
    atoms = []  # list of dicts for ATOM/HETATM lines
    matches1, matches2 = [], []
    for li, line in enumerate(lines):
        if line.startswith(("ATOM", "HETATM")):
            try:
                coords = _atom_coords(line)
            except ValueError:
                continue
            reskey = (line[21:22], line[22:27])  # chain, resseq+icode
            rec = {"li": li, "coords": coords, "reskey": reskey}
            atoms.append(rec)
            if _atom_matches(line, sel1):
                matches1.append(len(atoms) - 1)
            if _atom_matches(line, sel2):
                matches2.append(len(atoms) - 1)

    if not matches1 or not matches2:
        print(f"  Warning: rotate_bond — atom not found "
              f"({atom1!r}->{len(matches1)}, {atom2!r}->{len(matches2)}); nothing rotated")
        return content
    # A rotate must resolve each end to exactly one atom; an ambiguous
    # selection (e.g. "CA" with no chain/resnum) would otherwise silently
    # rotate about whichever match comes last.
    if len(matches1) > 1 or len(matches2) > 1:
        raise ValueError(
            f"rotate_bond: ambiguous selection — {atom1!r} matched {len(matches1)} atom(s), "
            f"{atom2!r} matched {len(matches2)} atom(s). Each end must resolve to exactly one "
            f"atom; qualify the selection with chain/resnum (e.g. 'A64.C64')."
        )
    idx1, idx2 = matches1[0], matches2[0]

    # Restrict connectivity to atom2's residue (intra-residue fragment).
    res = atoms[idx2]["reskey"]
    members = [i for i, a in enumerate(atoms) if a["reskey"] == res]

    def dist2(a, b):
        return sum((a[k] - b[k]) ** 2 for k in range(3))

    BOND2 = 1.9 ** 2
    # Build adjacency within the residue, excluding the atom1-atom2 bond itself.
    adj = {i: [] for i in members}
    for ii in range(len(members)):
        for jj in range(ii + 1, len(members)):
            a, b = members[ii], members[jj]
            if {a, b} == {idx1, idx2}:
                continue  # do not cross the rotatable bond
            if dist2(atoms[a]["coords"], atoms[b]["coords"]) <= BOND2:
                adj[a].append(b)
                adj[b].append(a)

    # BFS from atom2; everything reachable (without crossing the bond) moves.
    moving, stack = set(), [idx2]
    while stack:
        cur = stack.pop()
        if cur in moving:
            continue
        moving.add(cur)
        for nb in adj[cur]:
            if nb not in moving:
                stack.append(nb)

    if idx1 in moving:
        print(f"  Warning: rotate_bond — atom1 reachable from atom2 (ring/no clean "
              f"split across {atom1}-{atom2}); nothing rotated")
        return content

    # Rotation axis: unit vector atom1 -> atom2, pivot at atom2.
    p1, p2 = atoms[idx1]["coords"], atoms[idx2]["coords"]
    ax = [p2[k] - p1[k] for k in range(3)]
    norm = math.sqrt(sum(c * c for c in ax))
    if norm == 0:
        print("  Warning: rotate_bond — atom1 and atom2 coincide; nothing rotated")
        return content
    ux, uy, uz = (ax[0] / norm, ax[1] / norm, ax[2] / norm)
    theta = math.radians(angle)
    ct, st = math.cos(theta), math.sin(theta)

    def rotate(pt):
        # Rodrigues about axis u through pivot p2.
        vx, vy, vz = (pt[0] - p2[0], pt[1] - p2[1], pt[2] - p2[2])
        dot = ux * vx + uy * vy + uz * vz
        cx = uy * vz - uz * vy
        cy = uz * vx - ux * vz
        cz = ux * vy - uy * vx
        rx = vx * ct + cx * st + ux * dot * (1 - ct)
        ry = vy * ct + cy * st + uy * dot * (1 - ct)
        rz = vz * ct + cz * st + uz * dot * (1 - ct)
        return (rx + p2[0], ry + p2[1], rz + p2[2])

    for i in moving:
        if i == idx2:
            continue  # pivot stays fixed
        a = atoms[i]
        nx, ny, nz = rotate(a["coords"])
        # PDB coordinate columns are fixed 8-char %8.3f: anything outside
        # [-999.999, 9999.999] overflows the field and corrupts the record.
        for c in (nx, ny, nz):
            if not (-999.999 <= c <= 9999.999):
                raise ValueError(
                    f"rotate_bond: rotated coordinate {c:.3f} (atom serial near "
                    f"line {a['li'] + 1}) exceeds the PDB 8-char column range "
                    f"[-999.999, 9999.999]."
                )
        line = lines[a["li"]]
        lines[a["li"]] = (line[:30] + f"{nx:8.3f}{ny:8.3f}{nz:8.3f}" + line[54:])

    print(f"  rotate_bond: rotated {len(moving) - 1} atom(s) on {atom2}'s side of "
          f"{atom1}-{atom2} by {angle} deg")
    return '\n'.join(lines)


def remove_waters_from_content(content: str, format: str) -> str:
    """
    Remove water molecules from structure content.

    Args:
        content: Structure file content
        format: File format ("pdb" or "cif")

    Returns:
        Structure content with waters removed
    """
    if format != "pdb":
        # For CIF format, just return as-is for now
        return content

    lines = content.split('\n')
    filtered_lines = []

    for line in lines:
        # Skip water molecules (HOH) and other common solvent molecules
        if line.startswith(('ATOM', 'HETATM')):
            res_name = field_res_name(line)
            if res_name in ['HOH', 'WAT', 'H2O', 'SOL', 'TIP3', 'TIP4', 'SPC']:
                continue  # Skip water lines
        filtered_lines.append(line)

    return '\n'.join(filtered_lines)


def filter_chain_from_content(content: str, format: str, chain: str) -> str:
    """
    Filter structure content to keep only the specified chain.

    Args:
        content: Structure file content
        format: File format ("pdb" or "cif")
        chain: Chain identifier to keep (e.g., "A", "E")

    Returns:
        Structure content with only the specified chain
    """
    if format == "pdb":
        lines = content.split('\n')
        filtered_lines = []

        for line in lines:
            if line.startswith(('ATOM', 'HETATM')):
                if len(line) > 21 and line[21] != chain:
                    continue
            filtered_lines.append(line)

        return '\n'.join(filtered_lines)

    elif format == "cif":
        # CIF _atom_site lines are whitespace-delimited.
        # auth_asym_id (chain ID) column index is determined from the header.
        lines = content.split('\n')
        filtered_lines = []
        in_atom_site = False
        column_names = []
        chain_col_idx = None

        for line in lines:
            stripped = line.strip()

            # Detect start of _atom_site loop
            if stripped == 'loop_' and not in_atom_site:
                filtered_lines.append(line)
                in_atom_site = True
                column_names = []
                chain_col_idx = None
                continue

            if in_atom_site:
                if stripped.startswith('_atom_site.'):
                    column_names.append(stripped)
                    if stripped == '_atom_site.auth_asym_id':
                        chain_col_idx = len(column_names) - 1
                    filtered_lines.append(line)
                    continue

                # Data line within _atom_site block
                if column_names and (stripped.startswith('ATOM') or stripped.startswith('HETATM')):
                    if chain_col_idx is not None:
                        fields = stripped.split()
                        if len(fields) > chain_col_idx and fields[chain_col_idx] != chain:
                            continue
                    filtered_lines.append(line)
                    continue

                # End of _atom_site block
                in_atom_site = False
                column_names = []
                chain_col_idx = None

            filtered_lines.append(line)

        return '\n'.join(filtered_lines)

    return content


def extract_sequence_from_structure(content: str, format: str, chain: str = "longest") -> str:
    """
    Extract protein sequence from structure content.

    By default returns the sequence of the longest chain. Structures often contain
    multiple identical or non-identical chains (e.g. dimers, complexes); concatenating
    them would produce an artificially long sequence that inflates memory usage in
    downstream tools such as LigandMPNN.

    Args:
        content: Structure file content
        format: File format ("pdb" or "cif")
        chain: Which chain to extract. "longest" (default) picks the longest chain.
            Specify a chain letter (e.g. "A", "B") to select that chain.

    Returns:
        Protein sequence from the selected chain
    """
    import tempfile

    if format == "cif":
        return extract_sequence_from_cif(content, chain=chain)

    # PDB format
    try:
        # Write content to temporary file for parsing
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as tmp:
            tmp.write(content)
            tmp.flush()

            # Parse using the existing PDB parser
            atoms = parse_pdb_file(tmp.name)
            sequences_dict = get_protein_sequence(atoms)

            # Clean up temporary file
            os.unlink(tmp.name)

            if not sequences_dict:
                return ""

            # If a specific chain is requested, return it
            if _chain_filters_structure(chain):
                if chain in sequences_dict:
                    return sequences_dict[chain]
                else:
                    available = ', '.join(sorted(sequences_dict.keys()))
                    print(f"Warning: Chain '{chain}' not found in structure. Available chains: {available}. Falling back to longest chain.")

            # Return the longest chain sequence (sequences_dict already contains only
            # standard amino-acid chains; water/ions are excluded by get_protein_sequence)
            candidates = [seq for seq in sequences_dict.values() if len(seq) >= 10]
            if not candidates:
                return max(sequences_dict.values(), key=len)
            return max(candidates, key=len)

    except Exception as e:
        print(f"Warning: Could not extract sequence from PDB - {str(e)}")
        return ""


def extract_sequence_from_cif(content: str, chain: str = "longest") -> str:
    """
    Extract protein sequence from CIF content using BioPython.

    By default returns the sequence of the longest chain. Structures often contain
    multiple identical or non-identical chains (e.g. dimers, complexes); concatenating
    them would produce an artificially long sequence that inflates memory usage in
    downstream tools such as LigandMPNN.

    Args:
        content: CIF file content
        chain: Which chain to extract. "longest" (default) picks the longest chain.
            Specify a chain letter (e.g. "A", "B") to select that chain.

    Returns:
        Protein sequence from the selected chain
    """
    import tempfile

    # Standard amino acid mapping (3-letter to 1-letter)
    aa_map = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
    }

    try:
        from Bio.PDB import MMCIFParser

        # Write content to temporary file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.cif', delete=False) as tmp:
            tmp.write(content)
            tmp.flush()
            cif_path = tmp.name

        try:
            # Parse CIF file
            parser = MMCIFParser(QUIET=True)
            structure = parser.get_structure("structure", cif_path)

            # Extract sequences from all chains (use chain_obj to avoid shadowing the chain parameter)
            sequences = {}
            for model in structure:
                for chain_obj in model:
                    chain_id = chain_obj.id
                    residues = []
                    for residue in chain_obj:
                        res_name = residue.get_resname()
                        if res_name in aa_map:
                            res_id = residue.get_id()[1]  # Residue number
                            residues.append((res_id, aa_map[res_name]))

                    if residues:
                        # Sort by residue number and concatenate
                        residues.sort(key=lambda x: x[0])
                        sequences[chain_id] = ''.join([r[1] for r in residues])

            # Clean up
            os.unlink(cif_path)

            if not sequences:
                return ""

            # If a specific chain is requested, return it
            if _chain_filters_structure(chain):
                if chain in sequences:
                    print(f"  Extracted sequence from CIF: {len(sequences[chain])} residues (chain {chain})")
                    return sequences[chain]
                else:
                    available = ', '.join(sorted(sequences.keys()))
                    print(f"Warning: Chain '{chain}' not found in CIF. Available chains: {available}. Falling back to longest chain.")

            # Return the longest chain sequence (sequences already contains only
            # standard amino-acid chains; water/ions are excluded by the aa_map filter above)
            candidates = [seq for seq in sequences.values() if len(seq) >= 10]
            if not candidates:
                longest = max(sequences.values(), key=len)
            else:
                longest = max(candidates, key=len)
            print(f"  Extracted sequence from CIF: {len(longest)} residues (longest of {len(sequences)} chain(s))")
            return longest

        except Exception as e:
            # Clean up on error
            if os.path.exists(cif_path):
                os.unlink(cif_path)
            raise e

    except ImportError:
        print("Warning: BioPython not available for CIF sequence extraction")
        return ""
    except Exception as e:
        print(f"Warning: Could not extract sequence from CIF - {str(e)}")
        return ""


def extract_ligands_from_structure(content: str, format: str) -> List[str]:
    """
    Extract ligand identifiers (3-letter codes) from structure content.

    Args:
        content: Structure file content
        format: File format ("pdb" or "cif")

    Returns:
        List of unique ligand 3-letter codes (e.g., ['ATP', 'GDP'])
    """
    # Common non-ligand residues to exclude (solvents, ions, crystallization agents, caps)
    exclude_residues = {
        'HOH', 'WAT', 'H2O', 'SOL', 'TIP3', 'TIP4', 'SPC',
        'NA', 'CL', 'K', 'CA', 'MG', 'ZN', 'MN', 'FE', 'CU', 'NI', 'CO',
        'SO4', 'PO4', 'NO3',
        'GOL', 'EDO', 'PEG', 'PGE', 'PE4', 'PE3', 'P6G', 'PG4', '1PE',
        'ACT', 'ACE', 'ACY',
        'PYR', 'PYO',
        'DMS', 'BME', 'MPD', 'TRS', 'EPE',
        'NME', 'NH2'
    }

    if format == "cif":
        ligands = set()
        lines = content.split('\n')
        in_atom_site = False
        column_names = []
        group_col_idx = None   # _atom_site.group_PDB  (ATOM / HETATM)
        comp_col_idx = None    # _atom_site.label_comp_id  (residue/ligand name)

        for line in lines:
            stripped = line.strip()

            if stripped == 'loop_':
                in_atom_site = False
                column_names = []
                group_col_idx = None
                comp_col_idx = None
                # Don't append yet — wait to see if this is an _atom_site loop
                continue

            if stripped.startswith('_atom_site.'):
                in_atom_site = True
                column_names.append(stripped)
                if stripped == '_atom_site.group_PDB':
                    group_col_idx = len(column_names) - 1
                elif stripped == '_atom_site.label_comp_id':
                    comp_col_idx = len(column_names) - 1
                continue

            if in_atom_site and column_names:
                if stripped.startswith('ATOM') or stripped.startswith('HETATM'):
                    fields = stripped.split()
                    if (group_col_idx is not None and comp_col_idx is not None
                            and len(fields) > max(group_col_idx, comp_col_idx)):
                        if fields[group_col_idx] == 'HETATM':
                            res_name = fields[comp_col_idx]
                            if res_name and res_name not in exclude_residues:
                                ligands.add(res_name)
                    continue
                elif not stripped.startswith('_'):
                    in_atom_site = False

        return sorted(set(ligands))

    ligands = set()
    lines = content.split('\n')

    for line in lines:
        if line.startswith('HETATM'):
            res_name = field_res_name(line)
            # Only include if not in exclusion list (1-5 char extended CCD)
            if res_name and len(res_name) <= 5 and res_name not in exclude_residues:
                ligands.add(res_name)

    return sorted(list(ligands))


def fetch_ligand_smiles_from_rcsb(ligand_code: str) -> Optional[str]:
    """
    Fetch SMILES string for a ligand from RCSB REST API.

    Args:
        ligand_code: 3-letter ligand code (e.g., 'ATP', 'PVY')

    Returns:
        SMILES string or None if not found
    """
    try:
        import requests

        # RCSB REST API endpoint for ligand info
        url = f"https://data.rcsb.org/rest/v1/core/chemcomp/{ligand_code}"

        headers = {
            'User-Agent': 'BioPipelines-PDB/1.0 (https://github.com/locbp-uzh/biopipelines)'
        }

        response = requests.get(url, headers=headers, timeout=10)
        response.raise_for_status()

        data = response.json()

        # Extract SMILES from the response
        # RCSB provides SMILES in the 'chem_comp' section
        if 'chem_comp' in data:
            chem_comp = data['chem_comp']
            # Try different possible fields for SMILES
            for field in ['pdbx_smiles_canonical', 'smiles', 'smiles_canonical']:
                if field in chem_comp and chem_comp[field]:
                    print(f"  Found SMILES for {ligand_code}: {chem_comp[field][:50]}...")
                    return chem_comp[field]

        # Alternative location in descriptors
        if 'rcsb_chem_comp_descriptor' in data:
            descriptors = data['rcsb_chem_comp_descriptor']
            if isinstance(descriptors, dict):
                if 'smiles_canonical' in descriptors:
                    print(f"  Found SMILES for {ligand_code}: {descriptors['smiles_canonical'][:50]}...")
                    return descriptors['smiles_canonical']
                elif 'smiles' in descriptors:
                    print(f"  Found SMILES for {ligand_code}: {descriptors['smiles'][:50]}...")
                    return descriptors['smiles']

        print(f"  No SMILES found for {ligand_code} in RCSB response")
        return None

    except Exception as e:
        print(f"  Warning: Could not fetch SMILES for {ligand_code}: {str(e)}")
        return None


def find_local_structure(pdb_id: str, convert: Optional[str], local_folder: str,
                        repo_pdbs_folder: str) -> Optional[Tuple[str, str]]:
    """
    Find structure file locally.

    When convert is None: accepts any format found (pdb or cif), no conversion.
    When convert is "pdb": tries pdb first, then cif (will convert cif -> pdb)
    When convert is "cif": tries cif first, then pdb (will convert pdb -> cif)
    Search locations: local_folder (if given) -> repo_pdbs_folder

    Args:
        pdb_id: PDB identifier (with or without extension)
        convert: Target format to convert to ("pdb", "cif", or None to keep whatever is found)
        local_folder: Custom local folder (can be None)
        repo_pdbs_folder: Repository PDBs folder

    Returns:
        Tuple of (path to local file, actual format) or None if not found
    """
    # Remove extension if provided
    pdb_id_base = pdb_id.replace('.pdb', '').replace('.cif', '')

    search_locations = []
    if local_folder:
        search_locations.append(local_folder)
    search_locations.append(repo_pdbs_folder)

    if convert is None:
        # Accept any format: try pdb first, then cif - no conversion needed
        for ext, fmt in [(".pdb", "pdb"), (".cif", "cif")]:
            for location in search_locations:
                candidate = os.path.join(location, f"{pdb_id_base}{ext}")
                if os.path.exists(candidate):
                    print(f"Found {pdb_id_base} locally: {candidate}")
                    return candidate, fmt
        return None

    # Specific conversion target: primary then fallback with conversion
    if convert == "pdb":
        primary_ext, primary_fmt = ".pdb", "pdb"
        fallback_ext, fallback_fmt = ".cif", "cif"
    else:
        primary_ext, primary_fmt = ".cif", "cif"
        fallback_ext, fallback_fmt = ".pdb", "pdb"

    # Try primary format first (no conversion needed)
    for location in search_locations:
        candidate = os.path.join(location, f"{pdb_id_base}{primary_ext}")
        if os.path.exists(candidate):
            print(f"Found {pdb_id_base} locally: {candidate}")
            return candidate, primary_fmt

    # Try fallback format (will need conversion)
    for location in search_locations:
        candidate = os.path.join(location, f"{pdb_id_base}{fallback_ext}")
        if os.path.exists(candidate):
            print(f"Found {pdb_id_base} as {fallback_fmt.upper()} (will convert to {convert.upper()}): {candidate}")
            return candidate, fallback_fmt

    return None


def copy_local_structure(pdb_id: str, custom_id: str, source_path: str,
                        source_format: str, convert: Optional[str], remove_waters: bool,
                        output_folder: str, operations: List[Dict[str, Any]] = None,
                        chain: str = "longest",
                        fetch_compounds: bool = True) -> Tuple[bool, str, str, List[Dict[str, str]], Dict[str, Any]]:
    """
    Copy local structure file to output folder with optional format conversion.

    Args:
        pdb_id: PDB identifier
        custom_id: Custom ID for output filename
        source_path: Path to local structure file
        source_format: Format of source file ("pdb" or "cif")
        convert: Target format to convert to ("pdb", "cif", or None to keep source format as-is)
        remove_waters: Whether to remove water molecules
        output_folder: Directory to save the structure
        operations: List of operations to apply (e.g., rename)

    Returns:
        Tuple of (success: bool, file_path: str, sequence: str, ligands: List[Dict], metadata: dict)
    """
    if operations is None:
        operations = []

    # None means keep the source format as-is (no conversion)
    target_format = convert if convert is not None else source_format

    try:
        with open(source_path, 'r') as f:
            content = f.read()

        # Convert format if needed
        actual_format = source_format
        if source_format != target_format:
            if source_format == "cif" and target_format == "pdb":
                print(f"  Converting CIF to PDB format...")
                content = convert_cif_to_pdb(content)
                actual_format = "pdb"
            elif source_format == "pdb" and target_format == "cif":
                print(f"  Converting PDB to CIF format...")
                content = convert_pdb_to_cif(content)
                actual_format = "cif"

        if remove_waters:
            content = remove_waters_from_content(content, actual_format)

        # Filter to specific chain if requested
        if _chain_filters_structure(chain):
            content = filter_chain_from_content(content, actual_format, chain)

        # Extract ligands BEFORE applying rename operations (to get original CCD codes for SMILES lookup)
        original_ligand_codes = extract_ligands_from_structure(content, actual_format)
        rename_mapping = build_rename_mapping(operations) if operations else {}

        # Apply operations (e.g., rename, remove)
        if operations:
            content = apply_operations(content, actual_format, operations, structure_id=custom_id)

        extension = ".pdb" if actual_format == "pdb" else ".cif"
        output_path = os.path.join(output_folder, f"{custom_id}{extension}")

        with open(output_path, 'w') as f:
            f.write(content)

        file_size = os.path.getsize(output_path)
        sequence = extract_sequence_from_structure(content, actual_format, chain=chain)

        # Build ligands list using original codes for SMILES lookup, renamed codes for output
        ligands = []
        if original_ligand_codes:
            print(f"  Found {len(original_ligand_codes)} ligand(s) in structure: {', '.join(original_ligand_codes)}")
            for original_code in original_ligand_codes:
                # Use original code for SMILES fetch from RCSB
                smiles = fetch_ligand_smiles_from_rcsb(original_code) if fetch_compounds else None
                # Use renamed code (if any) for output
                output_code = rename_mapping.get(original_code, original_code)
                ligands.append({
                    'id': f"{custom_id}_{output_code}",
                    'code': output_code,
                    'format': 'smiles' if smiles else '',
                    'smiles': smiles if smiles else '',
                    'ccd': original_code  # Keep original CCD code for reference
                })

        source_info = f"local ({source_format.upper()})"
        if source_format != actual_format:
            source_info += f" -> converted to {actual_format.upper()}"

        metadata = {
            "file_size": file_size,
            "source": source_info,
            "source_path": source_path,
            "actual_format": actual_format
        }

        print(f"Successfully processed {pdb_id} as {custom_id}: {file_size} bytes ({source_info})")
        return True, output_path, sequence, ligands, metadata

    except Exception as e:
        error_msg = f"Error processing local file {pdb_id}: {str(e)}"
        print(f"Error: {error_msg}")
        metadata = {
            "error_message": error_msg,
            "source": "local_copy_failed",
            "attempted_path": source_path
        }
        return False, "", "", [], metadata


def download_from_rcsb(pdb_id: str, custom_id: str, convert: Optional[str], biological_assembly: bool,
                   remove_waters: bool, output_folder: str, repo_pdbs_folder: str,
                   operations: List[Dict[str, Any]] = None,
                   chain: str = "longest",
                   fetch_compounds: bool = True) -> Tuple[bool, str, str, List[Dict[str, str]], Dict[str, Any]]:
    """
    Download a single structure from RCSB PDB and save to both pdbs/ and output folder.

    If convert is None, tries PDB first and falls back to CIF if PDB is unavailable,
    keeping whatever format downloads without converting.
    If convert is "pdb", downloads PDB (falls back to CIF and converts if PDB unavailable).
    If convert is "cif", downloads CIF directly.

    Args:
        pdb_id: PDB identifier (4 characters)
        custom_id: Custom ID for renaming the structure
        convert: Target format ("pdb", "cif", or None to keep whatever is found on RCSB)
        biological_assembly: Whether to download biological assembly
        remove_waters: Whether to remove water molecules
        output_folder: Directory to save the structure
        repo_pdbs_folder: Repository PDBs folder for caching
        operations: List of operations to apply (e.g., rename)

    Returns:
        Tuple of (success: bool, file_path: str, sequence: str, ligands: List[Dict], metadata: dict)
    """
    if operations is None:
        operations = []

    pdb_id = pdb_id.upper()

    # Validate RCSB PDB ID format before attempting download
    if len(pdb_id) != 4 or not pdb_id.isalnum():
        error_msg = f"Invalid RCSB PDB ID format: {pdb_id}. RCSB IDs must be 4 alphanumeric characters (e.g., '4UFC'). File not found locally."
        print(f"Error: {error_msg}")
        metadata = {
            "error_message": error_msg,
            "source": "rcsb_invalid_id",
            "attempted_path": "N/A"
        }
        return False, "", "", [], metadata

    def attempt_download(download_fmt: str):
        """Helper to attempt download in specified format."""
        if download_fmt == "pdb":
            if biological_assembly:
                url = f"https://files.rcsb.org/download/{pdb_id}.pdb1.gz"
            else:
                url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        else:  # cif
            if biological_assembly:
                url = f"https://files.rcsb.org/download/{pdb_id}-assembly1.cif.gz"
            else:
                url = f"https://files.rcsb.org/download/{pdb_id}.cif"
        return url

    def do_download(url: str) -> str:
        """Execute HTTP download and return content string."""
        try:
            import requests
        except ImportError:
            raise RuntimeError("Cannot download from RCSB: 'requests' module not available. Please install with: pip install requests")
        headers = {'User-Agent': 'BioPipelines-PDB/1.0 (https://github.com/locbp-uzh/biopipelines)'}
        response = requests.get(url, headers=headers, timeout=30)
        response.raise_for_status()
        if url.endswith('.gz'):
            import gzip
            return gzip.decompress(response.content).decode('utf-8')
        return response.text

    def process_and_save(downloaded_content: str, download_fmt: str, actual_fmt: str) -> Tuple[bool, str, str, list, dict]:
        """Process downloaded content and save to output folder."""
        content = downloaded_content

        # Convert format if needed
        needs_conversion = (download_fmt != actual_fmt)
        if needs_conversion:
            if download_fmt == "cif" and actual_fmt == "pdb":
                print(f"  Converting CIF to PDB format...")
                content = convert_cif_to_pdb(content)
            elif download_fmt == "pdb" and actual_fmt == "cif":
                print(f"  Converting PDB to CIF format...")
                content = convert_pdb_to_cif(content)

        if remove_waters:
            content = remove_waters_from_content(content, actual_fmt)

        if _chain_filters_structure(chain):
            content = filter_chain_from_content(content, actual_fmt, chain)

        # Validate file content
        if actual_fmt == "pdb":
            if not (content.startswith("HEADER") or content.startswith("ATOM") or content.startswith("MODEL") or content.startswith("REMARK")):
                raise ValueError(f"Downloaded file does not appear to be valid PDB format")
        else:  # cif
            if not ("data_" in content or "_entry.id" in content):
                raise ValueError(f"Downloaded file does not appear to be valid CIF format")

        original_ligand_codes = extract_ligands_from_structure(content, actual_fmt)
        rename_mapping = build_rename_mapping(operations) if operations else {}

        if operations:
            content = apply_operations(content, actual_fmt, operations)

        output_extension = ".pdb" if actual_fmt == "pdb" else ".cif"
        output_path = os.path.join(output_folder, f"{custom_id}{output_extension}")
        with open(output_path, 'w') as f:
            f.write(content)

        file_size = os.path.getsize(output_path)
        sequence = extract_sequence_from_structure(content, actual_fmt, chain=chain)

        ligands = []
        if original_ligand_codes:
            print(f"  Found {len(original_ligand_codes)} ligand(s) in structure: {', '.join(original_ligand_codes)}")
            for original_code in original_ligand_codes:
                smiles = fetch_ligand_smiles_from_rcsb(original_code) if fetch_compounds else None
                output_code = rename_mapping.get(original_code, original_code)
                ligands.append({
                    'id': f"{custom_id}_{output_code}",
                    'code': output_code,
                    'format': 'smiles' if smiles else '',
                    'smiles': smiles if smiles else '',
                    'ccd': original_code
                })

        source_info = f"rcsb_download ({download_fmt.upper()})"
        if needs_conversion:
            source_info += f" -> converted to {actual_fmt.upper()}"

        metadata = {
            "file_size": file_size,
            "source": source_info,
            "url": url,
            "actual_format": actual_fmt
        }

        print(f"Successfully downloaded {pdb_id} as {custom_id}: {file_size} bytes ({source_info})")
        return True, output_path, sequence, ligands, metadata

    # Determine download order based on convert parameter
    if convert is None:
        # No conversion requested: try PDB first, fall back to CIF, keep whatever downloads
        formats_to_try = [("pdb", "pdb"), ("cif", "cif")]
    elif convert == "pdb":
        # Want PDB: try PDB directly, fall back to CIF then convert
        formats_to_try = [("pdb", "pdb"), ("cif", "pdb")]
    else:  # convert == "cif"
        # Want CIF: download CIF directly
        formats_to_try = [("cif", "cif")]

    last_error = None
    for download_fmt, actual_fmt in formats_to_try:
        url = attempt_download(download_fmt)
        try:
            print(f"Downloading {pdb_id} from RCSB: {url}")
            downloaded_content = do_download(url)

            # Cache to pdbs/ folder (in download format, before conversion)
            os.makedirs(repo_pdbs_folder, exist_ok=True)
            cache_ext = ".pdb" if download_fmt == "pdb" else ".cif"
            cache_path = os.path.join(repo_pdbs_folder, f"{pdb_id}{cache_ext}")
            with open(cache_path, 'w') as f:
                f.write(downloaded_content)
            print(f"Cached to pdbs/ folder: {cache_path}")

            return process_and_save(downloaded_content, download_fmt, actual_fmt)

        except Exception as e:
            last_error = e
            if download_fmt == "pdb" and len(formats_to_try) > 1:
                print(f"  PDB download failed ({str(e)}), trying CIF as fallback...")
            # else fall through to next format or fail

    # All attempts failed
    error_type = type(last_error).__name__
    http_status = 'unknown'
    if 'requests' in sys.modules and hasattr(last_error, 'response') and last_error.response:
        http_status = getattr(last_error.response, 'status_code', 'unknown')

    if 'RequestException' in error_type or 'HTTPError' in error_type:
        error_msg = f"HTTP error downloading {pdb_id}: {str(last_error)}"
        source = f"rcsb_download_failed_{http_status}"
    else:
        error_msg = f"Unexpected error downloading {pdb_id}: {str(last_error)}"
        source = "rcsb_processing_error"

    print(f"Error: {error_msg}")
    metadata = {
        "error_message": error_msg,
        "source": source,
        "attempted_path": url
    }
    return False, "", "", [], metadata


def resolve_upstream_file(pdb_id: str, upstream_files: List[str],
                         files_contain_wildcards: bool = False) -> Optional[str]:
    """
    Resolve the source file for a structure from upstream tool output.

    Args:
        pdb_id: The structure ID to resolve
        upstream_files: List of file paths from the upstream tool
        files_contain_wildcards: Deprecated, ignored. Wildcards detected from '*' in paths.

    Returns:
        Resolved file path, or None if not found
    """
    import glob as glob_module

    if not upstream_files:
        return None

    # Detect wildcards by checking for '*' in file paths
    has_wildcards = any('*' in f for f in upstream_files)
    if has_wildcards:
        # Try each pattern
        for pattern in upstream_files:
            expanded = glob_module.glob(pattern)
            if not expanded:
                continue
            # Try to match by ID
            for fp in expanded:
                basename = os.path.splitext(os.path.basename(fp))[0]
                if basename == pdb_id or basename.startswith(f"{pdb_id}_") or basename.startswith(f"{pdb_id}-"):
                    return fp
            # If single pattern and single expansion, use it
            if len(upstream_files) == 1 and len(expanded) == 1:
                return expanded[0]
        return None

    # Direct file list - match by index or by name
    if len(upstream_files) == 1:
        # Single file for all IDs
        if os.path.exists(upstream_files[0]):
            return upstream_files[0]
        return None

    # Try to match by name
    for fp in upstream_files:
        basename = os.path.splitext(os.path.basename(fp))[0]
        if basename == pdb_id:
            return fp

    return None


def fetch_structures(config_data: Dict[str, Any]) -> int:
    """
    Fetch multiple structures with priority-based lookup: local_folder -> pdbs/ -> RCSB download.
    Also handles upstream tool outputs where files already exist at known paths.

    Args:
        config_data: Configuration dictionary with fetch parameters

    Returns:
        Number of failed fetches
    """
    raw_pdb_ids = config_data['pdb_ids']
    pdb_ids = expand_ids(raw_pdb_ids) if any(contains_pattern(s) for s in raw_pdb_ids) else raw_pdb_ids
    raw_custom_ids = config_data.get('custom_ids', raw_pdb_ids)
    custom_ids = expand_ids(raw_custom_ids) if any(contains_pattern(s) for s in raw_custom_ids) else raw_custom_ids
    convert = config_data.get('convert')  # May be None (keep whatever format is found)
    local_folder = config_data.get('local_folder')
    repo_pdbs_folder = config_data['repo_pdbs_folder']
    biological_assembly = config_data.get('biological_assembly', False)
    remove_waters = config_data.get('remove_waters', True)
    output_folder = config_data['output_folder']
    structures_table = config_data['structures_table']
    sequences_table = config_data['sequences_table']
    failed_table = config_data['failed_table']
    missing_table = config_data.get('missing_table', os.path.join(output_folder, 'missing_structures.csv'))
    compounds_table = config_data.get('compounds_table', os.path.join(output_folder, 'compounds.csv'))
    fetch_compounds = config_data.get('fetch_compounds', True)
    operations = config_data.get('operations', [])
    chain = config_data.get('chain', 'auto')
    split_chains = bool(config_data.get('split_chains', False))
    from_upstream = config_data.get('from_upstream', False)
    raw_upstream_files = config_data.get('upstream_files', [])
    if raw_upstream_files and len(raw_upstream_files) == 1 and '<id>' in raw_upstream_files[0]:
        upstream_files = [expand_file_pattern(raw_upstream_files[0], eid) for eid in pdb_ids]
    elif raw_upstream_files and any(contains_pattern(s) for s in raw_upstream_files):
        upstream_files = expand_ids(raw_upstream_files)
    else:
        upstream_files = raw_upstream_files
    upstream_wildcards = config_data.get('upstream_files_contain_wildcards', False)  # Legacy, ignored

    # Resolve upstream ids to real files via the map table's file_path column (table-only stages
    # declare an <id>.pdb template that points at their own empty dir; the truth is in the map).
    upstream_id_to_file = {}
    _map = config_data.get('upstream_map_table')
    if from_upstream and _map and os.path.exists(_map):
        try:
            import pandas as _pd
            _df = _pd.read_csv(_map)
            if 'id' in _df.columns and 'file_path' in _df.columns:
                upstream_id_to_file = {
                    str(r['id']): str(r['file_path'])
                    for _, r in _df.iterrows()
                    if str(r['file_path']) and str(r['file_path']) != 'nan'
                }
        except Exception as e:
            print(f"  Warning: could not read upstream map table ({e}); falling back to file-list resolution")

    # Ids an upstream filter already dropped — their files are legitimately
    # absent and must NOT count as fetch failures (would otherwise hard-exit).
    # PDB owns tables/missing.csv: rows are re-keyed to THIS tool's id
    # (custom_id) so a rename via ids= still excuses the renamed output.
    from biopipelines.biopipelines_io import read_upstream_missing
    from biopipelines.id_map_utils import get_mapped_ids
    upstream_missing_rows = read_upstream_missing(config_data.get('upstream_missing'))
    upstream_missing_by_id = {str(r.get('id', '')).strip(): r for r in upstream_missing_rows}
    upstream_missing_ids = set(upstream_missing_by_id)
    missing_out_path = config_data.get('missing_out')
    propagated_missing = []  # rows re-keyed to custom_id, written to missing_out
    # Upstream missing.csv is keyed in the UPSTREAM id space (e.g. a ligand id
    # `nilotinib`), while this tool's ids are the combinatorial product
    # (`2HYY+nilotinib`). Map each output id to its upstream missing component
    # via provenance/suffix matching (same matcher the completion check uses),
    # so a product id is excused when any of its components was dropped upstream.
    excused_to_upstream = {}  # custom_id/pdb_id -> upstream missing id
    if upstream_missing_ids:
        print(f"Excusing {len(upstream_missing_ids)} ids dropped upstream: {sorted(upstream_missing_ids)}")
        candidate_ids = sorted(set(pdb_ids) | set(custom_ids))
        mapped = get_mapped_ids(candidate_ids, sorted(upstream_missing_ids), unique=True)
        for out_id, up_id in mapped.items():
            if up_id:
                excused_to_upstream[out_id] = up_id

    if from_upstream:
        print(f"Processing {len(pdb_ids)} structures from upstream tool")
    else:
        convert_display = f"convert to {convert.upper()}" if convert else "keep as-is (pdb|cif)"
        print(f"Fetching {len(pdb_ids)} structures ({convert_display})")
        print(f"Priority: {'local_folder -> ' if local_folder else ''}pdbs/ -> RCSB download")
    if biological_assembly:
        print("Including biological assemblies")
    if remove_waters:
        print("Water molecules will be removed")
    if operations:
        op_summaries = []
        for op in operations:
            if op.get("op") == "rename":
                op_summaries.append(f"Rename({op['old']} → {op['new']})")
            else:
                op_summaries.append(op.get("op", "unknown"))
        print(f"Operations: {', '.join(op_summaries)}")

    # Create output directory
    os.makedirs(output_folder, exist_ok=True)

    # Track results
    successful_downloads = []
    successful_sequences = []
    all_ligands = []
    failed_downloads = []
    missing_structures = []

    # Fetch each structure
    for i, (pdb_id, custom_id) in enumerate(zip(pdb_ids, custom_ids), 1):
        print(f"\n[{i}/{len(pdb_ids)}] Processing {pdb_id} -> {custom_id}")

        # Upstream missing.csv names the UPSTREAM stream id (== pdb_id when
        # from_upstream); custom_id may be a rename via ids=. Match either so a
        # renamed input is still excused, not re-failed, and re-key the row to
        # custom_id so completion excuses THIS tool's (possibly renamed) output.
        matched_key = pdb_id if pdb_id in upstream_missing_ids else (
            custom_id if custom_id in upstream_missing_ids else
            excused_to_upstream.get(pdb_id) or excused_to_upstream.get(custom_id))
        if matched_key is not None:
            src = upstream_missing_by_id[matched_key]
            propagated_missing.append({
                "id": custom_id,
                "removed_by": src.get("removed_by", ""),
                "kind": src.get("kind", "filter"),
                "cause": src.get("cause", ""),
            })
            print(f"  Skipping {pdb_id} -> {custom_id}: dropped by an upstream filter (excused)")
            continue

        if from_upstream:
            source_path = upstream_id_to_file.get(pdb_id)
            if not source_path:
                if len(upstream_files) == len(pdb_ids):
                    source_path = upstream_files[i - 1]
                else:
                    source_path = resolve_upstream_file(pdb_id, upstream_files, upstream_wildcards)

            if source_path and os.path.exists(source_path):
                # Detect source format from extension
                source_format = "cif" if source_path.endswith(".cif") else "pdb"
                success, file_path, sequence, ligands, metadata = copy_local_structure(
                    pdb_id, custom_id, source_path, source_format, convert, remove_waters, output_folder, operations, chain=chain, fetch_compounds=fetch_compounds
                )
            else:
                error_msg = f"Upstream file not found for '{pdb_id}': {source_path}"
                print(f"Error: {error_msg}")
                success = False
                file_path = ""
                sequence = ""
                ligands = []
                metadata = {
                    "error_message": error_msg,
                    "source": "upstream_not_found",
                    "attempted_path": str(source_path)
                }
        else:
            # Try to find locally first (returns tuple of (path, source_format) or None)
            local_result = find_local_structure(pdb_id, convert, local_folder, repo_pdbs_folder)

            if local_result:
                # Copy from local (may need conversion)
                local_path, source_format = local_result
                success, file_path, sequence, ligands, metadata = copy_local_structure(
                    pdb_id, custom_id, local_path, source_format, convert, remove_waters, output_folder, operations, chain=chain, fetch_compounds=fetch_compounds
                )
            else:
                # Download from RCSB (tries PDB first if convert=None, falls back to CIF)
                print(f"{pdb_id} not found locally, downloading from RCSB")
                success, file_path, sequence, ligands, metadata = download_from_rcsb(
                    pdb_id, custom_id, convert, biological_assembly, remove_waters,
                    output_folder, repo_pdbs_folder, operations, chain=chain, fetch_compounds=fetch_compounds
                )

        if success:
            actual_format = metadata.get('actual_format', convert or 'unknown')
            # When chain="all" + split_chains=True, replace the single
            # parent file with one file per chain letter present in the
            # parsed structure. The parent file is removed so the
            # structures stream cleanly maps each id to exactly one file.
            if split_chains and _chain_is_multi(chain) and actual_format in ("pdb", "cif"):
                with open(file_path, 'r') as fh:
                    parent_content = fh.read()
                atoms = parse_pdb_file(file_path) if actual_format == "pdb" else None
                if atoms is not None:
                    chains_in_file = sorted({a.chain for a in atoms if a.chain})
                else:
                    chains_in_file = []  # CIF: skip splitting (parser limitation)
                # Restrict to the explicitly-requested chain list when given.
                if isinstance(chain, list):
                    chain_letters = [c for c in chain if c in chains_in_file]
                else:
                    chain_letters = chains_in_file
                if not chain_letters:
                    print(f"  Warning: split_chains=True but no chains detected in {file_path}; keeping parent file")
                    successful_downloads.append({
                        'id': custom_id,
                        'pdb_id': pdb_id,
                        'file_path': file_path,
                        'format': actual_format,
                        'file_size': metadata['file_size'],
                        'source': metadata['source'],
                    })
                else:
                    extension = ".pdb" if actual_format == "pdb" else ".cif"
                    for ch_letter in chain_letters:
                        chain_content = filter_chain_from_content(parent_content, actual_format, ch_letter)
                        chain_path = os.path.join(output_folder, f"{custom_id}_{ch_letter}{extension}")
                        with open(chain_path, 'w') as fh:
                            fh.write(chain_content)
                        successful_downloads.append({
                            'id': f"{custom_id}_{ch_letter}",
                            'pdb_id': pdb_id,
                            'file_path': chain_path,
                            'format': actual_format,
                            'file_size': os.path.getsize(chain_path),
                            'source': metadata['source'] + f" (split chain {ch_letter})",
                        })
                    try:
                        os.remove(file_path)
                    except OSError:
                        pass
            else:
                successful_downloads.append({
                    'id': custom_id,
                    'pdb_id': pdb_id,
                    'file_path': file_path,
                    'format': actual_format,
                    'file_size': metadata['file_size'],
                    'source': metadata['source'],
                })
            # Sequence emission honours the chain parameter:
            #   chain == "auto"        -> single <custom_id>,<longest_seq> row
            #   chain == "all"         -> one <custom_id>_<chain_letter>,<seq>
            #                              row per chain in the structure
            #                              (no aggregate longest row)
            #   chain == ["A","C",...] -> one <custom_id>_<chain_letter>,<seq>
            #                              row per requested chain
            #   chain == "A"/"B"/...   -> single <custom_id>,<seq_for_that_chain>
            #                              row, queried from RCSB FASTA when
            #                              possible, else from the parsed
            #                              structure.
            sequence_rows = []
            if not from_upstream and _is_rcsb_pdb_code(pdb_id):
                try:
                    chain_pairs = fetch_rcsb_chain_sequences(pdb_id)
                except Exception as e:
                    print(f"  Warning: RCSB FASTA fetch failed for {pdb_id} ({e}); falling back to structure-extracted sequence")
                    chain_pairs = []

                if chain_pairs:
                    # Deduplicate: RCSB FASTA can repeat the same sequence
                    # for each chain letter listed in an entity record.
                    unique_letters = []
                    seen = set()
                    for ch_letter, seq in chain_pairs:
                        if ch_letter in seen:
                            continue
                        seen.add(ch_letter)
                        unique_letters.append((ch_letter, seq))

                    if isinstance(chain, list):
                        wanted = set(chain)
                        emitted = set()
                        for ch_letter, seq in unique_letters:
                            if ch_letter in wanted:
                                sequence_rows.append({
                                    'id': f"{custom_id}_{ch_letter}",
                                    'sequence': seq
                                })
                                emitted.add(ch_letter)
                        missing_letters = wanted - emitted
                        if missing_letters:
                            available = ', '.join(sorted({ch for ch, _ in unique_letters}))
                            print(f"  Warning: requested chain(s) {sorted(missing_letters)} not present in RCSB FASTA for {pdb_id}. Available: {available}")
                    elif chain == "all":
                        for ch_letter, seq in unique_letters:
                            sequence_rows.append({
                                'id': f"{custom_id}_{ch_letter}",
                                'sequence': seq
                            })
                    elif chain == "auto":
                        longest_seq = max((seq for _, seq in unique_letters), key=len)
                        sequence_rows.append({
                            'id': custom_id,
                            'sequence': longest_seq
                        })
                    else:
                        match = next(((ch, seq) for ch, seq in unique_letters if ch == chain), None)
                        if match is not None:
                            _, seq = match
                            sequence_rows.append({
                                'id': custom_id,
                                'sequence': seq
                            })
                        elif sequence:
                            sequence_rows.append({'id': custom_id, 'sequence': sequence})

            # Local file or upstream input — fall back to the structure-extracted sequence.
            if not sequence_rows and sequence:
                sequence_rows.append({'id': custom_id, 'sequence': sequence})

            successful_sequences.extend(sequence_rows)
            if ligands:  # Add ligands to the collection
                all_ligands.extend(ligands)
        else:
            failure_entry = {
                'pdb_id': pdb_id,
                'error_message': metadata['error_message'],
                'source': metadata['source'],
                'attempted_path': metadata['attempted_path']
            }
            failed_downloads.append(failure_entry)
            missing_structures.append(failure_entry)
    
    # Save successful downloads table
    if successful_downloads:
        df_success = pd.DataFrame(successful_downloads)
        df_success.to_csv(structures_table, index=False)
        print(f"\nSuccessful fetches saved: {structures_table} ({len(successful_downloads)} structures)")
    else:
        # Create empty table with proper columns
        empty_df = pd.DataFrame(columns=["id", "pdb_id", "file_path", "format", "file_size", "source"])
        empty_df.to_csv(structures_table, index=False)
        print(f"No successful fetches - created empty table: {structures_table}")

    # Save sequences table
    if successful_sequences:
        df_sequences = pd.DataFrame(successful_sequences)
        df_sequences.to_csv(sequences_table, index=False)
        print(f"Sequences saved: {sequences_table} ({len(successful_sequences)} sequences)")
    else:
        # Create empty sequences table with proper columns
        empty_seq_df = pd.DataFrame(columns=["id", "sequence"])
        empty_seq_df.to_csv(sequences_table, index=False)
        print(f"No sequences extracted - created empty table: {sequences_table}")

    # Save compounds table
    if all_ligands:
        df_compounds = pd.DataFrame(all_ligands)
        df_compounds.to_csv(compounds_table, index=False)
        print(f"Compounds saved: {compounds_table} ({len(all_ligands)} ligands)")
    else:
        # Create empty compounds table with proper columns
        empty_compounds_df = pd.DataFrame(columns=["id", "code", "format", "smiles", "ccd"])
        empty_compounds_df.to_csv(compounds_table, index=False)
        print(f"No ligands found - created empty table: {compounds_table}")

    # Save failed downloads table (always create, even if empty)
    if failed_downloads:
        df_failed = pd.DataFrame(failed_downloads)
        df_failed.to_csv(failed_table, index=False)
        print(f"Failed fetches saved: {failed_table} ({len(failed_downloads)} failures)")
    else:
        # Create empty failed downloads table
        empty_failed_df = pd.DataFrame(columns=["pdb_id", "error_message", "source", "attempted_path"])
        empty_failed_df.to_csv(failed_table, index=False)
        print("No failed fetches")

    # Save missing structures table (always create, even if empty)
    if missing_table:
        if missing_structures:
            df_missing = pd.DataFrame(missing_structures)
            df_missing.to_csv(missing_table, index=False)
            print(f"Missing structures saved: {missing_table} ({len(missing_structures)} missing)")
        else:
            empty_missing_df = pd.DataFrame(columns=["pdb_id", "error_message", "source", "attempted_path"])
            empty_missing_df.to_csv(missing_table, index=False)

    # Own the declared `missing` table (id|removed_by|kind|cause), re-keyed to
    # this tool's output ids, so completion excuses the (possibly renamed)
    # upstream-filtered outputs. Always written when declared.
    if missing_out_path:
        os.makedirs(os.path.dirname(missing_out_path), exist_ok=True)
        pd.DataFrame(propagated_missing,
                     columns=["id", "removed_by", "kind", "cause"]).to_csv(missing_out_path, index=False)
        print(f"Propagated missing.csv: {missing_out_path} ({len(propagated_missing)} excused)")

    # Summary
    print(f"\n=== FETCH SUMMARY ===")
    print(f"Requested: {len(pdb_ids)} structures")
    print(f"Successful: {len(successful_downloads)}")
    print(f"Failed: {len(failed_downloads)}")
    print(f"Success rate: {len(successful_downloads)/len(pdb_ids)*100:.1f}%")

    if successful_downloads:
        total_size = sum(item['file_size'] for item in successful_downloads)
        print(f"Total downloaded: {total_size:,} bytes ({total_size/1024/1024:.2f} MB)")

    print(f"Sequences extracted: {len(successful_sequences)}")
    print(f"Ligands extracted: {len(all_ligands)}")
    
    # Log any failed fetches
    if failed_downloads:
        print(f"\nFailed fetches:")
        for failure in failed_downloads:
            print(f"  - {failure['pdb_id']}: {failure['error_message']}")

    # Return the number of failures
    return len(failed_downloads)


def main():
    parser = argparse.ArgumentParser(description='Fetch protein structures with priority-based lookup')
    parser.add_argument('--config', required=True, help='JSON config file with fetch parameters')

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
    required_params = ['pdb_ids', 'repo_pdbs_folder', 'output_folder', 'structures_table', 'sequences_table', 'failed_table']
    for param in required_params:
        if param not in config_data:
            print(f"Error: Missing required parameter: {param}")
            sys.exit(1)

    try:
        failed_count = fetch_structures(config_data)

        # Don't hard-exit on partial failures: the tables (structures, failed,
        # missing) are written, and the completion check adjudicates — it flags
        # FAILED only for ids that are missing AND not excused by missing.csv.
        if failed_count > 0:
            print(f"\n{failed_count} structure(s) not produced — recorded in failed/missing tables")
        else:
            print("\nAll structures fetched successfully")

    except Exception as e:
        print(f"Error fetching structures: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()