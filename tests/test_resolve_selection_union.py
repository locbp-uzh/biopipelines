"""Tests for the '+' union branch of resolve_selection (atom-set references).

resolve_selection treats a selection containing BOTH '+' and '.' as a union of
independently-resolved terms (e.g. 'LIG.O3+LIG.Cl'). Each term is resolved on
its own, so chain/ligand context does NOT distribute across the '+'. These tests
pin that behavior, including the documented limitation that a bare residue
number after a '+' resolves globally (any chain)."""

import pytest
from biopipelines.pdb_parser import Atom, resolve_selection


def _atoms():
    # Two chains, a ligand copy on each chain.
    return [
        Atom(0.0, 0.0, 0.0, "CA", "ALA", 10, "A", "C"),
        Atom(1.0, 0.0, 0.0, "CB", "ALA", 10, "A", "C"),
        Atom(0.0, 1.0, 0.0, "CA", "ALA", 11, "A", "C"),
        Atom(0.0, 0.0, 1.0, "CA", "GLY", 11, "B", "C"),
        Atom(2.0, 0.0, 0.0, "O3", "LIG", 200, "A", "O"),
        Atom(3.0, 0.0, 0.0, "CL", "LIG", 200, "A", "CL"),
    ]


def test_union_of_ligand_atoms():
    """'LIG.O3+LIG.Cl' selects both named ligand atoms (the documented form)."""
    atoms = _atoms()
    out = resolve_selection("LIG.O3+LIG.CL", atoms)
    names = sorted(a.atom_name for a in out)
    assert names == ["CL", "O3"]


def test_union_dedupes_overlapping_terms():
    """Overlapping terms don't double-count the same atom."""
    atoms = _atoms()
    out = resolve_selection("LIG.O3+LIG.O3", atoms)
    assert len(out) == 1


def test_pure_number_plus_is_not_union():
    """'10+11' has no '.', so it stays on the residue-number branch (both residues,
    every chain) rather than the atom-set union branch."""
    atoms = _atoms()
    out = resolve_selection("10+11", atoms)
    resnums = sorted({a.res_num for a in out})
    assert resnums == [10, 11]


def test_chain_context_does_not_distribute_across_plus():
    """DOCUMENTED LIMITATION: in 'A10.CA+11', the chain qualifier on the first
    term does NOT carry to the second. The bare '11' resolves globally, so the
    chain-B residue 11 is included even though the user may have meant A11."""
    atoms = _atoms()
    out = resolve_selection("A10.CA+11", atoms)
    chains_resnums = sorted({(a.chain, a.res_num) for a in out})
    # A10.CA -> (A,10); 11 -> every chain's residue 11 = (A,11) and (B,11)
    assert ("A", 10) in chains_resnums
    assert ("A", 11) in chains_resnums
    assert ("B", 11) in chains_resnums  # global, not restricted to chain A
