"""Tests for plane_side_selection (and the chainless remap in its dispatch):

A plane passes through an anchor; its normal points normal_from -> normal_to.
keep="positive" keeps residues whose CA is toward normal_to. Two regressions
covered: (1) one_atom must raise on a multi-atom ref (plane must not depend on
atom order); (2) chainless residue inputs (("", n)) must be remapped to the
structure's real chain or they silently drop.
"""

import os
import sys
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "pipe_scripts"))
from pipe_selection import plane_side_selection, remap_chainless_residues  # noqa: E402


def _atom(serial, name, resname, chain, resnum, x, y, z, hetatm=False):
    rec = "HETATM" if hetatm else "ATOM  "
    elem = name.strip()[0]
    return (f"{rec}{serial:>5} {name:<4} {resname:<3} {chain}{resnum:>4}    "
            f"{x:>8.3f}{y:>8.3f}{z:>8.3f}  1.00  0.00          {elem:>2}\n")


@pytest.fixture
def pdb(tmp_path):
    """Chain A, residues 10-13 with CA along x (-2,-1,+1,+2); two marker atoms
    FROM/TO on chain B define a plane normal along +x through the origin; a 2-atom
    ligand LIG on chain B for the multi-atom-ref test."""
    lines = [
        # protein CA: residues split across x=0 plane
        _atom(1, "CA", "ALA", "A", 10, -2.0, 0.0, 0.0),
        _atom(2, "CA", "ALA", "A", 11, -1.0, 0.0, 0.0),
        _atom(3, "CA", "ALA", "A", 12, 1.0, 0.0, 0.0),
        _atom(4, "CA", "ALA", "A", 13, 2.0, 0.0, 0.0),
        # plane markers: from (-1,0,0) -> to (1,0,0); normal = +x; midpoint = origin
        _atom(5, "C1", "MRK", "B", 1, -1.0, 0.0, 0.0, hetatm=True),
        _atom(6, "C2", "MRK", "B", 2, 1.0, 0.0, 0.0, hetatm=True),
        # res 12 also gets a CB so the bare ref "12" resolves to 2 atoms
        _atom(7, "CB", "ALA", "A", 12, 1.5, 0.5, 0.0),
    ]
    p = tmp_path / "t.pdb"
    p.write_text("".join(lines))
    return str(p)


PLANE = {"normal_from": "B1.C1", "normal_to": "B2.C2",
         "anchor": "midpoint", "keep": "positive"}


def test_keeps_positive_side(pdb):
    # CA at x=+1,+2 (res 12,13) are toward normal_to; x=-1,-2 (res 10,11) are not.
    kept = plane_side_selection({("A", 10), ("A", 11), ("A", 12), ("A", 13)}, PLANE, pdb)
    assert kept == {("A", 12), ("A", 13)}


def test_keep_negative_side(pdb):
    plane = dict(PLANE, keep="negative")
    kept = plane_side_selection({("A", 10), ("A", 11), ("A", 12), ("A", 13)}, plane, pdb)
    assert kept == {("A", 10), ("A", 11)}


def test_multi_atom_ref_raises(pdb):
    # bare ref "12" resolves to res-12's CA+CB (2 atoms) -> the plane would depend
    # on atom order, so it must raise instead of silently taking sel[0].
    plane = dict(PLANE, normal_from="12", normal_to="B2.C2")
    with pytest.raises(ValueError, match="expected exactly one|selected 2 atoms"):
        plane_side_selection({("A", 13)}, plane, pdb)


def test_no_atom_ref_raises(pdb):
    plane = dict(PLANE, normal_from="B9.ZZ")
    with pytest.raises(ValueError, match="selected no atoms"):
        plane_side_selection({("A", 12)}, plane, pdb)


def test_chainless_remap_then_plane(pdb):
    # Chainless input ("", n) must remap to chain A before the CA lookup; without
    # the remap the lookup ca[("", n)] misses and everything drops.
    valid = {("A", 10), ("A", 11), ("A", 12), ("A", 13)}
    chainless = {("", 10), ("", 11), ("", 12), ("", 13)}
    remapped = remap_chainless_residues(chainless, valid)
    assert remapped == valid
    kept = plane_side_selection(remapped, PLANE, pdb)
    assert kept == {("A", 12), ("A", 13)}


def test_chainless_without_remap_drops(pdb):
    # Documents the bug the dispatch fix prevents: chainless keys don't match the
    # ("A", n)-keyed CA dict, so plane_side_selection alone keeps nothing.
    kept = plane_side_selection({("", 12), ("", 13)}, PLANE, pdb)
    assert kept == set()
