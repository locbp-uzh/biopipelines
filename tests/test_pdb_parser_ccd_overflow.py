"""A 5-character CCD code must not cost the parser its atoms.

The residue field is 3 wide, but the PDB issues 5-character codes (A1EI4 and
similar). Such a record shifts every later field right by 2, which 164ba7c set
out to handle by falling back to ``split()`` past column 17. Tokenizing works
only when every field is present and separated:

  * a **blank chain id** -- routine in this repo's own tool-generated ligand
    PDBs -- shifts the token list by one, so field_res_seq returned the x
    coordinate, int() raised, and parse_pdb_file dropped the atom;
  * an **insertion code** glues onto resSeq ("201A"), so int() raised;
  * two **adjacent full-width coordinates** ("-100.123-100.123") are one token.

All three silently discarded atoms, which is the same class of failure the
commit was fixing -- and there was no test for the fix at all. Parsing by
shifted column instead handles every case, and the shift is derived from the
code's width rather than assumed.

The reference layout here is a real RCSB record for 9RTM's tetramethylrhodamine,
the ligand that motivated the original commit, not a hand-built line.
"""

import os
import tempfile

import pytest

from biopipelines.pdb_parser import (
    _residue_shift, field_chain, field_coords, field_element, field_res_name,
    field_res_seq, parse_pdb_file,
)

# Real record, RCSB 9RTM. Columns measured, not assumed: the 5-char code occupies
# 17:22 and every later field sits exactly 2 to the right of canonical.
REAL_5CHAR = "HETATM 2422  C1  A1EI4 B 201      21.843  14.177  53.326  1.00 37.79           C  "
CANONICAL  = "ATOM   2425  CA  ALA A  10      10.000  11.000  12.000  1.00 20.00           C  "
BLANK_CHAIN = "HETATM 2423  C3  A1EI4   201      20.809  14.602  55.491  1.00 35.99           C  "
INSERT_CODE = "HETATM 2424  C4  A1EI4 B 201A     21.728  12.391  54.992  1.00 28.53           C  "


def test_shift_is_derived_from_the_code_width():
    assert _residue_shift(CANONICAL) == 0
    assert _residue_shift(REAL_5CHAR) == 2
    assert _residue_shift(BLANK_CHAIN) == 2
    assert _residue_shift(INSERT_CODE) == 2


@pytest.mark.parametrize("line,expected", [
    (REAL_5CHAR,  ("A1EI4", "B", "201", (21.843, 14.177, 53.326), "C")),
    (CANONICAL,   ("ALA",   "A", "10",  (10.0, 11.0, 12.0),       "C")),
    (BLANK_CHAIN, ("A1EI4", "",  "201", (20.809, 14.602, 55.491), "C")),
    (INSERT_CODE, ("A1EI4", "B", "201", (21.728, 12.391, 54.992), "C")),
])
def test_every_field_reads_correctly(line, expected):
    got = (field_res_name(line), field_chain(line), field_res_seq(line),
           field_coords(line), field_element(line))
    assert got == expected


def test_blank_chain_does_not_shift_the_coordinates():
    """The specific regression: tokenizing put the x coordinate in resSeq."""
    assert field_res_seq(BLANK_CHAIN) == "201", (
        "resSeq picked up a neighbouring field; int() would raise and the atom "
        "would be dropped")
    assert field_coords(BLANK_CHAIN) == (20.809, 14.602, 55.491)
    assert field_chain(BLANK_CHAIN) == "", "a blank chain must read as empty, not as a number"


def test_adjacent_full_width_coordinates_stay_separate():
    """%8.3f fills its field exactly at -100.123, so two in a row have no gap."""
    line = ("HETATM 2422  C1  A1EI4 B 201    -100.123-100.123-100.123  1.00 37.79           C  ")
    assert field_coords(line) == (-100.123, -100.123, -100.123)


def test_element_is_not_read_from_a_segid():
    """The element column moves with the shift; a segid must not stand in for it."""
    # element column blank, segid PROA present
    line = "HETATM 2422  C1  A1EI4 B 201      21.843  14.177  53.326  1.00 37.79 PROA         "
    assert field_element(line) == "C", "fell back to something other than the atom name"


def _write(tmpdir, lines):
    path = os.path.join(tmpdir, "t.pdb")
    with open(path, "w") as f:
        f.write("".join(l + "\n" for l in lines) + "END\n")
    return path


def test_no_atom_is_silently_dropped():
    """The headline case: every record in, every atom out."""
    with tempfile.TemporaryDirectory() as d:
        path = _write(d, [REAL_5CHAR, BLANK_CHAIN, INSERT_CODE, CANONICAL])
        atoms = parse_pdb_file(path)

    assert len(atoms) == 4, (
        f"{4 - len(atoms)} of 4 records were dropped. A partially parsed ligand "
        f"yields wrong distances, autobox centres and atom counts, with no error.")
    assert [a.res_name for a in atoms] == ["A1EI4", "A1EI4", "A1EI4", "ALA"]
    assert [a.chain for a in atoms] == ["B", "", "B", "A"]
    assert [a.res_num for a in atoms] == [201, 201, 201, 10]


def test_a_genuinely_unparseable_record_is_reported(capsys):
    """Dropping is sometimes right, but never silently."""
    broken = "HETATM 2422  C1  A1EI4 B 201      notanumber  14.177  53.326  1.00 37.79       C  "
    with tempfile.TemporaryDirectory() as d:
        path = _write(d, [REAL_5CHAR, broken])
        atoms = parse_pdb_file(path)

    assert len(atoms) == 1, "the good record should still parse"
    err = capsys.readouterr().err
    assert "skipped 1 unparseable" in err, (
        f"a dropped record produced no warning; stderr was {err!r}. Silence is "
        f"what kept the original CCD bug invisible.")
