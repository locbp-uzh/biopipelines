"""Tests for resolve_selection_in_sequence: the common selection grammar applied
to a raw sequence (used e.g. to resolve a Boltz2 covalent position by motif)."""

import pytest
from biopipelines.pdb_parser import resolve_selection_in_sequence as resolve

# One Cys, in the conserved SNAP ILIPCH motif near the C-terminus.
SEQ = ("SGSRVLDLSPELRALILERVGPEGFAAAGAKPGERVVLGGGYGTLVDATTGNTSGSELHHPVFQQESF"
       "TRQVLWKLLKVVKFGEVISYSHLAALAGNPAATAAVKTALSGNPVPILIPCHRVVQGDLDVGGYEGG"
       "LAVKEWLLAHEGHRLGKPGLG")


def test_sequence_context_finds_unique_cys():
    # ILIPCH motif -> the C in it; renumbering-independent.
    assert SEQ.count("C") == 1
    out = resolve("C in ILIPCH", SEQ)
    assert len(out) == 1
    assert SEQ[out[0] - 1] == "C"   # 1-indexed position points at the Cys


def test_simple_residue_number():
    assert resolve("10", SEQ) == [10]


def test_residue_range():
    assert resolve("10-12", SEQ) == [10, 11, 12]


def test_first_last_keywords():
    assert resolve("first", SEQ) == [1]
    assert resolve("last", SEQ) == [len(SEQ)]


def test_context_not_found_raises():
    with pytest.raises(ValueError):
        resolve("C in ZZZZZ", SEQ)
