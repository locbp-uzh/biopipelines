"""ProteinMPNN positions are 1-based ranks within a chain, not PDB residue numbers.

Upstream indexes a mask sized to the residue COUNT:

    fixed_position_mask = np.ones(chain_length)
    fixed_position_mask[np.array(fixed_pos_list)-1] = 0.0

so a residue number only works when the chain is numbered 1..N with no gaps -- which is
what Boltz2 and RFdiffusion3 emit, and why this stayed hidden. A trimmed chain numbered
35-174 (140 residues) sends 174, indexing 173 into a 140-long array.
"""


def _rank_map(residues):
    """The conversion the pipe script applies before writing fixed_pos.jsonl."""
    return {resnum: i + 1 for i, resnum in enumerate(residues)}


def _convert(residues, positions):
    rank = _rank_map(residues)
    return [rank[p] for p in positions if p in rank]


def test_offset_chain_positions_are_ranked():
    """A chain numbered 35-174 must map onto 1-140, not stay at 35-174."""
    residues = list(range(35, 175))          # 140 residues, contiguous
    assert len(residues) == 140

    converted = _convert(residues, [35, 100, 174])
    assert converted == [1, 66, 140]
    assert max(converted) <= len(residues), "must never index past the mask"


def test_gapped_chain_uses_order_not_arithmetic():
    """With a gap, subtracting a constant offset is wrong; rank by parsed order."""
    residues = [10, 11, 12, 50, 51]          # 5 residues, gap 13-49
    converted = _convert(residues, [10, 50, 51])
    assert converted == [1, 4, 5]
    # A naive `resnum - first + 1` would give 41 and 42 here, both out of bounds.
    assert max(converted) <= len(residues)


def test_one_based_chain_is_unchanged():
    """The case that always worked keeps working: number == rank."""
    residues = list(range(1, 166))           # Boltz2/RFD3 output, 1-165
    positions = [1, 83, 165]
    assert _convert(residues, positions) == positions


def test_absent_positions_are_dropped_not_indexed():
    """A position with no residue in the chain is dropped rather than mis-indexed."""
    residues = [1, 2, 3]
    assert _convert(residues, [2, 99]) == [2]
