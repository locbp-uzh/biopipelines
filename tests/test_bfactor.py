"""BFactor reads the per-residue B-factor a predictor wrote, so its selection parsing decides what is averaged.

1.0 averaged a chainless span across every chain of a complex, turned a typo into an empty selection reported as n=0, and merged insertion-code residues (100 and 100A) into one.
"""

import importlib.util
import os

import pytest

pytest.importorskip("yaml")

SCRIPT = os.path.join(os.path.dirname(__file__), "..", "pipe_scripts", "pipe_bfactor.py")


@pytest.fixture(scope="module")
def bf():
    spec = importlib.util.spec_from_file_location("pipe_bfactor", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _pdb(path, atoms):
    """atoms: [(chain, resi, icode, bfactor)], one CA each."""
    lines = [f"ATOM  {n:5d}  CA  ALA {c}{r:4d}{i or ' '}   {n:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00{b:6.2f}           C"
             for n, (c, r, i, b) in enumerate(atoms, 1)]
    path.write_text("\n".join(lines) + "\nEND\n")
    return str(path)


def test_a_typo_is_an_error_not_an_empty_selection(bf):
    with pytest.raises(ValueError, match="not a residue span"):
        bf.parse_selection("A75-77+A2x4")


def test_a_chainless_span_on_a_complex_is_refused(bf, tmp_path):
    residues = bf.read_residues(_pdb(tmp_path / "x.pdb", [("A", 75, "", 90.0), ("B", 75, "", 30.0)]))
    with pytest.raises(ValueError, match="chains A, B"):
        bf.check_chains("site", bf.parse_selection("75-77"), residues)


def test_a_ligand_or_water_chain_does_not_make_a_span_ambiguous(bf, tmp_path):
    """A predicted complex carries a ligand chain; refusing it would fail BFactor's main use case."""
    residues = bf.read_residues(_pdb(tmp_path / "x.pdb", [("A", 75, "", 90.0), ("B", 1, "", 50.0),
                                                          ("W", 1, "", 10.0)]))
    bf.check_chains("site", bf.parse_selection("75-77"), residues)


def test_a_chainless_span_on_a_monomer_is_fine(bf, tmp_path):
    residues = bf.read_residues(_pdb(tmp_path / "x.pdb", [("A", 75, "", 90.0)]))
    bf.check_chains("site", bf.parse_selection("75-77"), residues)


def test_insertion_codes_stay_separate_residues(bf, tmp_path):
    residues = bf.read_residues(_pdb(tmp_path / "x.pdb", [("A", 100, "", 80.0), ("A", 100, "A", 20.0)]))
    assert [(c, r, i, b) for c, r, i, b in residues] == [("A", 100, "", 80.0), ("A", 100, "A", 20.0)]


def test_the_wrapper_rejects_a_literal_typo_at_configuration(local_config, isolated_cwd):
    from biopipelines.bfactor import BFactor
    from biopipelines.mock import Mock
    from biopipelines.pipeline import Pipeline

    with Pipeline(project="TestSuite", job="bf", description="x", on_the_fly=False,
                  local_output=True, config="local"):
        source = Mock(ids=["s"], streams={"structures": {"format": "pdb", "file": "<id>.pdb"}})
        with pytest.raises(ValueError, match="not a residue span"):
            BFactor(structures=source.streams.structures, selections={"site": "A75-77+Ax"})
