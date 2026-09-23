"""RFdiffusion3's PDB element repair has to leave a spec-compliant file alone.

The repair was written against Biopython 1.86, whose PDBIO writes segID/element/charge one column early. Biopython 1.87 writes the spec layout, and the unconditional repair then read blanks, fell back to the atom name and turned every alpha carbon into calcium.
"""

import importlib.util
import os

import pytest

Bio = pytest.importorskip("Bio")

SCRIPT = os.path.join(os.path.dirname(__file__), "..", "pipe_scripts", "pipe_rfdiffusion3_postprocess.py")


@pytest.fixture(scope="module")
def post():
    spec = importlib.util.spec_from_file_location("pipe_rfdiffusion3_postprocess", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _written_by_pdbio(path):
    from Bio.PDB import PDBIO, StructureBuilder
    builder = StructureBuilder.StructureBuilder()
    builder.init_structure("x")
    builder.init_model(0)
    builder.init_chain("A")
    builder.init_seg("    ")
    builder.init_residue("ALA", " ", 1, " ")
    builder.init_atom("CA", [0.0, 0.0, 0.0], 0.0, 1.0, " ", " CA ", 1, "C")
    builder.init_residue("H_LIG", "H", 2, " ")
    builder.init_atom("SI1", [1.0, 0.0, 0.0], 0.0, 1.0, " ", "SI1 ", 2, "SI")
    writer = PDBIO()
    writer.set_structure(builder.get_structure())
    writer.save(str(path))


def _elements(path):
    return [line[76:78].strip() for line in open(path) if line.startswith(("ATOM", "HETATM"))]


def test_the_installed_pdbio_output_keeps_its_elements(post, tmp_path):
    target = tmp_path / "x.pdb"
    _written_by_pdbio(target)
    post._fix_element_columns(str(target))
    assert _elements(target) == ["C", "SI"]


def test_a_shifted_line_is_moved_into_place(post, tmp_path, monkeypatch):
    monkeypatch.setattr(post, "_pdbio_writes_shifted_columns", lambda: None)
    target = tmp_path / "shifted.pdb"
    shifted = "HETATM    2 SI1  LIG A   2       1.000   0.000   0.000  1.00  0.00         SI  \n"
    assert len(shifted.rstrip("\n")) == 79
    target.write_text(shifted)
    post._fix_element_columns(str(target))
    assert _elements(target) == ["SI"]


def test_the_repair_is_idempotent(post, tmp_path):
    target = tmp_path / "x.pdb"
    _written_by_pdbio(target)
    post._fix_element_columns(str(target))
    once = target.read_text()
    post._fix_element_columns(str(target))
    assert target.read_text() == once


@pytest.mark.parametrize("name, element", [(" CA ", "C"), ("SI1 ", "SI"), ("1HB ", "H"), ("FE  ", "FE")])
def test_a_blank_element_follows_the_atom_name_convention(post, name, element):
    assert post._element_from_atom_name(name) == element
