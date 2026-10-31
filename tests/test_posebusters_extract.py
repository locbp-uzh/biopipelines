"""PoseBusters: the protein written beside each extracted ligand keeps its ATOM records.

gemmi leaves the entity type of every residue unknown in a PDB without SEQRES, which is what a docking program writes. The protein was then written entirely as HETATM, PoseBusters read it as organic cofactors and ignored it, and every protein-ligand check came out False (`protein-ligand_maximum_distance` failed for twenty poses that sit in the pocket).

gemmi is in neither the CI environment nor a laptop one (only the posebusters environment has it), so the rule runs here against a stand-in gemmi, and the real parser is checked by the last test whenever gemmi is installed.
"""

import enum
import importlib.util
import os
import sys
import types

import pytest

pytest.importorskip("rdkit", reason="rdkit is not installed: the ligand SDF of the extraction is not written")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT = os.path.join(ROOT, "pipe_scripts", "pipe_posebusters.py")


@pytest.fixture(scope="module")
def pipe():
    spec = importlib.util.spec_from_file_location("pipe_posebusters", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class Container(list):
    """A gemmi structure, model, chain or residue: iterates its members and carries attributes."""

    def __init__(self, members=(), **attributes):
        super().__init__(members)
        self.__dict__.update(attributes)


@pytest.fixture
def stand_in_gemmi(monkeypatch):
    """The part of gemmi the extraction touches: `EntityType` and `read_structure`, with the structure set by the test."""
    module = types.ModuleType("gemmi")
    module.EntityType = enum.Enum("EntityType", "Unknown Polymer NonPolymer Branched Water")
    module.structure = None
    module.read_structure = lambda path: module.structure
    monkeypatch.setitem(sys.modules, "gemmi", module)
    return module


def stand_in_atom(serial, name, x, y, z, element):
    return types.SimpleNamespace(serial=serial, name=name, pos=types.SimpleNamespace(x=x, y=y, z=z), occ=1.0,
                                 b_iso=0.0, element=types.SimpleNamespace(name=element))


def stand_in_residue(gemmi, name, number, entity, het_flag, atoms):
    return Container(atoms, name=name, seqid=types.SimpleNamespace(num=number), het_flag=het_flag,
                     entity_type=getattr(gemmi.EntityType, entity))


@pytest.mark.parametrize("entity, het_flag, is_polymer", [
    ("Unknown", "A", True),
    ("Unknown", "H", False),
    ("Polymer", "A", True),
    ("Polymer", "H", True),
    ("NonPolymer", "H", False),
    ("NonPolymer", "A", False),
    ("Water", "H", False),
])
def test_a_residue_is_polymer_by_entity_or_by_its_atom_flag_when_there_is_no_entity(
    pipe, stand_in_gemmi, entity, het_flag, is_polymer,
):
    residue = stand_in_residue(stand_in_gemmi, "ALA", 1, entity, het_flag, [])
    assert pipe.is_polymer_residue(residue) is is_polymer


def test_a_protein_without_entity_information_is_written_as_atom_records_by_the_extraction(
    pipe, stand_in_gemmi, tmp_path,
):
    """The call site, not only the rule: the records of the protein file come from `is_polymer_residue`."""
    protein = Container([
        stand_in_residue(stand_in_gemmi, "ALA", 1, "Unknown", "A", [
            stand_in_atom(1, "N", 0.0, 0.0, 0.0, "N"), stand_in_atom(2, "CA", 1.4, 0.0, 0.0, "C")]),
        stand_in_residue(stand_in_gemmi, "HOH", 2, "Unknown", "H", [stand_in_atom(3, "O", 0.0, 5.0, 0.0, "O")]),
    ], name="A")
    ligand = Container([
        stand_in_residue(stand_in_gemmi, "LIG", 1, "Unknown", "H", [
            stand_in_atom(4, "C1", 2.0, 3.0, 0.0, "C"), stand_in_atom(5, "C2", 3.5, 3.0, 0.0, "C"),
            stand_in_atom(6, "O1", 4.2, 4.1, 0.0, "O")])
    ], name="Z")
    stand_in_gemmi.structure = Container([Container([protein, ligand])])

    results = pipe.extract_ligand_and_protein("docked.pdb", "LIG", str(tmp_path), "docked", ligand_smiles="CCO")

    assert len(results) == 1
    lines = [line for line in open(results[0][1]).read().splitlines() if line.strip() and not line.startswith("END")]
    assert [(line[:6].strip(), line[17:20]) for line in lines] == [("ATOM", "ALA"), ("ATOM", "ALA"), ("HETATM", "HOH")]


def atom(record, serial, name, residue, chain, number, xyz, element):
    return (f"{record:<6s}{serial:5d} {name:^4s} {residue:>3s} {chain}{number:4d}    "
            f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00  0.00          {element:>2s}\n")


def write_docked_complex(path):
    """Two residues of chain A as ATOM records, a water, and an ethanol named LIG on chain Z; no SEQRES, no entities."""
    lines = []
    serial = 1
    for number, x in ((1, 0.0), (2, 3.8)):
        for name, dx, element in (("N", 0.0, "N"), ("CA", 1.4, "C"), ("C", 2.5, "C")):
            lines.append(atom("ATOM", serial, name, "ALA", "A", number, (x + dx, 0.0, 0.0), element))
            serial += 1
    lines.append(atom("HETATM", serial, "O", "HOH", "A", 3, (0.0, 5.0, 0.0), "O"))
    serial += 1
    for name, xyz, element in (("C1", (2.0, 3.0, 0.0), "C"), ("C2", (3.5, 3.0, 0.0), "C"), ("O1", (4.2, 4.1, 0.0), "O")):
        lines.append(atom("HETATM", serial, name, "LIG", "Z", 1, xyz, element))
        serial += 1
    path.write_text("".join(lines) + "END\n")
    return str(path)


def test_a_protein_without_entity_information_is_written_as_atom_records(pipe, tmp_path):
    pytest.importorskip("gemmi", reason="gemmi is not installed in this environment: the real PDB parser is not exercised")
    complex_path = write_docked_complex(tmp_path / "docked.pdb")
    results = pipe.extract_ligand_and_protein(complex_path, "LIG", str(tmp_path), "docked", ligand_smiles="CCO")
    assert len(results) == 1
    protein = [line for line in open(results[0][1]) if line.strip() and not line.startswith("END")]

    records = {line[17:20]: line[:6].strip() for line in protein}
    assert records == {"ALA": "ATOM", "HOH": "HETATM"}
    assert sum(1 for line in protein if line.startswith("ATOM")) == 6
