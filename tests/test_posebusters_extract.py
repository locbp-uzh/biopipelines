"""PoseBusters: the protein written beside each extracted ligand keeps its ATOM records.

gemmi leaves the entity type of every residue unknown in a PDB without SEQRES, which is what a docking program writes. The protein was then written entirely as HETATM, PoseBusters read it as organic cofactors and ignored it, and every protein-ligand check came out False (`protein-ligand_maximum_distance` failed for twenty poses that sit in the pocket). Needs gemmi, which only the posebusters environment has.
"""

import importlib.util
import os

import pytest

pytest.importorskip("gemmi")
pytest.importorskip("rdkit")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT = os.path.join(ROOT, "pipe_scripts", "pipe_posebusters.py")


@pytest.fixture(scope="module")
def pipe():
    spec = importlib.util.spec_from_file_location("pipe_posebusters", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


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
    complex_path = write_docked_complex(tmp_path / "docked.pdb")
    results = pipe.extract_ligand_and_protein(complex_path, "LIG", str(tmp_path), "docked", ligand_smiles="CCO")
    assert len(results) == 1
    protein = [line for line in open(results[0][1]) if line.strip() and not line.startswith("END")]

    records = {line[17:20]: line[:6].strip() for line in protein}
    assert records == {"ALA": "ATOM", "HOH": "HETATM"}
    assert sum(1 for line in protein if line.startswith("ATOM")) == 6
