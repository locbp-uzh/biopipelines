# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Runtime conversion of ESMFold2 constraint specs into upstream inputs.

The wrapper takes 1-indexed positions (framework convention) and upstream's
``Modification`` and ``CovalentBond`` are 0-indexed, so
the whole surface hinges on one subtraction happening exactly once. These tests
pin it without the ``esm`` package: ``to_structure_input`` takes the input
builder as a parameter, so a stub stands in for it.

Not covered here: ``resolve_covalent_atom_indices``, which imports
``esm.models.esmfold2.prepare_input`` to enumerate a residue's CCD atoms and
therefore needs the real model environment.
"""

import importlib.util
import os
from dataclasses import dataclass, field
from typing import Any, List, Optional

import pytest


def _inference_module():
    repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    spec = importlib.util.spec_from_file_location(
        "pipe_esmfold2_inference",
        os.path.join(repo, "pipe_scripts", "pipe_esmfold2_inference.py"),
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ── stub input builder ────────────────────────────────────────────────────────

@dataclass
class _Modification:
    position: int
    ccd: str
    smiles: Optional[str] = None


@dataclass
class _ProteinInput:
    id: str
    sequence: str
    modifications: Optional[List[_Modification]] = None
    msa: Any = None


@dataclass
class _DNAInput:
    id: str
    sequence: str
    modifications: Optional[List[_Modification]] = None
    msa: Any = None


@dataclass
class _RNAInput:
    id: str
    sequence: str
    modifications: Optional[List[_Modification]] = None
    msa: Any = None


@dataclass
class _LigandInput:
    id: str
    smiles: Optional[str] = None
    ccd: Optional[List[str]] = None


@dataclass
class _CovalentBond:
    chain_id1: str
    res_idx1: int
    atom_idx1: int
    chain_id2: str
    res_idx2: int
    atom_idx2: int


@dataclass
class _StructurePredictionInput:
    sequences: List[Any]
    pocket: Any = None
    distogram_conditioning: Any = None
    covalent_bonds: Any = None


class _StubBuilder:
    Modification = _Modification
    ProteinInput = _ProteinInput
    DNAInput = _DNAInput
    RNAInput = _RNAInput
    LigandInput = _LigandInput
    CovalentBond = _CovalentBond
    StructurePredictionInput = _StructurePredictionInput


SEQ = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ"


def _chains():
    return [
        dict(entity_type="protein", id="A", sequence=SEQ, msa_path=None),
        dict(entity_type="protein", id="B", sequence=SEQ, msa_path=None),
        dict(entity_type="ligand", id="C", ccd="ATP", smiles=""),
    ]


def _build(constraints, chains=None):
    mod = _inference_module()
    return mod.to_structure_input(
        chains if chains is not None else _chains(),
        _StubBuilder, None, 1024, constraints,
    )


# ── position conversion ───────────────────────────────────────────────────────

def test_modification_position_becomes_zero_indexed():
    spi = _build({"modifications": [{"chain": "A", "position": 12, "ccd": "SEP"}]})
    chain_a = next(s for s in spi.sequences if s.id == "A")
    assert chain_a.modifications == [_Modification(position=11, ccd="SEP")]
    # An unnamed chain gets no modifications rather than an empty list.
    assert next(s for s in spi.sequences if s.id == "B").modifications is None


def test_modifications_group_per_chain():
    spi = _build({"modifications": [
        {"chain": "A", "position": 1, "ccd": "SEP"},
        {"chain": "B", "position": 3, "ccd": "TPO"},
        {"chain": "A", "position": len(SEQ), "ccd": "PTR"},
    ]})
    by_id = {s.id: s for s in spi.sequences}
    assert [m.position for m in by_id["A"].modifications] == [0, len(SEQ) - 1]
    assert [m.ccd for m in by_id["A"].modifications] == ["SEP", "PTR"]
    assert [m.position for m in by_id["B"].modifications] == [2]


def test_covalent_bond_residues_become_zero_indexed():
    mod = _inference_module()
    bonds = mod.build_covalent_bonds(
        _StubBuilder,
        {"covalent_bonds": [{"atom1": ["A", 12, "SG"], "atom2": ["C", 1, "C1"]}]},
        ["A", "B", "C"],
    )
    assert bonds == [_CovalentBond(chain_id1="A", res_idx1=11, atom_idx1=0,
                                  chain_id2="C", res_idx2=0, atom_idx2=0)]


def test_no_constraints_leaves_the_input_bare():
    spi = _build({})
    assert spi.covalent_bonds is None
    assert all(getattr(s, "modifications", None) is None
               for s in spi.sequences if hasattr(s, "modifications"))


# ── rejections ────────────────────────────────────────────────────────────────

def test_modification_on_absent_chain_raises():
    with pytest.raises(ValueError, match="chain 'Z' is not in this complex"):
        _build({"modifications": [{"chain": "Z", "position": 12, "ccd": "SEP"}]})


def test_modification_past_end_of_chain_raises():
    with pytest.raises(ValueError, match="past the end of chain A"):
        _build({"modifications": [{"chain": "A", "position": len(SEQ) + 1, "ccd": "SEP"}]})


def test_modification_on_ligand_chain_raises():
    with pytest.raises(ValueError, match="modifications apply to polymer chains only"):
        _build({"modifications": [{"chain": "C", "position": 1, "ccd": "SEP"}]})


def test_covalent_bond_on_absent_chain_raises():
    mod = _inference_module()
    with pytest.raises(ValueError, match=r"covalent_bonds\[0\].atom2: chain 'Z'"):
        mod.build_covalent_bonds(
            _StubBuilder,
            {"covalent_bonds": [{"atom1": ["A", 12, "SG"], "atom2": ["Z", 1, "C1"]}]},
            ["A", "B", "C"],
        )
