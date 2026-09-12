"""One ligand chemistry per stream, or none -- never a silent pick.

resolve_ligand_smiles returned the FIRST non-empty smiles in a compounds stream,
and both callers (OpenMM's ligand= parameterisation, PoseBusters' bond-order
template) apply that one answer to every structure they process. A stream holding
a small library therefore parameterised every structure with compound #1's
chemistry -- wrong bond orders, wrong charges, wrong formula -- and reported the
resulting energies as if they described the molecule asked for.

It also wrapped the whole lookup in `except Exception: pass`, so an unreadable
stream was indistinguishable from a code-only ligand that legitimately has no
SMILES.
"""

import json
import os
import tempfile

import pytest

pytest.importorskip("pandas")

import pandas as pd  # noqa: E402

from biopipelines.ligand_utils import resolve_ligand_smiles  # noqa: E402


def _stream(tmpdir, rows):
    """A value-based compounds stream: files=[], content in the map_table."""
    csv = os.path.join(tmpdir, "compounds_map.csv")
    pd.DataFrame(rows).to_csv(csv, index=False)
    js = os.path.join(tmpdir, "compounds.json")
    with open(js, "w") as f:
        json.dump({"name": "compounds", "ids": [r["id"] for r in rows],
                   "files": [], "map_table": csv, "format": "csv"}, f)
    return js


def test_single_compound_resolves():
    with tempfile.TemporaryDirectory() as d:
        js = _stream(d, [{"id": "lig1", "code": "LIG", "smiles": "CCO"}])
        assert resolve_ligand_smiles(js) == "CCO"


def test_repeated_identical_smiles_is_not_ambiguous():
    """The same molecule under several ids is still one chemistry."""
    with tempfile.TemporaryDirectory() as d:
        js = _stream(d, [{"id": "a", "code": "LIG", "smiles": "CCO"},
                         {"id": "b", "code": "LIG", "smiles": "CCO"}])
        assert resolve_ligand_smiles(js) == "CCO"


def test_code_only_ligand_has_no_smiles():
    """None is a legitimate answer -- Ligand(code=...) carries no chemistry."""
    with tempfile.TemporaryDirectory() as d:
        js = _stream(d, [{"id": "zit", "code": "ZIT", "smiles": ""}])
        assert resolve_ligand_smiles(js) is None


def test_multiple_distinct_smiles_raises():
    """The wrong-molecule case: no single answer exists, so do not invent one."""
    with tempfile.TemporaryDirectory() as d:
        js = _stream(d, [{"id": "a", "code": "LIG", "smiles": "CCO"},
                         {"id": "b", "code": "LIG", "smiles": "c1ccccc1"}])
        with pytest.raises(ValueError, match="distinct SMILES"):
            resolve_ligand_smiles(js)


def test_an_unreadable_stream_is_not_reported_as_no_chemistry():
    """A missing file must not masquerade as a code-only ligand."""
    with pytest.raises(Exception):
        resolve_ligand_smiles(os.path.join(tempfile.gettempdir(), "does_not_exist.json"))
