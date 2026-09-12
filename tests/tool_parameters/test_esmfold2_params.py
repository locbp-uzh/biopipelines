# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Parameter coverage for the ESMFold2 constraint parameters.

``modifications`` and ``covalent_bonds`` are the StructurePredictionInput
constraints the wrapper passes through. They are not
inlined into bash — they are serialized to
``_configuration/constraints_config.json`` — so the emitted-artifact grep covers
both the flag and the file. The rejection cases matter as much as the happy
path: upstream drops an unresolvable covalent bond silently, and its
``Modification.smiles`` field is an unimplemented stub, so both have to fail at
construction rather than at fold time.
"""

import json
from pathlib import Path

import pytest

from ._helpers import assert_substrings_in, read_all_emitted_artifacts


pytestmark = pytest.mark.tool_parameters


def _protein(ids="p1", sequence="MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ"):
    """A real Sequence entity: the combinatorics axis needs a map_table."""
    from biopipelines.sequence import Sequence
    return Sequence(seq=sequence, type="protein", ids=ids)


def _constraints_json(script_path):
    """The single constraints_config.json emitted under the run root."""
    run_root = Path(script_path).resolve().parent.parent
    hits = sorted(run_root.rglob("_configuration/constraints_config.json"))
    assert hits, f"constraints_config.json not emitted under {run_root}"
    return json.loads(hits[0].read_text())


# ── happy paths ───────────────────────────────────────────────────────────────

def test_modifications_reach_the_constraints_config(local_config, isolated_cwd, new_pipeline):
    from biopipelines.esmfold2 import ESMFold2
    pipeline = new_pipeline("esmfold2_modifications")
    with pipeline:
        q = _protein()
        ESMFold2(proteins=q,
                 modifications=[{"chain": "A", "position": 12, "ccd": "SEP"},
                                {"chain": "A", "position": 20, "ccd": "TPO"}])
        script_path = pipeline.save()
    content = read_all_emitted_artifacts(script_path)
    assert_substrings_in(content, ["--constraints", "constraints_config.json", "SEP", "TPO"])

    payload = _constraints_json(script_path)
    # Positions stay 1-indexed on this side of the boundary.
    assert payload["modifications"] == [
        {"chain": "A", "position": 12, "ccd": "SEP"},
        {"chain": "A", "position": 20, "ccd": "TPO"},
    ]


def test_covalent_bonds_reach_the_constraints_config(
        local_config, isolated_cwd, new_pipeline):
    from biopipelines.esmfold2 import ESMFold2
    pipeline = new_pipeline("esmfold2_bonds")
    with pipeline:
        q = _protein()
        ESMFold2(proteins=q,
                 covalent_bonds=[{"atom1": ["A", 12, "SG"], "atom2": ["B", 1, "C1"]}])
        script_path = pipeline.save()
    payload = _constraints_json(script_path)
    assert payload["covalent_bonds"] == [{"atom1": ["A", 12, "SG"], "atom2": ["B", 1, "C1"]}]


def test_pocket_is_not_a_parameter(local_config, isolated_cwd, new_pipeline):
    """Upstream accepts PocketConditioning and discards it (esm 3.3.0 zeroes
    pocket_feature unconditionally), so the wrapper must not offer one -- and
    BaseConfig would swallow the kwarg rather than reject it."""
    from biopipelines.esmfold2 import ESMFold2
    with new_pipeline("esmfold2_no_pocket"):
        q = _protein()
        with pytest.raises(ValueError, match="no pocket parameter"):
            ESMFold2(proteins=q, pocket={"binder": "B", "contacts": [["A", 12]]})


def test_no_constraints_config_when_unused(local_config, isolated_cwd, new_pipeline):
    from biopipelines.esmfold2 import ESMFold2
    pipeline = new_pipeline("esmfold2_no_constraints")
    with pipeline:
        q = _protein()
        ESMFold2(proteins=q)
        script_path = pipeline.save()
    run_root = Path(script_path).resolve().parent.parent
    assert not list(run_root.rglob("_configuration/constraints_config.json"))
    assert "--constraints" not in read_all_emitted_artifacts(script_path)


# ── rejections ────────────────────────────────────────────────────────────────

@pytest.mark.parametrize("modifications, expected", [
    ([{"chain": "A", "position": 12, "smiles": "OP(=O)(O)O"}], "smiles"),
    ([{"chain": "A", "position": 0, "ccd": "SEP"}], "1-indexed"),
    ([{"chain": "A", "position": True, "ccd": "SEP"}], "must be an int"),
    ([{"chain": "A", "position": 12}], "missing key"),
    ([{"chain": "A", "position": 12, "ccd": "SEP", "chian": "A"}], "unknown key"),
    ([{"chain": "", "position": 12, "ccd": "SEP"}], "non-empty chain id"),
    ([{"chain": "A", "position": 12, "ccd": ""}], "non-empty CCD code"),
    ([], "non-empty list"),
])
def test_bad_modifications_raise(local_config, isolated_cwd, new_pipeline,
                                 modifications, expected):
    from biopipelines.esmfold2 import ESMFold2
    with new_pipeline("esmfold2_bad_modifications"):
        q = _protein()
        with pytest.raises(ValueError, match=expected):
            ESMFold2(proteins=q, modifications=modifications)


@pytest.mark.parametrize("bonds, expected", [
    ([{"atom1": ["A", 12, "SG"]}], "missing key"),
    ([{"atom1": ["A", 12], "atom2": ["B", 1, "C1"]}],
     r"must be \[chain, position, atom_name\]"),
    ([{"atom1": ["A", 12, ""], "atom2": ["B", 1, "C1"]}], "atom_name must be"),
    ([{"atom1": ["A", 0, "SG"], "atom2": ["B", 1, "C1"]}], "1-indexed"),
])
def test_bad_covalent_bonds_raise(local_config, isolated_cwd, new_pipeline, bonds, expected):
    from biopipelines.esmfold2 import ESMFold2
    with new_pipeline("esmfold2_bad_bonds"):
        q = _protein()
        with pytest.raises(ValueError, match=expected):
            ESMFold2(proteins=q, covalent_bonds=bonds)
