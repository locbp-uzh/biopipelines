# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""What PLACER declares it will produce is what its driver produces.

`get_output_files` predicts ids through `predict_output_ids_with_provenance`, which honors `Bundle` and `Each`, while `pipe_placer_driver.py` composes them with a nested loop over the two axes. That pairing is what broke Mock: a hand-rolled product has no notion of bundle versus each, so a bundled axis iterated where it should have contributed one collapsed prefix, and the completion check marked the step FAILED.

PLACER cannot reach that state, and these cases pin the two independent reasons why. First, the constructor's type gate accepts only a `DataStream` or a `StandardizedOutput`, so a wrapper on either axis is refused before any prediction happens -- both axes are therefore always `each`, and for two iterated axes a flat product is exactly what the framework composer produces. Second, predict and produce are asserted equal on the reachable combinations, including the no-ligand degradation to a bare structure id.

The rejection assertions are the alarm, not a preference: if PLACER ever widens that gate to take `Bundle`, the flat product in the driver becomes wrong on the spot and these tests are what say so. The fix at that point is Mock's -- record each axis's mode in the config and call `predict_single_output_id` from `compose_output_id`.
"""

import importlib.util
import os

import pytest

from biopipelines.combinatorics import (
    Bundle,
    Each,
    predict_output_ids_with_provenance,
    predict_single_output_id,
)
from biopipelines.datastream import DataStream, create_map_table
from biopipelines.placer import PLACER


REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

AXIS_MODES = [
    ("each", "each", Each, Each),
    ("bundle", "each", Bundle, Each),
    ("each", "bundle", Each, Bundle),
    ("bundle", "bundle", Bundle, Bundle),
]


def _driver():
    spec = importlib.util.spec_from_file_location(
        "pipe_placer_driver_under_test",
        os.path.join(REPO_ROOT, "pipe_scripts", "pipe_placer_driver.py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _stream(tmp_path, name, ids, fmt, ext):
    table = tmp_path / f"{name}_map.csv"
    create_map_table(str(table), ids=list(ids))
    files = []
    if ext:
        for one in ids:
            (tmp_path / f"{one}.{ext}").write_text("ATOM\n")
        files = [str(tmp_path / f"<id>.{ext}")]
    return DataStream(name=name, ids=list(ids), map_table=str(table), format=fmt, files=files)


@pytest.fixture
def structures(tmp_path):
    return _stream(tmp_path, "structures", ["p1", "p2"], "pdb", "pdb")


@pytest.fixture
def ligands(tmp_path):
    return _stream(tmp_path, "compounds", ["LIG1", "LIG2", "LIG3"], "csv", None)


def _produce(structure_ids, ligand_ids):
    """Replay the driver's nested composition loop over the two axes."""
    compose = _driver().compose_output_id
    return [compose(sid, lid) for sid in structure_ids for lid in ligand_ids]


# ── why no axis mode can reach the driver ─────────────────────────────────────

@pytest.mark.parametrize("structures_mode, ligand_mode, wrap_structures, wrap_ligand", AXIS_MODES)
def test_a_wrapped_axis_is_refused_before_any_id_is_predicted(
    structures, ligands, structures_mode, ligand_mode, wrap_structures, wrap_ligand, record_case,
):
    """Neither axis accepts Bundle/Each, so the driver's flat product can never meet a collapsed prefix."""
    with pytest.raises(ValueError) as excinfo:
        PLACER(structures=wrap_structures(structures), ligand=wrap_ligand(ligands))

    record_case(input=f"structures={structures_mode}, ligand={ligand_mode}",
                expected="ValueError naming the accepted input types",
                actual=str(excinfo.value))
    assert "DataStream or StandardizedOutput" in str(excinfo.value)


def test_a_wrapped_structures_axis_is_refused_in_sidechain_mode_too(structures):
    with pytest.raises(ValueError, match="DataStream or StandardizedOutput"):
        PLACER(structures=Bundle(structures), target_res="A-149")


# ── predict equals produce on the reachable combinations ──────────────────────

def test_predict_equals_produce_for_the_two_iterated_axes(structures, ligands, record_case):
    predicted, _prov = predict_output_ids_with_provenance(
        structures=(structures, "structures"), compounds=(ligands, "compounds"))
    produced = _produce(structures.ids_expanded, ligands.ids_expanded)

    record_case(input="structures[2] x compounds[3], both iterated",
                expected=predicted, actual=produced)
    assert predicted == ["p1+LIG1", "p1+LIG2", "p1+LIG3", "p2+LIG1", "p2+LIG2", "p2+LIG3"]
    assert produced == predicted


def test_a_missing_ligand_axis_degrades_to_the_bare_structure_id(structures, record_case):
    """Sidechain/apo mode predicts the structure ids untouched, and the driver's `lig_id is None` branch must match."""
    predicted = list(structures.ids)
    produced = _produce(structures.ids_expanded, [None])

    record_case(input="structures[2], no ligand axis", expected=predicted, actual=produced)
    assert produced == ["p1", "p2"]
    assert produced == predicted


def test_the_shared_composer_agrees_on_the_single_axis_case(structures):
    """A one-axis call to the framework composer returns the bare id, so switching over would not change apo naming."""
    through_composer = [
        predict_single_output_id(structures=("each", list(structures.ids), i, [], False))
        for i in range(len(structures.ids))
    ]
    assert through_composer == _produce(structures.ids_expanded, [None])


def test_the_shared_composer_agrees_on_the_two_axis_case(structures, ligands):
    """The driver's flat product is the framework's product for two iterated axes, id for id and in order."""
    through_composer = [
        predict_single_output_id(
            structures=("each", list(structures.ids), i, [], False),
            compounds=("each", list(ligands.ids), j, [], False))
        for i in range(len(structures.ids))
        for j in range(len(ligands.ids))
    ]
    assert through_composer == _produce(structures.ids_expanded, ligands.ids_expanded)


# ── the bonds= path ───────────────────────────────────────────────────────────

def test_the_driver_can_read_a_bonds_json(tmp_path):
    """`bonds=` reaches `json.load` in the driver, which had no `json` import -- a NameError on a documented parameter."""
    import json

    bonds_json = tmp_path / "bonds.json"
    bonds_json.write_text(json.dumps([["A145.SG", "LIG.C12", 1.8]]))

    driver = _driver()
    with open(bonds_json) as handle:
        assert driver.json.load(handle) == [["A145.SG", "LIG.C12", 1.8]]
