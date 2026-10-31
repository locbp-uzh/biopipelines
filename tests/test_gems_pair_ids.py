# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""GEMS: a row of `affinity` is named `<protein>+<ligand>`, and the staging name is not an id.

The pair used to be staged as `<protein>__<ligand>` and that same string became the row id of `affinity` and of `missing`, a spelling no other tool uses for a two-axis id. The two axes are recorded in `structures.id` and `ligands.id`, so the id has to be the axis rule's `+` join. The files GEMS reads are named by a staging name that carries no id at all, and the predictions it writes are mapped back to the pair by exact lookup on that name.
"""

import importlib.util
import os

import pytest

pd = pytest.importorskip("pandas")

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# `a__b` + `c` and `a` + `b__c` both spell `a__b__c` under the old double-underscore join.
PROTEINS = ["a__b", "a"]
LIGANDS = ["c", "b__c"]


@pytest.fixture
def gems(monkeypatch):
    spec = importlib.util.spec_from_file_location("pipe_gems_under_test", os.path.join(REPO_ROOT, "pipe_scripts", "pipe_gems.py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    # the ligand conversion needs RDKit; the staging names are what is under test
    monkeypatch.setattr(module, "write_ligand_sdf", lambda src, dst, smiles: open(dst, "w").write("sdf\n"))
    return module


def _stage(gems, tmp_path):
    structures = []
    for pid in PROTEINS:
        path = tmp_path / f"{pid}.pdb"
        path.write_text("ATOM      1  N   ALA A   1       0.0   0.0   0.0\nHETATM    2  C1  LIG A   2       1.0   1.0   1.0\n")
        structures.append((pid, str(path)))
    ligands = [(lid, str(tmp_path / f"{lid}.sdf")) for lid in LIGANDS]
    work = tmp_path / "pairs"
    meta = gems.stage_pairs(structures, ligands, str(work), {})
    return meta, work


def test_every_pair_gets_the_axis_rule_id(gems, tmp_path):
    meta, _ = _stage(gems, tmp_path)

    assert [m["id"] for m in meta] == ["a__b+c", "a__b+b__c", "a+c", "a+b__c"]


def test_provenance_carries_the_parent_of_each_axis(gems, tmp_path):
    meta, _ = _stage(gems, tmp_path)

    assert [(m["structures.id"], m["ligands.id"]) for m in meta] == [
        ("a__b", "c"), ("a__b", "b__c"), ("a", "c"), ("a", "b__c")]


def test_two_distinct_pairs_never_share_a_staged_file(gems, tmp_path):
    meta, work = _stage(gems, tmp_path)

    names = [m["stage_name"] for m in meta]
    assert len(set(names)) == len(meta)
    assert sorted(os.listdir(work)) == sorted(f"{n}{ext}" for n in names for ext in (".pdb", ".sdf"))


def test_the_staging_name_spells_neither_axis(gems, tmp_path):
    meta, _ = _stage(gems, tmp_path)

    for m in meta:
        assert m["stage_name"] not in (m["id"],)
        assert "__" not in m["stage_name"] and "+" not in m["stage_name"]


def test_the_staged_protein_drops_its_hetatm_lines(gems, tmp_path):
    meta, work = _stage(gems, tmp_path)

    staged = (work / f"{meta[0]['stage_name']}.pdb").read_text()
    assert "HETATM" not in staged and "ATOM" in staged


def test_rows_are_mapped_back_from_the_predictions_by_exact_lookup(gems, tmp_path):
    meta, _ = _stage(gems, tmp_path)
    predictions = {m["stage_name"]: 5.0 + i for i, m in enumerate(meta)}

    rows, missing = gems.affinity_rows(meta, predictions, step_id="004")

    assert missing == []
    assert [(r["id"], r["structures.id"], r["ligands.id"], r["pkd_pred"]) for r in rows] == [
        ("a__b+c", "a__b", "c", 5.0), ("a__b+b__c", "a__b", "b__c", 6.0),
        ("a+c", "a", "c", 7.0), ("a+b__c", "a", "b__c", 8.0)]


def test_a_pair_without_a_prediction_is_a_failure_under_its_own_id(gems, tmp_path):
    meta, _ = _stage(gems, tmp_path)
    predictions = {meta[0]["stage_name"]: 5.0}

    rows, missing = gems.affinity_rows(meta, predictions, step_id="004")

    assert [r["id"] for r in rows] == ["a__b+c", "a__b+b__c", "a+c", "a+b__c"]
    assert [r["pkd_pred"] == r["pkd_pred"] for r in rows] == [True, False, False, False]
    assert [(m["id"], m["kind"]) for m in missing] == [
        ("a__b+b__c", "failure"), ("a+c", "failure"), ("a+b__c", "failure")]
    assert all(m["removed_by"] == "004" for m in missing)
