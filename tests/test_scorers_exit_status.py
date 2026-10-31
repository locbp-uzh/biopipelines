"""ProLIF and RTMScore: the step fails when nothing was attempted or every attempt raised, not when the table is empty.

Both long-format tables hold several rows per id, so "no rows" and "no id succeeded" are different statements. A complex
whose ligand makes no interaction ProLIF detects, or a pair RTMScore answers with no pose, is attempted and succeeds and
writes no row, and the step then exited non-zero and blocked the batch. The discriminator is the count of attempts and of
failures, as in the pipe template.
"""

import importlib.util
import os
import sys
import types

import pytest

pd = pytest.importorskip("pandas")

from biopipelines.datastream import DataStream

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _load(filename):
    spec = importlib.util.spec_from_file_location(f"{filename}_exit_status", os.path.join(REPO_ROOT, "pipe_scripts", filename))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _streams(tmp_path, structure_ids, with_files=True):
    if with_files:
        for sid in structure_ids:
            (tmp_path / f"{sid}.pdb").write_text("ATOM      1  N   ALA A   1       0.0   0.0   0.0\n")
    (tmp_path / "l1.sdf").write_text("sdf\n")
    (tmp_path / "compounds.csv").write_text("id,code,smiles\nl1,LIG,CCO\n")
    DataStream(name="structures", ids=list(structure_ids), files=[str(tmp_path / "<id>.pdb")], format="pdb",
               map_table=str(tmp_path / "structures_map.csv")).save_json(str(tmp_path / "structures.json"))
    DataStream(name="compounds", ids=["l1"], files=[str(tmp_path / "<id>.sdf")], format="sdf",
               map_table=str(tmp_path / "compounds.csv")).save_json(str(tmp_path / "compounds.json"))


def _exit_code(main):
    try:
        main()
    except SystemExit as stopped:
        return stopped.code
    return 0


# ── RTMScore ────────────────────────────────────────────────────────────────

def _run_rtmscore(monkeypatch, tmp_path, structure_ids, outcome):
    """`outcome[structure id]` is a list of scores, or an exception to raise."""
    module = _load("pipe_rtmscore.py")
    _streams(tmp_path, structure_ids)
    monkeypatch.setattr(module, "write_ligand_sdf", lambda src, dst, smiles: open(dst, "w").write("sdf\n"))

    def run(script, model, protein, ligand, cutoff, work, prefix):
        result = outcome[os.path.basename(protein).replace("_protein.pdb", "")]
        if isinstance(result, Exception):
            raise result
        path = os.path.join(work, "out.csv")
        pd.DataFrame({"id": [f"pose{i}" for i in range(len(result))], "score": result}).to_csv(path, index=False)
        return path

    monkeypatch.setattr(module, "run_rtmscore", run)
    scores, missing = tmp_path / "scores.csv", tmp_path / "tables" / "004_RTMScore" / "missing.csv"
    monkeypatch.setattr(sys, "argv", [
        "pipe_rtmscore.py", "--structures-json", str(tmp_path / "structures.json"),
        "--compounds-json", str(tmp_path / "compounds.json"), "--rtmscore-script", "rtmscore.py",
        "--model-path", "model.pth", "--scratch-dir", str(tmp_path / "scratch"),
        "--scores-csv", str(scores), "--missing-csv", str(missing)])
    return _exit_code(module.main), pd.read_csv(scores), pd.read_csv(missing)


def test_rtmscore_pairs_that_answered_with_no_pose_do_not_fail_the_step(monkeypatch, tmp_path):
    code, scores, missing = _run_rtmscore(monkeypatch, tmp_path, ["s1", "s2"], {"s1": [], "s2": []})

    assert code == 0
    assert scores.empty and missing.empty


def test_rtmscore_a_single_failed_pair_does_not_fail_the_step(monkeypatch, tmp_path):
    code, scores, missing = _run_rtmscore(monkeypatch, tmp_path, ["s1", "s2"], {"s1": [1.5], "s2": RuntimeError("boom")})

    assert code == 0
    assert list(scores["id"]) == ["s1+l1"]
    assert list(missing["id"]) == ["s2+l1"]


def test_rtmscore_every_pair_failing_fails_the_step(monkeypatch, tmp_path):
    code, scores, missing = _run_rtmscore(monkeypatch, tmp_path, ["s1", "s2"],
                                          {"s1": RuntimeError("boom"), "s2": RuntimeError("boom")})

    assert code == 1
    assert len(missing) == 2


# ── ProLIF ──────────────────────────────────────────────────────────────────

class _Atoms:
    atoms = [object()]

    def write(self, path):
        open(path, "w").write("pdb\n")


def _install_fake_prolif(monkeypatch, interactions):
    """Stand-ins for MDAnalysis and prolif; `interactions[structure id]` is a row count or an exception."""
    mda = types.ModuleType("MDAnalysis")
    mda.Universe = lambda path, guess_bonds=True: types.SimpleNamespace(
        select_atoms=lambda selection: _Atoms(), path=path)
    plf = types.ModuleType("prolif")
    state = {}

    class Molecule:
        def __init__(self, mol):
            self.mol = mol

        @staticmethod
        def from_mda(group):
            return Molecule(group)

    class Fingerprint:
        def run_from_iterable(self, ligands, protein, progress=False):
            state["current"] = ligands[0].mol

        def to_dataframe(self):
            result = interactions[state["current"]]
            if isinstance(result, Exception):
                raise result
            if not result:
                return pd.DataFrame(index=[0], columns=pd.MultiIndex.from_arrays([[], [], []]))
            columns = pd.MultiIndex.from_tuples([("l1", f"ALA{i + 1}.A", "Hydrophobic") for i in range(result)])
            return pd.DataFrame([[True] * result], columns=columns)

    plf.Molecule, plf.Fingerprint = Molecule, Fingerprint
    monkeypatch.setitem(sys.modules, "MDAnalysis", mda)
    monkeypatch.setitem(sys.modules, "prolif", plf)


def _run_prolif(monkeypatch, tmp_path, structure_ids, interactions):
    _install_fake_prolif(monkeypatch, interactions)
    module = _load("pipe_prolif.py")
    _streams(tmp_path, structure_ids)
    # The fake ligand carries the id of its structure through the fingerprint.
    monkeypatch.setattr(module, "templated_ligand_mol", lambda pdb, smiles: os.path.basename(pdb).replace("_lig.pdb", ""))
    fingerprints, missing = tmp_path / "fp.csv", tmp_path / "tables" / "004_ProLIF" / "missing.csv"
    monkeypatch.setattr(sys, "argv", [
        "pipe_prolif.py", "--structures-json", str(tmp_path / "structures.json"),
        "--ligand-json", str(tmp_path / "compounds.json"), "--extras-dir", str(tmp_path / "extras"),
        "--fingerprints-csv", str(fingerprints), "--missing-csv", str(missing)])
    code = _exit_code(module.main)
    return code, pd.read_csv(fingerprints), pd.read_csv(missing)


def test_prolif_complexes_with_no_interaction_do_not_fail_the_step(monkeypatch, tmp_path):
    code, fingerprints, missing = _run_prolif(monkeypatch, tmp_path, ["s1", "s2"], {"s1": 0, "s2": 0})

    assert code == 0
    assert fingerprints.empty and missing.empty


def test_prolif_a_single_failed_complex_does_not_fail_the_step(monkeypatch, tmp_path):
    code, fingerprints, missing = _run_prolif(monkeypatch, tmp_path, ["s1", "s2"], {"s1": 2, "s2": RuntimeError("boom")})

    assert code == 0
    assert set(fingerprints["id"]) == {"s1"}
    assert list(missing["id"]) == ["s2"]


def test_prolif_every_complex_failing_fails_the_step(monkeypatch, tmp_path):
    code, fingerprints, missing = _run_prolif(monkeypatch, tmp_path, ["s1", "s2"],
                                              {"s1": RuntimeError("boom"), "s2": RuntimeError("boom")})

    assert code == 1
    assert len(missing) == 2
