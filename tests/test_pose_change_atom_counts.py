"""PoseChange must refuse ligands whose atom counts differ.

`cmd.rms_cur(..., matchmaker=-1)` pairs atoms by their order in the selection. Given
selections of different sizes PyMOL writes `ExecutiveRMS-Error: Atom counts between
selections don't match` to stderr and **returns 0.0** — it does not raise. Without a
guard the tool therefore records a perfect 0.000 A ligand RMSD for a comparison it never
performed, which is the worst possible failure: silent, and indistinguishable from the
best possible result.

The case that surfaced it: comparing a reacted dye fragment grafted from a design (44
atoms, no benzylguanine) against a docked intact BG-dye (55 atoms).

`calculate_pose_change` takes the PyMOL command object as its first argument, so these
tests pass a stand-in rather than stubbing the pymol module.
"""

import importlib.util
import os
import sys

import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HELPER = os.path.join(REPO_ROOT, "pipe_scripts", "pipe_pose_change.py")


def _load_helper():
    spec = importlib.util.spec_from_file_location("pipe_pose_change", HELPER)
    module = importlib.util.module_from_spec(spec)
    sys.modules["pipe_pose_change"] = module
    spec.loader.exec_module(module)
    return module


class _FakeCmd:
    """Reproduces rms_cur's return-0.0-on-mismatch behaviour."""

    def __init__(self, counts):
        self._counts = counts
        self.rms_cur_called = False

    def load(self, path, name):
        pass

    def delete(self, name):
        pass

    def count_atoms(self, selection):
        for resn, count in self._counts.items():
            if resn in selection:
                return count
        return 0

    def align(self, mobile, target):
        return (0.5, 0, 0, 0, 0, 0, 0)

    def iterate(self, selection, expr, space):
        """Distinct per-ligand atom names, so RMSD pairs by name as it does live."""
        for i in range(self.count_atoms(selection)):
            space["names"].append(f"C{i}")

    def get_coords(self, selection):
        return [(float(i), 0.0, 0.0) for i in range(self.count_atoms(selection))]

    def centerofmass(self, selection):
        n = self.count_atoms(selection)
        return [sum(range(n)) / n if n else 0.0, 0.0, 0.0]

    def get_pdbstr(self, selection):
        return ""

    def rms_cur(self, mobile, target, matchmaker=-1):
        self.rms_cur_called = True
        return 0.0  # what PyMOL actually returns when the counts disagree


def _call(module, cmd, tmp_path):
    ref = tmp_path / "ref.pdb"
    tgt = tmp_path / "tgt.pdb"
    ref.write_text("END\n")
    tgt.write_text("END\n")
    return module.calculate_pose_change(
        cmd, str(ref), "UNL", str(tgt), "t1", "LIG",
        "not resn UNL", "not resn LIG",
    )


def test_mismatched_ligand_atom_counts_raise(tmp_path, record_case):
    """44 vs 55 atoms must raise, not silently report 0.000 A."""
    module = _load_helper()
    cmd = _FakeCmd({"UNL": 44, "LIG": 55})

    with pytest.raises(ValueError) as excinfo:
        _call(module, cmd, tmp_path)

    message = str(excinfo.value)
    # Pairing each count with its own selection catches a swap; bare "44"/"55" would not.
    expected_ref, expected_tgt = "reference 'UNL' has 44", "target 'LIG' has 55"
    record_case(
        input="reference UNL=44 atoms, target LIG=55 atoms",
        expected=(True, False),
        actual=(
            expected_ref in message and expected_tgt in message,
            cmd.rms_cur_called,
        ),
    )
    assert expected_ref in message, message
    assert expected_tgt in message, message
    assert not cmd.rms_cur_called, "rms_cur must not run on mismatched selections"


def test_matching_ligand_atom_counts_proceed(tmp_path, record_case):
    """Equal counts still reach the RMSD calculation."""
    module = _load_helper()
    cmd = _FakeCmd({"UNL": 44, "LIG": 44})

    result = _call(module, cmd, tmp_path)

    record_case(
        input="reference UNL=44 atoms, target LIG=44 atoms",
        expected=("atom-name", 44),
        actual=(result["rmsd_pairing"], result["num_ligand_atoms"]),
    )
    assert result["rmsd_pairing"] == "atom-name"
    assert result["ligand_rmsd"] == 0.0
    assert result["num_ligand_atoms"] == 44
