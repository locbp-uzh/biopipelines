"""Ligand pose RMSD must pair atoms by identity, not by position in the file.

Two programs writing the same molecule order its atoms differently -- Boltz2's
co-folded JF646 lists the same 44 atoms in a different order than the grafted
pose. Pairing by selection order then compares unrelated atoms. It happens not
to matter when the poses are far apart (every atom is far from every candidate),
which is exactly why it survived unnoticed: it corrupts the near-native poses
that a design campaign selects on.
"""

import importlib.util
import os
import sys

import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HELPER = os.path.join(REPO_ROOT, "pipe_scripts", "pipe_pose_change.py")

pytest.importorskip("pandas")


def _load_helper():
    spec = importlib.util.spec_from_file_location("pipe_pose_change", HELPER)
    module = importlib.util.module_from_spec(spec)
    sys.modules["pipe_pose_change"] = module
    spec.loader.exec_module(module)
    return module


class _FakeCmd:
    """Minimal stand-in for the PyMOL cmd object.

    Selections are "<obj>" or "(<obj>) and name A B C"; the latter is how the
    helper imposes a shared atom order before calling rms_cur.
    """

    def __init__(self, atoms):
        self.atoms = atoms  # {obj: [(name, (x, y, z)), ...]}

    def _resolve(self, sel):
        """Always file order: PyMOL ignores the order of a name list, verified
        against a real structure -- asking for reversed names returns the file
        order unchanged. A stub that honoured the request would hide that."""
        if ") and name " in sel:
            obj, _, names = sel.partition(") and name ")
            obj = obj.lstrip("(")
            wanted = {n.strip('"') for n in names.replace("+", " ").split()}
            return [(n, c) for n, c in self.atoms[obj] if n in wanted]
        return list(self.atoms[sel])

    def iterate(self, sel, expr, space):
        for name, _ in self._resolve(sel):
            space["names"].append(name)

    def get_coords(self, sel):
        return [c for _, c in self._resolve(sel)]

    def rms_cur(self, tgt_sel, ref_sel, matchmaker=-1):
        import math
        t = self._resolve(tgt_sel)
        r = self._resolve(ref_sel)
        assert len(t) == len(r)
        sq = sum(sum((a - b) ** 2 for a, b in zip(tc, rc)) for (_, tc), (_, rc) in zip(t, r))
        return math.sqrt(sq / len(t))

    def get_pdbstr(self, sel):
        return ""


def test_shuffled_atom_order_still_gives_zero_rmsd():
    """The same pose written in a different atom order is the same pose."""
    module = _load_helper()
    coords = {"C1": (0.0, 0.0, 0.0), "C2": (1.5, 0.0, 0.0),
              "N3": (3.0, 0.0, 0.0), "O4": (4.5, 0.0, 0.0)}
    ref = [(n, coords[n]) for n in ("C1", "C2", "N3", "O4")]
    tgt = [(n, coords[n]) for n in ("N3", "C1", "O4", "C2")]

    cmd = _FakeCmd({"ref": ref, "tgt": tgt})
    rmsd, pairing = module._ligand_rmsd(cmd, "ref", "tgt")

    assert pairing == "atom-name"
    assert rmsd == pytest.approx(0.0, abs=1e-9), (
        f"identical pose in shuffled order reported {rmsd:.3f} A")


def test_positional_pairing_would_have_failed_this():
    """Guard the guard: the shuffled order must be wrong if paired by position."""
    module = _load_helper()
    coords = {"C1": (0.0, 0.0, 0.0), "C2": (1.5, 0.0, 0.0),
              "N3": (3.0, 0.0, 0.0), "O4": (4.5, 0.0, 0.0)}
    ref = [(n, coords[n]) for n in ("C1", "C2", "N3", "O4")]
    tgt = [(n, coords[n]) for n in ("N3", "C1", "O4", "C2")]

    cmd = _FakeCmd({"ref": ref, "tgt": tgt})
    naive = cmd.rms_cur("tgt", "ref", matchmaker=-1)

    assert naive > 1.0, "test fixture no longer exercises the bug"


def test_real_displacement_is_still_measured():
    """Pairing by name must not mask a genuine pose difference."""
    module = _load_helper()
    ref = [("C1", (0.0, 0.0, 0.0)), ("C2", (1.5, 0.0, 0.0))]
    tgt = [("C2", (1.5, 0.0, 2.0)), ("C1", (0.0, 0.0, 2.0))]

    cmd = _FakeCmd({"ref": ref, "tgt": tgt})
    rmsd, pairing = module._ligand_rmsd(cmd, "ref", "tgt")

    assert pairing == "atom-name"
    assert rmsd == pytest.approx(2.0, abs=1e-9)


def test_uncorrespondable_ligands_raise_rather_than_guess():
    """No name correspondence and no RDKit match: refuse, do not pair blindly."""
    module = _load_helper()
    ref = [("C1", (0.0, 0.0, 0.0)), ("C2", (1.5, 0.0, 0.0))]
    tgt = [("XX", (0.0, 0.0, 0.0)), ("YY", (1.5, 0.0, 0.0))]

    cmd = _FakeCmd({"ref": ref, "tgt": tgt})
    with pytest.raises(ValueError, match="atom correspondence"):
        module._ligand_rmsd(cmd, "ref", "tgt")


def test_duplicate_atom_names_do_not_take_the_name_path():
    """Duplicate names cannot identify atoms; must not silently pair by order."""
    module = _load_helper()
    ref = [("C", (0.0, 0.0, 0.0)), ("C", (1.5, 0.0, 0.0))]
    tgt = [("C", (0.0, 0.0, 0.0)), ("C", (1.5, 0.0, 0.0))]

    cmd = _FakeCmd({"ref": ref, "tgt": tgt})
    # RDKit gets an empty PDB block from the stub, so this must raise rather
    # than fall through to positional pairing.
    with pytest.raises(ValueError, match="atom correspondence"):
        module._ligand_rmsd(cmd, "ref", "tgt")
