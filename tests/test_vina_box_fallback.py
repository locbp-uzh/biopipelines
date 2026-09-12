"""Vina must not dock inside a box that misses the protein.

With no explicit center+size and no crystal ligand in the receptor, box_args
falls back to the coordinates of the INPUT ligand -- the conformer about to be
docked, since nothing has been docked yet. Whether that is sensible depends
entirely on where the conformer came from:

  * carved from this receptor (Ligand(structures=..., codes=...)), its
    coordinates ARE the crystal pose in the receptor's frame, and boxing on it is
    the crystal-ligand autobox by another route;
  * embedded by RDKit from SMILES, it is centred on the origin while a PDB
    receptor sits wherever the crystallographers left it -- so the box lands tens
    of angstroms away and Vina searches empty space.

The second case used to run silently and return scores that look like any other
scores. It now raises.
"""

import importlib.util
import os
import pathlib
import tempfile

import pytest

pytest.importorskip("rdkit")


@pytest.fixture(scope="module")
def vina_backend():
    """Load the pipe script directly; pipe_scripts is not an importable package."""
    path = (pathlib.Path(__file__).resolve().parent.parent
            / "pipe_scripts" / "pipe_vina_backend.py")
    spec = importlib.util.spec_from_file_location("pipe_vina_backend", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# Receptor atoms around (21, 14, 53) -- roughly where 9RTM's chain sits.
RECEPTOR = "\n".join(
    f"ATOM  {i:5d}  CA  ALA A{i:4d}    "
    f"{20.0 + i:8.3f}{13.0 + i:8.3f}{52.0 + i:8.3f}  1.00 20.00           C"
    for i in range(1, 12)
) + "\nEND\n"


def _write(tmpdir, name, text):
    path = os.path.join(tmpdir, name)
    with open(path, "w") as f:
        f.write(text)
    return path


def _ligand_sdf(tmpdir, name, translate=(0.0, 0.0, 0.0)):
    """An aspirin conformer, optionally moved into the receptor's frame."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.AddHs(Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O"))
    AllChem.EmbedMolecule(mol, randomSeed=0xF00D)
    conf = mol.GetConformer()
    for i in range(mol.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        conf.SetAtomPosition(i, (p.x + translate[0], p.y + translate[1], p.z + translate[2]))
    path = os.path.join(tmpdir, name)
    w = Chem.SDWriter(path)
    w.write(mol)
    w.close()
    return path


def test_generated_conformer_is_refused(vina_backend):
    """The origin-centred case: a box nowhere near the protein must not run."""
    with tempfile.TemporaryDirectory() as d:
        receptor = _write(d, "rec.pdb", RECEPTOR)
        ligand = _ligand_sdf(d, "gen.sdf")            # untranslated -> origin
        with pytest.raises(RuntimeError, match="not positioned in the receptor"):
            vina_backend.box_args({}, ligand_sdf=ligand, receptor_file=receptor)


def test_carved_ligand_is_accepted(vina_backend):
    """The legitimate case: a ligand already in the receptor's frame boxes fine."""
    with tempfile.TemporaryDirectory() as d:
        receptor = _write(d, "rec.pdb", RECEPTOR)
        ligand = _ligand_sdf(d, "carved.sdf", translate=(25.0, 18.0, 57.0))
        args = vina_backend.box_args({}, ligand_sdf=ligand, receptor_file=receptor)

    flags = dict(zip(args[::2], args[1::2]))
    assert float(flags["--center_x"]) == pytest.approx(25.0, abs=3.0)
    assert float(flags["--center_z"]) == pytest.approx(57.0, abs=3.0)


def test_explicit_box_bypasses_the_check_entirely(vina_backend):
    """An explicit box is the user's call; do not second-guess it."""
    args = vina_backend.box_args({"center": "0,0,0", "size": 20})
    flags = dict(zip(args[::2], args[1::2]))
    assert flags["--center_x"] == "0" and flags["--size_x"] == "20"


def test_autobox_reference_is_not_treated_as_the_input_ligand(vina_backend):
    """A supplied reference is trusted; only the input-ligand fallback is checked."""
    with tempfile.TemporaryDirectory() as d:
        receptor = _write(d, "rec.pdb", RECEPTOR)
        reference = _ligand_sdf(d, "ref.sdf")   # at the origin, but explicitly given
        args = vina_backend.box_args({"autobox_ligand": reference},
                                     ligand_sdf=reference, receptor_file=receptor)
    assert args, "an explicit autobox_ligand must still produce a box"


def test_no_receptor_means_no_check(vina_backend):
    """Callers without a receptor path keep the old behaviour rather than crash."""
    with tempfile.TemporaryDirectory() as d:
        ligand = _ligand_sdf(d, "gen.sdf")
        args = vina_backend.box_args({}, ligand_sdf=ligand)
    assert args, "the box should still be built when no receptor is supplied"
