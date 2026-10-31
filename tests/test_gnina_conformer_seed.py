"""Gnina `seed` is the base seed of the conformer embedding as well as of the docking runs.

Both ETKDGv3 embeddings hard-coded `randomSeed = 42` while `seed` reached only `--seed seed + run`, so
`Gnina(generate_conformers=True, seed=7)` and `seed=99` docked the very same conformers. The default seed is 42,
so a run that never set `seed` embeds exactly what it embedded before.
"""

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from biopipelines.datastream import DataStream

from tests.test_gnina_ids import pipe_gnina  # noqa: F401 (a fixture)

FLEXIBLE = "CCCCCCCC(=O)NCCCC(=O)O"


def _compounds(tmp_path):
    (tmp_path / "compounds.csv").write_text(f"id,code,smiles\nlig,LIG,{FLEXIBLE}\n")
    return DataStream(name="compounds", ids=["lig"], files=[], map_table=str(tmp_path / "compounds.csv"),
                      format="csv")


def _coordinates(sdf):
    mol = next(iter(Chem.SDMolSupplier(str(sdf), removeHs=False)))
    return [float(c) for atom in mol.GetConformer().GetPositions() for c in atom]


def _conformers(pipe, tmp_path, seed, generate, tag):
    out = tmp_path / f"{tag}_seed{seed}"
    out.mkdir()
    config = {"seed": seed, "generate_conformers": generate, "num_conformers": 4, "energy_window": 50.0,
              "conformer_rmsd": 0.0, "conformer_energies_ref": None}
    return pipe.prepare_conformers(_compounds(tmp_path), str(out), config)


@pytest.mark.parametrize("generate", [False, True], ids=["single embedding", "generated conformers"])
def test_a_different_seed_embeds_different_coordinates(pipe_gnina, tmp_path, generate):
    first = _conformers(pipe_gnina, tmp_path, 7, generate, "a")
    second = _conformers(pipe_gnina, tmp_path, 99, generate, "b")
    assert first and second
    spread = max(abs(a - b) for a, b in zip(_coordinates(first[0]["sdf_file"]), _coordinates(second[0]["sdf_file"])))
    assert spread > 0.05


@pytest.mark.parametrize("generate", [False, True], ids=["single embedding", "generated conformers"])
def test_the_same_seed_embeds_the_same_coordinates(pipe_gnina, tmp_path, generate):
    first = _conformers(pipe_gnina, tmp_path, 7, generate, "a")
    second = _conformers(pipe_gnina, tmp_path, 7, generate, "b")
    assert len(first) == len(second)
    for one, other in zip(first, second):
        assert _coordinates(one["sdf_file"]) == pytest.approx(_coordinates(other["sdf_file"]), abs=1e-3)


def test_the_default_seed_embeds_what_the_fixed_seed_embedded(pipe_gnina, tmp_path):
    """Seed 42 is the default, so the coordinates a run got before `seed` reached the embedding are unchanged."""
    mol = Chem.AddHs(Chem.MolFromSmiles(FLEXIBLE))
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    AllChem.EmbedMolecule(mol, params)
    AllChem.MMFFOptimizeMolecule(mol)
    expected = [float(c) for atom in mol.GetConformer().GetPositions() for c in atom]

    (conformer,) = _conformers(pipe_gnina, tmp_path, 42, False, "a")
    assert _coordinates(conformer["sdf_file"]) == pytest.approx(expected, abs=1e-3)
