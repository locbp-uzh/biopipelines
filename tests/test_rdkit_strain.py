"""Conformer-strain tests for the RDKit tool.

The strain core (ligand_utils.conformer_strain) is exercised for real on
in-memory molecules — no cluster, no PDBs. The point under test is the REFERENCE
STATE: strain is measured against a torsion-restrained minimum, not against the
raw pose. A naive E(pose) - E(free minimum) charges every predicted/docked pose a
large constant for its bond-length/angle mismatch with force-field ideals, which
swamps the conformational signal; test_geometry_noise_does_not_inflate_strain is
the regression guard for that.

Config-time tests assert the tool's declared tables and its validation rules.
"""

import pytest

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolTransforms

from biopipelines.ligand_utils import conformer_strain, _rotatable_bonds, _torsion_for_bond

# A flexible molecule MMFF types cleanly: butyl + methoxy on a benzene ring.
FLEXIBLE_SMILES = "CCCCc1ccccc1OC"
# MMFF has no parameters for boron; UFF does.
BORONIC_SMILES = "OB(O)c1ccccc1"


def _embedded(smiles, optimize=True, seed=0xf00d):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    assert AllChem.EmbedMolecule(mol, randomSeed=seed) == 0
    if optimize:
        AllChem.MMFFOptimizeMolecule(mol)
    return mol


def _heavy_torsion(mol):
    """A rotatable torsion whose four atoms are all heavy (so eclipsing it bites)."""
    for i, j in _rotatable_bonds(mol):
        torsion = _torsion_for_bond(mol, i, j)
        if torsion and all(mol.GetAtomWithIdx(a).GetAtomicNum() > 1 for a in torsion):
            return torsion
    raise AssertionError("no all-heavy rotatable torsion in test molecule")


def test_unstrained_conformer_has_zero_strain(record_case):
    mol = _embedded(FLEXIBLE_SMILES)

    e_pose, e_relaxed, strain, engine = conformer_strain(mol)
    record_case(input="MMFF-optimised conformer fed back in",
                expected="strain ~= 0 (|strain| <= 1 kcal/mol)",
                actual=f"strain={strain:.2f} ({engine})")

    # The restrained reference can give a small negative; that is expected.
    assert abs(strain) <= 1.0
    assert engine == "MMFF"


def test_eclipsed_torsion_is_clearly_strained(record_case):
    mol = _embedded(FLEXIBLE_SMILES)
    _, _, relaxed_strain, _ = conformer_strain(mol)

    eclipsed = Chem.Mol(mol)
    rdMolTransforms.SetDihedralDeg(eclipsed.GetConformer(), *_heavy_torsion(eclipsed), 0.0)

    e_pose, e_relaxed, strain, engine = conformer_strain(eclipsed)
    record_case(input="same molecule, one rotatable torsion set to 0 deg (eclipsed)",
                expected="strain clearly positive and >> the unstrained case",
                actual=f"strain={strain:.2f} vs unstrained {relaxed_strain:.2f} ({engine})")

    assert strain > 2.0
    assert strain > relaxed_strain + 2.0
    assert e_pose > e_relaxed


def test_geometry_noise_does_not_inflate_strain(record_case):
    """The reference-state guard: a pose whose bond lengths/angles are perturbed
    but whose CONFORMATION is unchanged must not be charged for the mismatch.
    Scoring the raw pose (the naive reference) charges it tens of kcal/mol."""
    import numpy as np

    mol = _embedded(FLEXIBLE_SMILES)
    noisy = Chem.Mol(mol)
    conf = noisy.GetConformer()
    rng = np.random.default_rng(1)
    for idx in range(noisy.GetNumAtoms()):
        p = conf.GetAtomPosition(idx)
        conf.SetAtomPosition(idx, (p.x + rng.normal(0, 0.03),
                                   p.y + rng.normal(0, 0.03),
                                   p.z + rng.normal(0, 0.03)))

    _, _, strain, _ = conformer_strain(noisy)

    # What the naive reference state would have reported for the same pose.
    props = AllChem.MMFFGetMoleculeProperties(noisy)
    e_raw_pose = AllChem.MMFFGetMoleculeForceField(noisy, props).CalcEnergy()
    free = Chem.Mol(noisy)
    field = AllChem.MMFFGetMoleculeForceField(free, AllChem.MMFFGetMoleculeProperties(free))
    field.Minimize(maxIts=2000)
    naive_strain = e_raw_pose - field.CalcEnergy()

    record_case(input="unchanged conformation, coordinates perturbed by 0.03 A noise",
                expected="strain ~= 0; the naive E(pose)-E(free min) reference inflates it",
                actual=f"strain={strain:.2f} vs naive={naive_strain:.2f}")

    assert abs(strain) <= 1.0
    assert naive_strain > 5 * max(abs(strain), 1.0)


def test_non_convergence_raises_rather_than_scoring(record_case):
    """A minimisation that hits max_iters has not reached a minimum, so the difference
    of the two energies is not a strain. It must raise, not emit a plausible number."""
    mol = _embedded(FLEXIBLE_SMILES, optimize=False)  # unoptimised embed: far from a minimum

    with pytest.raises(RuntimeError, match="did not converge") as exc:
        conformer_strain(mol, max_iters=1)

    converged = conformer_strain(mol)[2]
    record_case(input="max_iters=1 (cannot converge) vs the 2000 default",
                expected="RuntimeError naming non-convergence; the default still converges",
                actual=f"{exc.value} | default strain={converged:.2f}")


def test_uff_fallback_when_mmff_cannot_type(record_case):
    mol = _embedded(BORONIC_SMILES, optimize=False)
    assert AllChem.MMFFGetMoleculeProperties(mol) is None, "test molecule must be MMFF-untypable"

    e_pose, e_relaxed, strain, engine = conformer_strain(mol, ff="auto")
    record_case(input=f"ff='auto' on an MMFF-untypable molecule ({BORONIC_SMILES})",
                expected="falls back to UFF, ff_engine == 'UFF'",
                actual=f"engine={engine}, strain={strain:.2f}")

    assert engine == "UFF"
    assert strain == pytest.approx(e_pose - e_relaxed)

    with pytest.raises(RuntimeError, match="MMFF cannot type"):
        conformer_strain(mol, ff="mmff")


def test_mmff_types_silicon(record_case):
    """Si does not force UFF — MMFF types a Si-aryl (as it does a Si-rhodamine)."""
    mol = _embedded("[SiH3]c1ccccc1", optimize=False, seed=7)
    engine = conformer_strain(mol)[3]
    record_case(input="Si-containing molecule, ff='auto'",
                expected="MMFF (Si is typable; do not assume it forces UFF)",
                actual=engine)
    assert engine == "MMFF"


def _compounds(ids=("aspirin",)):
    from biopipelines.datastream import DataStream
    return DataStream(name="compounds", ids=list(ids), files=[],
                      map_table="compounds.csv", format="csv")


def _poses(ids=("pose_1",)):
    from biopipelines.datastream import DataStream
    return DataStream(name="structures", ids=list(ids), files=["<id>.pdb"],
                      map_table="structures_map.csv", format="pdb")


def test_structures_requires_a_bond_order_template(record_case):
    """Coordinate-only perception mis-assigns conjugated/charged systems, so a
    template is mandatory — no silent fallback."""
    from biopipelines.rdkit_descriptors import RDKit

    with pytest.raises(ValueError, match="bond-order template") as exc:
        RDKit(structures=_poses())
    record_case(input="RDKit(structures=poses) with no smiles= and no compounds=",
                expected="ValueError naming the missing bond-order template",
                actual=str(exc.value))

    # A literal SMILES satisfies it, and so does a compounds stream's smiles column.
    RDKit(structures=_poses(), smiles=FLEXIBLE_SMILES)
    RDKit(structures=_poses(), compounds=_compounds())


def test_strain_table_declared_only_with_structures(record_case):
    from biopipelines.rdkit_descriptors import RDKit, STRAIN_COLUMNS

    tool = RDKit(structures=_poses(), smiles=FLEXIBLE_SMILES, output_folder="out")
    tables = tool.get_output_files()["tables"]
    record_case(input="RDKit(structures=..., smiles=...)",
                expected=("strain table with the documented columns", "no descriptors table"),
                actual=(list(tables.get("strain").info.columns), "descriptors" in tables))

    assert list(tables["strain"].info.columns) == STRAIN_COLUMNS
    assert "descriptors" not in tables


def test_compounds_only_behaviour_unchanged(record_case):
    """Regression: the SMILES-only descriptors path keeps its table and schema."""
    from biopipelines.rdkit_descriptors import RDKit, DEFAULT_DESCRIPTORS

    tool = RDKit(compounds=_compounds(), output_folder="out")
    tables = tool.get_output_files()["tables"]
    columns = list(tables["descriptors"].info.columns)
    record_case(input="RDKit(compounds=library)",
                expected=("descriptors table id|smiles|<defaults>", "no strain table"),
                actual=(columns[:2] + ["..."], "strain" in tables))

    assert columns == ["id", "smiles"] + list(DEFAULT_DESCRIPTORS)
    assert "strain" not in tables


def test_requires_at_least_one_input():
    from biopipelines.rdkit_descriptors import RDKit

    with pytest.raises(ValueError, match="compounds.*and/or structures"):
        RDKit()


def _load_pipe():
    import importlib.util
    import os

    repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    spec = importlib.util.spec_from_file_location(
        "pipe_rdkit_descriptors",
        os.path.join(repo, "pipe_scripts", "pipe_rdkit_descriptors.py"),
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _strain_cfg(tmp_path, structure_ids, compound_rows, provenance):
    """Write the on-disk streams a docking/co-folding step would leave behind:
    a compounds map (id|smiles) and a structures map whose rows carry compounds.id."""
    import json
    import pandas as pd
    from biopipelines.datastream import DataStream

    compounds_map = tmp_path / "compounds.csv"
    pd.DataFrame(compound_rows).to_csv(compounds_map, index=False)
    compounds = DataStream(name="compounds", ids=[r["id"] for r in compound_rows],
                           files=[], map_table=str(compounds_map), format="csv")

    structures_map = tmp_path / "structures_map.csv"
    pd.DataFrame([
        {"id": sid, "file": str(tmp_path / f"{sid}.pdb"), "compounds.id": provenance[sid]}
        for sid in structure_ids
    ]).to_csv(structures_map, index=False)
    structures = DataStream(name="structures", ids=list(structure_ids),
                            files=[str(tmp_path / "<id>.pdb")],
                            map_table=str(structures_map), format="pdb")

    compounds_json = tmp_path / "compounds.json"
    structures_json = tmp_path / "structures.json"
    compounds.save_json(str(compounds_json))
    structures.save_json(str(structures_json))

    cfg = {"compounds_json": str(compounds_json), "structures_json": str(structures_json)}
    return cfg, json.loads(structures_json.read_text())


def test_smiles_resolves_through_compounds_provenance(tmp_path, record_case):
    """A pose id is not a compound id: a co-folded complex is `prot+lig` while the
    compound it carries is `lig`. Resolution goes through the structures map's
    compounds.id provenance (as Gnina does), not an exact-id dict lookup."""
    from biopipelines.biopipelines_io import load_datastream

    mod = _load_pipe()
    cfg, _ = _strain_cfg(
        tmp_path,
        structure_ids=["prot+ligA", "prot+ligB"],
        compound_rows=[{"id": "ligA", "smiles": "CCO"}, {"id": "ligB", "smiles": "CCN"}],
        provenance={"prot+ligA": "ligA", "prot+ligB": "ligB"},
    )
    resolve = mod._smiles_resolver(cfg, load_datastream(cfg["structures_json"]))
    got = {sid: resolve(sid) for sid in ("prot+ligA", "prot+ligB")}

    record_case(input="structures ids prot+ligA/prot+ligB, compounds ids ligA/ligB",
                expected={"prot+ligA": "CCO", "prot+ligB": "CCN"},
                actual=got)
    assert got == {"prot+ligA": "CCO", "prot+ligB": "CCN"}


def test_single_compound_broadcasts_to_every_pose(tmp_path, record_case):
    from biopipelines.biopipelines_io import load_datastream

    mod = _load_pipe()
    cfg, _ = _strain_cfg(
        tmp_path,
        structure_ids=["pose_1", "pose_2"],
        compound_rows=[{"id": "dye", "smiles": "CCO"}],
        provenance={"pose_1": "dye", "pose_2": "dye"},
    )
    resolve = mod._smiles_resolver(cfg, load_datastream(cfg["structures_json"]))
    got = [resolve("pose_1"), resolve("pose_2")]

    record_case(input="one compound, two poses", expected=["CCO", "CCO"], actual=got)
    assert got == ["CCO", "CCO"]
