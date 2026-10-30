"""Gnina ids: axes join with `+`, conformer, run and pose are multiplier suffixes, no two rows collide.

The wrapper and the pipe script each joined protein and ligand with `_`, so a ligand id containing `_`
could merge two pairs into one id (`a_b`+`c` and `a`+`b_c` both became `a_b_c`) and they shared one
scratch SDF and one `<id>_best.pdb`. Two conformers of one pair also shared the summary id and the best
file, and their pose ids repeated because the conformer was not part of them. The declared ids and the
ids the runtime writes now both come from the shared axis composition, and the conformer, the run and
the pose each append one suffix to their parent.
"""

import importlib.util
import os
import pathlib
import subprocess

import pandas as pd
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

REPO = pathlib.Path(__file__).resolve().parent.parent

PROTEIN_LINE = "ATOM      1  CA  GLY A   1       0.000   0.000   0.000  1.00  0.00           C\n"


@pytest.fixture(scope="module")
def pipe_gnina():
    """Load the pipe script directly; pipe_scripts is not an importable package."""
    spec = importlib.util.spec_from_file_location("pipe_gnina", REPO / "pipe_scripts" / "pipe_gnina.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _write_sdf(path, affinities):
    """One pose per affinity, ethanol each, carrying the score properties gnina writes."""
    mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
    AllChem.EmbedMolecule(mol, randomSeed=7)
    writer = Chem.SDWriter(str(path))
    for affinity in affinities:
        mol.SetProp("minimizedAffinity", str(affinity))
        mol.SetProp("CNNscore", "0.5")
        mol.SetProp("CNNaffinity", "5.0")
        writer.write(mol)
    writer.close()


class _FakeGnina:
    """Stands in for the gnina binary: one SDF per call, scored by (protein id, conformer file stem)."""

    def __init__(self, scores):
        self.scores = scores
        self.outputs = []

    def __call__(self, cmd, **kwargs):
        receptor = pathlib.Path(cmd[cmd.index("-r") + 1]).stem
        ligand = pathlib.Path(cmd[cmd.index("-l") + 1]).stem
        out = cmd[cmd.index("-o") + 1]
        self.outputs.append(out)
        _write_sdf(out, self.scores[(receptor, ligand)])
        return subprocess.CompletedProcess(cmd, 0, "", "")


def _dock(pipe_gnina, monkeypatch, tmp_path, proteins, conformers, scores, num_runs=1, best_dir=None):
    """Run the docking, consistency and aggregation stages of the pipe script on synthetic inputs.

    ``conformers`` is a list of (ligand_id, conformer_id, energy); ``scores`` maps
    (protein_id, "<ligand_id>_conf<conformer_id>") to the poses' affinities.
    """
    prep = tmp_path / "prepared"
    prep.mkdir()
    prepared = {}
    for protein_id in proteins:
        pdb = prep / f"{protein_id}.pdb"
        pdb.write_text(PROTEIN_LINE + "END\n")
        prepared[protein_id] = {"protein": str(pdb), "crystal_ligand": None}

    conformer_records = [
        {"ligand_id": ligand_id, "conformer_id": conformer_id,
         "sdf_file": str(tmp_path / f"{ligand_id}_conf{conformer_id}.sdf"), "energy": energy}
        for ligand_id, conformer_id, energy in conformers
    ]

    best_dir = str(best_dir or tmp_path / "structures")
    config = {
        "output_folder": str(tmp_path / "execution"), "backend": "gnina", "gnina_binary": "gnina",
        "exhaustiveness": 8, "num_modes": 2, "num_runs": num_runs, "seed": 1,
        "cnn_scoring": "rescore", "cnn_score_threshold": 0.0, "scoring": "vina",
        "box": {"center": "0,0,0", "size": 20, "autobox_add": 4.0}, "extra_args": [],
        "dock_timeout": 0, "rmsd_threshold": 2.0, "best_poses_dir": best_dir,
        "docking_results_csv": str(tmp_path / "docking_results.csv"),
        "docking_summary_csv": str(tmp_path / "docking_summary.csv"),
        "structures_map": str(tmp_path / "structures_map.csv"),
    }

    fake = _FakeGnina(scores)
    monkeypatch.setattr(pipe_gnina.subprocess, "run", fake)
    all_poses, _failures = pipe_gnina.run_docking(prepared, conformer_records, config)
    consistency = pipe_gnina.analyze_pose_consistency(all_poses, config)
    pipe_gnina.aggregate_results(all_poses, conformer_records, consistency, config, prepared)

    return {
        "results": pd.read_csv(config["docking_results_csv"], dtype={"id": str}),
        "summary": pd.read_csv(config["docking_summary_csv"], dtype={"id": str}),
        "map": pd.read_csv(config["structures_map"], dtype={"id": str}),
        "outputs": fake.outputs,
        "best_dir": best_dir,
    }


def _declared_stream(new_pipeline, protein_ids, ligand_ids, cls="Gnina", **kwargs):
    """The structures stream a Gnina or Vina step declares at configuration time."""
    from biopipelines import gnina as gnina_module
    from biopipelines.mock import Mock

    pipeline = new_pipeline("gnina_ids")
    with pipeline:
        struct = Mock(ids=list(protein_ids), streams={"structures": {"format": "pdb", "file": "<id>.pdb"}},
                      map_table_strategy="config")
        comp = Mock(ids=list(ligand_ids), streams={"compounds": {"format": "csv", "file": "<id>.csv"}},
                    map_table_strategy="config")
        out = getattr(gnina_module, cls)(structures=struct.streams.structures,
                                         compounds=comp.streams.compounds, **kwargs)
    return out.streams.structures, struct.streams.structures, comp.streams.compounds


# ── configuration time ──────────────────────────────────────────────────────

@pytest.mark.parametrize("cls", ["Gnina", "Vina"])
def test_declared_ids_are_the_axes_joined_with_plus(local_config, isolated_cwd, new_pipeline, cls, record_case):
    from biopipelines.combinatorics import predict_output_ids

    declared, structures, compounds = _declared_stream(new_pipeline, ["p_1", "p2"], ["l_1", "c2"], cls=cls)
    shared = predict_output_ids(structures=(structures, "structures"), compounds=(compounds, "compounds"))
    expected = ["p_1+l_1", "p_1+c2", "p2+l_1", "p2+c2"]
    record_case(input="structures p_1,p2 x compounds l_1,c2", expected=expected, actual=list(declared.ids))

    assert list(declared.ids) == expected
    assert list(declared.ids) == shared, "the wrapper must declare what the shared composition produces"


@pytest.mark.parametrize("mode", ["score", "minimize"])
def test_no_search_modes_keep_the_complex_id(local_config, isolated_cwd, new_pipeline, mode):
    """One row per complex: the stream is 1:1 with its input, so it keeps the input's id."""
    declared, _, _ = _declared_stream(new_pipeline, ["s_1", "s2"], ["l_1"], mode=mode)

    assert list(declared.ids) == ["s_1", "s2"]
    assert declared.files[0].endswith(f"<id>_{mode}.pdb")


# ── runtime ─────────────────────────────────────────────────────────────────

def test_pairs_that_underscore_joining_would_merge_stay_distinct(pipe_gnina, monkeypatch, tmp_path, record_case):
    proteins = ["a_b", "a"]
    ligands = ["c", "b_c"]
    conformers = [(lig, 0, "") for lig in ligands]
    scores = {(prot, f"{lig}_conf0"): [-7.0 - 0.5 * i]
              for i, (prot, lig) in enumerate((p, l) for p in proteins for l in ligands)}

    out = _dock(pipe_gnina, monkeypatch, tmp_path, proteins, conformers, scores)

    expected = ["a_b+c", "a_b+b_c", "a+c", "a+b_c"]
    record_case(input="proteins a_b,a x ligands c,b_c", expected=expected, actual=list(out["map"]["id"]))
    assert list(out["map"]["id"]) == expected
    assert len(set(out["map"]["file"])) == 4, "two pairs wrote the same best pose file"
    assert all(os.path.exists(path) for path in out["map"]["file"])
    assert len(set(out["outputs"])) == len(out["outputs"]) == 4, "two pairs shared one docking scratch file"
    assert out["results"]["id"].is_unique
    assert out["summary"]["id"].is_unique and len(out["summary"]) == 4


@pytest.mark.parametrize("energies, winner", [
    ((None, None), 0),
    ((3.0, 0.0), 1),
])
def test_conformers_get_distinct_ids_and_the_pair_keeps_one_file(
        pipe_gnina, monkeypatch, tmp_path, energies, winner, record_case):
    """Conformer 0 has the better Vina score; strain flips the ranking when energies are known."""
    conformers = [("l_1", 0, "" if energies[0] is None else energies[0]),
                  ("l_1", 1, "" if energies[1] is None else energies[1])]
    scores = {("p", "l_1_conf0"): [-8.0, -6.0], ("p", "l_1_conf1"): [-7.0, -5.0]}

    out = _dock(pipe_gnina, monkeypatch, tmp_path, ["p"], conformers, scores, num_runs=2)

    pose_ids = sorted(f"p+l_1_conf{c}_r{r}_p{k}" for c in (0, 1) for r in (0, 1) for k in (0, 1))
    record_case(input="one pair, conformers 0 and 1, 2 runs, 2 poses", expected=pose_ids,
                actual=sorted(out["results"]["id"]))
    assert sorted(out["results"]["id"]) == pose_ids
    assert sorted(out["summary"]["id"]) == ["p+l_1_conf0", "p+l_1_conf1"]
    assert list(out["summary"]["structures.id"]) == ["p", "p"]
    assert list(out["summary"]["compounds.id"]) == ["l_1", "l_1"]

    stream_file = os.path.join(out["best_dir"], "p+l_1_best.pdb")
    assert list(out["map"]["id"]) == ["p+l_1"], "the stream holds one row per pair"
    assert list(out["map"]["file"]) == [stream_file]
    by_conformer = out["summary"].set_index("conformer_id")["best_pose_file"]
    assert by_conformer[winner] == stream_file, "the pair's file must hold the winning conformer's pose"
    loser_file = by_conformer[1 - winner]
    assert loser_file != stream_file and os.path.exists(loser_file), \
        "the other conformer's pose must keep a file of its own"


def test_declared_ids_equal_the_ids_the_runtime_writes(
        pipe_gnina, monkeypatch, tmp_path, local_config, isolated_cwd, new_pipeline, record_case):
    proteins = ["p_1", "p2"]
    ligands = ["l_1", "c2"]
    from biopipelines.combinatorics import predict_output_ids

    declared, structures, compounds = _declared_stream(new_pipeline, proteins, ligands)
    shared = predict_output_ids(structures=(structures, "structures"), compounds=(compounds, "compounds"))

    conformers = [(lig, 0, "") for lig in ligands]
    scores = {(prot, f"{lig}_conf0"): [-7.0] for prot in proteins for lig in ligands}
    out = _dock(pipe_gnina, monkeypatch, tmp_path, proteins, conformers, scores,
                best_dir=os.path.dirname(declared.files[0]))

    record_case(input="declared at configuration vs written at runtime",
                expected=list(declared.ids), actual=list(out["map"]["id"]))
    assert list(declared.ids) == shared
    assert list(out["map"]["id"]) == list(declared.ids)
    for declared_id in declared.ids:
        assert os.path.exists(declared.files[0].replace("<id>", declared_id)), \
            f"declared file for {declared_id} was not written"


def test_a_pose_id_is_its_conformer_run_and_pose_appended_to_the_pair(pipe_gnina, tmp_path):
    """Provenance columns name the parts; the id is the pair with one suffix per multiplier."""
    sdf = tmp_path / "run.sdf"
    _write_sdf(sdf, [-7.0, -6.0])

    poses = pipe_gnina._parse_gnina_output(str(sdf), "prot", "lig_2", 3, 4, 0.0)

    assert [p["id"] for p in poses] == ["prot+lig_2_conf3_r4_p0", "prot+lig_2_conf3_r4_p1"]
    assert {(p["structures.id"], p["compounds.id"], p["conformer_id"], p["run"]) for p in poses} == \
        {("prot", "lig_2", 3, 4)}
