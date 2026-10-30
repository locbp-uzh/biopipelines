# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Gnina and Vina: which conformer stands for a (protein, ligand) pair, and what `top_only` surfaces.

The structures stream holds one file per pair, so among a pair's conformers one has to stand for it. That
choice was a fixed rule (lowest `pseudo_binding_energy` when every conformer has one, else lowest
`best_vina`); `rank_conformers_by` now names the metric, with the direction that metric needs, and
`top_only=False` surfaces every conformer under its `docking_summary` id instead of choosing.

Nothing here runs gnina or vina; the binary is replaced by a function that writes the scores a test asks for.
"""

import os
import pathlib
import subprocess

import pandas as pd
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from biopipelines import id_patterns
from biopipelines.datastream import DataStream, create_map_table
from biopipelines.gnina import CONFORMER_RANK_METRICS, Gnina, Vina
from tests.test_gnina_ids import PROTEIN_LINE, _declared_stream, pipe_gnina  # noqa: F401 (a fixture)

SEED = 1

# One pair, three conformers, two runs. Per conformer and run, the poses are (vina, cnn_score, cnn_affinity).
# best_vina / mean_vina / best_cnn_score / mean_cnn_affinity / pseudo_binding_energy of each conformer:
#   0: -8.0 / -7.0 / 0.20 / 4.0 / -2.0   1: -7.5 / -7.5 / 0.95 / 5.0 / -4.5   2: -5.0 / -5.0 / 0.60 / 9.0 / -5.0
POSES = {
    0: [[(-8.0, 0.20, 4.0)], [(-6.0, 0.20, 4.0)]],
    1: [[(-7.5, 0.95, 5.0)], [(-7.5, 0.95, 5.0)]],
    2: [[(-5.0, 0.60, 9.0)], [(-5.0, 0.60, 9.0)]],
}
ENERGIES = {0: 6.0, 1: 3.0, 2: 0.0}

# The conformer each metric picks, in the direction that metric needs: no two directions agree on all of them.
WINNER = {"auto": 2, "pseudo_binding_energy": 2, "best_vina": 0, "mean_vina": 1,
          "best_cnn_score": 1, "mean_cnn_affinity": 2}


def _write_sdf(path, poses):
    mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
    AllChem.EmbedMolecule(mol, randomSeed=7)
    writer = Chem.SDWriter(str(path))
    for vina, cnn_score, cnn_affinity in poses:
        mol.SetProp("minimizedAffinity", str(vina))
        mol.SetProp("CNNscore", str(cnn_score))
        mol.SetProp("CNNaffinity", str(cnn_affinity))
        writer.write(mol)
    writer.close()


class _FakeGnina:
    """Stands in for the binary: the poses of run r of conformer c are POSES[c][r], whichever protein it docks into."""

    def __call__(self, cmd, **kwargs):
        conformer = int(pathlib.Path(cmd[cmd.index("-l") + 1]).stem.rsplit("_conf", 1)[1])
        run = int(cmd[cmd.index("--seed") + 1]) - SEED
        _write_sdf(cmd[cmd.index("-o") + 1], POSES[conformer][run])
        return subprocess.CompletedProcess(cmd, 0, "", "")


def _dock(pipe, monkeypatch, tmp_path, proteins, ligands, energies=None, **config_extra):
    """Dock `proteins` x `ligands` (one id each, three conformers each) and return what aggregation needs.

    `energies` maps (ligand id, conformer) to its energy and defaults to ENERGIES for every ligand.
    """
    prep = tmp_path / "prepared"
    prep.mkdir()
    prepared = {}
    for protein_id in proteins:
        pdb = prep / f"{protein_id}.pdb"
        pdb.write_text(PROTEIN_LINE + "END\n")
        prepared[protein_id] = {"protein": str(pdb), "crystal_ligand": None}
    energies = energies or {(ligand, k): ENERGIES[k] for ligand in ligands for k in POSES}
    conformers = [{"ligand_id": ligand, "conformer_id": k, "sdf_file": str(tmp_path / f"{ligand}_conf{k}.sdf"),
                   "energy": energies[(ligand, k)]} for ligand in ligands for k in POSES]
    config = {
        "output_folder": str(tmp_path / "execution"), "backend": "gnina", "gnina_binary": "gnina",
        "exhaustiveness": 8, "num_modes": 2, "num_runs": 2, "seed": SEED, "cnn_scoring": "rescore",
        "cnn_score_threshold": 0.0, "scoring": "vina", "box": {"center": "0,0,0", "size": 20, "autobox_add": 4.0},
        "extra_args": [], "dock_timeout": 0, "rmsd_threshold": 2.0,
        "best_poses_dir": str(tmp_path / "structures"),
        "docking_results_csv": str(tmp_path / "tables" / "docking_results.csv"),
        "docking_summary_csv": str(tmp_path / "tables" / "docking_summary.csv"),
        "structures_map": str(tmp_path / "structures_map.csv"),
        "top_only": True, "rank_conformers_by": "auto", **config_extra,
    }
    (tmp_path / "tables").mkdir()
    monkeypatch.setattr(pipe.subprocess, "run", _FakeGnina())
    all_poses, _failures = pipe.run_docking(prepared, conformers, config)
    consistency = pipe.analyze_pose_consistency(all_poses, config)
    return all_poses, conformers, consistency, config, prepared


def _aggregate(pipe, args):
    pipe.aggregate_results(*args)
    config = args[3]
    return {"summary": pd.read_csv(config["docking_summary_csv"], dtype={"id": str}),
            "map": pd.read_csv(config["structures_map"], dtype={"id": str})}


# ── which conformer stands for the pair ───────────────────────────────────────

@pytest.mark.parametrize("metric", sorted(WINNER))
def test_each_metric_picks_its_conformer(pipe_gnina, monkeypatch, tmp_path, metric, record_case):
    out = _aggregate(pipe_gnina, _dock(pipe_gnina, monkeypatch, tmp_path, ["p"], ["l"], rank_conformers_by=metric))
    stream_file = os.path.join(str(tmp_path / "structures"), "p+l_best.pdb")
    by_conformer = out["summary"].set_index("conformer_id")["best_pose_file"]
    record_case(input=f"rank_conformers_by={metric!r}", expected=f"conformer {WINNER[metric]}",
                actual=[k for k in POSES if by_conformer[k] == stream_file])
    assert list(out["map"]["id"]) == ["p+l"]
    assert by_conformer[WINNER[metric]] == stream_file
    assert all(by_conformer[k] != stream_file for k in POSES if k != WINNER[metric])


def test_the_directions_are_lower_for_energies_and_higher_for_cnn_scores():
    assert {name: CONFORMER_RANK_METRICS[name] for name in WINNER if name != "auto"} == {
        "pseudo_binding_energy": "lower", "best_vina": "lower", "mean_vina": "lower",
        "best_cnn_score": "higher", "mean_cnn_affinity": "higher"}


def test_a_higher_is_better_metric_does_not_pick_the_lowest(pipe_gnina):
    rows = [{"id": f"p+l_conf{k}", "structures.id": "p", "compounds.id": "l", "best_cnn_score": score,
             "best_vina": -5.0, "pseudo_binding_energy": ""} for k, score in enumerate((0.2, 0.95, 0.6))]
    winners, unrankable = pipe_gnina._pair_winners(rows, {row["id"]: {} for row in rows}, "best_cnn_score")
    assert winners == {"p+l": "p+l_conf1"} and not unrankable


def test_a_tie_goes_to_the_first_conformer(pipe_gnina):
    rows = [{"id": f"p+l_conf{k}", "structures.id": "p", "compounds.id": "l", "best_vina": -6.0,
             "pseudo_binding_energy": ""} for k in range(3)]
    winners, _ = pipe_gnina._pair_winners(rows, {row["id"]: {} for row in rows}, "best_vina")
    assert winners == {"p+l": "p+l_conf0"}


def test_a_metric_missing_for_a_conformer_fails_loudly_and_keeps_the_rest(pipe_gnina, monkeypatch, tmp_path):
    """Ranking by an energy that one conformer lacks must not fall back to another metric for that pair."""
    energies = {(ligand, k): ENERGIES[k] for ligand in ("l_a", "l_b") for k in POSES}
    energies[("l_b", 1)] = ""
    args = _dock(pipe_gnina, monkeypatch, tmp_path, ["p"], ["l_a", "l_b"], energies=energies,
                 rank_conformers_by="pseudo_binding_energy")
    with pytest.raises(ValueError, match=r"pseudo_binding_energy.*p\+l_b.*p\+l_b_conf1"):
        pipe_gnina.aggregate_results(*args)
    config = args[3]
    assert list(pd.read_csv(config["structures_map"], dtype={"id": str})["id"]) == ["p+l_a"]
    assert len(pd.read_csv(config["docking_summary_csv"])) == 6, "the summary of every conformer is kept"
    assert os.path.exists(tmp_path / "execution" / "conformer_poses" / "p+l_b_conf1_best.pdb")


def test_auto_falls_back_to_best_vina_where_an_energy_is_missing(pipe_gnina, monkeypatch, tmp_path):
    energies = {("l", k): ENERGIES[k] for k in POSES}
    energies[("l", 1)] = ""
    out = _aggregate(pipe_gnina, _dock(pipe_gnina, monkeypatch, tmp_path, ["p"], ["l"], energies=energies))
    stream_file = os.path.join(str(tmp_path / "structures"), "p+l_best.pdb")
    assert out["summary"].set_index("conformer_id")["best_pose_file"][0] == stream_file


# ── top_only ──────────────────────────────────────────────────────────────────

def test_top_only_false_declares_every_conformer_and_writes_it(
        pipe_gnina, monkeypatch, tmp_path, local_config, isolated_cwd, new_pipeline, record_case):
    proteins, ligands = ["p_1", "p2"], ["l_1"]
    declared, _, _ = _declared_stream(new_pipeline, proteins, ligands, top_only=False)
    assert list(declared.ids) == ["p_1+l_1[_conf<#>]", "p2+l_1[_conf<#>]"]

    args = _dock(pipe_gnina, monkeypatch, tmp_path, proteins, ligands, top_only=False,
                 best_poses_dir=os.path.dirname(declared.files[0]))
    out = _aggregate(pipe_gnina, args)
    written = list(out["map"]["id"])
    record_case(input="2 proteins x 1 ligand x 3 conformers, top_only=False",
                expected=list(declared.ids), actual=written)
    assert written == [f"{p}+l_1_conf{k}" for p in proteins for k in POSES]
    assert len(set(written)) == len(written)
    assert id_patterns.select_ids(list(declared.ids), written) == written, "a written id is not declared"
    assert list(out["summary"]["id"].sort_values()) == sorted(written), "stream and summary ids must be one set"
    for stream_id in written:
        assert os.path.exists(declared.files[0].replace("<id>", stream_id)), stream_id
    assert not os.path.exists(tmp_path / "execution" / "conformer_poses")
    assert set(out["summary"]["best_pose_file"]) == {declared.files[0].replace("<id>", i) for i in written}


def test_top_only_true_still_declares_the_bare_pair(local_config, isolated_cwd, new_pipeline):
    declared, _, _ = _declared_stream(new_pipeline, ["p_1", "p2"], ["l_1"], top_only=True)
    assert list(declared.ids) == ["p_1+l_1", "p2+l_1"]


# ── construction ──────────────────────────────────────────────────────────────

def _streams(tmp_path):
    s_map, c_map = tmp_path / "s.csv", tmp_path / "c.csv"
    create_map_table(str(s_map), ids=["p"], files=["/x/<id>.pdb"])
    create_map_table(str(c_map), ids=["l"])
    return (DataStream(name="structures", ids=["p"], files=["/x/<id>.pdb"], map_table=str(s_map), format="pdb"),
            DataStream(name="compounds", ids=["l"], files=[], map_table=str(c_map), format="csv"))


def _tool(tmp_path, cls=Gnina, **kwargs):
    structures, compounds = _streams(tmp_path)
    tool = cls(structures=structures, compounds=compounds, **kwargs)
    tool.validate_params()
    return tool


@pytest.fixture
def construct(tmp_path, local_config, isolated_cwd):
    return lambda cls=Gnina, **kwargs: _tool(tmp_path, cls, **kwargs)


def test_top_only_must_be_a_bool(construct):
    with pytest.raises(ValueError, match="top_only must be a bool"):
        construct(top_only="yes")


@pytest.mark.parametrize("kwargs, message", [
    ({"rank_conformers_by": "best_vinna"}, "rank_conformers_by must be 'auto' or one of"),
    ({"rank_conformers_by": "conformer_energy"}, "rank_conformers_by must be 'auto' or one of"),
    ({"rank_conformers_by": "best_vina", "top_only": False}, "only applies with top_only=True"),
    ({"rank_conformers_by": "best_cnn_score", "cnn_scoring": "none"}, "needs CNN scores"),
    ({"rank_conformers_by": "pseudo_binding_energy"}, "needs conformer energies"),
    ({"mode": "score", "top_only": False}, "docking-only"),
    ({"mode": "minimize", "rank_conformers_by": "best_vina"}, "docking-only"),
])
def test_a_ranking_that_cannot_be_honored_is_refused(construct, kwargs, message):
    with pytest.raises(ValueError, match=message):
        construct(**kwargs)


def test_vina_has_no_cnn_column_to_rank_by(construct):
    with pytest.raises(ValueError, match="Vina has no CNN"):
        construct(Vina, rank_conformers_by="mean_cnn_affinity")


@pytest.mark.parametrize("cls, metric", [(Gnina, "best_cnn_score"), (Vina, "best_vina")])
def test_a_named_metric_that_exists_is_accepted(construct, cls, metric):
    assert construct(cls, rank_conformers_by=metric).rank_conformers_by == metric


def test_pseudo_binding_energy_is_accepted_once_energies_exist(construct):
    assert construct(rank_conformers_by="pseudo_binding_energy", generate_conformers=True)


def test_the_choices_reach_the_pipe_script_through_its_config(local_config, isolated_cwd, new_pipeline):
    import json
    _declared_stream(new_pipeline, ["p"], ["l"], top_only=False)
    (config_path,) = pathlib.Path.cwd().rglob("gnina_config.json")
    with open(config_path) as handle:
        written = json.load(handle)
    assert written["top_only"] is False and written["rank_conformers_by"] == "auto"
