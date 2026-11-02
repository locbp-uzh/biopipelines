# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Gnina and Vina docking: a pair that fails keeps a NaN row of `docking_summary`, under its pair id.

`docking_summary` was keyed on conformers, so a (protein, ligand) pair whose receptor could not be prepared, whose
ligand gave no conformer, that posed nothing or whose `rank_conformers_by` metric was absent had no row there and was
recorded only in `missing`; a `Panda` merge of the table against the one that fed it came back shorter. Every pair that
entered now has a row. A pair that failed has one NaN row: its id is the pair id, `structures.id` and `compounds.id`
name the parents, `conformer_id` and every score are empty. The NaN row is written after the conformers are ranked, so
it can never stand for a pair.

Nothing here runs gnina or vina: the binary is a function and the pipe script runs on the files it writes.
"""

import json
import pathlib
import re
import subprocess
import sys

import pandas as pd
import pytest

from biopipelines import gnina as gnina_module
from biopipelines.mock import Mock
from tests._helpers import (
    GninaBinary, configured_gnina_step, fake_the_gnina_tools, load_pipe_script, materialize_gnina_inputs,
    write_pose_sdf,
)

PAIRS = {("good", "lig"), ("bad", "lig")}
KEY = ["structures.id", "compounds.id"]


def _empty_binary(cmd, **kwargs):
    pathlib.Path(cmd[cmd.index("-o") + 1]).write_text("")
    return subprocess.CompletedProcess(cmd, 0, "", "")


def _vina_only_binary(cmd, **kwargs):
    """A binary whose poses carry a Vina score and no CNN score, for the receptor named `bad`; every other pose is whole."""
    path = pathlib.Path(cmd[cmd.index("-o") + 1])
    write_pose_sdf(path, [(-7.0, 0.8, 5.0)])
    if pathlib.Path(cmd[cmd.index("-r") + 1]).stem.startswith("bad"):
        path.write_text(re.sub(r">  <CNN\w+>[^\n]*\n[^\n]*\n\n", "", path.read_text()))
    return subprocess.CompletedProcess(cmd, 0, "", "")


def _dock(pipe, monkeypatch, new_pipeline, cls_name="Gnina", binary=None, smiles=None, **config_edits):
    """Run the docking step of two complexes, `good` and `bad`, against one ligand; return (summary, exit code, config)."""
    config_path, config, _script, _step = configured_gnina_step(new_pipeline, cls_name, "docking", False)
    materialize_gnina_inputs(config, False)
    if smiles is not None:
        with open(config["compounds_json"]) as handle:
            pathlib.Path(json.load(handle)["map_table"]).write_text(f"id,smiles,code\nlig,{smiles},LIG\n")
    if config_edits:
        config.update(config_edits)
        config_path.write_text(json.dumps(config))
    fake_the_gnina_tools(pipe, monkeypatch, cls_name, "docking")
    if binary is not None:
        monkeypatch.setattr(pipe.subprocess, "run", binary)
    monkeypatch.setattr(sys, "argv", ["pipe_gnina.py", str(config_path)])
    code = 0
    try:
        pipe.main()
    except SystemExit as stopped:
        code = stopped.code
    return pd.read_csv(config["docking_summary_csv"], dtype={"id": str}), code, config


def _nan_rows(summary):
    """The rows that stand for a failed pair: no conformer."""
    return summary[summary["conformer_id"].isna()]


def _score_columns(summary):
    return [c for c in summary.columns if c not in ("id", "structures.id", "compounds.id", "conformer_id")]


@pytest.mark.parametrize("cls_name", ["Gnina", "Vina"])
def test_a_pair_whose_receptor_fails_keeps_a_nan_row_beside_the_pair_that_docked(
        pipe_gnina, monkeypatch, local_config, isolated_cwd, new_pipeline, cls_name):
    summary, code, config = _dock(pipe_gnina, monkeypatch, new_pipeline, cls_name)

    assert code == 0
    assert set(zip(summary["structures.id"], summary["compounds.id"])) == PAIRS, "every pair that entered has a row"
    assert list(summary["id"]) == ["good+lig_conf0", "bad+lig"]
    bad = _nan_rows(summary)
    assert list(bad["id"]) == ["bad+lig"], "the id is the pair id the failure row of missing carries"
    assert (bad.iloc[0]["structures.id"], bad.iloc[0]["compounds.id"]) == ("bad", "lig")
    assert bad[_score_columns(summary)].isna().all(axis=None), "every score is NaN, best_pose_file included"
    assert summary.loc[summary["id"] == "good+lig_conf0", "best_vina"].notna().all()
    local = pd.read_csv(config["local_missing_csv"], dtype=str)
    assert list(local["id"]) == ["bad+lig"] and list(local["kind"]) == ["failure"], "the failure is still a row of missing"
    assert list(pd.read_csv(config["structures_map"])["id"]) == ["good+lig"], "the NaN row adds nothing to the stream"


def test_a_ligand_that_cannot_be_prepared_keeps_a_nan_row_of_every_pair_though_the_step_fails(
        pipe_gnina, monkeypatch, local_config, isolated_cwd, new_pipeline):
    summary, code, _config = _dock(pipe_gnina, monkeypatch, new_pipeline, binary=GninaBinary(), smiles="not a smiles")

    assert code == 1, "no pair gave a structure"
    assert list(summary["id"]) == ["good+lig", "bad+lig"]
    assert len(_nan_rows(summary)) == 2 and summary[_score_columns(summary)].isna().all(axis=None)


def test_a_pair_that_posed_nothing_has_a_header_and_a_nan_row(
        pipe_gnina, monkeypatch, local_config, isolated_cwd, new_pipeline):
    summary, code, config = _dock(pipe_gnina, monkeypatch, new_pipeline, binary=_empty_binary)

    assert code == 1
    assert list(summary["id"]) == ["good+lig", "bad+lig"]
    assert list(summary.columns[:4]) == ["id", "structures.id", "compounds.id", "conformer_id"]
    assert pd.read_csv(config["docking_results_csv"]).empty


def test_a_pair_whose_rank_metric_is_absent_keeps_its_conformer_rows_and_gets_a_nan_row(
        pipe_gnina, monkeypatch, local_config, isolated_cwd, new_pipeline):
    summary, code, config = _dock(pipe_gnina, monkeypatch, new_pipeline, binary=_vina_only_binary,
                                  rank_conformers_by="best_cnn_score")

    assert code == 0
    assert list(summary["id"]) == ["good+lig_conf0", "bad+lig_conf0", "bad+lig"]
    assert list(_nan_rows(summary)["id"]) == ["bad+lig"], "the posed conformer keeps its own row beside the pair's NaN row"
    assert summary.loc[summary["id"] == "bad+lig_conf0", "best_vina"].notna().all()
    assert list(pd.read_csv(config["structures_map"])["id"]) == ["good+lig"], "no conformer stands for bad+lig"
    local = pd.read_csv(config["local_missing_csv"], dtype=str)
    assert list(local["id"]) == ["bad+lig"] and "best_cnn_score" in local.iloc[0]["cause"]


def test_top_only_false_adds_the_nan_row_only_for_the_pair_that_failed(
        pipe_gnina, monkeypatch, local_config, isolated_cwd, new_pipeline):
    summary, code, config = _dock(pipe_gnina, monkeypatch, new_pipeline, top_only=False)

    assert code == 0
    assert list(summary["id"]) == ["good+lig_conf0", "bad+lig"]
    assert list(pd.read_csv(config["structures_map"])["id"]) == ["good+lig_conf0"]


@pytest.mark.parametrize("cls_name", ["Gnina", "Vina"])
def test_the_table_has_the_columns_the_wrapper_declares(
        pipe_gnina, monkeypatch, local_config, isolated_cwd, new_pipeline, cls_name):
    summary, _code, _config = _dock(pipe_gnina, monkeypatch, new_pipeline, cls_name)
    with new_pipeline("gnina_columns"):
        structures = Mock(ids=["good"], streams={"structures": {"format": "pdb", "file": "<id>.pdb"}},
                          map_table_strategy="config")
        compounds = Mock(ids=["lig"], streams={"compounds": {"format": "csv", "file": "<id>.csv"}},
                         map_table_strategy="config")
        out = getattr(gnina_module, cls_name)(structures=structures.streams.structures,
                                              compounds=compounds.streams.compounds)
    declared = out.tables._tables["docking_summary"].info.columns

    assert list(summary.columns) == list(declared)
    assert ("best_cnn_score" in declared) == (cls_name == "Gnina")


def test_the_nan_row_is_never_the_best_row_of_a_consumer(
        pipe_gnina, monkeypatch, local_config, isolated_cwd, new_pipeline):
    summary, _code, config = _dock(pipe_gnina, monkeypatch, new_pipeline)
    assert len(_nan_rows(summary)) == 1
    panda = load_pipe_script("pipe_panda_summary_rows", "pipe_panda.py")
    table = panda.read_table_csv(config["docking_summary_csv"])

    for column in ("best_vina", "mean_vina", "best_cnn_score"):
        for ascending in (True, False):
            top = panda.execute_head(panda.execute_sort(table, {"by": column, "ascending": ascending}), {"n": 1})
            assert list(top["id"]) == ["good+lig_conf0"], (column, ascending)
    assert list(panda.execute_filter(table, {"expr": "best_vina < 0"})["id"]) == ["good+lig_conf0"]
    rows = [{**row, "pseudo_binding_energy": ""} for row in summary.fillna("").to_dict("records")]
    winners, unrankable = pipe_gnina._pair_winners(rows, {"good+lig_conf0": {}}, "best_vina")
    assert winners == {"good+lig": "good+lig_conf0"} and not unrankable
