# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Gnina and Vina: a pair or complex that fails here is a `failure` row of this step's `missing.csv`.

The pipe script wrote its failures to `tables/missing.csv` itself, and only when an upstream `missing` table existed.
When one did, the propagation step rewrote that same file from the upstream manifests alone and erased them; when none
did, the table was not declared and nothing was written. The pipe now writes `execution/local_missing.csv`, the wrapper
always declares `missing` and passes the local file to the propagation, and every pair that has no entry in the
structures stream, whatever the reason, is one row under its own id.

Nothing here runs gnina or vina: the binary is a function, and the propagation and the completion check run on the
files the pipe wrote.
"""

import json
import pathlib
import shlex
import subprocess
import sys

import pandas as pd
import pytest

from biopipelines import gnina as gnina_module
from biopipelines.mock import Mock
from tests._helpers import (
    GninaBinary, configured_gnina_step, fake_the_gnina_tools, load_pipe_script, materialize_gnina_inputs,
)


def _run_step_and_propagate(pipe, monkeypatch, new_pipeline, cls_name, mode, upstream):
    config_path, config, script, step_id = configured_gnina_step(new_pipeline, cls_name, mode, upstream)
    materialize_gnina_inputs(config, upstream)
    fake_the_gnina_tools(pipe, monkeypatch, cls_name, mode)
    monkeypatch.setattr(sys, "argv", ["pipe_gnina.py", str(config_path)])
    try:
        pipe.main()
    except SystemExit:
        pass
    # What the generated script does next: merge the upstream manifests and the local file into tables/missing.csv.
    script = script.replace("\r\n", "\n")
    block = script[script.rindex("python", 0, script.index("pipe_propagate_missing.py")):].split("\n\n")[0]
    command = shlex.split(block.replace("\\\n", " "))
    arguments = {}
    for token in command[2:]:
        if token.startswith("--"):
            arguments[token] = []
        else:
            arguments[list(arguments)[-1]].append(token)
    propagate = load_pipe_script("pipe_propagate_missing_under_test", "pipe_propagate_missing.py")
    propagate.propagate_missing_table(
        upstream_folders=arguments["--upstream-folders"], output_folder=arguments["--output-folder"][0],
        missing_csv_path=arguments["--missing-csv"][0],
        local_missing=(arguments.get("--local-missing") or [None])[0])
    return pd.read_csv(arguments["--missing-csv"][0], dtype=str, keep_default_na=False), step_id


CASES = [("Gnina", "docking"), ("Vina", "docking"), ("Gnina", "score"), ("Gnina", "minimize")]
# What each case's failure row says: a binary that exits 1, a receptor that cannot be prepared, a complex without the ligand.
CAUSES = {("Gnina", "docking"): "GNINA failed", ("Vina", "docking"): "receptor PDBQT conversion failed",
          ("Gnina", "score"): "no HETATM records", ("Gnina", "minimize"): "no HETATM records"}


@pytest.mark.parametrize("upstream", [False, True], ids=["no upstream table", "upstream table"])
@pytest.mark.parametrize("cls_name, mode", CASES)
def test_a_failed_pair_or_complex_is_a_failure_row_of_the_final_missing_table(
        pipe_gnina, monkeypatch, local_config, isolated_cwd, new_pipeline, cls_name, mode, upstream):
    missing, step_id = _run_step_and_propagate(pipe_gnina, monkeypatch, new_pipeline, cls_name, mode, upstream)
    failed_id = "bad+lig" if mode == "docking" else "bad"
    own = missing[missing["removed_by"] == step_id]
    assert [(row["id"], row["kind"]) for _, row in own.iterrows()] == [(failed_id, "failure")]
    assert CAUSES[(cls_name, mode)] in own.iloc[0]["cause"], own.iloc[0]["cause"]
    if upstream:
        assert [(row["id"], row["removed_by"], row["kind"]) for _, row in missing.iterrows()
                if row["id"] == "up"] == [("up", "001_Mock", "filter")]
    else:
        assert list(missing["id"]) == [failed_id]


@pytest.mark.parametrize("cls_name, mode", [("Gnina", "docking"), ("Gnina", "score")])
def test_the_completion_check_excuses_the_upstream_row_and_not_the_failure(
        pipe_gnina, monkeypatch, local_config, isolated_cwd, new_pipeline, tmp_path, cls_name, mode):
    checker = load_pipe_script("pipe_check_completion_under_test", "pipe_check_completion.py")
    missing, step_id = _run_step_and_propagate(pipe_gnina, monkeypatch, new_pipeline, cls_name, mode, True)
    path = tmp_path / "missing.csv"
    missing.to_csv(path, index=False)
    excused = checker._load_expected_missing_ids({"tables": {"missing": {"path": str(path)}}}, step_id, cls_name)
    assert excused == ["up"]


@pytest.mark.parametrize("cls_name, mode", CASES)
def test_the_missing_table_is_declared_whether_or_not_an_upstream_one_exists(
        local_config, isolated_cwd, new_pipeline, cls_name, mode):
    pipeline = new_pipeline("gnina_declares")
    with pipeline:
        structures = Mock(ids=["good"], streams={"structures": {"format": "pdb", "file": "<id>.pdb"}},
                          map_table_strategy="config")
        compounds = Mock(ids=["lig"], streams={"compounds": {"format": "csv", "file": "<id>.csv"}},
                         map_table_strategy="config")
        out = getattr(gnina_module, cls_name)(structures=structures.streams.structures,
                                              compounds=compounds.streams.compounds, mode=mode)
        script_path = pipeline.save()
    assert "missing" in out.tables._tables
    body = next(pathlib.Path(script_path).parent.glob(f"*_{cls_name}.sh")).read_text(encoding="utf-8", errors="replace")
    assert "--local-missing" in body


# ── a step in which no pair gave a structure fails, and says why for each pair ──

def _docking_step(pipe, monkeypatch, new_pipeline, binary, smiles="CCO"):
    config_path, config, _script, step_id = configured_gnina_step(new_pipeline, "Gnina", "docking", False)
    materialize_gnina_inputs(config, False)
    with open(config["compounds_json"]) as handle:
        pathlib.Path(json.load(handle)["map_table"]).write_text(f"id,smiles,code\nlig,{smiles},LIG\n")
    monkeypatch.setattr(pipe.subprocess, "run", binary)
    monkeypatch.setattr(sys, "argv", ["pipe_gnina.py", str(config_path)])
    with pytest.raises(SystemExit) as stopped:
        pipe.main()
    assert stopped.value.code == 1
    local = pd.read_csv(config["local_missing_csv"], dtype=str, keep_default_na=False)
    return {row["id"]: (row["removed_by"], row["kind"], row["cause"]) for _, row in local.iterrows()}, step_id


def test_a_ligand_that_cannot_be_prepared_fails_every_pair_of_it(pipe_gnina, monkeypatch, local_config,
                                                                   isolated_cwd, new_pipeline):
    rows, step_id = _docking_step(pipe_gnina, monkeypatch, new_pipeline, GninaBinary(), smiles="not a smiles")
    assert set(rows) == {"good+lig", "bad+lig"}
    for removed_by, kind, cause in rows.values():
        assert (removed_by, kind) == (step_id, "failure")
        assert "no conformer could be prepared for the ligand" in cause


def test_a_pair_without_an_accepted_pose_is_a_failure_row(pipe_gnina, monkeypatch, local_config, isolated_cwd,
                                                          new_pipeline):
    def empty(cmd, **kwargs):
        pathlib.Path(cmd[cmd.index("-o") + 1]).write_text("")
        return subprocess.CompletedProcess(cmd, 0, "", "")

    rows, step_id = _docking_step(pipe_gnina, monkeypatch, new_pipeline, empty)
    assert rows == {"good+lig": (step_id, "failure", "no accepted pose"),
                    "bad+lig": (step_id, "failure", "no accepted pose")}
