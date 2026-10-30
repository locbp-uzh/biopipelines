"""Gnina and Vina: a receptor that cannot be read is a `failure` row, and a failed complex keeps a NaN row of `scores`.

Open Babel reads an ATOM line whose coordinates are not numbers as an atom at the origin. Vina boxed on the origin-centred
conformer, found the box inside that one-point "receptor", and docked into it, so the pair scored about -1 kcal/mol
instead of failing; GNINA stopped for lack of a box. The pipe script now refuses such a receptor before it is protonated,
converted or boxed, and says why in `missing`. Score and minimize wrote a `missing` row for a complex that failed and no
row of `scores`, although a result table keeps one row per id that entered, with NaN where nothing was measured.
"""

import json
import pathlib
import subprocess
import sys

import pandas as pd
import pytest

from tests.test_gnina_conformer_selection import _write_sdf
from tests.test_gnina_ids import PROTEIN_LINE, pipe_gnina  # noqa: F401 (a fixture)
from tests.test_gnina_missing_table import (
    _configured_step, _fake_the_tools, _materialize,
)

UNPARSEABLE_LINE = "ATOM      1  CA  GLY A   1         x.xxx   y.yyy   z.zzz  1.00  0.00           C\n"


def test_a_receptor_with_unparseable_coordinates_is_not_prepared(pipe_gnina, monkeypatch, tmp_path):
    files = {}
    for name, line in (("good", PROTEIN_LINE), ("bad", PROTEIN_LINE + UNPARSEABLE_LINE)):
        files[name] = tmp_path / f"{name}.pdb"
        files[name].write_text(line + "END\n")
    monkeypatch.setattr(pipe_gnina, "iterate_files", lambda pairs: [(i, str(f)) for i, f in pairs])
    failed = {}

    prepared = pipe_gnina.prepare_proteins(list(files.items()), str(tmp_path / "out"), False, 7.4, failed=failed)

    assert list(prepared) == ["good"]
    assert list(failed) == ["bad"] and "unparseable coordinates" in failed["bad"]


def test_a_pair_of_an_unprepared_protein_carries_the_cause(pipe_gnina):
    causes = pipe_gnina.unprepared_causes(
        ["good", "bad"], ["lig"], {"good": {}}, [{"ligand_id": "lig"}],
        failed={"bad": "the receptor has unparseable coordinates"})

    assert causes == {"bad+lig": ["the receptor has unparseable coordinates"]}


@pytest.mark.parametrize("cls_name", ["Gnina", "Vina"])
def test_an_unparseable_receptor_fails_its_pair_and_the_rest_of_the_batch_goes_on(
        pipe_gnina, monkeypatch, local_config, isolated_cwd, new_pipeline, cls_name):
    config_path, config, _script, step_id = _configured_step(new_pipeline, cls_name, "docking", False)
    _materialize(config, False)
    with open(config["structures_json"]) as handle:
        template = json.load(handle)["files"][0]
    pathlib.Path(template.replace("<id>", "bad")).write_text(PROTEIN_LINE + UNPARSEABLE_LINE + "END\n")
    _fake_the_tools(pipe_gnina, monkeypatch, cls_name, "docking")
    # Binaries that accept every receptor, so only the pipe script itself can refuse this one.
    monkeypatch.setattr(pipe_gnina.subprocess, "run", lambda cmd, **kw: (
        _write_sdf(cmd[cmd.index("-o") + 1], [(-7.0, 0.8, 5.0)]) or subprocess.CompletedProcess(cmd, 0, "", "")))
    monkeypatch.setattr(pipe_gnina.vina_backend, "receptor_pdbqt", lambda src, dst: dst)
    monkeypatch.setattr(sys, "argv", ["pipe_gnina.py", str(config_path)])

    pipe_gnina.main()

    local = pd.read_csv(config["local_missing_csv"], dtype=str, keep_default_na=False)
    assert list(local["id"]) == ["bad+lig"]
    assert (local.iloc[0]["removed_by"], local.iloc[0]["kind"]) == (step_id, "failure")
    assert "unparseable coordinates" in local.iloc[0]["cause"]
    assert list(pd.read_csv(config["structures_map"])["id"]) == ["good+lig"]


@pytest.mark.parametrize("cls_name, mode", [("Gnina", "score"), ("Gnina", "minimize"), ("Vina", "score")])
def test_a_failed_complex_keeps_a_nan_row_of_scores(
        pipe_gnina, monkeypatch, local_config, isolated_cwd, new_pipeline, cls_name, mode):
    config_path, config, _script, _step = _configured_step(new_pipeline, cls_name, mode, False)
    _materialize(config, False)
    _fake_the_tools(pipe_gnina, monkeypatch, cls_name, mode)
    monkeypatch.setattr(pipe_gnina.vina_backend, "run", lambda *a, **k: subprocess.CompletedProcess(
        ["vina"], 0, "Estimated Free Energy of Binding   :  -7.0 (kcal/mol)\n", ""))
    monkeypatch.setattr(sys, "argv", ["pipe_gnina.py", str(config_path)])

    pipe_gnina.main()

    scores = pd.read_csv(config["scores_csv"], dtype={"id": str})
    assert list(scores["id"]) == ["good", "bad"], "one row per complex that entered"
    bad = scores[scores["id"] == "bad"].iloc[0]
    assert (bad["structures.id"], bad["compounds.id"]) == ("bad", "lig")
    numeric = [column for column in scores.columns if column not in ("id", "structures.id", "compounds.id")]
    assert scores.loc[scores["id"] == "bad", numeric].isna().all(axis=None)
    assert scores.loc[scores["id"] == "good", "vina_affinity"].notna().all()
    local = pd.read_csv(config["local_missing_csv"], dtype=str)
    assert list(local["id"]) == ["bad"], "the failure is still a row of missing"
