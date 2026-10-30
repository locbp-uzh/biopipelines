# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Fuse: a declared fusion id it could not build is a `failure` row of this step's `missing.csv`, naming the slot sequence it lacked.

A slot sequence that the producer never wrote (the producer made fewer than it declared, or an upstream filter dropped it) leaves its row out of the slot's map table. Fuse declared every combination of the slot ids, built the ones it had sequences for and logged "skipping combo" for the rest, so over a short MutationComposer it declared 100 ids, wrote 36 and left the other 64 on no record. The pipe now writes `execution/local_missing.csv` with one row per declared id it did not build, the wrapper always declares `missing` and merges that file with the `missing` table of every slot's producer, and the step fails only when it built no fusion at all.
"""

import json
import pathlib
import subprocess
import sys

import pandas as pd
import pytest

from biopipelines.fuse import Fuse
from biopipelines.mock import Mock
from tests.test_mutation_composer_missing import _load, _propagate

REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent
PIPE = REPO_ROOT / "pipe_scripts" / "pipe_fuse_queries.py"
UPSTREAM_ROWS = "id,removed_by,kind,cause\na2,001_Mock,filter,dropped upstream\n"
DECLARED = ["a1+2+b1", "a1+2+b2", "a2+2+b1", "a2+2+b2", "a3+2+b1", "a3+2+b2"]


def _configured_step(new_pipeline, upstream):
    """Save a pipeline with one Fuse of two slots and return (config path, config, script text, step id)."""
    pipeline = new_pipeline("fuse_missing")
    with pipeline:
        tables = {"missing": {"columns": ["id", "removed_by", "kind", "cause"]}} if upstream else None
        first = Mock(ids=["a1", "a2", "a3"], tables=tables,
                     streams={"sequences": {"format": "csv", "values": ["AC", "DE", "FG"]}}, map_table_strategy="config")
        second = Mock(ids=["b1", "b2"], streams={"sequences": {"format": "csv", "values": ["HI", "KL"]}},
                      map_table_strategy="config")
        Fuse(sequences=[first, second], linker="GS", linker_lengths=["2"])
        pipeline.save()
    root = pathlib.Path.cwd()
    (config_path,) = root.rglob("fuse_config.json")
    (script_path,) = root.rglob("*_Fuse.sh")
    with open(config_path) as handle:
        config = json.load(handle)
    return config_path, config, script_path.read_text(encoding="utf-8", errors="replace"), script_path.stem


def _write_slots(config, present):
    """The map table of each slot, holding only the ids in `present` (a list of ids per slot)."""
    for slot, ids in zip(config["slots"], present):
        path = pathlib.Path(slot["map_table"])
        path.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame({"id": ids, "sequence": ["ACDE"[: 2 + i % 3] for i in range(len(ids))]}).to_csv(path, index=False)


def _run(config_path):
    return subprocess.run([sys.executable, str(PIPE), "--config", str(config_path)], capture_output=True, text=True)


@pytest.mark.parametrize("upstream", [False, True], ids=["no upstream table", "upstream table"])
def test_a_fusion_id_that_was_not_built_is_a_failure_row_naming_the_missing_slot_sequence(
        local_config, isolated_cwd, new_pipeline, upstream):
    config_path, config, script, step_id = _configured_step(new_pipeline, upstream)
    _write_slots(config, [["a1", "a3"], ["b1", "b2"]])

    done = _run(config_path)

    assert done.returncode == 0, done.stdout + done.stderr
    written = list(pd.read_csv(config["output_csv"])["id"])
    assert written == ["a1+2+b1", "a1+2+b2", "a3+2+b1", "a3+2+b2"]
    missing = _propagate(script, config, UPSTREAM_ROWS)
    own = missing[missing["removed_by"] == step_id]
    assert sorted(own["id"]) == sorted(set(DECLARED) - set(written)) == ["a2+2+b1", "a2+2+b2"]
    assert set(own["kind"]) == {"failure"}
    assert all("slot 1 sequence 'a2'" in cause for cause in own["cause"]), list(own["cause"])
    if upstream:
        assert [(row["removed_by"], row["kind"]) for _, row in missing.iterrows() if row["id"] == "a2"] == [
            ("001_Mock", "filter")]
    else:
        assert set(missing["id"]) == set(own["id"])


def test_a_combination_lacking_two_slots_names_both(local_config, isolated_cwd, new_pipeline):
    config_path, config, script, step_id = _configured_step(new_pipeline, False)
    _write_slots(config, [["a1", "a3"], ["b1"]])

    assert _run(config_path).returncode == 0

    local = pd.read_csv(config["local_missing_csv"], dtype=str, keep_default_na=False)
    both = local[local["id"] == "a2+2+b2"]
    assert len(both) == 1
    assert "slot 1 sequence 'a2'" in both.iloc[0]["cause"] and "slot 2 sequence 'b2'" in both.iloc[0]["cause"]


def test_the_step_fails_only_when_no_fusion_was_built(local_config, isolated_cwd, new_pipeline):
    config_path, config, script, step_id = _configured_step(new_pipeline, False)
    _write_slots(config, [["x"], ["b1", "b2"]])

    done = _run(config_path)

    assert done.returncode != 0
    local = pd.read_csv(config["local_missing_csv"], dtype=str, keep_default_na=False)
    assert sorted(local["id"]) == sorted(DECLARED)
    assert set(local["removed_by"]) == {step_id} and set(local["kind"]) == {"failure"}


def test_a_complete_slot_leaves_no_failure_row(local_config, isolated_cwd, new_pipeline):
    config_path, config, script, step_id = _configured_step(new_pipeline, False)
    _write_slots(config, [["a1", "a2", "a3"], ["b1", "b2"]])

    assert _run(config_path).returncode == 0

    assert list(pd.read_csv(config["output_csv"])["id"]) == DECLARED
    assert len(pd.read_csv(config["local_missing_csv"])) == 0


def test_the_missing_table_is_declared_whether_or_not_an_upstream_one_exists(local_config, isolated_cwd, new_pipeline):
    pipeline = new_pipeline("fuse_declares_missing")
    with pipeline:
        first = Mock(ids=["a1"], streams={"sequences": {"format": "csv", "values": ["AC"]}}, map_table_strategy="config")
        second = Mock(ids=["b1"], streams={"sequences": {"format": "csv", "values": ["HI"]}}, map_table_strategy="config")
        out = Fuse(sequences=[first, second])
        script_path = pipeline.save()
    assert "missing" in out.tables._tables
    body = next(pathlib.Path(script_path).parent.glob("*_Fuse.sh")).read_text(encoding="utf-8", errors="replace")
    assert "--local-missing" in body
