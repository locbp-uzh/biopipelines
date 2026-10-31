# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""MutationComposer: a declared sequence id it could not generate is a `failure` row of this step's `missing.csv`.

The wrapper declares `num_sequences` ids, and a short frequency table gives fewer sequences than that (a mode that draws each mutation once runs out of mutations, a stack finds every position taken). The step still exited 0 with the shorter file, so the ids it declared and never wrote were nowhere on record, and a Fuse downstream skipped the combinations built on them with a warning. The pipe now writes `execution/local_missing.csv` with one row for each declared id that has no sequence, the wrapper always declares `missing` and passes the local file to the propagation, and the step exits non-zero only when it wrote no sequence at all.

Nothing here needs the cluster: the pipe, the propagation and the completion check run on the files the wrapper planned.
"""

import json
import os
import pathlib
import subprocess
import sys

import pandas as pd
import pytest

from biopipelines import combinatorics
from biopipelines.mock import Mock
from tests._helpers import load_pipe_script, propagate_missing

REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent
SCRIPTS = REPO_ROOT / "pipe_scripts"
PIPE = SCRIPTS / "pipe_mutation_composer.py"
AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")
REFERENCE = "MKTA"
NUM_SEQUENCES = 10


def _frequency_table(path, bias=0):
    """One substitution at 0.4 per position, the rest at 0: four mutations in all, so ten sequences cannot all be drawn."""
    rows = []
    for position, original in enumerate(REFERENCE, start=1):
        row = {"position": position, "original": original}
        for aa in AMINO_ACIDS:
            row[aa] = 0.0
        row[original] = 0.6
        row[[a for a in AMINO_ACIDS if a != original][(position + bias) % 19]] = 0.4
        rows.append(row)
    pd.DataFrame(rows).to_csv(path, index=False)
    return str(path)


def _configured_step(new_pipeline, mode, strategy, tables, upstream, min_frequency=0.01):
    """Save a pipeline holding one MutationComposer and return (config path, config, script text, step id)."""
    from biopipelines.mutation_composer import MutationComposer

    pipeline = new_pipeline(f"composer_missing_{mode}_{strategy}")
    with pipeline:
        inputs = []
        for index in range(tables):
            spec = {"frequencies": {"columns": ["position", "original"] + AMINO_ACIDS}}
            if upstream and index == 0:
                spec["missing"] = {"columns": ["id", "removed_by", "kind", "cause"]}
            inputs.append(Mock(ids=["profile"], tables=spec, map_table_strategy="config"))
        MutationComposer(frequencies=inputs if tables > 1 else inputs[0], num_sequences=NUM_SEQUENCES,
                         mode=mode, combination_strategy=strategy, hotspot_count=3, random_seed=0,
                         min_frequency=min_frequency)
        pipeline.save()
    root = pathlib.Path.cwd()
    (config_path,) = root.rglob("mutation_composer_config.json")
    (script_path,) = root.rglob("*_MutationComposer.sh")
    with open(config_path) as handle:
        config = json.load(handle)
    step_id = script_path.stem
    for index, path in enumerate(config["frequency_tables"]):
        pathlib.Path(path).parent.mkdir(parents=True, exist_ok=True)
        _frequency_table(path, bias=index)
    return config_path, config, script_path.read_text(encoding="utf-8", errors="replace"), step_id


def _run(config_path):
    return subprocess.run([sys.executable, str(PIPE), "--config", str(config_path)], capture_output=True, text=True)


UPSTREAM_ROWS = "id,removed_by,kind,cause\nup,001_Mock,filter,dropped upstream\n"
# (mode, strategy, tables): every case here has fewer mutations than the ten sequences asked for.
SHORT_CASES = [
    ("single_point", "average", 1),
    ("single_point", "stack", 2),
    ("single_point", "round_robin", 2),
    ("hotspot_focused", "stack", 2),
]


@pytest.mark.parametrize("upstream", [False, True], ids=["no upstream table", "upstream table"])
@pytest.mark.parametrize("mode, strategy, tables", SHORT_CASES)
def test_a_declared_id_without_a_sequence_is_a_failure_row_of_the_final_missing_table(
        local_config, isolated_cwd, new_pipeline, mode, strategy, tables, upstream):
    config_path, config, script, step_id = _configured_step(new_pipeline, mode, strategy, tables, upstream)

    done = _run(config_path)
    assert done.returncode == 0, done.stdout + done.stderr
    written = list(pd.read_csv(config["sequences_output"])["id"])
    declared = [combinatorics.numbered_id(config["prefix"], i) for i in range(1, NUM_SEQUENCES + 1)]
    assert 0 < len(written) < NUM_SEQUENCES

    missing = propagate_missing(script, UPSTREAM_ROWS)
    own = missing[missing["removed_by"] == step_id]
    assert sorted(own["id"]) == sorted(set(declared) - set(written))
    assert set(own["kind"]) == {"failure"}
    assert all(cause for cause in own["cause"])
    if upstream:
        assert [(row["removed_by"], row["kind"]) for _, row in missing.iterrows() if row["id"] == "up"] == [
            ("001_Mock", "filter")]
    else:
        assert set(missing["id"]) == set(own["id"])


def test_the_completion_check_excuses_the_upstream_row_and_not_the_failure(
        local_config, isolated_cwd, new_pipeline, tmp_path):
    checker = load_pipe_script("pipe_check_completion_under_test", "pipe_check_completion.py")
    config_path, config, script, step_id = _configured_step(new_pipeline, "single_point", "average", 1, True)
    assert _run(config_path).returncode == 0
    missing = propagate_missing(script, UPSTREAM_ROWS)
    path = tmp_path / "missing.csv"
    missing.to_csv(path, index=False)
    excused = checker._load_expected_missing_ids({"tables": {"missing": {"path": str(path)}}}, step_id, "MutationComposer")
    assert excused == ["up"]


def test_the_step_fails_only_when_no_sequence_was_written(local_config, isolated_cwd, new_pipeline):
    config_path, config, script, step_id = _configured_step(
        new_pipeline, "single_point", "average", 1, False, min_frequency=0.9)

    done = _run(config_path)

    assert done.returncode != 0
    local = pd.read_csv(config["local_missing_csv"], dtype=str, keep_default_na=False)
    assert list(local["id"]) == [combinatorics.numbered_id(config["prefix"], i) for i in range(1, NUM_SEQUENCES + 1)]
    assert set(local["kind"]) == {"failure"} and set(local["removed_by"]) == {step_id}


def test_a_table_that_gives_every_sequence_leaves_no_failure_row(local_config, isolated_cwd, new_pipeline):
    config_path, config, script, step_id = _configured_step(new_pipeline, "weighted_random", "average", 1, False)

    assert _run(config_path).returncode == 0

    assert len(pd.read_csv(config["sequences_output"])) == NUM_SEQUENCES
    assert len(pd.read_csv(config["local_missing_csv"])) == 0


def test_the_missing_table_is_declared_whether_or_not_an_upstream_one_exists(local_config, isolated_cwd, new_pipeline):
    from biopipelines.mutation_composer import MutationComposer

    pipeline = new_pipeline("composer_declares_missing")
    with pipeline:
        out = MutationComposer(frequencies=_frequency_table(pathlib.Path("freq.csv")), num_sequences=3)
        script_path = pipeline.save()
    assert "missing" in out.tables._tables
    body = next(pathlib.Path(script_path).parent.glob("*_MutationComposer.sh")).read_text(encoding="utf-8", errors="replace")
    assert "--local-missing" in body
