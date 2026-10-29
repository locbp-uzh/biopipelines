"""Panda keeps numeric and boolean columns usable when a failed id left a blank cell.

A result table keeps one row per id that entered and writes NaN as an empty cell, so any table with a failed id has blanks. Panda reads its inputs with ``keep_default_na=False`` (an id named "nan" or "None" must stay a string), which made every column holding a blank a column of text: ``score > 0.5`` raised TypeError, ``passed == True`` kept no row, and ``sort`` ordered "10.5" before "9.5".
"""

import csv
import json
import os
import subprocess
import sys

import pytest

ROWS = [
    ["a", 9.5, True, 3],
    ["b", None, None, None],
    ["c", 10.5, False, 12],
    ["d", -2.0, True, 1],
]


def _repo_root():
    import biopipelines.panda as _p
    return os.path.dirname(os.path.dirname(os.path.abspath(_p.__file__)))


def _run_pipe(name, config_path, config_flag=True):
    script = os.path.join(_repo_root(), "pipe_scripts", f"pipe_{name}.py")
    args = [sys.executable, script] + (["--config", str(config_path)] if config_flag else [str(config_path)])
    subprocess.run(args, check=True)


def _panda(new_pipeline, name, operations):
    from biopipelines.mock import Mock
    from biopipelines.panda import Panda

    pipeline = new_pipeline(name)
    with pipeline:
        m = Mock(
            ids=[r[0] for r in ROWS],
            streams={"structures": {"format": "pdb", "file": "<id>.pdb"}},
            tables={"scores": {"columns": ["id", "score", "passed", "count"], "rows": ROWS}},
        )
        pan = Panda(tables=m.tables.scores, operations=operations, pool=m)
        pipeline.save()
    _run_pipe("mock", os.path.join(m.output_folder, "_configuration", "mock_config.json"), config_flag=False)
    config = os.path.join(pan.output_folder, "_configuration", "panda_config.json")
    _run_pipe("panda", config)
    with open(config) as f:
        output_csv = json.load(f)["output_csv"]
    with open(output_csv, newline="") as f:
        return list(csv.DictReader(f))


def test_filter_on_a_numeric_column_with_a_blank(local_config, isolated_cwd, new_pipeline):
    from biopipelines.panda import Panda
    rows = _panda(new_pipeline, "panda_blank_numeric", [Panda.filter("score > 0.5")])
    assert [r["id"] for r in rows] == ["a", "c"]


def test_filter_on_a_boolean_column_with_a_blank(local_config, isolated_cwd, new_pipeline):
    from biopipelines.panda import Panda
    rows = _panda(new_pipeline, "panda_blank_boolean", [Panda.filter("passed == True")])
    assert [r["id"] for r in rows] == ["a", "d"]


def test_filter_on_an_integer_column_with_a_blank_keeps_integers(local_config, isolated_cwd, new_pipeline):
    from biopipelines.panda import Panda
    rows = _panda(new_pipeline, "panda_blank_integer", [Panda.filter("count >= 3")])
    assert [(r["id"], r["count"]) for r in rows] == [("a", "3"), ("c", "12")]


def test_sort_on_a_numeric_column_with_a_blank_is_numeric(local_config, isolated_cwd, new_pipeline):
    from biopipelines.panda import Panda
    rows = _panda(new_pipeline, "panda_blank_sort", [Panda.sort("score", ascending=True)])
    assert [r["score"] for r in rows] == ["-2.0", "9.5", "10.5", ""]


@pytest.mark.parametrize("text, expected", [
    ("id,x\nnan,1\nnull,\nNone,3\n", ["nan", "null", "None"]),
])
def test_ids_that_look_like_missing_values_stay_text(tmp_path, text, expected):
    import importlib.util
    path = os.path.join(_repo_root(), "pipe_scripts", "pipe_panda.py")
    spec = importlib.util.spec_from_file_location("_pipe_panda_read", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    csv_path = tmp_path / "t.csv"
    csv_path.write_text(text)
    df = mod.read_table_csv(str(csv_path))
    assert df["id"].tolist() == expected
