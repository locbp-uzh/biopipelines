"""End-to-end runtime test for the Scripting runner (pipe_scripting.py).

Drives the runner against a real script + manifests and asserts the handle API materializes each output kind: file stream (``file``), standalone table (``row``), value stream (``add``), ``dataframe``, adopt-CSV (``outputs[name] = path``), and ``drop`` → scripting_dropped.csv with the framework-canonical ``kind='filter'`` (or ``'unmatched'``, which the completion check does not excuse). ``max_dropped`` fails the step when too many ids were dropped.
"""

import json
import os
import subprocess
import sys

import pandas as pd
import pytest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RUNNER = os.path.join(REPO, "pipe_scripts", "pipe_scripting.py")

SCRIPT = '''
def configuration(inputs):
    from biopipelines.scripting_api import Stream, Table
    ids = inputs["structures"].ids
    return {"kept": Stream("pdb", ids),
            "vals": Stream("csv", ids),
            "count": Table(columns=["id", "n_atoms"]),
            "adopted": Table(columns=["id", "x"])}

def execution(inputs, outputs):
    import pandas as pd
    for sid, path in inputs["structures"].iterate():
        n = sum(1 for L in open(path) if L.startswith("ATOM"))
        if n < 2:
            outputs.drop(sid, cause=f"n {n} < 2")
            continue
        out = outputs["kept"].file(sid, f"{sid}.pdb")
        open(out, "w").write(open(path).read())
        outputs["vals"].add(id=sid, value=n)
        outputs["count"].row({"id": sid, "n_atoms": n})
    outputs["adopted"].dataframe(pd.DataFrame([{"id": "z", "x": 1}]))
'''


@pytest.fixture
def scenario(tmp_path):
    """A 2-id structures stream (p1 keeps, p2 drops) + script + manifests."""
    sf = tmp_path / "sf"
    sf.mkdir()
    (sf / "p1.pdb").write_text("ATOM\nATOM\nATOM\n")   # 3 atoms -> kept
    (sf / "p2.pdb").write_text("ATOM\n")                # 1 atom  -> dropped

    map_csv = sf / "structures_map.csv"
    pd.DataFrame({"id": ["p1", "p2"],
                  "file": [str(sf / "p1.pdb"), str(sf / "p2.pdb")]}).to_csv(map_csv, index=False)
    ds_json = sf / "structures.json"
    ds_json.write_text(json.dumps({
        "name": "structures", "ids": ["p1", "p2"],
        "files": [str(sf / "p1.pdb"), str(sf / "p2.pdb")],
        "map_table": str(map_csv), "format": "pdb",
    }))

    script = sf / "user.py"
    script.write_text(SCRIPT)

    inputs_json = sf / "inputs.json"
    inputs_json.write_text(json.dumps({
        "structures": {"kind": "stream", "payload": str(ds_json)},
    }))

    dropped = sf / "scripting_dropped.csv"
    outputs_json = sf / "outputs.json"
    outputs_json.write_text(json.dumps({
        "streams": {
            "kept": {"format": "pdb", "map_table": str(sf / "kept_map.csv"), "folder": str(sf / "kept")},
            "vals": {"format": "csv", "map_table": str(sf / "vals_map.csv"), "folder": str(sf / "vals")},
        },
        "tables": {
            "count": {"columns": ["id", "n_atoms"], "path": str(sf / "count.csv")},
            "adopted": {"columns": ["id", "x"], "path": str(sf / "adopted.csv")},
        },
        "tool_name": "Scripting",
        "missing_csv": str(dropped),
    }))

    return {"sf": sf, "script": script, "inputs": inputs_json,
            "outputs": outputs_json, "dropped": dropped}


def _run(scenario):
    subprocess.run(
        [sys.executable, RUNNER,
         "--script", str(scenario["script"]),
         "--inputs", str(scenario["inputs"]),
         "--outputs", str(scenario["outputs"])],
        check=True, cwd=REPO,
    )


def test_file_stream_writes_map_and_files(scenario):
    _run(scenario)
    sf = scenario["sf"]
    kept = pd.read_csv(sf / "kept_map.csv")
    assert list(kept["id"]) == ["p1"]                      # p2 was dropped
    assert (sf / "kept" / "p1.pdb").exists()


def test_table_row(scenario):
    _run(scenario)
    count = pd.read_csv(scenario["sf"] / "count.csv")
    assert list(count["id"]) == ["p1"]
    assert int(count.loc[0, "n_atoms"]) == 3


def test_value_stream_add(scenario):
    _run(scenario)
    vals = pd.read_csv(scenario["sf"] / "vals_map.csv")
    assert list(vals["id"]) == ["p1"]
    assert list(vals["value"]) == [3]


def test_dataframe_handover(scenario):
    _run(scenario)
    adopted = pd.read_csv(scenario["sf"] / "adopted.csv")
    assert list(adopted["id"]) == ["z"]
    assert list(adopted["x"]) == [1]


def test_drop_writes_filter_kind(scenario):
    _run(scenario)
    dropped = pd.read_csv(scenario["dropped"])
    assert list(dropped["id"]) == ["p2"]
    # Must be the framework-canonical value the completion check excuses.
    assert list(dropped["kind"]) == ["filter"]
    assert list(dropped["removed_by"]) == ["Scripting"]


def test_whole_table_input(tmp_path):
    """A whole-table input (kind='table_full') exposes rows/columns/value(id, col)."""
    sf = tmp_path / "sf"
    sf.mkdir()
    (sf / "p1.pdb").write_text("ATOM\n")
    map_csv = sf / "m.csv"
    pd.DataFrame({"id": ["p1"], "file": [str(sf / "p1.pdb")]}).to_csv(map_csv, index=False)
    ds_json = sf / "s.json"
    ds_json.write_text(json.dumps({"name": "structures", "ids": ["p1"],
                                   "files": [str(sf / "p1.pdb")], "map_table": str(map_csv), "format": "pdb"}))
    scores = sf / "scores.csv"
    pd.DataFrame([{"id": "p1", "plddt": 80, "rmsd": 1.2}]).to_csv(scores, index=False)

    script = sf / "u.py"
    script.write_text(
        "def configuration(inputs):\n"
        "    from biopipelines.scripting_api import Table\n"
        "    return {'out': Table(columns=['id', 'plddt', 'rmsd', 'ncols'])}\n"
        "def execution(inputs, outputs):\n"
        "    t = inputs['scores']\n"
        "    for sid, path in inputs['structures'].iterate():\n"
        "        r = t.row(sid)\n"
        "        outputs['out'].row({'id': sid, 'plddt': r['plddt'],\n"
        "                            'rmsd': t.value(sid, 'rmsd'), 'ncols': len(t.columns)})\n"
    )
    inputs_json = sf / "i.json"
    inputs_json.write_text(json.dumps({
        "structures": {"kind": "stream", "payload": str(ds_json)},
        "scores": {"kind": "table_full", "payload": str(scores)},
    }))
    target = sf / "out.csv"
    outputs_json = sf / "o.json"
    outputs_json.write_text(json.dumps({
        "streams": {}, "tables": {"out": {"columns": ["id", "plddt", "rmsd", "ncols"], "path": str(target)}},
        "tool_name": "Scripting", "missing_csv": str(sf / "d.csv"),
    }))

    subprocess.run([sys.executable, RUNNER, "--script", str(script),
                    "--inputs", str(inputs_json), "--outputs", str(outputs_json)],
                   check=True, cwd=REPO)
    out = pd.read_csv(target)
    assert list(out["id"]) == ["p1"]
    assert int(out.loc[0, "plddt"]) == 80
    assert float(out.loc[0, "rmsd"]) == 1.2
    assert int(out.loc[0, "ncols"]) == 3   # id, plddt, rmsd


def test_adopt_existing_csv(tmp_path):
    """outputs[name] = path adopts a CSV wholesale into the managed target."""
    sf = tmp_path / "sf"
    sf.mkdir()
    (sf / "p1.pdb").write_text("ATOM\n")
    map_csv = sf / "m.csv"
    pd.DataFrame({"id": ["p1"], "file": [str(sf / "p1.pdb")]}).to_csv(map_csv, index=False)
    ds_json = sf / "s.json"
    ds_json.write_text(json.dumps({"name": "structures", "ids": ["p1"],
                                   "files": [str(sf / "p1.pdb")], "map_table": str(map_csv), "format": "pdb"}))
    src = sf / "external.csv"
    pd.DataFrame([{"id": "a", "n_atoms": 9}]).to_csv(src, index=False)

    script = sf / "u.py"
    script.write_text(
        "def configuration(inputs):\n"
        "    from biopipelines.scripting_api import Table\n"
        "    return {'count': Table(columns=['id', 'n_atoms'])}\n"
        "def execution(inputs, outputs):\n"
        f"    outputs['count'] = {str(src)!r}\n"
    )
    inputs_json = sf / "i.json"
    inputs_json.write_text(json.dumps({"structures": {"kind": "stream", "payload": str(ds_json)}}))
    target = sf / "count.csv"
    outputs_json = sf / "o.json"
    outputs_json.write_text(json.dumps({
        "streams": {}, "tables": {"count": {"columns": ["id", "n_atoms"], "path": str(target)}},
        "tool_name": "Scripting", "missing_csv": str(sf / "d.csv"),
    }))

    subprocess.run([sys.executable, RUNNER, "--script", str(script),
                    "--inputs", str(inputs_json), "--outputs", str(outputs_json)],
                   check=True, cwd=REPO)
    out = pd.read_csv(target)
    assert list(out["id"]) == ["a"]
    assert int(out.loc[0, "n_atoms"]) == 9


# ── drop(kind=...) ───────────────────────────────────────────────────────────
#
# Every drop used to be recorded as kind="filter", which the completion check
# excuses, so a lookup that found no partner for most ids passed as a filter.

DROP_CALL = 'outputs.drop(sid, cause=f"n {n} < 2")'


def _use_script(scenario, source):
    scenario["script"].write_text(source)


def _set_manifest(scenario, **extra):
    manifest = json.loads(scenario["outputs"].read_text())
    manifest.update(extra)
    scenario["outputs"].write_text(json.dumps(manifest))


def _run_raw(scenario):
    return subprocess.run(
        [sys.executable, RUNNER,
         "--script", str(scenario["script"]),
         "--inputs", str(scenario["inputs"]),
         "--outputs", str(scenario["outputs"])],
        cwd=REPO, capture_output=True, text=True,
    )


def test_drop_unmatched_is_recorded_with_its_kind(scenario):
    assert DROP_CALL in SCRIPT
    _use_script(scenario, SCRIPT.replace(DROP_CALL, DROP_CALL[:-1] + ', kind="unmatched")'))
    _run(scenario)
    dropped = pd.read_csv(scenario["dropped"])
    assert list(dropped["id"]) == ["p2"]
    assert list(dropped["kind"]) == ["unmatched"]
    assert list(dropped["removed_by"]) == ["Scripting"]


@pytest.mark.parametrize("kind", ["failure", "filtered", "", None])
def test_drop_rejects_a_kind_outside_the_controlled_set(scenario, kind):
    _use_script(scenario, SCRIPT.replace(DROP_CALL, DROP_CALL[:-1] + f", kind={kind!r})"))
    proc = _run_raw(scenario)
    assert proc.returncode != 0
    assert "kind must be one of" in proc.stderr


def _completion(scenario, tmp_path, kind):
    """Run the step with drops of ``kind``, merge its missing rows, run the completion check."""
    import pipe_scripts.pipe_check_completion as completion
    import pipe_scripts.pipe_propagate_missing as propagate

    _use_script(scenario, SCRIPT.replace(DROP_CALL, DROP_CALL[:-1] + f', kind="{kind}")'))
    _run(scenario)
    sf = scenario["sf"]
    step = tmp_path / "007_Scripting"
    missing = tmp_path / "tables" / "missing.csv"
    propagate.propagate_missing_table([], str(step), str(missing), str(scenario["dropped"]))
    expected = {
        "output_folder": str(step),
        "kept": {"name": "kept", "format": "pdb", "ids": ["p1", "p2"],
                 "files": [str(sf / "kept" / "<id>.pdb")],
                 "map_table": str(sf / "kept_map.csv")},
        "tables": {"missing": {"path": str(missing)}},
    }
    return completion.check_expected_outputs(expected, "Scripting", str(step))


def test_a_filter_drop_is_excused_by_the_completion_check(scenario, tmp_path):
    ok, missing = _completion(scenario, tmp_path, "filter")
    assert ok, missing


def test_an_unmatched_drop_is_flagged_by_the_completion_check(scenario, tmp_path):
    ok, missing = _completion(scenario, tmp_path, "unmatched")
    assert not ok
    assert [os.path.basename(p) for p in missing["kept"]] == ["p2.pdb"]


# ── max_dropped ──────────────────────────────────────────────────────────────

def test_max_dropped_exceeded_fails_the_step_with_both_counts(scenario):
    _set_manifest(scenario, max_dropped=0.4)   # 1 of 2 ids is dropped: 50% > 40%
    proc = _run_raw(scenario)
    assert proc.returncode != 0
    assert "max_dropped=0.4 exceeded" in proc.stderr
    assert "1 of 2 ids" in proc.stderr
    # outputs and the dropped ids are written before the check, so they stay inspectable
    assert list(pd.read_csv(scenario["dropped"])["id"]) == ["p2"]
    assert list(pd.read_csv(scenario["sf"] / "kept_map.csv")["id"]) == ["p1"]


def test_max_dropped_at_the_share_passes(scenario):
    _set_manifest(scenario, max_dropped=0.5)   # exactly 1 of 2: not "exceeds"
    assert _run_raw(scenario).returncode == 0


def test_max_dropped_of_one_never_fails(scenario):
    _set_manifest(scenario, max_dropped=1)
    assert _run_raw(scenario).returncode == 0


def test_no_max_dropped_keeps_todays_behavior(scenario):
    _set_manifest(scenario, max_dropped=None)
    assert _run_raw(scenario).returncode == 0


def test_max_dropped_counts_distinct_ids_not_calls(scenario):
    twice = DROP_CALL + "\n            " + DROP_CALL
    _use_script(scenario, SCRIPT.replace(DROP_CALL, twice))
    _set_manifest(scenario, max_dropped=0.5)   # p2 dropped twice is still 1 of 2
    assert _run_raw(scenario).returncode == 0
    assert len(pd.read_csv(scenario["dropped"])) == 2   # both calls stay on record


def test_max_dropped_entered_ids_are_the_union_across_inputs(scenario, tmp_path):
    """A table input adds ids to those that entered: 1 dropped of {p1, p2, t1, t2} is 25%."""
    table = tmp_path / "extra.csv"
    pd.DataFrame({"id": ["t1", "t2"], "x": [1, 2]}).to_csv(table, index=False)
    manifest = json.loads(scenario["inputs"].read_text())
    manifest["extra"] = {"kind": "table_full", "payload": str(table)}
    scenario["inputs"].write_text(json.dumps(manifest))

    _set_manifest(scenario, max_dropped=0.2)
    proc = _run_raw(scenario)
    assert proc.returncode != 0 and "1 of 4 ids" in proc.stderr

    _set_manifest(scenario, max_dropped=0.25)
    assert _run_raw(scenario).returncode == 0


# ── wrapper: max_dropped validation and manifest ─────────────────────────────

HELPER = '''
def configuration(inputs):
    from biopipelines.scripting_api import Table
    return {"out": Table(columns=["id", "n"])}

def execution(inputs, outputs):
    pass
'''


def _scripting(tmp_path, new_pipeline, **options):
    from biopipelines.mock import Mock
    from biopipelines.scripting import Scripting

    helper = tmp_path / "helper.py"
    helper.write_text(HELPER, encoding="utf-8")
    pipeline = new_pipeline("scripting_max_dropped")
    with pipeline:
        m = Mock(ids=["a", "b"], streams={"structures": {"format": "pdb", "file": "<id>.pdb"}})
        Scripting(script=str(helper), inputs={"structures": m}, **options)
        return pipeline.save()


@pytest.mark.parametrize("bad", [0, 0.0, -0.1, 1.5, True, "0.5", float("nan")])
def test_max_dropped_must_be_a_fraction_in_zero_one(tmp_path, local_config, isolated_cwd,
                                                    new_pipeline, bad):
    with pytest.raises(ValueError, match="max_dropped"):
        _scripting(tmp_path, new_pipeline, max_dropped=bad)


@pytest.mark.parametrize("value", [None, 0.25, 1, 1.0])
def test_max_dropped_reaches_the_runner_through_the_output_manifest(tmp_path, local_config,
                                                                    isolated_cwd, new_pipeline,
                                                                    value):
    import glob

    script_path = _scripting(tmp_path, new_pipeline, max_dropped=value)
    run_root = os.path.dirname(os.path.dirname(script_path))
    manifests = glob.glob(os.path.join(run_root, "**", "scripting_outputs.json"), recursive=True)
    assert manifests, f"no scripting_outputs.json under {run_root}"
    assert json.load(open(manifests[0], encoding="utf-8"))["max_dropped"] == value
