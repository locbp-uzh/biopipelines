"""End-to-end runtime test for the Scripting runner (pipe_scripting.py).

Drives the runner against a real script + manifests and asserts the handle API materializes each output kind: file stream (``file``), standalone table (``row``), value stream (``add``), ``dataframe``, adopt-CSV (``outputs[name] = path``), and ``drop`` → scripting_dropped.csv with the framework-canonical ``kind='filter'``.
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
