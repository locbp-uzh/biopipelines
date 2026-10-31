"""A Scripting step keeps, and runs, a copy of the script it was given.

`my_scripts/` is gitignored, so without a copy in the run directory nothing records
which version of a helper produced the outputs. The step snapshots the resolved file
into its `_configuration/` folder and points the runner at that copy, so an edit
between submission and execution cannot change what runs.
"""

import glob
import json
import os
import subprocess
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RUNNER = os.path.join(REPO, "pipe_scripts", "pipe_scripting.py")

HELPER = '''
def configuration(inputs):
    from biopipelines.scripting_api import Table
    return {"out": Table(columns=["id", "n"])}

def execution(inputs, outputs):
    for sid in inputs["structures"].ids:
        outputs["out"].row({"id": sid, "n": 1})
'''


def _build(tmp_path, new_pipeline):
    from biopipelines.mock import Mock
    from biopipelines.scripting import Scripting

    helper = tmp_path / "helper_step.py"
    helper.write_text(HELPER, encoding="utf-8")

    pipeline = new_pipeline("scripting_snapshot")
    with pipeline:
        m = Mock(ids=["a", "b"],
                 streams={"structures": {"format": "pdb", "file": "<id>.pdb"}})
        Scripting(script=str(helper), inputs={"structures": m})
        script_path = pipeline.save()
    return helper, script_path


def test_snapshot_is_written_next_to_the_manifests(tmp_path, local_config,
                                                   isolated_cwd, new_pipeline):
    helper, script_path = _build(tmp_path, new_pipeline)

    run_root = os.path.dirname(os.path.dirname(script_path))  # job folder, above RunTime/
    copies = glob.glob(os.path.join(run_root, "**", "_configuration", "helper_step.py"),
                       recursive=True)
    assert copies, f"no snapshot under {run_root}"
    assert open(copies[0], encoding="utf-8").read() == helper.read_text(encoding="utf-8")


def test_generated_bash_runs_the_snapshot(tmp_path, local_config,
                                          isolated_cwd, new_pipeline):
    helper, script_path = _build(tmp_path, new_pipeline)

    step_sh = [p for p in glob.glob(os.path.join(os.path.dirname(script_path), "*.sh"))
               if "Scripting" in os.path.basename(p)]
    assert step_sh, f"no Scripting step script in {os.path.dirname(script_path)}"
    content = open(step_sh[0], encoding="utf-8", errors="replace").read()

    run_line = [L for L in content.splitlines() if "--script " in L]
    assert run_line, "no --script line in the generated bash"
    assert "_configuration" in run_line[0], run_line[0]
    assert f'--script-origin-dir "{os.path.dirname(str(helper))}"' in content


def test_origin_dir_keeps_a_sibling_import_working(tmp_path):
    """The snapshot runs from elsewhere, so the script's own folder goes on sys.path."""
    home = tmp_path / "home"
    home.mkdir()
    (home / "sidekick.py").write_text("VALUE = 7\n", encoding="utf-8")
    (home / "step.py").write_text(
        "def configuration(inputs):\n"
        "    from biopipelines.scripting_api import Table\n"
        "    return {'out': Table(columns=['id', 'n'])}\n\n"
        "def execution(inputs, outputs):\n"
        "    import sidekick\n"
        "    outputs['out'].row({'id': 'x', 'n': sidekick.VALUE})\n",
        encoding="utf-8")

    elsewhere = tmp_path / "run" / "_configuration"
    elsewhere.mkdir(parents=True)
    snapshot = elsewhere / "step.py"
    snapshot.write_text((home / "step.py").read_text(encoding="utf-8"), encoding="utf-8")

    out_csv = tmp_path / "run" / "out.csv"
    inputs_json = tmp_path / "run" / "in.json"
    outputs_json = tmp_path / "run" / "out.json"
    inputs_json.write_text(json.dumps({}), encoding="utf-8")
    outputs_json.write_text(json.dumps({
        "streams": {},
        "tables": {"out": {"columns": ["id", "n"], "path": str(out_csv)}},
        "tool_name": "Scripting",
        "missing_csv": str(tmp_path / "run" / "dropped.csv"),
    }), encoding="utf-8")

    def run(extra):
        return subprocess.run(
            [sys.executable, RUNNER, "--script", str(snapshot),
             "--inputs", str(inputs_json), "--outputs", str(outputs_json)] + extra,
            capture_output=True, text=True)

    without = run([])
    assert "sidekick" in (without.stderr + without.stdout), without.stderr[-400:]

    with_dir = run(["--script-origin-dir", str(home)])
    assert with_dir.returncode == 0, with_dir.stderr[-400:]
    assert "7" in open(out_csv, encoding="utf-8").read()


def test_output_streams_carry_their_producers_manifest(tmp_path, local_config,
                                                       isolated_cwd, new_pipeline):
    """A consumer must be able to find where its input's producer records drops."""
    from biopipelines.mock import Mock

    pipeline = new_pipeline("manifest_stamp")
    with pipeline:
        m = Mock(ids=["a"], streams={"structures": {"format": "pdb", "file": "<id>.pdb"}})

    stamped = m.streams.structures.metadata.get("missing_table")
    assert stamped, "no missing_table stamped on the output stream"
    assert stamped.endswith(os.path.join("tables", "missing.csv")), stamped


def test_the_step_manifest_names_the_step_for_removed_by(tmp_path, local_config,
                                                         isolated_cwd, new_pipeline):
    """`removed_by` is the step id `<NNN>_Scripting`, which tells a chain of two Scripting steps apart."""
    _helper, script_path = _build(tmp_path, new_pipeline)

    run_root = os.path.dirname(os.path.dirname(script_path))
    manifests = glob.glob(os.path.join(run_root, "**", "_configuration", "scripting_outputs.json"),
                          recursive=True)
    assert manifests, f"no output manifest under {run_root}"
    step_folder = os.path.dirname(os.path.dirname(manifests[0]))
    with open(manifests[0], encoding="utf-8") as fh:
        assert json.load(fh)["step_id"] == os.path.basename(step_folder)
