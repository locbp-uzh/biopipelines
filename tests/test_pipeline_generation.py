"""Integration smoke test: Pipeline(...).save() produces a bash script.

This verifies the full config-time path (tool registration, stream wiring,
bash emission) end-to-end without executing any ML model. The pipeline uses
only the Sequence tool, which builds CSV/FASTA files from an inline sequence
and has no GPU / external-weight dependencies.
"""

import os


def test_pipeline_save_emits_bash_script(local_config, isolated_cwd):
    """Build a single-step Sequence pipeline and save it to disk."""
    from biopipelines.pipeline import Pipeline
    from biopipelines.sequence import Sequence

    # No Resources() call: the fixture config uses scheduler=none, so resource
    # specs (which are only meaningful for SLURM) are auto-initialized.
    pipeline = Pipeline(
        project="TestSuite",
        job="smoke",
        description="Smoke test from the automated suite",
        on_the_fly=False,
        local_output=True,
        config="local",
    )
    with pipeline:
        Sequence(seq="MKTAYIAKQRQISFVKSHFSRQLEERLGL", type="protein", ids="demo")
        script_path = pipeline.save()

    # Assertions: a bash script was written, looks like bash, mentions the
    # pipeline identity.
    assert os.path.isfile(script_path), f"pipeline.sh not found at {script_path}"
    assert script_path.endswith("pipeline.sh")

    content = open(script_path, encoding="utf-8").read()
    assert content.startswith("#!/bin/bash"), "pipeline.sh missing shebang"
    assert "TestSuite" in content, "project name missing from generated script"
    assert "smoke" in content, "job name missing from generated script"
    assert len(content) > 200, "pipeline.sh suspiciously short"


# ── Multi-tool chain ──────────────────────────────────────────────────────────

def test_pipeline_three_tool_chain_via_mock(
    local_config, isolated_cwd, new_pipeline, assert_valid_script, record_case,
):
    """Mock → ReMap → Panda — verify a 3-tool chain emits a valid script with
    tools appearing in the right order."""
    from biopipelines.mock import Mock
    from biopipelines.remap import ReMap
    from biopipelines.panda import Panda

    pipeline = new_pipeline("chain_three")
    with pipeline:
        m = Mock(
            ids=["a", "b", "c"],
            streams={"structures": {"format": "pdb", "file": "<id>.pdb"}},
            tables={"scores": {"columns": ["score"], "fill": 1.0}},
            map_table_strategy="config",
        )
        ReMap(source=m, onto="design")
        Panda(
            tables=m.tables.scores,
            operations=[Panda.sort("id", ascending=True)],
        )
        script_path = pipeline.save()

    content = open(script_path, encoding="utf-8").read()
    mock_pos = content.find("Mock")
    remap_pos = content.find("ReMap")
    panda_pos = content.find("Panda")
    record_case(
        input="Mock → ReMap → Panda chain",
        expected="all 3 markers present, Mock before downstream",
        actual=f"Mock@{mock_pos}, ReMap@{remap_pos}, Panda@{panda_pos}",
    )
    assert_valid_script(script_path, "Mock", "ReMap", "Panda")
    assert mock_pos >= 0 and mock_pos < remap_pos
    assert mock_pos < panda_pos


def test_pipeline_missing_propagates_through_downstream(
    local_config, isolated_cwd, new_pipeline, assert_valid_script, record_case,
):
    """Mock with `missing` feeding Panda — downstream wiring must not crash
    when upstream drops IDs."""
    from biopipelines.mock import Mock
    from biopipelines.panda import Panda

    pipeline = new_pipeline("chain_missing")
    with pipeline:
        m = Mock(
            ids=["a", "b", "c"],
            streams={"s": {"format": "pdb", "file": "<id>.pdb"}},
            tables={"scores": {"columns": ["score"], "fill": 0.5}},
            missing=["b"],
            map_table_strategy="config",
        )
        Panda(
            tables=m.tables.scores,
            operations=[Panda.sort("id")],
        )
        script_path = pipeline.save()

    record_case(
        input="Mock(missing=['b']) → Panda",
        expected="script emitted, both tools wired",
        actual=os.path.basename(script_path),
    )
    assert_valid_script(script_path, "Mock", "Panda", "chain_missing")


# ── Pipeline output folder structure ──────────────────────────────────────────

def test_pipeline_save_creates_runtime_folder_layout(
    local_config, isolated_cwd, new_pipeline, record_case,
):
    """pipeline.save() must create RunTime/ alongside pipeline.sh and a per-tool
    output_folder for each registered tool."""
    from biopipelines.mock import Mock

    pipeline = new_pipeline("layout_check")
    with pipeline:
        m = Mock(
            ids=["a"],
            streams={"s": {"format": "pdb", "file": "<id>.pdb"}},
        )
        script_path = pipeline.save()

    runtime_dir = os.path.dirname(script_path)
    tool_out = m.output_folder
    record_case(
        input="pipeline.save() layout",
        expected=("RunTime dir exists, tool output_folder exists", True, True),
        actual=("RunTime dir, tool output_folder",
                os.path.isdir(runtime_dir), os.path.isdir(tool_out)),
    )
    assert os.path.isdir(runtime_dir)
    assert os.path.isdir(tool_out)
    assert os.path.basename(runtime_dir) == "RunTime"


# ── End-to-end: execute the generated pipeline.sh ────────────────────────────

def test_generated_pipeline_sh_executes_end_to_end(
    local_config, isolated_cwd, new_pipeline, record_case,
):
    """Execute the generated ``pipeline.sh`` itself (not the individual
    pipe_*.py scripts) for a 3-stage Mock-only chain and verify each stage
    materialized its expected stub output files.

    This is the integration test the A2 response promises: the master bash
    script emitted by ``Pipeline.save()`` is chmod'd and invoked, exercising
    the config-time wiring, the per-tool ``.sh`` files, the master script's
    orchestration, and the runtime ``pipe_mock.py`` execution in one shot.
    Requires ``bash`` on PATH; pip-mode fixture avoids any conda dependency.
    """
    import shutil
    import subprocess
    import sys

    import pytest

    if sys.platform.startswith("win"):
        pytest.skip(
            "generated pipeline.sh embeds Windows paths with backslashes "
            "that MSYS/Git-bash reinterprets as escapes; the script is "
            "designed for POSIX hosts (Linux CI, HPC). See tests.yml."
        )
    bash = shutil.which("bash")
    if bash is None:
        pytest.skip("bash not available on PATH")

    from biopipelines.mock import Mock

    pipeline = new_pipeline("e2e_pipeline_sh")
    with pipeline:
        m1 = Mock(
            ids=["s1", "s2"],
            streams={"structures": {"format": "pdb", "file": "<id>.pdb"}},
        )
        m2 = Mock(
            source=m1.streams.structures,
            children="<1..2>",
            streams={"designs": {"format": "pdb", "file": "<id>.pdb"}},
        )
        m3 = Mock(
            source=m2.streams.designs,
            streams={"refined": {"format": "pdb", "file": "<id>.pdb"}},
        )
        script_path = pipeline.save()

    # Execute pipeline.sh end-to-end.
    result = subprocess.run(
        [bash, script_path], capture_output=True, text=True, timeout=120,
    )
    stdout_tail = result.stdout[-600:]
    stderr_tail = result.stderr[-600:]

    # Files each stage should have produced (one .pdb per output ID).
    # Mock writes per-stream files into <output_folder>/<stream_name>/
    # under the canonical layout (stream_path helper on BaseConfig).
    s1 = os.path.join(m1.output_folder, "structures", "s1.pdb")
    s2 = os.path.join(m1.output_folder, "structures", "s2.pdb")
    d11 = os.path.join(m2.output_folder, "designs", "s1_1.pdb")
    d22 = os.path.join(m2.output_folder, "designs", "s2_2.pdb")
    r11 = os.path.join(m3.output_folder, "refined", "s1_1.pdb")

    produced = {
        "m1.s1": os.path.exists(s1),
        "m1.s2": os.path.exists(s2),
        "m2.s1_1": os.path.exists(d11),
        "m2.s2_2": os.path.exists(d22),
        "m3.s1_1": os.path.exists(r11),
    }

    record_case(
        input="bash pipeline.sh for Mock → Mock(<1..2>) → Mock",
        expected=("rc=0", {k: True for k in produced}),
        actual=(f"rc={result.returncode}", produced),
    )
    assert result.returncode == 0, (
        f"pipeline.sh failed (rc={result.returncode})\n"
        f"stdout: {stdout_tail}\nstderr: {stderr_tail}"
    )
    assert all(produced.values()), (
        f"missing stub outputs: {produced}\nstdout: {stdout_tail}"
    )
    # Master-script orchestration markers: every stage's echo header fired.
    assert result.stdout.count("Mock") >= 3, (
        f"expected 3 Mock stage headers in stdout, got:\n{stdout_tail}"
    )
    assert "Job done" in result.stdout
