"""A step's own failure must survive the blocks that run after it.

A generated tool script is a sequence: the tool's own command, then optionally
the missing-manifest propagation, then the completion check. There is no
`set -e`, so the script's exit status is the LAST command's — and both trailing
blocks succeed on their own. A tool that failed therefore exited 0, and
`_step_failure_guard` read `${PIPESTATUS[0]}` as clean and wrote no marker.

Found on Daint: a `Scripting` step whose `execution()` raised was stamped
COMPLETED. The traceback was in its log, the propagation block ran after it,
and the step exited 0.

The second half is the completion check. A value-based stream (`files == []`)
keeps all its rows in its map_table, so it contributes no per-id path and every
existence check is vacuous for it. That step declared one such stream, wrote
nothing, and still reported "1 of 1 declared outputs found" -- the one output
verified was the `missing` table the trailing block had just created.

Either fix alone would have caught that run; both are needed because they are
different holes.
"""
from __future__ import annotations

import json
import os
import shutil
import subprocess
from pathlib import Path

import pytest

BASH = shutil.which("bash")

pytestmark = pytest.mark.skipif(
    BASH is None, reason="drives generated bash and needs a real bash",
)


# ── the exit-status half ──────────────────────────────────────────────────────

def _run_script(tmp_path: Path, body: str) -> subprocess.CompletedProcess:
    script = tmp_path / "step.sh"
    script.write_text("#!/bin/bash\n" + body, encoding="utf-8")
    script.chmod(0o755)
    return subprocess.run([BASH, str(script)], capture_output=True, text=True)


def test_a_failing_main_command_sets_the_scripts_exit_status(tmp_path):
    """The shape the emitters produce: main, then two succeeding blocks."""
    body = (
        "false\n"                      # the tool's own command, failing
        "BP_MAIN_RC=$?\n"              # propagation block's first line
        "true\n"                       # the rest of the propagation block
        "true\n"                       # the completion check
        'if [ "${BP_MAIN_RC:-0}" -ne 0 ]; then\n'
        '    echo "ERROR: exited ${BP_MAIN_RC}"\n'
        '    exit "${BP_MAIN_RC}"\n'
        "fi\n"
    )
    proc = _run_script(tmp_path, body)

    assert proc.returncode != 0, "a failed step still exited 0"
    assert "ERROR: exited 1" in proc.stdout


def test_a_succeeding_main_command_still_exits_zero(tmp_path):
    body = (
        "true\n"
        "BP_MAIN_RC=$?\n"
        "true\n"
        'if [ "${BP_MAIN_RC:-0}" -ne 0 ]; then exit "${BP_MAIN_RC}"; fi\n'
    )
    assert _run_script(tmp_path, body).returncode == 0


def test_the_footer_captures_the_status_when_no_propagation_block_ran(tmp_path):
    """Most tools emit no propagation block, so the footer must capture it."""
    body = (
        "false\n"                              # the tool's own command
        "BP_MAIN_RC=${BP_MAIN_RC:-$?}\n"       # footer's first line
        "true\n"
        'if [ "${BP_MAIN_RC:-0}" -ne 0 ]; then exit "${BP_MAIN_RC}"; fi\n'
    )
    assert _run_script(tmp_path, body).returncode != 0


def test_the_footer_does_not_overwrite_an_earlier_capture(tmp_path):
    """With a propagation block the real status is already saved; `:-` keeps it."""
    body = (
        "false\n"
        "BP_MAIN_RC=$?\n"                      # propagation captured the truth
        "true\n"                               # ... then something succeeded
        "BP_MAIN_RC=${BP_MAIN_RC:-$?}\n"       # footer must NOT clobber it
        'if [ "${BP_MAIN_RC:-0}" -ne 0 ]; then exit "${BP_MAIN_RC}"; fi\n'
    )
    assert _run_script(tmp_path, body).returncode != 0


def test_both_emitters_carry_the_capture():
    """Pin the emitted text, so a refactor cannot silently drop it."""
    import inspect
    from biopipelines.base_config import BaseConfig

    footer = inspect.getsource(BaseConfig.generate_completion_check_footer)
    assert "BP_MAIN_RC=${{BP_MAIN_RC:-$?}}" in footer
    assert 'exit "${{BP_MAIN_RC}}"' in footer

    propagation = inspect.getsource(BaseConfig.generate_missing_propagation)
    assert "BP_MAIN_RC=$?" in propagation


# ── the vacuous-check half ────────────────────────────────────────────────────

def _check(tmp_path, expected_outputs) -> subprocess.CompletedProcess:
    """Drive pipe_check_completion over a declaration written to disk."""
    folder = tmp_path / "007_Thing"
    folder.mkdir(parents=True, exist_ok=True)
    manifest = folder / ".expected_outputs.json"
    manifest.write_text(json.dumps({
        "tool_name": "Thing", "tool_class": "Thing",
        "output_structure": expected_outputs,
    }), encoding="utf-8")
    checker = Path(__file__).resolve().parents[1] / "pipe_scripts" / "pipe_check_completion.py"
    import sys
    return subprocess.run(
        [sys.executable, str(checker), str(folder), "Thing", str(manifest)],
        capture_output=True, text=True,
    )


def _value_stream(folder: Path, name: str):
    return {
        "name": name, "ids": ["a", "b"], "files": [],
        "map_table": str(folder / name / f"{name}_map.csv"),
        "format": "csv", "metadata": {},
    }


def test_a_value_based_stream_with_no_map_table_fails(tmp_path):
    """The whole content of such a stream is the map; absent means nothing ran."""
    folder = tmp_path / "007_Thing"
    folder.mkdir(parents=True, exist_ok=True)
    out = {"output_folder": str(folder), "out": _value_stream(folder, "out")}

    proc = _check(tmp_path, out)

    assert proc.returncode != 0, "a value-based stream that wrote nothing passed"
    assert "out" in proc.stdout


def test_a_value_based_stream_with_its_map_table_passes(tmp_path):
    folder = tmp_path / "007_Thing"
    (folder / "out").mkdir(parents=True, exist_ok=True)
    (folder / "out" / "out_map.csv").write_text("id,value\na,1\nb,2\n", encoding="utf-8")
    out = {"output_folder": str(folder), "out": _value_stream(folder, "out")}

    proc = _check(tmp_path, out)

    assert proc.returncode == 0, proc.stdout + proc.stderr


def test_a_file_based_stream_is_unaffected(tmp_path):
    """A stream that declares files keeps the per-id check; its map is not required."""
    folder = tmp_path / "007_Thing"
    (folder / "s").mkdir(parents=True, exist_ok=True)
    for i in ("a", "b"):
        (folder / "s" / f"{i}.pdb").write_text("ATOM\nEND\n", encoding="utf-8")
    out = {
        "output_folder": str(folder),
        "s": {
            "name": "s", "ids": ["a", "b"],
            "files": [str(folder / "s" / "<id>.pdb")],
            # Deliberately absent on disk: a file stream is judged by its files.
            "map_table": str(folder / "s" / "s_map.csv"),
            "format": "pdb", "metadata": {},
        },
    }

    proc = _check(tmp_path, out)

    assert proc.returncode == 0, proc.stdout + proc.stderr
