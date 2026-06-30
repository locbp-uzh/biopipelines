"""Helpers for the tool_parameters test suite.

Each test in this suite drives a single-tool ``Pipeline`` at config time and
asserts that user-facing kwargs land in the emitted artifacts. Most wrappers
inline CLI flags into generated shell scripts, while others serialize runtime
helper arguments into configuration files. The helpers below support grepping
either surface depending on what a wrapper emits.
"""

import os
from pathlib import Path
from typing import Iterable


def read_pipeline_sh(script_path: str) -> str:
    """Read pipeline.sh and every per-step ``RunTime/NNN_<Tool>.sh`` script in
    its sibling directory, returning their concatenated contents.

    Pipeline.save() emits a master ``pipeline.sh`` plus one ``NNN_<Tool>.sh``
    per step under the same ``RunTime/`` folder. Tool-specific CLI flags land
    in the per-step scripts, not the master, so the tool_parameters tests
    grep across all of them.
    """
    assert os.path.isfile(script_path), f"pipeline.sh missing: {script_path}"
    runtime_dir = Path(script_path).parent
    parts = [Path(script_path).read_text(encoding="utf-8", errors="replace")]
    for sub in sorted(runtime_dir.glob("*.sh")):
        if sub.name == "pipeline.sh":
            continue
        parts.append(sub.read_text(encoding="utf-8", errors="replace"))
    return "\n".join(parts)


def assert_substrings_in(content: str, substrings: Iterable[str], label: str = ""):
    """Assert every substring appears in content; raise with the missing ones listed."""
    missing = [s for s in substrings if s not in content]
    assert not missing, (
        f"{label or 'pipeline.sh'} missing expected substrings: {missing}"
    )


def read_all_emitted_artifacts(script_path: str) -> str:
    """Read pipeline.sh, every per-step ``RunTime/*.sh``, and every config
    file (JSON, YAML, CSV, JSONL, A3M) emitted under the run's tool output
    folders. Returns their concatenated contents.

    Some wrappers (Gnina, RFdiffusion3) write parameter values into a JSON
    config file at ``configuration/`` rather than inlining them in bash. The
    tool_parameters tests need to grep across both surfaces.
    """
    base = read_pipeline_sh(script_path)
    parts = [base]

    # The run root is the parent of RunTime/.
    run_root = Path(script_path).resolve().parent.parent
    if run_root.exists():
        for cfg in run_root.rglob("_configuration/*"):
            if cfg.is_file() and cfg.suffix.lower() in (".json", ".yaml", ".yml", ".csv", ".jsonl", ".a3m", ".txt"):
                try:
                    parts.append(cfg.read_text(encoding="utf-8", errors="replace"))
                except OSError:
                    continue
    return "\n".join(parts)


def assert_kwarg_emitted(content: str, kwarg_name: str, value, *, flag: str = None):
    """Assert that kwarg=value reached the emitted bash.

    ``flag`` is the upstream CLI flag string the value should appear next to;
    when ``flag`` is omitted we only check for the value.
    """
    str_value = str(value)
    assert str_value in content, (
        f"value {str_value!r} for kwarg {kwarg_name!r} not found in pipeline.sh"
    )
    if flag is not None:
        assert flag in content, (
            f"flag {flag!r} for kwarg {kwarg_name!r} not found in pipeline.sh"
        )
