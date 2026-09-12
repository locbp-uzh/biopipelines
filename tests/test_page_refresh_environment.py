# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""The run page is written at configuration time, when nothing has run, so the copy a user opens is only a plan until the job refreshes it at the end. That refresh is a `python -c` in the job shell, where the only thing applied is the scheduler's module load -- its bare `python` has no pandas, so an unqualified interpreter imports nothing and the page silently stays a plan. These tests hold the emitted bash to activating the framework environment first and to reporting a failure rather than discarding it.

Under `pip` there is no environment to activate and the ambient interpreter is the right one, so the manager is what decides; the tests say so for both.
"""

import re

import pytest

from biopipelines.config_manager import ConfigManager


def _saved_pipeline(new_pipeline, name):
    from biopipelines.mock import Mock

    pipeline = new_pipeline(name)
    with pipeline:
        Mock(ids=["a"], streams={"structures": {"format": "pdb", "file": "<id>.pdb"}})
        pipeline.save()
    return pipeline


def _refresh_block(new_pipeline, name, monkeypatch, env_manager="mamba"):
    monkeypatch.setattr(ConfigManager, "get_env_manager", lambda self: env_manager)
    return "\n".join(_saved_pipeline(new_pipeline, name)._page_refresh_lines())


def test_the_refresh_activates_the_framework_environment(local_config, isolated_cwd, new_pipeline, monkeypatch):
    block = _refresh_block(new_pipeline, "refresh_env", monkeypatch)
    assert "python -c" in block
    activate = re.search(r"^\s*(mamba|conda|micromamba) activate (\S+)", block, re.M)
    assert activate, f"the refresh runs an unqualified python:\n{block}"
    assert activate.group(2) == "biopipelines"
    assert block.index(activate.group(0)) < block.index("python -c"), "activation must precede the interpreter"


def test_under_pip_the_ambient_interpreter_is_the_right_one(local_config, isolated_cwd, new_pipeline, monkeypatch):
    block = _refresh_block(new_pipeline, "refresh_pip", monkeypatch, env_manager="pip")
    assert "python -c" in block
    assert re.search(r"^\s*(mamba|conda|micromamba) activate", block, re.M) is None


@pytest.mark.parametrize("env_manager", ["mamba", "pip"])
def test_a_failed_refresh_is_reported_not_discarded(local_config, isolated_cwd, new_pipeline, monkeypatch, env_manager):
    """A refresh whose output went to /dev/null made an unrefreshed page indistinguishable from a run that produced nothing, which is how this stayed broken."""
    block = _refresh_block(new_pipeline, f"refresh_loud_{env_manager}", monkeypatch, env_manager)
    assert ">/dev/null" not in block
    assert "WARNING" in block


@pytest.mark.parametrize("env_manager", ["mamba", "pip"])
def test_the_refresh_cannot_leak_its_environment_into_the_run(local_config, isolated_cwd, new_pipeline, monkeypatch, env_manager):
    """Activation happens in a subshell, so nothing the refresh sources or activates outlives it."""
    block = _refresh_block(new_pipeline, f"refresh_subshell_{env_manager}", monkeypatch, env_manager)
    lines = [ln.strip() for ln in block.splitlines()]
    assert "(" in lines and ")" in lines
    assert lines.index("(") < lines.index(")")


def test_no_tools_means_no_refresh(local_config, isolated_cwd, new_pipeline):
    """The activation snippet comes from a tool, so a pipeline with none has nothing to emit -- and nothing to refresh either."""
    assert new_pipeline("refresh_empty")._page_refresh_lines() == []
