# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Constructing a Pipeline without saying `on_the_fly` or `local_output`.

Every other test in the suite passes both explicitly, so the branch that decides them was never executed -- 1722 tests passed while `_detect_notebook` had lost its `@staticmethod` to a neighbouring extraction and raised `TypeError` on the first real submission. A default that only the real entry point reaches is a default nothing checks.
"""

import pytest

from biopipelines.pipeline import Pipeline


def test_detect_notebook_is_callable_from_an_instance(local_config, isolated_cwd):
    """It is a staticmethod, so `self._detect_notebook()` must pass nothing. Losing the decorator turns that into `takes 0 positional arguments but 1 was given` -- but only on the path where on_the_fly is left unset."""
    pipeline = Pipeline(project="TestSuite", job="defaults_detect",
                        description="detect", config="local", on_the_fly=False,
                        local_output=True)
    assert pipeline._detect_notebook() in (True, False)


def test_a_pipeline_can_be_built_without_saying_on_the_fly(local_config, isolated_cwd):
    """The real entry point: `./submit` runs a script whose Pipeline(...) names neither flag."""
    pipeline = Pipeline(project="TestSuite", job="defaults_unset",
                        description="no on_the_fly, no local_output", config="local")
    assert pipeline.on_the_fly in (True, False)
    assert pipeline.local_output in (True, False)


@pytest.mark.parametrize("override, expected", [("1", True), ("0", False)])
def test_the_environment_override_decides_local_output(monkeypatch, override, expected):
    """A single-node container backend repoints biopipelines_output at a persistent mount and sets this so results are not diverted to a cwd that dies with the container."""
    monkeypatch.setenv("BIOPIPELINES_LOCAL_OUTPUT", override)
    assert Pipeline._default_local_output(on_the_fly=True) is expected
    assert Pipeline._default_local_output(on_the_fly=False) is expected


def test_without_an_override_an_interactive_run_defaults_to_local(monkeypatch, local_config):
    monkeypatch.delenv("BIOPIPELINES_LOCAL_OUTPUT", raising=False)
    assert Pipeline._default_local_output(on_the_fly=True) is True
    assert Pipeline._default_local_output(on_the_fly=False) is False
