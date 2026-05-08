"""Parameter coverage for biopipelines.rfdiffusion3.RFdiffusion3."""

import pytest

from ._helpers import assert_substrings_in, read_all_emitted_artifacts


pytestmark = pytest.mark.tool_parameters


def _build(local_config, isolated_cwd, new_pipeline, **kwargs):
    from biopipelines.rfdiffusion3 import RFdiffusion3

    kwargs.setdefault("length", 80)

    pipeline = new_pipeline("rfd3_params")
    with pipeline:
        RFdiffusion3(**kwargs)
        script_path = pipeline.save()
    return read_all_emitted_artifacts(script_path)


def test_length(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline, length=120)
    assert "120" in content


def test_num_designs(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline, num_designs=8)
    assert "8" in content


def test_num_models(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline, num_models=2)
    assert "2" in content


def test_design_startnum(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline, design_startnum=99)
    assert "99" in content


def test_smoke_all_params(local_config, isolated_cwd, new_pipeline):
    content = _build(
        local_config, isolated_cwd, new_pipeline,
        length=100,
        num_designs=4,
        num_models=3,
        design_startnum=10,
    )
    assert_substrings_in(content, ["100", "10"])
