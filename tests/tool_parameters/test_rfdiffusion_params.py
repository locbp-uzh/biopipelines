"""Parameter coverage for biopipelines.rfdiffusion.RFdiffusion.

Verifies the kwargs CURRENTLY exposed today (Track 2 — no implementation
gaps closed for this tool, but the regression suite locks down what IS
exposed).
"""

import pytest

from ._helpers import assert_substrings_in, read_pipeline_sh


pytestmark = pytest.mark.tool_parameters


def _build(local_config, isolated_cwd, new_pipeline, contigs="100-200", **kwargs):
    from biopipelines.rfdiffusion import RFdiffusion

    pipeline = new_pipeline("rfd_params")
    with pipeline:
        RFdiffusion(contigs=contigs, **kwargs)
        script_path = pipeline.save()
    return read_pipeline_sh(script_path)


def test_contigs(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline, contigs="50-100")
    assert "contigmap.contigs=['50-100']" in content


def test_num_designs(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline, num_designs=8)
    assert "inference.num_designs=8" in content


def test_steps(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline, steps=25)
    assert "diffuser.T=25" in content


def test_reproducible(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline, reproducible=True)
    assert "inference.deterministic=True" in content


def test_design_startnum(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline, design_startnum=42)
    assert "inference.design_startnum=42" in content


def test_active_site(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline, active_site=True)
    assert "ActiveSite_ckpt" in content


def test_inpaint_seq(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline, inpaint="A10-20")
    assert "contigmap.inpaint_seq=['A10-20']" in content


def test_inpaint_str(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline, inpaint_str="A30-40")
    assert "contigmap.inpaint_str=['A30-40']" in content


def test_smoke_all_params(local_config, isolated_cwd, new_pipeline):
    content = _build(
        local_config, isolated_cwd, new_pipeline,
        contigs="A1-50/30-50",
        num_designs=4,
        steps=30,
        reproducible=True,
        design_startnum=10,
        active_site=False,
    )
    assert_substrings_in(content, [
        "contigmap.contigs=['A1-50/30-50']",
        "inference.num_designs=4",
        "diffuser.T=30",
        "inference.deterministic=True",
        "inference.design_startnum=10",
    ])
