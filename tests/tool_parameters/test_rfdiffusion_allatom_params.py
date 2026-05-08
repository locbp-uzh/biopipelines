"""Parameter coverage for biopipelines.rfdiffusion_allatom.RFdiffusionAllAtom."""

import pytest

from ._helpers import assert_substrings_in, read_pipeline_sh


pytestmark = pytest.mark.tool_parameters


def _build(local_config, isolated_cwd, new_pipeline, ligand="ATP", **kwargs):
    from biopipelines.rfdiffusion_allatom import RFdiffusionAllAtom

    pipeline = new_pipeline("rfdaa_params")
    with pipeline:
        RFdiffusionAllAtom(ligand=ligand, contigs=kwargs.pop("contigs", "100-200"), **kwargs)
        script_path = pipeline.save()
    return read_pipeline_sh(script_path)


def test_ligand(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline, ligand="HEM")
    assert "HEM" in content


def test_num_designs(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline, num_designs=6)
    assert "num_designs=6" in content


def test_steps(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline, steps=80)
    assert "diffuser.T=80" in content


def test_num_recycles(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline, num_recycles=3)
    assert "num_recycles=3" in content or "recycles" in content


def test_smoke_all_params(local_config, isolated_cwd, new_pipeline):
    content = _build(
        local_config, isolated_cwd, new_pipeline,
        contigs="100-150",
        num_designs=4,
        steps=100,
        num_recycles=2,
        deterministic=True,
        design_startnum=5,
    )
    assert_substrings_in(content, [
        "num_designs=4",
        "diffuser.T=100",
        "design_startnum=5",
    ])
