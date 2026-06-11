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


def test_symmetry_string(local_config, isolated_cwd, new_pipeline):
    # String shorthand -> {"id": "C3"} in the JSON template, and the symmetry
    # sampler is auto-selected on the run line.
    content = _build(local_config, isolated_cwd, new_pipeline, length=110, symmetry="C3")
    assert_substrings_in(content, ['"symmetry"', '"id"', "C3", "inference_sampler.kind=symmetry"])


def test_cfg_and_sampler_knobs(local_config, isolated_cwd, new_pipeline):
    # Diffused-ligand binder paper settings map to inference_sampler.* overrides.
    content = _build(
        local_config, isolated_cwd, new_pipeline,
        length=100, cfg=True, cfg_scale=2.0,
        step_scale=1.5, noise_scale=0.6, num_steps=200,
    )
    assert_substrings_in(content, [
        "inference_sampler.use_classifier_free_guidance=True",
        "inference_sampler.cfg_scale=2.0",
        "inference_sampler.step_scale=1.5",
        "inference_sampler.gamma_0=0.6",
        "inference_sampler.num_timesteps=200",
    ])


def test_unindex_motif(local_config, isolated_cwd, new_pipeline):
    # Atomic-motif enzyme path: unindexed motif key in the JSON template.
    content = _build(local_config, isolated_cwd, new_pipeline, contig="80-150", unindex="A50-52")
    assert_substrings_in(content, ['"unindex"', "A50-52"])


def test_partial_diffusion_rejects_length(local_config, isolated_cwd, new_pipeline):
    import pytest as _pytest
    with _pytest.raises(Exception):
        _build(local_config, isolated_cwd, new_pipeline, length=100, partial_t=10.0)


def test_coord_selector_without_input_rejected(local_config, isolated_cwd, new_pipeline):
    # select_buried needs coordinates; de-novo (no pdb, bare-code ligand) has none.
    import pytest as _pytest
    with _pytest.raises(Exception):
        _build(local_config, isolated_cwd, new_pipeline, length=100, select_buried="B1")


def test_bare_ligand_without_input_rejected(local_config, isolated_cwd, new_pipeline):
    # foundry appends a ligand to an input atom array; a bare code with
    # length-only design has nothing to bind to.
    from biopipelines.ligand import Ligand
    import pytest as _pytest
    with _pytest.raises(Exception):
        _build(local_config, isolated_cwd, new_pipeline, length=100, ligand=Ligand(code="SAM"))
