"""Parameter coverage for biopipelines.boltzgen.BoltzGen."""

import pytest

from ._helpers import assert_substrings_in, read_all_emitted_artifacts


pytestmark = pytest.mark.tool_parameters


def _build(local_config, isolated_cwd, new_pipeline, **kwargs):
    from biopipelines.mock import Mock
    from biopipelines.boltzgen import BoltzGen

    pipeline = new_pipeline("bg_params")
    with pipeline:
        target = Mock(
            ids=["t1"],
            streams={"structures": {"format": "pdb", "file": "<id>.pdb"}},
            map_table_strategy="config",
        )
        BoltzGen(
            target_structure=target.streams.structures,
            binder_spec="80",
            **kwargs,
        )
        script_path = pipeline.save()
    return read_all_emitted_artifacts(script_path)


def test_protocol(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline, protocol="protein-anything")
    assert "--protocol protein-anything" in content


def test_num_designs(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline, num_designs=500)
    assert "--num_designs 500" in content


def test_budget(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline, budget=42)
    assert "--budget 42" in content


def test_step_scale(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline, step_scale=1.25)
    assert "--step_scale 1.25" in content


def test_noise_scale(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline, noise_scale=0.6)
    assert "--noise_scale 0.6" in content


def test_diffusion_batch_size(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline, diffusion_batch_size=8)
    assert "--diffusion_batch_size 8" in content


def test_inverse_fold_num_sequences(local_config, isolated_cwd, new_pipeline):
    content = _build(
        local_config, isolated_cwd, new_pipeline, inverse_fold_num_sequences=4
    )
    assert "--inverse_fold_num_sequences 4" in content


def test_skip_inverse_folding(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline, skip_inverse_folding=True)
    assert "--skip_inverse_folding" in content


def test_alpha(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline, alpha=0.7)
    assert "--alpha 0.7" in content


def test_refolding_rmsd_threshold(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline, refolding_rmsd_threshold=2.5)
    assert "--refolding_rmsd_threshold 2.5" in content


def test_filter_biased(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline, filter_biased=False)
    assert "--filter_biased false" in content


def test_smoke_all_params(local_config, isolated_cwd, new_pipeline):
    content = _build(
        local_config, isolated_cwd, new_pipeline,
        protocol="protein-anything",
        num_designs=200,
        budget=20,
        step_scale=1.0,
        noise_scale=0.5,
        diffusion_batch_size=4,
        inverse_fold_num_sequences=2,
        alpha=0.5,
        refolding_rmsd_threshold=2.0,
    )
    assert_substrings_in(content, [
        "--protocol protein-anything",
        "--num_designs 200",
        "--budget 20",
        "--step_scale 1.0",
        "--noise_scale 0.5",
        "--diffusion_batch_size 4",
        "--inverse_fold_num_sequences 2",
        "--alpha 0.5",
        "--refolding_rmsd_threshold 2.0",
    ])
