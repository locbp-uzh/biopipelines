"""Parameter coverage for biopipelines.alphafold.AlphaFold."""

import pytest

from ._helpers import assert_substrings_in, read_pipeline_sh


pytestmark = pytest.mark.tool_parameters


def _build(local_config, isolated_cwd, new_pipeline, **af_kwargs):
    from biopipelines.sequence import Sequence
    from biopipelines.alphafold import AlphaFold

    pipeline = new_pipeline("af_params")
    with pipeline:
        s = Sequence(seq="MKTAYIAKQRQISFVKSHFSRQLEERLGL", type="protein", ids="p1")
        AlphaFold(proteins=s, **af_kwargs)
        script_path = pipeline.save()
    return read_pipeline_sh(script_path)


def test_num_relax(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline, num_relax=2)
    assert "--num-relax 2" in content
    assert "--amber" in content


def test_num_recycle(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline, num_recycle=5)
    assert "--num-recycle 5" in content


def test_rand_seed(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline, rand_seed=2026)
    assert "--random-seed 2026" in content


def test_smoke_all_params(local_config, isolated_cwd, new_pipeline):
    content = _build(
        local_config, isolated_cwd, new_pipeline,
        num_relax=1,
        num_recycle=4,
        rand_seed=99,
    )
    assert_substrings_in(content, [
        "--num-relax 1",
        "--num-recycle 4",
        "--random-seed 99",
        "--amber",
    ])
