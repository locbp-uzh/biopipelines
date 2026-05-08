"""Parameter coverage for biopipelines.mmseqs2.MMseqs2."""

import pytest

from ._helpers import assert_substrings_in, read_pipeline_sh


pytestmark = pytest.mark.tool_parameters


def _build(local_config, isolated_cwd, new_pipeline, **kwargs):
    from biopipelines.sequence import Sequence
    from biopipelines.mmseqs2 import MMseqs2

    pipeline = new_pipeline("mmseqs_params")
    with pipeline:
        s = Sequence(seq="MKTAYIAKQRQISFVKSHFSRQLEERLGL", type="protein", ids="p1")
        MMseqs2(sequences=s, **kwargs)
        script_path = pipeline.save()
    return read_pipeline_sh(script_path)


def test_output_format_csv(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline, output_format="csv")
    assert "--output_format csv" in content


def test_output_format_a3m(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline, output_format="a3m")
    assert "--output_format a3m" in content


def test_mask(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline, mask="10-20+30-40")
    assert "10-20" in content or "mask" in content.lower()


def test_smoke_all_params(local_config, isolated_cwd, new_pipeline):
    content = _build(
        local_config, isolated_cwd, new_pipeline,
        output_format="a3m",
        mask="5-15",
    )
    assert_substrings_in(content, ["--output_format a3m"])
