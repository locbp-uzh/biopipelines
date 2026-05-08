"""Parameter coverage for biopipelines.rbs_designer.RBSDesigner."""

import pytest

from ._helpers import assert_substrings_in, read_all_emitted_artifacts


pytestmark = pytest.mark.tool_parameters


def _build(local_config, isolated_cwd, new_pipeline, **kwargs):
    from biopipelines.sequence import Sequence
    from biopipelines.rbs_designer import RBSDesigner

    pipeline = new_pipeline("rbs_params")
    with pipeline:
        s = Sequence(seq="ATGAAAGCATCAGCT", type="dna", ids="g1")
        RBSDesigner(sequences=s, **kwargs)
        script_path = pipeline.save()
    return read_all_emitted_artifacts(script_path)


def test_tir_preset_low(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline, tir="low")
    assert '"tir": 100' in content


def test_tir_preset_high(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline, tir="high")
    assert '"tir": 10000' in content


def test_tir_numeric(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline, tir=2500)
    assert '"tir": 2500' in content


def test_pre_sequence(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline, pre_sequence="AAGGAG")
    assert '"pre_sequence": "AAGGAG"' in content


def test_add_start_codon(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline, add_start_codon=True)
    assert '"add_start_codon": true' in content


def test_smoke_all_params(local_config, isolated_cwd, new_pipeline):
    content = _build(
        local_config, isolated_cwd, new_pipeline,
        tir="medium",
        pre_sequence="GAATTC",
        add_start_codon=True,
    )
    assert_substrings_in(content, [
        '"tir": 1000',
        '"pre_sequence": "GAATTC"',
        '"add_start_codon": true',
    ])
