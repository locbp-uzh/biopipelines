"""Parameter coverage for biopipelines.dna_encoder.DNAEncoder."""

import pytest

from ._helpers import read_all_emitted_artifacts


pytestmark = pytest.mark.tool_parameters


def _build(local_config, isolated_cwd, new_pipeline, **kwargs):
    from biopipelines.sequence import Sequence
    from biopipelines.dna_encoder import DNAEncoder

    pipeline = new_pipeline("dna_params")
    with pipeline:
        s = Sequence(seq="MKTAYIAKQ", type="protein", ids="p1")
        DNAEncoder(sequences=s, **kwargs)
        script_path = pipeline.save()
    return read_all_emitted_artifacts(script_path)


def test_organism_EC(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline, organism="EC")
    assert '"organism": "EC"' in content


def test_organism_HS(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline, organism="HS")
    assert '"organism": "HS"' in content


def test_organism_combined(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline, organism="EC&HS")
    assert '"organism": "EC&HS"' in content


def test_smoke_all_params(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline, organism="SC")
    assert '"organism": "SC"' in content
