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


def test_exclude_sites_reaches_the_runtime_config(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline,
                     exclude_sites=["GAATTC", "GGATCC"])
    assert '"exclude_sites"' in content
    assert "GAATTC" in content and "GGATCC" in content


def test_exclude_sites_defaults_to_empty(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline)
    assert '"exclude_sites": []' in content


def _reject(local_config, isolated_cwd, new_pipeline, sites):
    """Construct inside a pipeline, which is where Sequence yields a stream.

    Outside a pipeline context Sequence(...) returns the raw tool instance rather
    than a StandardizedOutput, so DNAEncoder rejects it on type before any
    exclude_sites validation runs -- which is why these two assertions passed
    vacuously on an error message about the wrong thing.
    """
    from biopipelines.sequence import Sequence
    from biopipelines.dna_encoder import DNAEncoder

    pipeline = new_pipeline("dna_params_reject")
    with pipeline:
        s = Sequence(seq="MKTAYIAKQ", type="protein", ids="p1")
        DNAEncoder(sequences=s, exclude_sites=sites)


def test_non_iupac_site_is_rejected(local_config, isolated_cwd, new_pipeline):
    """Catch a typo at config time, not after the cluster job starts."""
    with pytest.raises(ValueError, match="non-IUPAC"):
        _reject(local_config, isolated_cwd, new_pipeline, ["GAATTZ"])


def test_too_short_site_is_rejected(local_config, isolated_cwd, new_pipeline):
    """A 3 bp 'site' occurs every ~64 bp and would make encoding impossible."""
    with pytest.raises(ValueError, match="only 3 bases"):
        _reject(local_config, isolated_cwd, new_pipeline, ["GAA"])
