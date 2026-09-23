"""Pre-computed MSAs cannot reach a multi-chain ColabFold query, so asking for them is refused.

ColabFold folds a complex with its own paired+unpaired pipeline and never reads a per-chain a3m, so `msas=` behind a `Bundle` was already refused. A `Grouped` complex was not, and its MSAs were dropped with only a warning.
"""

import pytest


@pytest.fixture
def af_inputs(local_config, isolated_cwd):
    from biopipelines.mock import Mock
    from biopipelines.pipeline import Pipeline
    from biopipelines.protein_mpnn import ProteinMPNN

    pipeline = Pipeline(project="TestSuite", job="af_msas", description="x",
                        on_the_fly=False, local_output=True, config="local")
    pipeline.__enter__()
    source = Mock(ids=["bb"], streams={"structures": {"format": "pdb", "file": "<id>.pdb"}})
    mpnn = ProteinMPNN(structures=source.streams.structures, num_sequences=2, chains=["A", "B"])
    msas = Mock(ids=["bb_1_A"], streams={"msas": {"format": "a3m", "file": "<id>.a3m"}})
    yield mpnn, msas
    pipeline.__exit__(None, None, None)


def test_msas_with_a_grouped_complex_is_refused(af_inputs):
    from biopipelines.alphafold import AlphaFold
    from biopipelines.combinatorics import Grouped
    mpnn, msas = af_inputs
    with pytest.raises(ValueError, match="multi-chain complex"):
        AlphaFold(proteins=Grouped(mpnn), msas=msas)


def test_a_grouped_member_of_a_bundle_is_detected():
    from biopipelines.alphafold import AlphaFold
    from biopipelines.combinatorics import Bundle, Each, Grouped
    assert AlphaFold._contains_grouped(Bundle(Each(Grouped("x")), "y"))
    assert not AlphaFold._contains_grouped(Each("x"))


def test_a_bare_grouped_complex_is_accepted_without_msas(af_inputs):
    from biopipelines.alphafold import AlphaFold
    from biopipelines.combinatorics import Grouped
    mpnn, _msas = af_inputs
    AlphaFold(proteins=Grouped(mpnn))
