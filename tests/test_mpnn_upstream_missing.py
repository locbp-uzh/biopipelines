"""The MPNN tools must find their upstream missing manifest.

_get_upstream_missing_table_path reads `.tables`, which a DataStream does not have -- only
the StandardizedOutput it came from does. Both MPNN tools resolved `structures` to a stream
and then handed that stream to the helper, so it returned None and no manifest propagated.

LigandMPNN is where that showed: its `fasta` stream declares one .fa per DECLARED id, so the
completion check demanded a file for every design an upstream filter had dropped -- 1200
expected against 126 produced from 127 survivors, reported as FAILED although the run had
succeeded. ProteinMPNN shares the defect but hides it, having no per-id file stream.
"""

import pytest

from biopipelines.ligand_mpnn import LigandMPNN
from biopipelines.protein_mpnn import ProteinMPNN


@pytest.mark.parametrize("tool", [ProteinMPNN, LigandMPNN])
def test_original_input_is_retained(tool):
    """Whatever was passed as `structures` must survive on the instance.

    Keeping only the resolved DataStream loses `.tables`, and with it any route to the
    upstream missing manifest.
    """
    import inspect
    body = inspect.getsource(tool.__init__)
    assert "self.structures_input = structures" in body, (
        f"{tool.__name__} must keep the original `structures` argument; a DataStream "
        f"alone has no .tables and the missing manifest cannot be found from it"
    )


@pytest.mark.parametrize("mod", ["biopipelines.ligand_mpnn", "biopipelines.protein_mpnn"])
def test_helper_is_called_with_the_original_input(mod):
    """The manifest lookup must use the retained input, not the resolved stream."""
    import importlib, inspect
    src = inspect.getsource(importlib.import_module(mod))
    i = src.find("_get_upstream_missing_table_path(")
    assert i != -1, "expected an upstream-missing lookup"
    call = src[i:i + 200]
    assert "self.structures_input" in call, (
        "the lookup must receive the StandardizedOutput; passing structures_stream "
        "returns None because a DataStream has no .tables"
    )
    assert "self.structures_stream" not in call.split(")")[0], (
        "passing the resolved stream reintroduces the bug"
    )
