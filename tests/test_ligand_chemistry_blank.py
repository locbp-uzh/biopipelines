"""`ligand_chemistry` reads every serialized spelling of a blank cell as empty.

`resolve_ligand_smiles` already took `<NA>` (the spelling of a pandas nullable string) as blank while `_clean`, which
`ligand_chemistry` is built on, did not, so a cell holding it counted as chemistry for every consumer of the helper.
"""

import pytest

from biopipelines.ligand_utils import ligand_chemistry


@pytest.mark.parametrize("cell", ["NaN", "None", "<NA>", "<na>", "  ", "", float("nan")])
def test_every_spelling_of_blank_is_empty(cell):
    assert ligand_chemistry({"smiles": cell, "ccd": cell}) == {"ccd": "", "smiles": ""}


def test_a_real_smiles_and_code_are_kept_and_stripped():
    assert ligand_chemistry({"smiles": " CCO ", "ccd": "LIG"}) == {"ccd": "LIG", "smiles": "CCO"}
