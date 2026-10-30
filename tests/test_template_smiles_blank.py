# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""PocketGen and ProLIF: a blank SMILES cell is no template, it is not the template "nan".

`load_template_smiles` read the `smiles` cell of the row whose `code` matches with `str(...)`. A blank cell is stored empty and read back as NaN, and `str(NaN)` is the non-empty string "nan", so the "has no SMILES" check passed and "nan" went on as the bond-order template; RDKit then failed with "failed to parse template SMILES: 'nan'" instead of the message that names the missing SMILES (GEMS and RTMScore had the same read, fixed in their own tests).
"""

import importlib.util
import os
import sys
from unittest import mock

import pytest

pytest.importorskip("pandas")

from biopipelines.datastream import DataStream, create_map_table

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _script(monkeypatch, filename):
    # pipe_prolif imports MDAnalysis and prolif at module level; only the SMILES read is under test
    for name in ("MDAnalysis", "prolif"):
        monkeypatch.setitem(sys.modules, name, mock.MagicMock())
    spec = importlib.util.spec_from_file_location(f"{filename}_under_test", os.path.join(REPO_ROOT, "pipe_scripts", filename))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _compounds(tmp_path, codes, smiles):
    ids = [f"lig_{n}" for n in range(len(codes))]
    table = tmp_path / "compounds.csv"
    create_map_table(str(table), ids=ids, additional_columns={"code": codes, "smiles": smiles})
    DataStream(name="compounds", ids=ids, map_table=str(table), format="csv", files=[]).save_json(str(tmp_path / "compounds.json"))
    return str(tmp_path / "compounds.json")


@pytest.mark.parametrize("filename", ["pipe_pocketgen_stage.py", "pipe_prolif.py"])
class TestTheTemplateOfTheLigandWithACode:

    def test_the_smiles_of_the_matching_row_is_the_template(self, filename, monkeypatch, tmp_path):
        script = _script(monkeypatch, filename)

        assert script.load_template_smiles(_compounds(tmp_path, ["AAA", "GG2"], ["CCN", "CCO"]), "gg2") == "CCO"

    def test_a_blank_cell_is_reported_as_a_missing_smiles(self, filename, monkeypatch, tmp_path):
        script = _script(monkeypatch, filename)

        with pytest.raises(ValueError, match="no SMILES"):
            script.load_template_smiles(_compounds(tmp_path, ["GG2"], [""]), "GG2")

    def test_a_code_no_row_has_is_reported(self, filename, monkeypatch, tmp_path):
        script = _script(monkeypatch, filename)

        with pytest.raises(ValueError, match="no row with code"):
            script.load_template_smiles(_compounds(tmp_path, ["AAA"], ["CCN"]), "GG2")
