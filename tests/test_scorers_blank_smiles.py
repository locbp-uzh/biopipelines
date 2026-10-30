# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""GEMS and RTMScore: a ligand without a SMILES has no template, it does not have the template "nan".

`Ligand(codes=..., structures=...)` without `template_smiles` carves a ligand whose `smiles` cell is blank, which the compounds table stores empty and pandas reads back as NaN. The scorers turned that cell into the string "nan" and handed it to RDKit as a bond-order template, so the whole step died on the first carved ligand ("RDKit failed to parse template SMILES: 'nan'") instead of perceiving the bonds from the coordinates, which is what a ligand without chemistry is documented to get.
"""

import importlib.util
import os

import pytest

pytest.importorskip("pandas")

from biopipelines.datastream import DataStream, create_map_table

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _script(filename):
    spec = importlib.util.spec_from_file_location(f"{filename}_under_test", os.path.join(REPO_ROOT, "pipe_scripts", filename))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _compounds(tmp_path, smiles):
    ids = [f"lig_{n}" for n in range(len(smiles))]
    table = tmp_path / "compounds.csv"
    create_map_table(str(table), ids=ids, additional_columns={"smiles": smiles})
    stream = DataStream(name="compounds", ids=ids, map_table=str(table), format="csv", files=[])
    stream.save_json(str(tmp_path / "compounds.json"))
    return str(tmp_path / "compounds.json")


@pytest.mark.parametrize("filename", ["pipe_gems.py", "pipe_rtmscore.py"])
class TestAScorerReadsTheTemplateOfALigand:

    def test_a_blank_cell_is_no_template(self, filename, tmp_path):
        lookup = _script(filename).load_smiles_lookup(_compounds(tmp_path, ["", "CCO"]))

        assert lookup == {"lig_1": "CCO"}

    def test_a_ligand_with_every_cell_blank_has_no_lookup_at_all(self, filename, tmp_path):
        assert _script(filename).load_smiles_lookup(_compounds(tmp_path, ["", ""])) == {}

    def test_no_stream_is_no_lookup(self, filename):
        assert _script(filename).load_smiles_lookup("") == {}
