# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""DiffDock, NeuralPLexer and DynamicBind: a ligand with a blank SMILES cell is named and refused, it is not handed to the model as "nan".

For a ligand stream without coordinate files each script reads the `smiles` cell and refuses a ligand that has none ("Compound ... has neither a file nor a 'smiles' value"). A blank cell comes back from pandas as NaN, which is truthy, so the refusal never fired: DiffDock and DynamicBind wrote the string "nan" as the ligand into the table the model reads, and NeuralPLexer passed a float to RDKit and died with a TypeError that names no ligand.
"""

import importlib.util
import os

import pytest

pytest.importorskip("pandas")

from biopipelines.datastream import DataStream, create_map_table

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

SCRIPTS = ["pipe_diffdock_build_csv.py", "pipe_neuralplexer_stage.py", "pipe_dynamicbind_build_csv.py"]


def _script(filename):
    spec = importlib.util.spec_from_file_location(f"{filename}_under_test", os.path.join(REPO_ROOT, "pipe_scripts", filename))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _compounds(tmp_path, smiles):
    ids = [f"lig_{n}" for n in range(len(smiles))]
    table = tmp_path / "compounds.csv"
    create_map_table(str(table), ids=ids, additional_columns={"smiles": smiles})
    DataStream(name="compounds", ids=ids, map_table=str(table), format="csv", files=[]).save_json(str(tmp_path / "compounds.json"))
    return str(tmp_path / "compounds.json")


def _smiles_of(collected):
    """The SMILES per ligand id from what collect_ligands returns, whichever tuple shape the script uses."""
    return {one[0]: next(part for part in one[1:] if isinstance(part, str)) for one in collected}


@pytest.mark.parametrize("filename", SCRIPTS)
class TestALigandWithoutASmiles:

    def test_every_ligand_with_a_smiles_is_read(self, filename, tmp_path):
        collected = _script(filename).collect_ligands(_compounds(tmp_path, ["CCO", "CCN"]))

        assert _smiles_of(collected) == {"lig_0": "CCO", "lig_1": "CCN"}

    def test_a_blank_cell_is_refused_and_names_the_ligand(self, filename, tmp_path):
        with pytest.raises(ValueError, match="lig_1"):
            _script(filename).collect_ligands(_compounds(tmp_path, ["CCO", ""]))

    def test_a_ligand_with_every_cell_blank_is_refused_too(self, filename, tmp_path):
        with pytest.raises(ValueError, match="lig_0"):
            _script(filename).collect_ligands(_compounds(tmp_path, ["", ""]))
