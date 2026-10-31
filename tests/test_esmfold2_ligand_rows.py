# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""ESMFold2 reads a ligand out of the compounds table the way the table is written.

A SMILES compound carries a residue `code` (LIG unless told otherwise) and a blank `ccd`.
Reading `code` as the CCD entry folded every SMILES ligand as the CCD entry LIG, with pLDDT and ipTM reported for the wrong molecule.
Reading the blank cells with pandas' default NA handling turned them into NaN, which is truthy, and a NaN CCD code crashed inside esm.
These tests pin the routing and the reading of blank cells through `build_complexes` and `to_structure_input`, with a stub in place of the esm input builder.

Not covered: that esm accepts the SMILES or the CCD code it is handed, which needs the model environment.
"""

import importlib.util
import math
import os
from dataclasses import dataclass
from typing import List, Optional

import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
INFERENCE = os.path.join(REPO_ROOT, "pipe_scripts", "pipe_esmfold2_inference.py")

SEQ = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ"
COLUMNS = "id,format,code,lookup,source,ccd,cid,cas,smiles,name,formula,file_path"


@dataclass
class _LigandInput:
    id: str
    smiles: Optional[str] = None
    ccd: Optional[List[str]] = None


@dataclass
class _ProteinInput:
    id: str
    sequence: str
    msa: object = None
    modifications: object = None


class _StubBuilder:
    ProteinInput = _ProteinInput
    LigandInput = _LigandInput
    CovalentBond = None
    StructurePredictionInput = lambda sequences, covalent_bonds=None: sequences  # noqa: E731


@pytest.fixture
def inference():
    spec = importlib.util.spec_from_file_location("pipe_esmfold2_inference", INFERENCE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _config(tmp_path, ligand_rows):
    proteins = tmp_path / "proteins.csv"
    proteins.write_text(f"id,sequence\np1,{SEQ}\n")
    ligands = tmp_path / "compounds.csv"
    ligands.write_text(COLUMNS + "\n" + "".join(row + "\n" for row in ligand_rows))
    return {"axes": {
        "proteins": {"mode": "each", "entity_type": "protein", "sources": [str(proteins)]},
        "ligands": {"mode": "each", "entity_type": "ligand", "sources": [str(ligands)]},
    }}


def _ligand_chains(inference, config):
    (_, chains), = inference.build_complexes(config, None)
    return [c for c in chains if c["entity_type"] == "ligand"]


SMILES_ROW = "lig,smiles,LIG,,,,,,CC(=O)Oc1ccccc1C(=O)O,,,"
CCD_ROW = "atp,ccd,ATP,ATP,ccd,ATP,,,,,,"


class TestTheCompoundsRowDecidesTheRoute:

    def test_a_smiles_ligand_keeps_its_smiles_and_is_not_looked_up_as_the_code_it_carries(self, inference, tmp_path):
        (chain,) = _ligand_chains(inference, _config(tmp_path, [SMILES_ROW]))
        assert chain["ccd"] == [] and chain["smiles"] == "CC(=O)Oc1ccccc1C(=O)O"

    def test_a_ccd_ligand_is_folded_by_its_ccd_code(self, inference, tmp_path):
        (chain,) = _ligand_chains(inference, _config(tmp_path, [CCD_ROW]))
        assert chain["ccd"] == ["ATP"] and chain["smiles"] == ""

    def test_the_stub_builder_receives_the_smiles_for_one_and_the_code_for_the_other(self, inference, tmp_path):
        for row, expected in ((SMILES_ROW, dict(smiles="CC(=O)Oc1ccccc1C(=O)O", ccd=None)),
                              (CCD_ROW, dict(smiles=None, ccd=["ATP"]))):
            (_, chains), = inference.build_complexes(_config(tmp_path, [row]), None)
            built = inference.to_structure_input(chains, _StubBuilder, None, 1024, {})
            (ligand,) = [s for s in built if isinstance(s, _LigandInput)]
            assert dict(smiles=ligand.smiles, ccd=ligand.ccd) == expected

    def test_a_compound_with_neither_names_what_it_lacks(self, inference, tmp_path):
        (_, chains), = inference.build_complexes(_config(tmp_path, ["bare,smiles,LIG,,,,,,,,,"]), None)
        with pytest.raises(ValueError, match="neither ccd nor smiles"):
            inference.to_structure_input(chains, _StubBuilder, None, 1024, {})


class TestBlankCellsStayBlank:

    def test_a_blank_ccd_cell_reads_back_as_an_empty_string_not_nan(self, inference, tmp_path):
        config = _config(tmp_path, [SMILES_ROW])
        (record,), _, _ = inference.load_axis_records(config["axes"]["ligands"])
        assert record["ccd"] == "" and not (isinstance(record["ccd"], float) and math.isnan(record["ccd"]))

    def test_an_id_that_looks_like_a_number_stays_a_string(self, inference, tmp_path):
        config = _config(tmp_path, ["001,smiles,LIG,,,,,,CCO,,,"])
        (record,), _, _ = inference.load_axis_records(config["axes"]["ligands"])
        assert record["id"] == "001"

    def test_a_ccd_row_with_no_smiles_does_not_hand_esm_the_string_nan(self, inference, tmp_path):
        (chain,) = _ligand_chains(inference, _config(tmp_path, [CCD_ROW]))
        assert chain["smiles"] == ""
