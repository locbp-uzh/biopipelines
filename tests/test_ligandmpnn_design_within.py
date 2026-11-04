"""`design_within` is the scoring cutoff of the ligand context, not a way to choose residues, and the tool says so.

`docs/tool/sequence_design.md` has always said `design_within` is for analysis and does not control the design,
and the code passes it to `--ligand_mpnn_cutoff_for_score`. The class docstring and the position helper's warning
said the opposite ("Distance from ligand to redesign", "Using design_within=5.0 around BTN"), so a user who gave a
ligand and no positions read it as a 5 A shell and got every residue redesigned. Found on a Colab run of 1STP with biotin:
121 residues redesigned after a warning that named 5 A.
"""

import json
import os
import subprocess
import sys

from biopipelines.ligand_mpnn import LigandMPNN
from tests.test_ligandmpnn_failures import GOOD_PDB, POSITIONS


def _ligand_source(tmp_path):
    pdb = tmp_path / "s1.pdb"
    pdb.write_bytes(GOOD_PDB if isinstance(GOOD_PDB, bytes) else GOOD_PDB.encode())
    ds = tmp_path / "ds.json"
    json.dump({"name": "structures", "ids": ["s1"], "files": [str(pdb)], "map_table": "", "format": "pdb"}, open(ds, "w"))
    compounds = tmp_path / "compounds.csv"
    compounds.write_text("id,code\nbiotin,BTN\n", encoding="utf-8")
    ligand = tmp_path / "ligand.json"
    json.dump({"name": "compounds", "ids": ["biotin"], "files": [], "map_table": str(compounds), "format": "csv"}, open(ligand, "w"))
    cfg, out = tmp_path / "cfg.json", tmp_path / "out.json"
    json.dump({"structures_json": str(ds), "input_source": "ligand", "input_table": "-", "fixed_positions": "-",
               "designed_positions": "-", "ligand": str(ligand), "design_within": "5.0", "output_file": str(out),
               "default_chain": "auto"}, open(cfg, "w"))
    result = subprocess.run([sys.executable, POSITIONS, str(cfg)], capture_output=True, text=True)
    return result, json.load(open(out))


def test_the_warning_says_every_position_is_redesigned_and_what_design_within_is(tmp_path):
    result, written = _ligand_source(tmp_path)
    assert result.returncode == 0, result.stderr
    assert "every position" in result.stdout and "does not select residues" in result.stdout, result.stdout
    assert "Using design_within" not in result.stdout


def test_a_ligand_without_positions_selects_no_residue(tmp_path):
    _, written = _ligand_source(tmp_path)
    assert written["s1"] == {"fixed_option": "", "redesigned_option": ""}


def test_the_docstring_does_not_promise_a_shell():
    doc = LigandMPNN.__init__.__doc__
    line = next(l for l in doc.splitlines() if "design_within:" in l)
    assert "redesign" not in line.lower().replace("redesigned", ""), line
    assert "scor" in line.lower(), line
