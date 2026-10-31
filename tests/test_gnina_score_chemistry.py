"""Gnina score and minimize mode: a compound without chemistry fails its own complexes, not the step.

`_load_score_chemistry` ran before the loop over the complexes and raised on the first compound with no SMILES, so
one such compound in a screen left the step with no `scores.csv` and no `local_missing.csv`, and the record of every
other complex was lost. The compound now fails only the complexes paired with it (a NaN row of `scores` and a
`failure` row of `missing`), and every other complex is scored. Every SMILES the pipe script reads goes through
`ligand_utils.ligand_chemistry`, so a cell holding `NaN`, `None`, `<NA>` or only spaces is no chemistry.
"""

import math
import pathlib
import subprocess

import pandas as pd
import pytest

from biopipelines.datastream import DataStream

from tests.test_gnina_ids import PROTEIN_LINE, pipe_gnina  # noqa: F401 (a fixture)
from tests.test_gnina_conformer_selection import _write_sdf
from tests.test_gnina_missing_table import LIGAND_PDB

SCORE_TERMS = ["vina_affinity", "cnn_score", "cnn_affinity", "cnn_vs", "cnn_affinity_variance"]


def _streams(tmp_path, compounds, complexes):
    """(structures stream, compounds stream) on disk; `compounds` is {id: (code, smiles)}, `complexes` {id: compound id}."""
    tmp_path.joinpath("compounds.csv").write_text(
        "id,code,smiles\n" + "".join(f"{cid},{code},{smiles}\n" for cid, (code, smiles) in compounds.items()))
    lines = ["id,file,compounds.id"]
    for structure_id, compound_id in complexes.items():
        pdb = tmp_path / f"{structure_id}.pdb"
        pdb.write_text(PROTEIN_LINE + LIGAND_PDB + "END\n")
        lines.append(f"{structure_id},{pdb},{compound_id}")
    tmp_path.joinpath("structures_map.csv").write_text("\n".join(lines) + "\n")
    structures = DataStream(name="structures", ids=list(complexes), files=[str(tmp_path / "<id>.pdb")],
                            map_table=str(tmp_path / "structures_map.csv"), format="pdb")
    compounds_ds = DataStream(name="compounds", ids=list(compounds), files=[],
                              map_table=str(tmp_path / "compounds.csv"), format="csv")
    return structures, compounds_ds


def _config(tmp_path):
    return {"mode": "score", "backend": "gnina", "gnina_binary": "gnina", "cnn_scoring": "rescore", "scoring": "vina",
            "output_folder": str(tmp_path / "execution"), "best_poses_dir": str(tmp_path / "structures"),
            "scores_csv": str(tmp_path / "tables" / "001_Gnina" / "scores.csv"),
            "structures_map": str(tmp_path / "structures" / "structures_map.csv")}


def _fake_the_tools(pipe, monkeypatch, templates):
    """gnina and the ligand extraction are functions; `templates` records the SMILES each ligand was built from."""
    def run(cmd, **kwargs):
        _write_sdf(cmd[cmd.index("-o") + 1], [(-7.0, 0.8, 5.0)])
        return subprocess.CompletedProcess(cmd, 0, "", "")

    def extract(complex_pdb, code, out_pdb):
        pathlib.Path(out_pdb).write_text(LIGAND_PDB)
        return out_pdb

    def write_sdf(lig_pdb, dst, smiles=None):
        templates.append(smiles)
        _write_sdf(dst, [(-7.0, 0.8, 5.0)])

    monkeypatch.setattr(pipe.subprocess, "run", run)
    monkeypatch.setattr(pipe, "_extract_ligand_hetatm", extract)
    monkeypatch.setattr(pipe, "write_ligand_sdf", write_sdf)


@pytest.mark.parametrize("blank", ["", "nan", "NaN", "NAN", "None", "none", "<NA>", "   "])
def test_a_compound_without_smiles_fails_its_own_complexes_and_the_rest_are_scored(
        pipe_gnina, monkeypatch, tmp_path, blank):
    structures, compounds = _streams(
        tmp_path, {"good": ("LIG", "CCO"), "blank": ("LIG", blank)},
        {"s_good1": "good", "s_blank": "blank", "s_good2": "good"})
    templates = []
    _fake_the_tools(pipe_gnina, monkeypatch, templates)

    failed = pipe_gnina.run_score_mode(structures, compounds, _config(tmp_path))

    scores = pd.read_csv(_config(tmp_path)["scores_csv"], dtype={"id": str})
    assert list(scores["id"]) == ["s_good1", "s_blank", "s_good2"], "every complex that entered keeps a row"
    assert scores.set_index("id").loc["s_blank", SCORE_TERMS].isna().all()
    for scored in ("s_good1", "s_good2"):
        assert not math.isnan(scores.set_index("id").loc[scored, "vina_affinity"])
    assert templates == ["CCO", "CCO"], "the blank compound's chemistry never reached the ligand builder"
    assert [(row["id"], row["kind"]) for row in failed] == [("s_blank", "failure")]
    assert "no SMILES" in failed[0]["cause"] and "blank" in failed[0]["cause"]


def test_a_compound_without_a_residue_code_fails_its_own_complexes(pipe_gnina, monkeypatch, tmp_path):
    structures, compounds = _streams(
        tmp_path, {"good": ("LIG", "CCO"), "nocode": ("", "CCN")}, {"s_good": "good", "s_nocode": "nocode"})
    _fake_the_tools(pipe_gnina, monkeypatch, [])

    failed = pipe_gnina.run_score_mode(structures, compounds, _config(tmp_path))

    scores = pd.read_csv(_config(tmp_path)["scores_csv"], dtype={"id": str}).set_index("id")
    assert list(scores.index) == ["s_good", "s_nocode"]
    assert not math.isnan(scores.loc["s_good", "vina_affinity"])
    assert [(row["id"], row["kind"]) for row in failed] == [("s_nocode", "failure")]
    assert "residue code" in failed[0]["cause"]


def test_the_step_writes_local_missing_when_one_compound_has_no_smiles(pipe_gnina, monkeypatch, tmp_path):
    """main() reaches write_local_missing, which the raise before the loop used to skip."""
    import json
    import sys

    structures, compounds = _streams(
        tmp_path, {"good": ("LIG", "CCO"), "blank": ("LIG", "")}, {"s_good": "good", "s_blank": "blank"})
    structures.save_json(str(tmp_path / "structures.json"))
    compounds.save_json(str(tmp_path / "compounds.json"))
    config = {**_config(tmp_path), "structures_json": str(tmp_path / "structures.json"),
              "compounds_json": str(tmp_path / "compounds.json"), "upstream_missing": None,
              "local_missing_csv": str(tmp_path / "execution" / "local_missing.csv")}
    (tmp_path / "gnina_config.json").write_text(json.dumps(config))
    _fake_the_tools(pipe_gnina, monkeypatch, [])
    monkeypatch.setattr(sys, "argv", ["pipe_gnina.py", str(tmp_path / "gnina_config.json")])

    pipe_gnina.main()

    local = pd.read_csv(config["local_missing_csv"], dtype=str, keep_default_na=False)
    assert list(local["id"]) == ["s_blank"] and list(local["kind"]) == ["failure"]
    assert pathlib.Path(config["scores_csv"]).exists()


@pytest.mark.parametrize("blank", ["  ", "NAN", "Nan", "none", ""])
def test_the_docking_ligand_loader_reads_a_blank_smiles_as_no_chemistry(pipe_gnina, monkeypatch, tmp_path, blank):
    """`_load_compound_mols` built a ligand from any blank that pandas did not already read as NaN: only the exact string `nan` counted."""
    (tmp_path / "compounds.csv").write_text(f'id,code,smiles\nlig,LIG,"{blank}"\n')
    compounds = DataStream(name="compounds", ids=["lig"], files=[], map_table=str(tmp_path / "compounds.csv"),
                           format="csv")
    seen = []
    original = pipe_gnina.Chem.MolFromSmiles
    monkeypatch.setattr(pipe_gnina.Chem, "MolFromSmiles", lambda text, *a, **k: (seen.append(text), original(text, *a, **k))[1])
    loaded = pipe_gnina._load_compound_mols(compounds, str(tmp_path), 42)
    assert loaded == []
    assert seen == [], f"a blank cell reached RDKit as chemistry: {seen}"
