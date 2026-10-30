# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""RFDAA_PrepareLigand: a prepared ligand keeps the id of the ligand it was made from.

The step turns each ligand structure into a ligand-plus-dummy-peptide PDB, one output per input, yet named its single output `prepared_ligand` whatever went in. A stream that is 1:1 with its input keeps the input's id, so the designs RFdiffusionAllAtom derives from it (`<pdb id>_<n>`) carry the ligand's name and not a constant, and it handled only the first of several ligands without saying so.
"""

import csv
import os
import subprocess
import sys

import pytest

from biopipelines.datastream import DataStream, create_map_table

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PIPE = os.path.join(REPO_ROOT, "pipe_scripts", "pipe_rfdaa_prepare_ligand.py")

LIGANDS = ["ZIT", "ATP"]

LIGAND_PDB = "HETATM    1  C1  {code} A   1       0.000   0.000   0.000  1.00  0.00           C\nEND\n"
PEPTIDE_PDB = "".join(
    f"ATOM  {i:5d}  CA  GLY A{i:4d}    {float(i):8.3f}   0.000   0.000  1.00  0.00           C\n" for i in (1, 2, 3))


def _ligand_stream(tmp_path, ids):
    folder = tmp_path / "ligands"
    folder.mkdir(exist_ok=True)
    for lid in ids:
        (folder / f"{lid}.pdb").write_text(LIGAND_PDB.format(code=lid))
    template = str(folder / "<id>.pdb")
    table = str(tmp_path / "ligand_map.csv")
    create_map_table(table, ids=list(ids), files=[template])
    return DataStream(name="structures", ids=list(ids), files=[template], map_table=table, format="pdb")


def _prepare(new_pipeline, tmp_path, ids):
    from biopipelines.rfdiffusion_allatom import RFDAA_PrepareLigand

    pipeline = new_pipeline("prepare_ligand_ids")
    with pipeline:
        prepared = RFDAA_PrepareLigand(ligand=_ligand_stream(tmp_path, ids))
        script_path = pipeline.save()
    return prepared, script_path


def _run(config, tmp_path):
    pdbs = tmp_path / "pdbs"
    pdbs.mkdir(exist_ok=True)
    (pdbs / "1PEF.pdb").write_text(PEPTIDE_PDB)
    done = subprocess.run(
        [sys.executable, PIPE,
         "--ligand_json", config.ligand_ds_json,
         "--structures_folder", config.stream_folder("structures"),
         "--output_csv", config.structures_csv,
         "--structures_map", config.structures_map,
         "--pdbs_folder", str(pdbs)],
        capture_output=True, text=True)
    assert done.returncode == 0, done.stdout + done.stderr


def _rows(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def test_the_declared_ids_are_the_ligand_ids(local_config, isolated_cwd, new_pipeline, tmp_path):
    prepared, _ = _prepare(new_pipeline, tmp_path, LIGANDS)

    assert list(prepared.streams.structures.ids) == LIGANDS


def test_every_ligand_is_prepared_under_its_own_id(local_config, isolated_cwd, new_pipeline, tmp_path):
    prepared, _ = _prepare(new_pipeline, tmp_path, LIGANDS)
    config = prepared.producer

    _run(config, tmp_path)

    assert sorted(os.listdir(config.stream_folder("structures"))) == [f"{lid}.pdb" for lid in sorted(LIGANDS)] + ["structures_map.csv"]
    mapped = _rows(config.structures_map)
    assert [r["id"] for r in mapped] == LIGANDS
    for row in mapped:
        assert os.path.basename(row["file"]) == f"{row['id']}.pdb" and os.path.exists(row["file"])
    assert [r["id"] for r in _rows(config.structures_csv)] == LIGANDS


def test_each_prepared_file_holds_its_own_ligand(local_config, isolated_cwd, new_pipeline, tmp_path):
    prepared, _ = _prepare(new_pipeline, tmp_path, LIGANDS)
    config = prepared.producer

    _run(config, tmp_path)

    for lid in LIGANDS:
        text = open(os.path.join(config.stream_folder("structures"), f"{lid}.pdb")).read()
        assert f"HETATM    1  C1  {lid} A" in text and "ATOM" in text


def test_the_script_hands_the_pipe_the_stream_not_one_resolved_file(local_config, isolated_cwd, new_pipeline, tmp_path):
    _, script_path = _prepare(new_pipeline, tmp_path, LIGANDS)
    runtime = os.path.dirname(script_path)
    step = next(f for f in sorted(os.listdir(runtime)) if f.endswith("_RFDAA_PrepareLigand.sh"))
    with open(os.path.join(runtime, step), encoding="utf-8", errors="replace") as handle:
        text = handle.read()

    assert "--ligand_json" in text and "--structures_folder" in text
    assert "--ligand_pdb" not in text and "LIGAND_ID" not in text
