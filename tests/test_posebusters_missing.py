"""PoseBusters declares a `missing` table and records the structures it could not validate.

Both inputs (`structures` and the `ligand` stream) can be filtered upstream, yet the tool declared no `missing` table
and wrote no local failures: a structure with no ligand to extract, or one on which PoseBusters raised, got an
`all_pass=False` row and no trace in `missing`. The wrapper now declares the table and merges its local failures with
the upstream manifests, as the other tools that consume streams do.

Nothing here needs `posebusters` or `gemmi`: the validator and the extraction are stand-ins.
"""

import importlib.util
import json
import os
import sys
import types

import pytest

pd = pytest.importorskip("pandas")

from biopipelines.datastream import DataStream
from biopipelines.mock import Mock

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _build(new_pipeline, ligand="HEM", **kwargs):
    from biopipelines.ligand import Ligand
    from biopipelines.posebusters import PoseBusters

    pipeline = new_pipeline("pb_missing")
    with pipeline:
        structures = Mock(ids=["s1"], streams={"structures": {"format": "pdb", "file": "<id>.pdb"}},
                          map_table_strategy="config")
        out = PoseBusters(structures=structures.streams.structures, ligand=Ligand(codes=ligand), **kwargs)
        script_path = pipeline.save()
    return out, script_path


def test_the_missing_table_is_declared(local_config, isolated_cwd, new_pipeline):
    out, _ = _build(new_pipeline)

    assert "missing" in out.tables.keys()
    assert out.tables.missing.info.columns == ["id", "removed_by", "kind", "cause"]


def test_the_script_merges_the_local_failures_into_the_missing_table(local_config, isolated_cwd, new_pipeline):
    import pathlib

    _, script_path = _build(new_pipeline)
    script = next(pathlib.Path(script_path).parent.glob("*_PoseBusters.sh")).read_text(encoding="utf-8", errors="replace")
    config = json.loads(next(pathlib.Path(script_path).parent.parent.rglob("posebusters_config.json")).read_text())

    assert "--local-missing" in script and "pipe_propagate_missing.py" in script
    assert config["local_missing_csv"].endswith("local_missing.csv")


def _load_runner(monkeypatch, outcome):
    """pipe_posebusters with a stand-in `posebusters`; `outcome[structure id]` is True/False or an exception."""
    module_spec = importlib.util.spec_from_file_location(
        "pipe_posebusters_missing", os.path.join(REPO_ROOT, "pipe_scripts", "pipe_posebusters.py"))
    module = importlib.util.module_from_spec(module_spec)
    module_spec.loader.exec_module(module)

    state = {}

    class PoseBusters:
        def __init__(self, config):
            pass

        def bust(self, mol_pred, mol_cond, **kwargs):
            result = outcome[state["current"]]
            if isinstance(result, Exception):
                raise result
            return pd.DataFrame([{"mol_pred_loaded": result, "bond_lengths": result}])

    fake = types.ModuleType("posebusters")
    fake.PoseBusters = PoseBusters
    monkeypatch.setitem(sys.modules, "posebusters", fake)

    def extract(structure_path, ligand_name, output_dir, structure_id, ligand_smiles=None):
        state["current"] = structure_id
        if structure_id == "nolig":
            return []
        return [(os.path.join(output_dir, "lig.sdf"), os.path.join(output_dir, "prot.pdb"), "")]

    monkeypatch.setattr(module, "extract_ligand_and_protein", extract)
    return module


def _run(monkeypatch, tmp_path, ids, outcome):
    for sid in ids:
        (tmp_path / f"{sid}.pdb").write_text("ATOM      1  N   ALA A   1       0.0   0.0   0.0\n")
    (tmp_path / "compounds.csv").write_text("id,code,smiles\nlig,LIG,\n")
    DataStream(name="structures", ids=list(ids), files=[str(tmp_path / "<id>.pdb")], format="pdb",
               map_table=str(tmp_path / "structures_map.csv")).save_json(str(tmp_path / "structures.json"))
    DataStream(name="compounds", ids=["lig"], files=[], format="csv",
               map_table=str(tmp_path / "compounds.csv")).save_json(str(tmp_path / "ligand.json"))
    output = tmp_path / "Job_001" / "007_PoseBusters"
    config = {"structures_json": str(tmp_path / "structures.json"), "ligand_json": str(tmp_path / "ligand.json"),
              "mode": "dock", "reference_pdb": None, "output_csv": str(output / "tables" / "posebusters.csv"),
              "execution_dir": str(output / "execution"), "check_columns": ["mol_pred_loaded", "bond_lengths"],
              "local_missing_csv": str(output / "execution" / "local_missing.csv")}
    os.makedirs(output / "tables")
    os.makedirs(output / "execution")
    module = _load_runner(monkeypatch, outcome)
    module.run_posebusters(config)
    return (pd.read_csv(config["output_csv"]),
            pd.read_csv(config["local_missing_csv"], dtype=str, keep_default_na=False))


def test_a_structure_without_the_ligand_and_a_failed_validation_are_failure_rows(monkeypatch, tmp_path):
    table, missing = _run(monkeypatch, tmp_path, ["good", "nolig", "boom"],
                          {"good": True, "boom": RuntimeError("PoseBusters could not read the molecule")})

    assert list(table["id"]) == ["good", "nolig", "boom"], "every structure keeps its row"
    assert table["all_pass"].tolist() == [True, False, False]
    assert list(missing["id"]) == ["nolig", "boom"]
    assert set(missing["kind"]) == {"failure"} and set(missing["removed_by"]) == {"007_PoseBusters"}
    assert "no ligand" in missing.loc[0, "cause"] and "could not read the molecule" in missing.loc[1, "cause"]


def test_a_step_where_every_structure_validated_has_an_empty_local_missing_file(monkeypatch, tmp_path):
    _, missing = _run(monkeypatch, tmp_path, ["a", "b"], {"a": True, "b": False})

    assert list(missing.columns) == ["id", "removed_by", "kind", "cause"]
    assert missing.empty, "a pose that fails a check is a result, not a failure of the step"
