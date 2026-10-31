"""PPIDdG: wrapper validation, generated script, result table shape and the shared scoring helpers.

Everything here runs without PyRosetta. The tests that need it are at the bottom and skip when it is not installed.
"""

import glob
import importlib.util
import json
import os
import re
import sys

import pytest

pytest.importorskip("pandas")
pytest.importorskip("yaml")

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HELPER = os.path.join(REPO_ROOT, "pipe_scripts", "pipe_ppi_ddg.py")


@pytest.fixture(scope="module")
def pipe():
    spec = importlib.util.spec_from_file_location("pipe_ppi_ddg", HELPER)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def step_script(pipeline, tool):
    """The step's bash script, as written by pipeline.save()."""
    found = glob.glob(os.path.join(pipeline.folders["output"], "**", f"*_{tool}.sh"), recursive=True)
    assert found, f"no {tool} step script under {pipeline.folders['output']}"
    # the framework writes step scripts in the platform encoding, so an em dash in a comment is not valid UTF-8 on Windows
    return open(found[0], encoding="utf-8", errors="replace").read()


def build(new_pipeline, name, fmt="pdb", **kwargs):
    from biopipelines.mock import Mock
    from biopipelines.ppi_ddg import PPIDdG

    pipeline = new_pipeline(name)
    with pipeline:
        source = Mock(ids=["s1", "s2"], streams={"structures": {"format": fmt, "file": "<id>.pdb"}})
        PPIDdG(structures=source.streams.structures, **kwargs)
        pipeline.save()
    return pipeline, pipeline.tools[-1]


# -- wrapper validation ---------------------------------------------------------------------------

@pytest.mark.parametrize("kwargs,match", [
    ({"binder_chain": "AB"}, "single-character"),
    ({"target_chain": ""}, "single-character"),
    ({"binder_chain": "C", "target_chain": "C"}, "different chains"),
    ({"relax_rounds": 0}, "positive integer"),
    ({"relax_rounds": True}, "positive integer"),
    ({"replicates": 0}, "positive integer"),
    ({"replicates": 1.5}, "positive integer"),
    ({"seed": 0}, "positive integer"),
    ({"seed": True}, "positive integer"),
    ({"constrain_to_start_coords": "yes"}, "must be a bool"),
    ({"constraint_weight": -1.0}, "number >= 0"),
])
def test_bad_parameters_are_refused_before_a_job_is_queued(local_config, isolated_cwd, new_pipeline, kwargs, match):
    with pytest.raises(ValueError, match=match):
        build(new_pipeline, "ppi_bad", **kwargs)


def test_cif_structures_are_refused(local_config, isolated_cwd, new_pipeline):
    with pytest.raises(ValueError, match="reads PDB files"):
        build(new_pipeline, "ppi_cif", fmt="cif")


def test_defaults_follow_the_chaperone_block(local_config, isolated_cwd, new_pipeline):
    _, tool = build(new_pipeline, "ppi_defaults")
    assert (tool.binder_chain, tool.target_chain, tool.relax_rounds, tool.replicates,
            tool.constrain_to_start_coords, tool.constraint_weight) == ("A", "B", 1, 1, True, 1.0)


# -- generated script and declared outputs ----------------------------------------------------------

def test_script_flags_match_the_pipe_script(local_config, isolated_cwd, new_pipeline):
    pipeline, _ = build(new_pipeline, "ppi_flags", replicates=3)
    script = step_script(pipeline, "PPIDdG")
    sent = set(re.findall(r"--[a-z][a-z-]*", script.split("pipe_ppi_ddg.py", 1)[1].split("Propagating", 1)[0]))
    accepted = set(re.findall(r'add_argument\("(--[a-z-]+)"', open(HELPER, encoding="utf-8").read()))
    assert sent == accepted, (sent, accepted)
    assert "Propagating upstream missing manifest(s)" in script


def test_parameters_reach_the_pipe_through_the_config_yaml(local_config, isolated_cwd, new_pipeline):
    import yaml

    pipeline, tool = build(new_pipeline, "ppi_yaml", binder_chain="H", target_chain="L", relax_rounds=3,
                           replicates=2, constrain_to_start_coords=False, constraint_weight=0.5)
    config = yaml.safe_load(open(tool.config_yaml))
    assert config == {"binder_chain": "H", "target_chain": "L", "relax_rounds": 3, "replicates": 2,
                      "constrain_to_start_coords": False, "constraint_weight": 0.5, "seed": 42}


def test_the_result_table_and_the_missing_table_are_declared(local_config, isolated_cwd, new_pipeline):
    _, tool = build(new_pipeline, "ppi_outputs")
    tables = tool.get_output_files()["tables"]
    assert list(tables["ppi_ddg"].info.columns) == ["id", "ddg", "ddg_std", "n_replicates", "relax_score"]
    assert list(tables["missing"].info.columns) == ["id", "removed_by", "kind", "cause"]


# -- install ------------------------------------------------------------------------------------------

def test_install_builds_the_env_and_fetches_molfile_to_params(local_config, monkeypatch):
    from biopipelines.ppi_ddg import PPIDdG

    monkeypatch.setattr(PPIDdG, "_install_env", classmethod(lambda cls, manager, name=None: name or "PyRosetta"))
    script = PPIDdG._install_script({"biopipelines": "/repo"}, env_manager="mamba")
    assert "mamba env create -f" in script and "PyRosetta" in script
    assert "git clone --depth 1 --filter=blob:none --sparse https://github.com/RosettaCommons/rosetta.git" in script
    assert "sparse-checkout set source/scripts/python/public" in script
    assert "share/rosetta_public/source/scripts/python/public/molfile_to_params.py" in script
    assert "license" in script.lower()
    assert 'touch "$INSTALL_SUCCESS"' in script


def test_the_env_files_carry_the_pyrosetta_index_and_the_chemistry_stack():
    pip = open(os.path.join(REPO_ROOT, "environments", "PyRosetta.pip.txt"), encoding="utf-8").read()
    assert "--find-links https://west.rosettacommons.org/pyrosetta/quarterly/release.cxx11thread.serialization" in pip
    requirements = [line.strip() for line in pip.splitlines() if line.strip() and not line.startswith(("#", "-"))]
    assert {"pyrosetta", "rdkit", "dimorphite-dl", "pandas", "pyyaml"} <= set(requirements)
    # unpinned on purpose: the run's env export is the reproducibility record
    assert not [r for r in requirements if "=" in r]


# -- shared scoring helpers ---------------------------------------------------------------------------

def test_mean_std_is_the_population_deviation():
    from biopipelines.rosetta_scoring import mean_std

    assert mean_std([1.0, 3.0]) == (2.0, 1.0)
    assert mean_std([5.0]) == (5.0, 0.0)


@pytest.mark.parametrize("values,expected", [
    ([-3.0], (-3.0, 0.0, 1)),
    ([-3.0, -1.0], (-2.0, 1.0, 2)),
    ([-9.0, -3.0, -1.0, 40.0], (-2.0, 1.0, 2)),
    ([1.0, 2.0, 3.0], (2.0, 0.0, 1)),
])
def test_trimmed_mean_drops_best_and_worst_from_three_replicates_on(values, expected):
    from biopipelines.rosetta_scoring import trimmed_mean_std

    assert trimmed_mean_std(values) == expected


def test_each_replicate_of_each_structure_gets_its_own_fixed_seed():
    """PyRosetta seeds itself from the system, so a rerun gave other numbers; the seeds are derived, not drawn."""
    from biopipelines.rosetta_scoring import MAX_SEED, replicate_seed

    seeds = {(sid, r): replicate_seed(42, sid, r) for sid in ("wt", "mutant") for r in range(5)}
    assert seeds == {(sid, r): replicate_seed(42, sid, r) for sid in ("wt", "mutant") for r in range(5)}
    assert len(set(seeds.values())) == 10
    assert replicate_seed(43, "wt", 0) != seeds[("wt", 0)]
    assert all(1 <= value <= MAX_SEED for value in seeds.values())


def test_the_pipe_seeds_rosetta_before_every_replicate(pipe, monkeypatch):
    import types

    from biopipelines.rosetta_scoring import replicate_seed

    seeded, relaxed = [], []
    monkeypatch.setitem(sys.modules, "pyrosetta", types.SimpleNamespace(
        pose_from_file=lambda path: "pose"))
    monkeypatch.setattr(pipe, "chains_of", lambda pose: {"A", "B"})
    monkeypatch.setattr(pipe, "seed_rosetta", lambda seed: seeded.append(seed))
    monkeypatch.setattr(pipe, "relax_and_ddg", lambda *a: relaxed.append(len(seeded)) or (-10.0, -300.0))
    config = {"binder_chain": "A", "target_chain": "B", "relax_rounds": 1, "replicates": 3,
              "constrain_to_start_coords": True, "constraint_weight": 1.0, "seed": 7}
    pipe.score_structure("s1.pdb", config, "s1")
    assert seeded == [replicate_seed(7, "s1", r) for r in range(3)]
    assert relaxed == [1, 2, 3]


def test_a_failed_id_keeps_its_row_with_nan_and_is_recorded_as_a_failure(pipe):
    from biopipelines.rosetta_scoring import collect_rows

    def score(sid, path):
        if sid == "bad":
            raise RuntimeError("pose failed to load")
        return {"ddg": -12.5, "ddg_std": 0.0, "n_replicates": 1, "relax_score": -300.0}

    rows, missing, attempted, failures = collect_rows(
        [("good", "g.pdb"), ("bad", "b.pdb")], score, pipe.COLUMNS, "007_PPIDdG")
    assert [r["id"] for r in rows] == ["good", "bad"]
    assert rows[0]["ddg"] == -12.5
    assert all(rows[1][c] != rows[1][c] for c in pipe.COLUMNS if c != "id")
    assert missing == [{"id": "bad", "removed_by": "007_PPIDdG", "kind": "failure", "cause": "pose failed to load"}]
    assert (attempted, failures) == (2, 1)


def test_the_step_fails_only_when_nothing_ran_or_everything_raised():
    from biopipelines.rosetta_scoring import step_exit_status

    assert step_exit_status(0, 0) is not None
    assert step_exit_status(3, 3) is not None
    assert step_exit_status(3, 2) is None
    assert step_exit_status(3, 0) is None


def _write_structures(tmp_path, ids):
    import pandas as pd

    files = [str(tmp_path / f"{i}.pdb") for i in ids]
    for f in files:
        open(f, "w").write("END\n")
    map_csv = tmp_path / "structures_map.csv"
    pd.DataFrame({"id": ids, "file": files}).to_csv(map_csv, index=False)
    path = tmp_path / "structures.json"
    json.dump({"name": "structures", "ids": ids, "files": files, "map_table": str(map_csv),
               "format": "pdb"}, open(path, "w"))
    return str(path)


def test_main_writes_one_row_per_id_a_local_missing_file_and_scores_with_the_yaml_parameters(pipe, tmp_path, monkeypatch):
    import pandas as pd
    import yaml

    structures = _write_structures(tmp_path, ["a", "b", "c"])
    config = tmp_path / "config.yaml"
    yaml.safe_dump({"binder_chain": "A", "target_chain": "B", "relax_rounds": 2, "replicates": 3,
                    "constrain_to_start_coords": True, "constraint_weight": 1.0}, open(config, "w"))
    seen = []

    def fake_score(path, cfg, sid):
        seen.append((os.path.basename(path), cfg["replicates"], cfg["relax_rounds"]))
        if path.endswith("b.pdb"):
            raise ValueError("chain 'B' is not in the structure")
        return {"ddg": -10.0, "ddg_std": 0.5, "n_replicates": 3, "relax_score": -250.0}

    monkeypatch.setattr(pipe, "score_structure", fake_score)
    monkeypatch.setattr(pipe, "init_pyrosetta", lambda options: None)
    table = tmp_path / "tables" / "ppi_ddg.csv"
    local = tmp_path / "_execution" / "local_missing.csv"
    pipe.main(["--config-yaml", str(config), "--structures-json", structures,
               "--table-csv", str(table), "--local-missing-csv", str(local)])

    assert seen == [("a.pdb", 3, 2), ("b.pdb", 3, 2), ("c.pdb", 3, 2)]
    frame = pd.read_csv(table)
    assert list(frame.columns) == pipe.COLUMNS
    assert list(frame["id"]) == ["a", "b", "c"]
    assert frame.loc[frame["id"] == "b", "ddg"].isna().all()
    assert frame.loc[frame["id"] == "a", "ddg"].item() == -10.0
    failed = pd.read_csv(local)
    assert list(failed["id"]) == ["b"] and list(failed["kind"]) == ["failure"]
    assert "not in the structure" in failed["cause"].item()


def test_main_exits_non_zero_when_every_id_raised(pipe, tmp_path, monkeypatch):
    import yaml

    structures = _write_structures(tmp_path, ["a", "b"])
    config = tmp_path / "config.yaml"
    yaml.safe_dump({"binder_chain": "A", "target_chain": "B"}, open(config, "w"))
    monkeypatch.setattr(pipe, "score_structure", lambda path, cfg, sid: (_ for _ in ()).throw(RuntimeError("boom")))
    monkeypatch.setattr(pipe, "init_pyrosetta", lambda options: None)
    with pytest.raises(SystemExit) as stop:
        pipe.main(["--config-yaml", str(config), "--structures-json", structures,
                   "--table-csv", str(tmp_path / "t.csv"), "--local-missing-csv", str(tmp_path / "m.csv")])
    assert stop.value.code == 1
    assert os.path.isfile(tmp_path / "t.csv")


# -- PyRosetta-dependent ------------------------------------------------------------------------------

def test_the_pyrosetta_api_the_pipes_call_exists():
    pyr = pytest.importorskip("pyrosetta", reason="pyrosetta is not installed in this environment (it needs its own licensed one): the Rosetta API the pipes call is not checked")
    ros = pyr.rosetta
    assert callable(pyr.create_score_function) and callable(pyr.pose_from_file)
    assert callable(pyr.generate_nonstandard_residue_set)
    assert hasattr(ros.core.scoring, "coordinate_constraint")
    assert callable(ros.protocols.relax.FastRelax)
    assert hasattr(ros.protocols.relax.FastRelax, "constrain_relax_to_start_coords")
    assert callable(ros.protocols.toolbox.pose_manipulation.rigid_body_move)
    for selector in ("ChainSelector", "NeighborhoodResidueSelector", "NotResidueSelector", "OrResidueSelector"):
        assert hasattr(ros.core.select.residue_selector, selector)
    for operation in ("RestrictToRepacking", "OperateOnResidueSubset", "PreventRepackingRLT", "IncludeCurrent"):
        assert hasattr(ros.core.pack.task.operation, operation)
    assert hasattr(ros.protocols.minimization_packing, "PackRotamersMover")
    assert hasattr(ros.protocols.rosetta_scripts, "XmlObjects")
