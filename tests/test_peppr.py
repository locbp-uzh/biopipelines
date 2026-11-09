"""Peppr scores predictions against references with the peppr package.

The wrapper, the id pairing and the failure handling are tested without the library. The metric values are tested on a synthetic complex (two helices and a ligand) and skipped when peppr is not installed.
"""

import glob
import importlib.util
import json
import math
import os
import re
import sys

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("yaml")

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT = os.path.join(REPO_ROOT, "pipe_scripts", "pipe_peppr.py")
ALL_METRICS = ["rmsd", "tm_score", "lddt", "lddt_pli", "fnat", "irmsd", "lrmsd", "dockq"]
PAIR_METRICS = ["fnat", "irmsd", "lrmsd", "dockq"]


@pytest.fixture(scope="module")
def pipe():
    spec = importlib.util.spec_from_file_location("pipe_peppr", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# -- the wrapper --------------------------------------------------------------

def _structures_mock(ids):
    from biopipelines.mock import Mock
    return Mock(ids=ids, streams={"structures": {"format": "pdb", "file": "<id>.pdb"}},
                map_table_strategy="config")


@pytest.fixture
def make_peppr(local_config, isolated_cwd, new_pipeline):
    """Build a Peppr inside a live pipeline, where validate_params runs at registration."""
    from biopipelines.peppr import Peppr

    with new_pipeline("peppr"):
        predictions, references = _structures_mock(["a", "b"]), _structures_mock(["r"])

        def _make(**kwargs):
            return Peppr(structures=predictions.streams.structures,
                         references=references, **kwargs)
        yield _make


def _columns(table):
    return table._info.columns


def test_the_table_is_id_plus_every_metric_by_default(make_peppr):
    out = make_peppr()
    assert _columns(out.tables.peppr) == ["id"] + ALL_METRICS


def test_the_metrics_parameter_picks_and_orders_the_columns(make_peppr):
    out = make_peppr(metrics=["dockq", "rmsd"])
    assert _columns(out.tables.peppr) == ["id", "dockq", "rmsd"]


def test_a_missing_table_is_declared(make_peppr):
    out = make_peppr()
    assert _columns(out.tables.missing) == ["id", "removed_by", "kind", "cause"]


@pytest.mark.parametrize("metrics,message", [
    (["rmsd", "gdt"], "unknown metric"),
    (["rmsd", "rmsd"], "more than once"),
    ([], "non-empty list"),
    ("dockq", "non-empty list"),
])
def test_bad_metrics_are_refused_at_configuration(make_peppr, metrics, message):
    with pytest.raises(ValueError, match=message):
        make_peppr(metrics=metrics)


def test_inputs_must_be_streams(local_config, isolated_cwd, new_pipeline):
    from biopipelines.peppr import Peppr

    with new_pipeline("peppr_inputs"):
        good = _structures_mock(["a"]).streams.structures
        with pytest.raises(ValueError, match="structures must be"):
            Peppr(structures="a.pdb", references=good)
        with pytest.raises(ValueError, match="references must be"):
            Peppr(structures=good, references="r.pdb")


def test_the_wrapper_and_the_pipe_agree_on_the_metric_keys(pipe):
    from biopipelines.peppr import PEPPR_METRICS
    assert list(PEPPR_METRICS) == ALL_METRICS
    assert set(pipe.METRIC_FACTORIES) == set(PEPPR_METRICS)


def test_the_generated_script_runs_the_helper_and_propagates_missing(
        local_config, isolated_cwd, new_pipeline):
    from biopipelines.peppr import Peppr

    with new_pipeline("peppr_script"):
        predictions, references = _structures_mock(["a", "b"]), _structures_mock(["r"])
        Peppr(structures=predictions, references=references, metrics=["rmsd", "dockq"])

    [step] = glob.glob(str(isolated_cwd / "outputs" / "**" / "*_Peppr.sh"), recursive=True)
    script = open(step, encoding="utf-8", errors="replace").read()
    assert "pipe_peppr.py" in script
    assert "pipe_propagate_missing.py" in script and "local_missing.csv" in script

    # Every flag the wrapper sends is one the helper declares.
    invocation = script.split("pipe_peppr.py")[1].split("\n\n")[0]
    sent = set(re.findall(r"^\s+(--[a-z-]+) ", invocation, re.M))
    declared = set(re.findall(r'add_argument\("(--[a-z-]+)"', open(SCRIPT, encoding="utf-8").read()))
    assert sent and sent <= declared

    [config_dir] = glob.glob(str(isolated_cwd / "outputs" / "**" / "003_Peppr" / "_configuration"),
                             recursive=True)
    import yaml
    assert yaml.safe_load(open(os.path.join(config_dir, "config.yaml"))) == {"metrics": ["rmsd", "dockq"]}
    assert json.load(open(os.path.join(config_dir, "references.json")))["ids"] == ["r"]
    assert json.load(open(os.path.join(config_dir, "structures.json")))["ids"] == ["a", "b"]


# -- pairing predictions with references --------------------------------------

def test_one_reference_scores_every_prediction(pipe):
    pairing = pipe.pair_references(["a", "b", "c"], {"r": "r.pdb"}, n_declared=1)
    assert pairing == {"a": "r", "b": "r", "c": "r"}


def test_several_references_are_matched_by_id(pipe):
    pairing = pipe.pair_references(["d_1", "d_2"], {"d_1": "1.pdb", "d_2": "2.pdb"}, n_declared=2)
    assert pairing == {"d_1": "d_1", "d_2": "d_2"}


def test_a_refold_finds_the_design_it_derives_from(pipe):
    pairing = pipe.pair_references(["d_1_1", "d_2_1"], {"d_1": "1.pdb", "d_2": "2.pdb"}, n_declared=2)
    assert pairing == {"d_1_1": "d_1", "d_2_1": "d_2"}


def test_a_prediction_without_a_reference_pairs_to_none(pipe):
    pairing = pipe.pair_references(["d_1", "other"], {"d_1": "1.pdb", "d_2": "2.pdb"}, n_declared=2)
    assert pairing == {"d_1": "d_1", "other": None}


def test_a_single_declared_reference_whose_file_is_absent_pairs_to_none(pipe):
    assert pipe.pair_references(["a"], {}, n_declared=1) == {"a": None}


def test_a_sibling_reference_is_not_a_match(pipe):
    """A filtered reference stream still declares every id, so the declared count decides between one-for-all and by-id.

    get_mapped_ids would answer d_1 for d_2 through its sibling tier; scoring d_2 against another molecule's structure gives a number that looks fine and is wrong.
    """
    pairing = pipe.pair_references(["d_1", "d_2"], {"d_1": "1.pdb"}, n_declared=2)
    assert pairing == {"d_1": "d_1", "d_2": None}


def test_a_renamed_pose_finds_its_reference_through_the_provenance_column(pipe, tmp_path):
    """Panda renames design_A to pose_1; no string relationship survives, so only the
    poses' map_table can say which reference the pose derives from."""
    poses_map = tmp_path / "structures_map.csv"
    poses_map.write_text("id,file,structures.id\n"
                         "pose_1,pose_1.pdb,design_A\n"
                         "pose_2,pose_2.pdb,design_B\n")
    references = {"design_A": "a.pdb", "design_B": "b.pdb"}

    blind = pipe.pair_references(["pose_1", "pose_2"], references, n_declared=2)
    assert blind == {"pose_1": None, "pose_2": None}, "without the map there is nothing to match on"

    seeing = pipe.pair_references(["pose_1", "pose_2"], references, n_declared=2,
                                  map_table_paths=[str(poses_map)])
    assert seeing == {"pose_1": "design_A", "pose_2": "design_B"}


def test_the_provenance_tier_does_not_resurrect_a_sibling_match(pipe, tmp_path):
    """Adding the map must not loosen the sibling refusal: a common ancestor is still not a match."""
    poses_map = tmp_path / "structures_map.csv"
    poses_map.write_text("id,file,structures.id\n"
                         "d_1,d_1.pdb,d_1\n"
                         "d_2,d_2.pdb,d_2\n")
    pairing = pipe.pair_references(["d_1", "d_2"], {"d_1": "1.pdb"}, n_declared=2,
                                   map_table_paths=[str(poses_map)])
    assert pairing == {"d_1": "d_1", "d_2": None}


def test_score_all_hands_the_map_down_to_the_pairing(pipe, tmp_path):
    """score_all is the only caller in the pipe, so a map it drops is a map nobody uses."""
    poses_map = tmp_path / "structures_map.csv"
    poses_map.write_text("id,file,structures.id\npose_1,pose_1.pdb,design_A\n")
    seen = {}

    def evaluate(sid, reference_path, pose_path, keys):
        seen[sid] = reference_path
        return {key: 1.0 for key in keys}

    rows, missing = pipe.score_all({"pose_1": "p.pdb"}, {"design_A": "a.pdb", "design_B": "b.pdb"},
                                   2, ["rmsd"], "003_Peppr", evaluate,
                                   map_table_paths=[str(poses_map)])
    assert seen == {"pose_1": "a.pdb"}
    assert missing == []
    assert rows[0]["rmsd"] == 1.0


# -- one row per id, NaN for what could not be evaluated ----------------------

def _fake(outcomes):
    """An evaluate() that returns or raises what `outcomes` says for each system id."""
    def evaluate(system_id, reference_path, pose_path, keys):
        outcome = outcomes[system_id]
        if isinstance(outcome, Exception):
            raise outcome
        return {key: outcome[key] for key in keys}
    return evaluate


def test_every_id_gets_a_row_and_a_failure_is_nan_with_a_missing_entry(pipe):
    poses = {"good": "good.pdb", "bad": "bad.pdb", "orphan": "orphan.pdb"}
    references = {"good": "ref_good.pdb", "bad": "ref_bad.pdb"}
    evaluate = _fake({"good": {"rmsd": 0.5, "dockq": 0.9},
                      "bad": ValueError("Reference and pose have different entities")})
    rows, missing = pipe.score_all(poses, references, 2, ["rmsd", "dockq"], "003_Peppr", evaluate)

    assert [row["id"] for row in rows] == ["good", "bad", "orphan"]
    assert (rows[0]["rmsd"], rows[0]["dockq"]) == (0.5, 0.9)
    for row in rows[1:]:
        assert math.isnan(row["rmsd"]) and math.isnan(row["dockq"])
    assert [(m["id"], m["kind"], m["removed_by"]) for m in missing] == [
        ("bad", "failure", "003_Peppr"), ("orphan", "unmatched", "003_Peppr")]
    assert "different entities" in missing[0]["cause"]
    assert "no reference" in missing[1]["cause"]


def test_a_prediction_with_only_a_sibling_reference_is_unmatched_not_scored(pipe):
    """`des_3_1` shares only the ancestor `des` with the reference `des_1`; it must not be scored against it."""
    seen = []

    def evaluate(system_id, reference_path, pose_path, keys):
        seen.append(reference_path)
        return {key: 1.0 for key in keys}

    rows, missing = pipe.score_all({"des_1_1": "p1.pdb", "des_3_1": "p3.pdb"},
                                   {"des_1": "r1.pdb", "des_2": "r2.pdb"}, 2, ["rmsd"], "003_Peppr", evaluate)
    assert seen == ["r1.pdb"]
    assert [(m["id"], m["kind"]) for m in missing] == [("des_3_1", "unmatched")]
    assert math.isnan(rows[1]["rmsd"])


def test_the_single_reference_removed_upstream_is_not_reported_as_an_id_mismatch(pipe):
    """One declared reference broadcasts to every pose; with no file the cause is the upstream removal, not a mismatch."""
    rows, missing = pipe.score_all({"a": "a.pdb", "b": "b.pdb"}, {}, 1, ["rmsd"], "003_Peppr", _fake({}))
    assert [(m["id"], m["kind"]) for m in missing] == [("a", "unmatched"), ("b", "unmatched")]
    for entry in missing:
        assert "removed upstream" in entry["cause"]
    assert all(math.isnan(row["rmsd"]) for row in rows)


def _fake_peppr(warnings_emitted, results):
    """A stand-in peppr module: `warnings_emitted` is [(category name, text)], `results` the per-metric values."""
    import types

    module = types.ModuleType("peppr")

    class MatchWarning(UserWarning):
        pass

    class EvaluationWarning(UserWarning):
        pass

    class Evaluator:
        def __init__(self, metrics, tolerate_exceptions=False):
            self.metrics = metrics

        def feed(self, system_id, reference, pose):
            import warnings as warn
            for category, text in warnings_emitted:
                warn.warn(text, {"match": MatchWarning, "evaluation": EvaluationWarning}[category])

        def get_results(self):
            return [[[value]] for value in results]

    module.MatchWarning = MatchWarning
    module.EvaluationWarning = EvaluationWarning
    module.Evaluator = Evaluator
    for name in ("MonomerRMSD", "MonomerTMScore", "MonomerLDDTScore", "LDDTPLIScore",
                 "ContactFraction", "InterfaceRMSD", "LigandRMSD", "DockQScore"):
        setattr(module, name, lambda *a, **k: object())
    return module


def test_several_warnings_from_one_metric_do_not_discard_the_metrics_that_computed(pipe, monkeypatch):
    """peppr can warn more than once per metric, so counting warnings against the metric count threw away good values."""
    monkeypatch.setitem(sys.modules, "peppr",
                        _fake_peppr([("evaluation", "interface absent"), ("evaluation", "interface absent"),
                                     ("evaluation", "interface absent")],
                                    [0.5, math.nan, 0.8]))
    monkeypatch.setattr(pipe, "load_system", lambda path: object())
    values = pipe.evaluate_pair("x", "ref.pdb", "pose.pdb", ["rmsd", "dockq", "lddt"])
    assert values["rmsd"] == 0.5 and values["lddt"] == 0.8 and math.isnan(values["dockq"])


def test_every_metric_failing_still_raises(pipe, monkeypatch):
    monkeypatch.setitem(sys.modules, "peppr",
                        _fake_peppr([("evaluation", "nothing to compare")], [math.nan, math.nan]))
    monkeypatch.setattr(pipe, "load_system", lambda path: object())
    with pytest.raises(ValueError, match="nothing to compare"):
        pipe.evaluate_pair("x", "ref.pdb", "pose.pdb", ["rmsd", "dockq"])


def test_a_metric_that_is_nan_for_one_id_is_not_a_failure(pipe):
    """A monomer has no interface; its dockq is NaN and the id is still evaluated."""
    evaluate = _fake({"a": {"rmsd": 0.1, "dockq": math.nan}})
    rows, missing = pipe.score_all({"a": "a.pdb"}, {"r": "r.pdb"}, 1, ["rmsd", "dockq"], "003_Peppr", evaluate)
    assert rows[0]["rmsd"] == 0.1 and math.isnan(rows[0]["dockq"])
    assert missing == []


def test_the_cause_is_trimmed_to_a_cell(pipe):
    rows, missing = pipe.score_all({"a": "a.pdb"}, {"r": "r.pdb"}, 1, ["rmsd"], "003_Peppr",
                                   _fake({"a": ValueError("x" * 500)}))
    assert len(missing[0]["cause"]) == 200


def _stream_json(tmp_path, name, paths_by_id):
    path = tmp_path / f"{name}.json"
    path.write_text(json.dumps({"name": "structures", "ids": list(paths_by_id),
                                "files": list(paths_by_id.values()),
                                "map_table": "", "format": "pdb"}))
    return str(path)


def _prepare(monkeypatch, tmp_path, poses, references, metrics):
    """Point sys.argv at stream files for poses and references; return the step folder main() writes to."""
    for path in list(poses.values()) + list(references.values()):
        if not os.path.exists(path):
            open(path, "w").close()
    (tmp_path / "config.yaml").write_text(json.dumps({"metrics": metrics}))
    step = tmp_path / "003_Peppr"
    monkeypatch.setattr(sys, "argv", [
        "pipe_peppr.py",
        "--config-yaml", str(tmp_path / "config.yaml"),
        "--structures-json", _stream_json(tmp_path, "poses", poses),
        "--references-json", _stream_json(tmp_path, "references", references),
        "--peppr-csv", str(step / "tables" / "peppr.csv"),
        "--local-missing-csv", str(step / "_execution" / "local_missing.csv")])
    return step


def _exit_code(pipe):
    try:
        pipe.main()
    except SystemExit as stop:
        return stop.code
    return 0


def _results(step):
    return (pd.read_csv(step / "tables" / "peppr.csv"),
            pd.read_csv(step / "_execution" / "local_missing.csv"))


def test_main_writes_a_row_per_id_and_exits_zero_when_some_were_evaluated(pipe, monkeypatch, tmp_path):
    step = _prepare(monkeypatch, tmp_path, {"a": str(tmp_path / "a.pdb"), "b": str(tmp_path / "b.pdb")},
                    {"r": str(tmp_path / "r.pdb")}, ["rmsd"])
    monkeypatch.setattr(pipe, "evaluate_pair",
                        _fake({"a": {"rmsd": 1.5}, "b": RuntimeError("could not read b")}))
    assert _exit_code(pipe) == 0

    table, missing = _results(step)
    assert list(table["id"]) == ["a", "b"] and list(table.columns) == ["id", "rmsd"]
    assert table.loc[0, "rmsd"] == 1.5 and math.isnan(table.loc[1, "rmsd"])
    assert list(missing["id"]) == ["b"] and missing.loc[0, "removed_by"] == "003_Peppr"


def test_main_exits_non_zero_when_every_id_failed_but_still_writes_the_table(pipe, monkeypatch, tmp_path):
    step = _prepare(monkeypatch, tmp_path, {"a": str(tmp_path / "a.pdb")}, {"r": str(tmp_path / "r.pdb")}, ["rmsd"])
    monkeypatch.setattr(pipe, "evaluate_pair", _fake({"a": RuntimeError("boom")}))
    assert _exit_code(pipe) == 1

    table, missing = _results(step)
    assert list(table["id"]) == ["a"] and math.isnan(table.loc[0, "rmsd"])
    assert list(missing["id"]) == ["a"]


def test_main_exits_non_zero_when_no_structure_file_exists(pipe, monkeypatch, tmp_path):
    step = _prepare(monkeypatch, tmp_path, {"a": str(tmp_path / "a.pdb")}, {"r": str(tmp_path / "r.pdb")}, ["rmsd"])
    os.remove(tmp_path / "a.pdb")
    assert _exit_code(pipe) == 1
    table, _ = _results(step)
    assert table.empty


def test_main_refuses_an_unknown_metric_key(pipe, monkeypatch, tmp_path):
    _prepare(monkeypatch, tmp_path, {"a": str(tmp_path / "a.pdb")}, {"r": str(tmp_path / "r.pdb")}, ["gdt"])
    assert _exit_code(pipe) == 1


# -- the metrics, against the real library ------------------------------------

def _place(a, b, c, length, angle, torsion):
    """Place the next atom of a chain from the previous three (NeRF)."""
    bc = (c - b) / np.linalg.norm(c - b)
    normal = np.cross(b - a, bc)
    normal /= np.linalg.norm(normal)
    offset = np.array([-length * math.cos(angle),
                       length * math.sin(angle) * math.cos(torsion),
                       length * math.sin(angle) * math.sin(torsion)])
    return c + offset[0] * bc + offset[1] * np.cross(normal, bc) + offset[2] * normal


def _helix(n_residues, phi=-57.0, psi=-47.0):
    """Backbone and CB of an alpha helix of alanines: a list of {atom name: xyz}."""
    rad = math.radians
    n_atom, ca = np.zeros(3), np.array([1.458, 0.0, 0.0])
    c = ca + 1.525 * np.array([math.cos(rad(69.0)), math.sin(rad(69.0)), 0.0])
    residues = []
    for _ in range(n_residues):
        residues.append({"N": n_atom, "CA": ca, "C": c,
                         "O": _place(n_atom, ca, c, 1.231, rad(120.8), rad(psi + 180)),
                         "CB": _place(c, n_atom, ca, 1.53, rad(110.5), rad(-122.6))})
        next_n = _place(n_atom, ca, c, 1.329, rad(116.2), rad(psi))
        next_ca = _place(ca, c, next_n, 1.458, rad(121.7), rad(180))
        c = _place(c, next_n, next_ca, 1.525, rad(111.0), rad(phi))
        n_atom, ca = next_n, next_ca
    return residues


def _record(record, serial, name, resname, chain, resi, xyz):
    return (f"{record:<6s}{serial:5d} {name:^4s} {resname:>3s} {chain}{resi:4d}    "
            f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00  0.00          {name[0]:>2s}")


def _complex_pdb(chain_b_shift=0.0, noise=0.0, seed=0, ligand=True, chain_b=True):
    """Two helices A and B, 11 A apart, and an ethylene glycol ligand C on chain A, with CONECT records.

    Peppr reads ligand bonds from CONECT in a PDB file, so the ligand carries them.
    """
    rng = np.random.default_rng(seed)
    helix = _helix(12)
    lines, serial = [], 1

    def jitter():
        return rng.normal(0.0, noise, 3) if noise else np.zeros(3)

    offsets = {"A": np.zeros(3)}
    if chain_b:
        offsets["B"] = np.array([0.0, 11.0 + chain_b_shift, 0.0])
    for chain, offset in offsets.items():
        for resi, atoms in enumerate(helix, 1):
            for name, xyz in atoms.items():
                lines.append(_record("ATOM", serial, name, "ALA", chain, resi, xyz + offset + jitter()))
                serial += 1
    if ligand:
        origin = helix[5]["CB"] + np.array([0.0, -3.6, 0.0])
        atoms = {"C1": origin, "O1": origin + [1.4, 0.3, 0.1],
                 "C2": origin + [-0.6, 1.2, 0.8], "O2": origin + [-1.5, 2.2, 0.3]}
        first = serial
        for name, xyz in atoms.items():
            lines.append(_record("HETATM", serial, name, "EDO", "C", 1, xyz + jitter()))
            serial += 1
        lines += [f"CONECT{first:5d}{first + 1:5d}{first + 2:5d}", f"CONECT{first + 2:5d}{first + 3:5d}"]
    return "\n".join(lines + ["END"]) + "\n"


@pytest.fixture(scope="module")
def complexes(tmp_path_factory):
    folder = tmp_path_factory.mktemp("peppr_complexes")
    texts = {
        "reference": _complex_pdb(),
        "pose": _complex_pdb(chain_b_shift=1.0, noise=0.3, seed=1),
        "monomer": _complex_pdb(ligand=False, chain_b=False),
    }
    paths = {}
    for name, text in texts.items():
        paths[name] = str(folder / f"{name}.pdb")
        open(paths[name], "w").write(text)
    return paths


@pytest.fixture(scope="module")
def library():
    package = pytest.importorskip("peppr", reason="peppr is not installed; the metric values are not checked")
    # With biopipelines/ on sys.path, a missing library resolves to the wrapper module of the same name.
    if not hasattr(package, "Evaluator"):
        pytest.skip("peppr is not installed; the metric values are not checked")
    pytest.importorskip("biotite", reason="biotite is not installed in this environment: the metric values are not checked")


@pytest.fixture(scope="module")
def scored(pipe, library, complexes):
    """Each case is scored once; the first call also pays for biotite loading the component dictionary."""
    ref = complexes["reference"]
    return {
        "identical": pipe.evaluate_pair("identical", ref, ref, ALL_METRICS),
        "pose": pipe.evaluate_pair("pose", ref, complexes["pose"], ALL_METRICS),
        "monomer": pipe.evaluate_pair("monomer", complexes["monomer"], complexes["monomer"], ALL_METRICS),
    }


def test_an_identical_complex_is_a_perfect_prediction(scored):
    perfect = scored["identical"]
    assert list(perfect) == ALL_METRICS
    assert perfect["rmsd"] == pytest.approx(0.0, abs=1e-3)
    for key in ("tm_score", "lddt", "lddt_pli", "fnat", "dockq"):
        assert perfect[key] == pytest.approx(1.0, abs=1e-3), key
    assert perfect["irmsd"] == pytest.approx(0.0, abs=1e-3)
    assert perfect["lrmsd"] == pytest.approx(0.0, abs=1e-3)


def test_a_displaced_pose_scores_worse_on_every_metric_but_is_computed(scored):
    pose, perfect = scored["pose"], scored["identical"]
    assert all(math.isfinite(value) for value in pose.values())
    assert 0.0 < pose["rmsd"] < 2.0
    assert pose["irmsd"] > perfect["irmsd"] and pose["lrmsd"] > perfect["lrmsd"]
    for key in ("tm_score", "lddt", "lddt_pli", "dockq"):
        assert pose[key] < perfect[key], key
    assert 0.0 <= pose["fnat"] <= 1.0


def test_a_monomer_has_no_interface_so_the_pair_metrics_are_nan(scored):
    monomer = scored["monomer"]
    assert math.isfinite(monomer["rmsd"]) and monomer["lddt"] == pytest.approx(1.0, abs=1e-3)
    for key in PAIR_METRICS + ["lddt_pli"]:
        assert math.isnan(monomer[key]), key


def test_an_mmcif_file_is_read_like_a_pdb_file(pipe, library, complexes, tmp_path):
    import biotite.structure.io.pdb as pdb
    import biotite.structure.io.pdbx as pdbx

    atoms = pdb.get_structure(pdb.PDBFile.read(complexes["reference"]), model=1, include_bonds=True)
    cif = pdbx.CIFFile()
    pdbx.set_structure(cif, atoms, include_bonds=True)
    path = str(tmp_path / "reference.cif")
    cif.write(path)

    values = pipe.evaluate_pair("x", path, path, ALL_METRICS)
    assert values["dockq"] == pytest.approx(1.0, abs=1e-3)
    assert values["lddt_pli"] == pytest.approx(1.0, abs=1e-3)


def test_a_metric_subset_returns_only_those_keys(pipe, library, complexes):
    values = pipe.evaluate_pair("x", complexes["reference"], complexes["pose"], ["dockq", "rmsd"])
    assert list(values) == ["dockq", "rmsd"]


def test_entities_that_cannot_be_matched_raise(pipe, library, complexes):
    with pytest.raises(ValueError, match="different entities"):
        pipe.evaluate_pair("x", complexes["reference"], complexes["monomer"], ALL_METRICS)


def _without_conect(text):
    return "".join(line for line in text.splitlines(keepends=True) if not line.startswith("CONECT"))


def test_a_ligand_without_bonds_is_named_in_the_failure(pipe, library, complexes, tmp_path):
    """A Boltz2 SMILES-ligand mmCIF has no chem_comp_bond and an unknown residue name, so the ligand has no bonds and peppr cannot match it."""
    bondless = tmp_path / "bondless.pdb"
    bondless.write_text(_without_conect(open(complexes["reference"]).read()))
    with pytest.raises(ValueError) as raised:
        pipe.evaluate_pair("x", complexes["reference"], str(bondless), ALL_METRICS)
    cause = str(raised.value)
    assert cause.startswith("pose ligand EDO has no bonds") and "give it a CCD code" in cause
    assert "give it a CCD code" in cause[:200]


def test_ligands_are_described_by_what_their_bonds_lack(pipe, library, complexes, tmp_path):
    import biotite.structure.io.pdb as pdb

    def ligands(path):
        return pipe.ligands_without_typed_bonds(pdb.get_structure(pdb.PDBFile.read(path), model=1, include_bonds=True))

    bondless = tmp_path / "bondless.pdb"
    bondless.write_text(_without_conect(open(complexes["reference"]).read()))
    assert ligands(complexes["reference"]) == {"EDO": "bonds of unknown order"}
    assert ligands(str(bondless)) == {"EDO": "no bonds"}
    assert ligands(complexes["monomer"]) == {}


def test_an_unsupported_file_type_raises(pipe, library, complexes, tmp_path):
    odd = tmp_path / "model.xyz"
    odd.write_text("")
    with pytest.raises(ValueError, match="unsupported structure file"):
        pipe.evaluate_pair("x", complexes["reference"], str(odd), ALL_METRICS)


def test_the_real_run_keeps_an_unmatched_id_as_a_nan_row(pipe, library, complexes, monkeypatch, tmp_path):
    """One prediction is fine, one has a different composition from the reference: both get a row."""
    step = _prepare(monkeypatch, tmp_path, {"fine": complexes["pose"], "different": complexes["monomer"]},
                    {"r": complexes["reference"]}, ALL_METRICS)
    assert _exit_code(pipe) == 0

    table, missing = _results(step)
    assert list(table["id"]) == ["fine", "different"]
    assert table.loc[0, ALL_METRICS].notna().all()
    assert table.loc[1, ALL_METRICS].isna().all()
    assert list(missing["id"]) == ["different"] and missing.loc[0, "kind"] == "failure"
