# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""FLOWR: parameter validation, the emitted bash, the output declaration, and the runtime half (staging, upstream command lines, SDF splitting, the tables) with upstream replaced by a fake that writes the SDF it would write."""

import glob
import importlib.util
import json
import math
import os
import re
import sys

import pandas as pd
import pytest

pytest.importorskip("yaml")
pytest.importorskip("rdkit")
import yaml  # noqa: E402
from rdkit import Chem  # noqa: E402
from rdkit.Chem import AllChem  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT = os.path.join(ROOT, "pipe_scripts", "pipe_flowr.py")

CFG = {
    "num_molecules": 5, "integration_steps": 100, "conditioning": None, "substructure": None,
    "sample_sizes": True, "num_heavy_atoms": None, "pocket_cutoff": 6.0, "protonate_pocket": False,
    "protonate_ligands": False, "batch_cost": 100, "seed": 42, "pocket_residues": None,
}

SELECTION_CFG = dict(CFG, sample_sizes=False, num_heavy_atoms=12, pocket_residues="A2-5+B6")


@pytest.fixture(scope="module")
def pipe():
    spec = importlib.util.spec_from_file_location("pipe_flowr", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# ── fixtures: a complex, a fake upstream, the stream JSONs ───────────────────

def _atom(record, serial, name, res, chain, seq, xyz, element):
    return (f"{record:<6s}{serial:5d} {name:^4s} {res:>3s} {chain}{seq:4d}    "
            f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00  0.00          {element:>2s}")


def write_complex(path, ligand=True):
    """Chain A residues 1-8, chain B residues 5-8, a three-atom ethanol named LIG and one water."""
    lines, serial = [], 1
    for chain, first in (("A", 1), ("B", 5)):
        for seq in range(first, 9):
            for name, dx in (("N", 0.0), ("CA", 1.4), ("C", 2.5)):
                lines.append(_atom("ATOM", serial, name, "ALA", chain, seq, (3.8 * seq + dx, 0.0 if chain == "A" else 8.0, 0.0), name[0]))
                serial += 1
    if ligand:
        for name, xyz, el in (("C1", (10.0, 4.0, 0.0), "C"), ("C2", (11.52, 4.0, 0.0), "C"), ("O1", (12.25, 5.16, 0.0), "O")):
            lines.append(_atom("HETATM", serial, name, "LIG", "A", 20, xyz, el))
            serial += 1
    lines.append(_atom("HETATM", serial, "O", "HOH", "A", 30, (0.0, 3.0, 0.0), "O"))
    path.write_text("\n".join(lines) + "\nEND\n")
    return str(path)


def sdf_with(path, smiles_list, bad_at=()):
    """An SDF of embedded 3D molecules; positions in `bad_at` hold a pentavalent carbon RDKit refuses to read."""
    blocks = []
    for i, smi in enumerate(smiles_list):
        if i in bad_at:
            mol = Chem.MolFromSmiles("C(C)(C)(C)(C)C", sanitize=False)
            blocks.append(Chem.MolToMolBlock(mol, kekulize=False) + "$$$$\n")
            continue
        mol = Chem.AddHs(Chem.MolFromSmiles(smi))
        AllChem.EmbedMolecule(mol, randomSeed=7)
        blocks.append(Chem.MolToMolBlock(Chem.RemoveHs(mol)) + "$$$$\n")
    path.write_text("".join(blocks))
    return str(path)


def write_stream(path, name, ids, files, fmt, map_table=""):
    path.write_text(json.dumps({"name": name, "ids": ids, "files": files, "map_table": map_table, "format": fmt}))
    return str(path)


def ligand_stream(tmp_path, code="LIG", smiles="CCO"):
    table = tmp_path / "ligand_map.csv"
    pd.DataFrame([{"id": "L1", "code": code, "smiles": smiles}]).to_csv(table, index=False)
    return write_stream(tmp_path / "ligand.json", "compounds", ["L1"], [], "csv", str(table))


class FakeUpstream:
    """Stands in for subprocess.run: records the command and writes the SDF upstream would."""

    def __init__(self, smiles_by_stem, returncode=0, log="Traceback (most recent call last):\nRuntimeError: no ligands found"):
        self.smiles_by_stem, self.returncode, self.log = smiles_by_stem, returncode, log
        self.commands = []

    def __call__(self, cmd, cwd=None, env=None, stdout=None, stderr=None):
        self.commands.append((cmd, cwd, env))
        arg = lambda flag: cmd[cmd.index(flag) + 1]  # noqa: E731
        stem = os.path.splitext(os.path.basename(arg("--pdb_file")))[0]
        if self.returncode:
            stdout.write(self.log + "\n")
        elif stem in self.smiles_by_stem:
            save_dir = arg("--save_dir")
            os.makedirs(save_dir, exist_ok=True)
            suffix = "_protonated" if "--protonate_generated_ligands" in cmd else ""
            sdf_with(__import__("pathlib").Path(save_dir) / f"samples_{stem}{suffix}.sdf", self.smiles_by_stem[stem])
        return type("Completed", (), {"returncode": self.returncode})()


def run_main(pipe, monkeypatch, tmp_path, fake, cfg=CFG, ids=("s1", "s2"), with_ligand=True):
    structures = [write_complex(tmp_path / f"{i}.pdb") for i in ids]
    checkpoint = tmp_path / "flowr_noHs.ckpt"
    checkpoint.write_text("weights")
    config = tmp_path / "config.yaml"
    config.write_text(yaml.safe_dump(cfg))
    out = {"compounds": tmp_path / "compounds" / "compounds.csv", "map": tmp_path / "structures" / "structures_map.csv",
           "pockets": tmp_path / "tables" / "pockets.csv", "local": tmp_path / "_execution" / "local_missing.csv"}
    argv = ["pipe_flowr.py", "--config-yaml", str(config),
            "--structures-json", write_stream(tmp_path / "s.json", "structures", list(ids), structures, "pdb"),
            "--flowr-repo", str(tmp_path / "repo"), "--checkpoint", str(checkpoint),
            "--runs-folder", str(tmp_path / "_execution" / "runs"),
            "--compounds-csv", str(out["compounds"]), "--structures-folder", str(tmp_path / "structures"),
            "--structures-map-csv", str(out["map"]), "--pockets-csv", str(out["pockets"]),
            "--local-missing-csv", str(out["local"])]
    if with_ligand:
        argv += ["--ligand-json", ligand_stream(tmp_path)]
    monkeypatch.setattr(pipe.subprocess, "run", fake)
    monkeypatch.setattr(sys, "argv", argv)
    pipe.main()
    return {k: pd.read_csv(v, keep_default_na=False) if k == "compounds" else pd.read_csv(v) for k, v in out.items()}


# ── upstream command lines ───────────────────────────────────────────────────

def command(pipe, cfg=CFG, ligand_file="ref.sdf", prefix=(), workers=4):
    return pipe.build_command(cfg, "pocket.pdb", ligand_file, "out", "w.ckpt", list(prefix), workers)


def value(cmd, flag):
    return cmd[cmd.index(flag) + 1]


def test_the_de_novo_command_carries_the_documented_defaults(pipe):
    cmd = command(pipe)
    assert cmd[:3] == ["python", "-m", "flowr.gen.generate_from_pdb"]
    assert value(cmd, "--sample_n_molecules_per_target") == "5" and value(cmd, "--integration_steps") == "100"
    assert value(cmd, "--arch") == "pocket" and value(cmd, "--pocket_type") == "holo" and value(cmd, "--gpus") == "1"
    assert value(cmd, "--ligand_file") == "ref.sdf" and "--cut_pocket" in cmd and value(cmd, "--pocket_cutoff") == "6.0"
    assert "--sample_mol_sizes" in cmd and "--filter_valid_unique" in cmd
    for absent in ("--num_heavy_atoms", "--protonate_pocket", "--compute_interactions", "--protonate_generated_ligands"):
        assert absent not in cmd
    assert not any(flag.endswith("_inpainting") for flag in cmd)


@pytest.mark.parametrize("task, flag", [
    ("interaction", "--interaction_inpainting"), ("functional_group", "--func_group_inpainting"),
    ("scaffold", "--scaffold_inpainting"), ("linker", "--linker_inpainting"),
    ("substructure", "--substructure_inpainting"),
])
def test_each_conditioning_task_sets_exactly_its_inpainting_flag(pipe, task, flag):
    cmd = command(pipe, dict(CFG, conditioning=task, sample_sizes=False, substructure="c1ccccc1"))
    assert [f for f in cmd if f.endswith("_inpainting")] == [flag]
    assert "--sample_mol_sizes" not in cmd


def test_only_interaction_conditioning_computes_interactions_and_only_substructure_takes_a_query(pipe):
    assert "--compute_interactions" in command(pipe, dict(CFG, conditioning="interaction", sample_sizes=False))
    assert "--compute_interactions" not in command(pipe, dict(CFG, conditioning="scaffold", sample_sizes=False))
    cmd = command(pipe, dict(CFG, conditioning="substructure", sample_sizes=False, substructure="C(=O)N"))
    assert value(cmd, "--substructure") == "C(=O)N"
    assert "--substructure" not in command(pipe, dict(CFG, conditioning="linker", sample_sizes=False))


def test_the_conditioning_names_of_the_wrapper_and_the_pipe_agree(pipe):
    from biopipelines.flowr import CONDITIONING
    assert set(CONDITIONING) == set(pipe.CONDITIONING_FLAGS)


def test_a_residue_pocket_passes_a_size_and_no_reference_ligand(pipe):
    cmd = command(pipe, SELECTION_CFG, ligand_file=None)
    assert value(cmd, "--num_heavy_atoms") == "12"
    for absent in ("--ligand_file", "--cut_pocket", "--pocket_cutoff", "--sample_mol_sizes"):
        assert absent not in cmd


def test_protonating_the_molecules_passes_the_worker_count(pipe):
    cmd = command(pipe, dict(CFG, protonate_ligands=True, protonate_pocket=True), workers=7)
    assert "--protonate_generated_ligands" in cmd and value(cmd, "--num_workers") == "7" and "--protonate_pocket" in cmd
    assert "--num_workers" not in command(pipe)


def test_a_container_prefix_leads_the_command(pipe):
    cmd = command(pipe, prefix=["apptainer", "exec", "img.sif"])
    assert cmd[:6] == ["apptainer", "exec", "img.sif", "python", "-m", "flowr.gen.generate_from_pdb"]


# ── staging ──────────────────────────────────────────────────────────────────

def test_the_reference_ligand_keeps_its_coordinates_and_takes_bond_orders_from_the_smiles(pipe, tmp_path):
    pdb = write_complex(tmp_path / "c.pdb")
    sdf = pipe.stage_reference_ligand(pdb, "lig", "CCO", str(tmp_path), "c")
    mol = Chem.MolFromMolFile(sdf)
    assert Chem.MolToSmiles(mol) == "CCO" and mol.GetNumAtoms() == 3
    assert [round(v, 2) for v in mol.GetConformer().GetAtomPosition(1)] == [11.52, 4.0, 0.0]


def test_a_missing_ligand_residue_is_named_in_the_error(pipe, tmp_path):
    pdb = write_complex(tmp_path / "c.pdb", ligand=False)
    with pytest.raises(ValueError, match="'ZIT'"):
        pipe.extract_ligand_block(pdb, "ZIT")


def test_only_the_first_copy_of_the_ligand_is_extracted(pipe, tmp_path):
    pdb = tmp_path / "c.pdb"
    pdb.write_text(_atom("HETATM", 1, "C1", "LIG", "A", 20, (0, 0, 0), "C") + "\n"
                   + _atom("HETATM", 2, "C1", "LIG", "B", 21, (9, 9, 9), "C") + "\n")
    assert pipe.extract_ligand_block(str(pdb), "LIG").count("HETATM") == 1


def test_a_residue_pocket_is_chain_aware(pipe, tmp_path):
    pdb = write_complex(tmp_path / "c.pdb")
    out = tmp_path / "pocket.pdb"
    assert pipe.stage_pocket_pdb(pdb, "A2-5+B6", str(out)) == 5
    kept = {(l[21], int(l[22:26])) for l in out.read_text().splitlines() if l.startswith("ATOM")}
    assert kept == {("A", 2), ("A", 3), ("A", 4), ("A", 5), ("B", 6)}
    assert "HETATM" not in out.read_text()


def test_a_bare_chain_letter_takes_the_whole_chain(pipe, tmp_path):
    out = tmp_path / "pocket.pdb"
    assert pipe.stage_pocket_pdb(write_complex(tmp_path / "c.pdb"), "B", str(out)) == 4


def test_a_selection_that_matches_nothing_is_an_error(pipe, tmp_path):
    with pytest.raises(ValueError, match="lies in"):
        pipe.stage_pocket_pdb(write_complex(tmp_path / "c.pdb"), "C1-5", str(tmp_path / "p.pdb"))


def test_an_unparseable_substructure_is_refused_before_any_run(pipe):
    pipe.check_config(dict(CFG, substructure="c1ccccc1"))
    with pytest.raises(ValueError, match="substructure"):
        pipe.check_config(dict(CFG, substructure="not a smiles ((("))


# ── collecting what upstream wrote ───────────────────────────────────────────

def test_molecules_become_compounds_rows_and_one_pose_each(pipe, tmp_path):
    sdf = sdf_with(tmp_path / "samples.sdf", ["CCO", "CCN", "c1ccccc1"])
    (tmp_path / "poses").mkdir()
    compounds, poses, missing = pipe.collect_molecules("s1", sdf, 3, str(tmp_path / "poses"), "004_FLOWR")
    assert [r["id"] for r in compounds] == ["s1_1", "s1_2", "s1_3"]
    assert [r["smiles"] for r in compounds] == ["CCO", "CCN", "c1ccccc1"]
    assert {r["structures.id"] for r in compounds} == {"s1"} and [r["sample"] for r in compounds] == [1, 2, 3]
    assert {(r["format"], r["code"], r["ccd"]) for r in compounds} == {("smiles", "LIG", "")}
    assert missing == []
    pose = Chem.MolFromMolFile(poses[0]["file"])
    assert poses[0]["id"] == "s1_1" and pose.GetProp("_Name") == "s1_1"
    original = Chem.SDMolSupplier(sdf)[0].GetConformer().GetPositions()
    assert abs(pose.GetConformer().GetPositions() - original).max() < 1e-3


def test_ids_upstream_did_not_return_are_filter_rows(pipe, tmp_path):
    sdf = sdf_with(tmp_path / "samples.sdf", ["CCO", "CCN"])
    (tmp_path / "poses").mkdir()
    _, _, missing = pipe.collect_molecules("s1", sdf, 5, str(tmp_path / "poses"), "004_FLOWR")
    assert [(r["id"], r["kind"], r["removed_by"]) for r in missing] == [(f"s1_{k}", "filter", "004_FLOWR") for k in (3, 4, 5)]


def test_an_unreadable_record_is_a_failure_row_and_keeps_the_numbering(pipe, tmp_path):
    sdf = sdf_with(tmp_path / "samples.sdf", ["CCO", "CCN", "c1ccccc1"], bad_at=(1,))
    (tmp_path / "poses").mkdir()
    compounds, _, missing = pipe.collect_molecules("s1", sdf, 3, str(tmp_path / "poses"), "004_FLOWR")
    assert [r["id"] for r in compounds] == ["s1_1", "s1_3"]
    assert [(r["id"], r["kind"]) for r in missing] == [("s1_2", "failure")]


def test_hydrogens_stay_in_the_pose_and_out_of_the_smiles(pipe, tmp_path):
    mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
    AllChem.EmbedMolecule(mol, randomSeed=3)
    sdf = tmp_path / "samples.sdf"
    writer = Chem.SDWriter(str(sdf))
    writer.write(mol)
    writer.close()
    (tmp_path / "poses").mkdir()
    compounds, poses, _ = pipe.collect_molecules("s1", str(sdf), 1, str(tmp_path / "poses"), "004_FLOWR")
    assert compounds[0]["smiles"] == "CCO"
    assert Chem.MolFromMolFile(poses[0]["file"], removeHs=False).GetNumAtoms() == 9


# ── one pocket, and the whole step ───────────────────────────────────────────

def paths(tmp_path):
    return {"runs_folder": str(tmp_path / "runs"), "structures_folder": str(tmp_path / "poses"), "checkpoint": "w.ckpt",
            "flowr_repo": str(tmp_path / "repo"), "num_workers": 2}


def test_a_failed_upstream_run_reports_the_last_line_of_its_log(pipe, monkeypatch, tmp_path):
    monkeypatch.setattr(pipe.subprocess, "run", FakeUpstream({}, returncode=1))
    pdb = write_complex(tmp_path / "s1.pdb")
    with pytest.raises(RuntimeError, match="status 1: RuntimeError: no ligands found"):
        pipe.run_pocket("s1", pdb, CFG, {"code": "LIG", "smiles": "CCO"}, paths(tmp_path), [], "004_FLOWR")


def test_the_reference_run_reads_a_copy_of_the_structure_without_conect_records(pipe, monkeypatch, tmp_path):
    """Upstream's reader raised IndexError on the CONECT records of a crystal (3RY2); everything else must be kept."""
    fake = FakeUpstream({"s1": ["CCO"]})
    monkeypatch.setattr(pipe.subprocess, "run", fake)
    pdb = write_complex(tmp_path / "s1.pdb")
    with open(pdb, "a") as handle:
        handle.write("CONECT   49   50\nCONECT   50   49   51\n")
    pipe.run_pocket("s1", pdb, CFG, {"code": "LIG", "smiles": "CCO"}, paths(tmp_path), [], "004_FLOWR")

    staged = value(fake.commands[0][0], "--pdb_file")
    text = open(staged).read()
    assert staged != pdb and "CONECT" not in text
    assert text.count("ATOM") + text.count("HETATM") == open(pdb).read().count("ATOM") + open(pdb).read().count("HETATM")
    assert "CONECT" in open(pdb).read()


def test_an_upstream_run_without_its_sdf_is_an_error(pipe, monkeypatch, tmp_path):
    monkeypatch.setattr(pipe.subprocess, "run", FakeUpstream({}))
    pdb = write_complex(tmp_path / "s1.pdb")
    with pytest.raises(RuntimeError, match="samples_s1.sdf"):
        pipe.run_pocket("s1", pdb, CFG, {"code": "LIG", "smiles": "CCO"}, paths(tmp_path), [], "004_FLOWR")


def test_the_step_writes_every_table_and_pose(pipe, monkeypatch, tmp_path):
    fake = FakeUpstream({"s1": ["CCO", "CCN", "c1ccccc1"], "s2": ["CCO"]})
    out = run_main(pipe, monkeypatch, tmp_path, fake)
    assert list(out["compounds"].columns) == pipe.COMPOUND_COLUMNS
    assert list(out["compounds"]["id"]) == ["s1_1", "s1_2", "s1_3", "s2_1"]
    assert list(out["map"]["id"]) == ["s1_1", "s1_2", "s1_3", "s2_1"] and list(out["map"].columns) == ["id", "file", "structures.id"]
    assert all(os.path.isfile(p) for p in out["map"]["file"])
    assert out["pockets"].to_dict("records") == [
        {"id": "s1", "pocket": "ligand:LIG", "n_requested": 5, "n_generated": 3},
        {"id": "s2", "pocket": "ligand:LIG", "n_requested": 5, "n_generated": 1}]
    local = out["local"]
    assert list(local.columns) == ["id", "removed_by", "kind", "cause"]
    assert set(local["kind"]) == {"filter"} and len(local) == 2 + 4


def test_every_declared_id_is_either_produced_or_excused(pipe, monkeypatch, tmp_path):
    from biopipelines.datastream import DataStream
    from biopipelines.combinatorics import generate_multiplied_ids_pattern

    out = run_main(pipe, monkeypatch, tmp_path, FakeUpstream({"s1": ["CCO", "CCN", "c1ccccc1"], "s2": ["CCO"]}))
    declared = DataStream(name="compounds", ids=generate_multiplied_ids_pattern(["s1", "s2"], "<1..5>")).ids_expanded
    accounted = set(out["compounds"]["id"]) | set(out["local"]["id"])
    assert accounted == set(declared) and len(declared) == 10


def test_each_pocket_runs_upstream_once_in_the_repository_with_it_on_the_path(pipe, monkeypatch, tmp_path):
    fake = FakeUpstream({"s1": ["CCO"], "s2": ["CCO"]})
    run_main(pipe, monkeypatch, tmp_path, fake)
    assert len(fake.commands) == 2
    cmd, cwd, env = fake.commands[0]
    assert cwd == str(tmp_path / "repo") and env["PYTHONPATH"].split(os.pathsep)[0] == str(tmp_path / "repo")
    assert value(cmd, "--ckpt_path") == str(tmp_path / "flowr_noHs.ckpt")


def test_a_pocket_that_fails_keeps_its_table_row_as_nan_and_a_failure_row(pipe, monkeypatch, tmp_path):
    class Half(FakeUpstream):
        def __call__(self, cmd, cwd=None, env=None, stdout=None, stderr=None):
            self.returncode = 1 if "s2" in cmd[cmd.index("--pdb_file") + 1] else 0
            return super().__call__(cmd, cwd, env, stdout, stderr)
    out = run_main(pipe, monkeypatch, tmp_path, Half({"s1": ["CCO", "CCN"]}))
    pockets = out["pockets"].set_index("id")
    assert pockets.loc["s1", "n_generated"] == 2 and math.isnan(pockets.loc["s2", "n_generated"])
    failures = out["local"][out["local"]["kind"] == "failure"]
    assert list(failures["id"]) == ["s2"] and "no ligands found" in failures["cause"].iloc[0]
    assert list(out["compounds"]["id"]) == ["s1_1", "s1_2"]


def test_a_step_where_every_pocket_fails_exits_non_zero(pipe, monkeypatch, tmp_path):
    with pytest.raises(SystemExit) as exit_info:
        run_main(pipe, monkeypatch, tmp_path, FakeUpstream({}, returncode=1))
    assert exit_info.value.code == 1
    assert (tmp_path / "tables" / "pockets.csv").exists()


def test_a_residue_pocket_runs_on_the_staged_pdb_without_a_reference(pipe, monkeypatch, tmp_path):
    fake = FakeUpstream({"s1": ["CCO"], "s2": ["CCN"]})
    out = run_main(pipe, monkeypatch, tmp_path, fake, cfg=SELECTION_CFG, with_ligand=False)
    cmd = fake.commands[0][0]
    staged = value(cmd, "--pdb_file")
    assert staged == str(tmp_path / "_execution" / "runs" / "s1" / "s1.pdb")
    assert "--ligand_file" not in cmd and value(cmd, "--num_heavy_atoms") == "12"
    assert set(out["pockets"]["pocket"]) == {"residues:A2-5+B6"}
    assert "HETATM" not in open(staged).read()


def test_a_missing_checkpoint_stops_the_step_before_any_run(pipe, monkeypatch, tmp_path):
    fake = FakeUpstream({"s1": ["CCO"]})
    monkeypatch.setattr(pipe.subprocess, "run", fake)
    config = tmp_path / "config.yaml"
    config.write_text(yaml.safe_dump(CFG))
    argv = ["pipe_flowr.py", "--config-yaml", str(config), "--structures-json", "s.json", "--flowr-repo", "r",
            "--checkpoint", str(tmp_path / "absent.ckpt"), "--runs-folder", "r", "--compounds-csv", "c",
            "--structures-folder", "s", "--structures-map-csv", "m", "--pockets-csv", "p", "--local-missing-csv", "l"]
    monkeypatch.setattr(sys, "argv", argv)
    with pytest.raises(SystemExit) as exit_info:
        pipe.main()
    assert exit_info.value.code == 1 and fake.commands == []


# ── wrapper ──────────────────────────────────────────────────────────────────

def _structures():
    from biopipelines.mock import Mock
    return Mock(ids=["s1", "s2"], streams={"structures": {"format": "pdb", "file": "<id>.pdb"}})


def _step_text(script_path):
    step = sorted(glob.glob(os.path.join(os.path.dirname(script_path), "*_FLOWR.sh")))[0]
    return open(step, encoding="utf-8", errors="replace").read()


def test_the_outputs_are_a_compounds_stream_and_one_sdf_pose_per_molecule(local_config, isolated_cwd, new_pipeline):
    from biopipelines.flowr import FLOWR

    with new_pipeline("flowr_outputs"):
        out = FLOWR(structures=_structures(), ligand="ZIT", num_molecules=3)
    compounds, structures = out.streams.compounds, out.streams.structures
    assert compounds.format == "csv" and compounds.files == [] and structures.format == "sdf"
    assert compounds.ids == structures.ids == ["s1_<1..3>", "s2_<1..3>"]
    assert list(out.tables.compounds.info.columns) == ["id", "structures.id", "sample", "format", "code", "ccd", "smiles"]
    assert list(out.tables.pockets.info.columns) == ["id", "pocket", "n_requested", "n_generated"]
    assert list(out.tables.missing.info.columns) == ["id", "removed_by", "kind", "cause"]
    assert out.tables.compounds.info.path == compounds.map_table


def test_exactly_one_pocket_definition_is_required(local_config, isolated_cwd, new_pipeline):
    from biopipelines.flowr import FLOWR

    with new_pipeline("flowr_pocket"):
        with pytest.raises(ValueError, match="exactly one"):
            FLOWR(structures=_structures())
        with pytest.raises(ValueError, match="exactly one"):
            FLOWR(structures=_structures(), ligand="ZIT", pocket_residues="A1-5", num_heavy_atoms=20)


@pytest.mark.parametrize("kwargs, message", [
    ({"num_molecules": 0}, "num_molecules"),
    ({"num_molecules": 2.5}, "num_molecules"),
    ({"integration_steps": 0}, "integration_steps"),
    ({"batch_cost": -1}, "batch_cost"),
    ({"seed": -1}, "seed"),
    ({"conditioning": "fragment"}, "conditioning"),
    ({"conditioning": "substructure"}, "substructure"),
    ({"substructure": "c1ccccc1"}, "only used with"),
    ({"conditioning": "scaffold", "sample_sizes": True}, "sample_sizes"),
    ({"num_heavy_atoms": 20}, "num_heavy_atoms"),
    ({"pocket_cutoff": 0}, "pocket_cutoff"),
])
def test_bad_parameters_with_a_reference_ligand_are_refused_at_configuration(local_config, isolated_cwd, new_pipeline, kwargs, message):
    from biopipelines.flowr import FLOWR

    with new_pipeline("flowr_bad"):
        with pytest.raises(ValueError, match=message):
            FLOWR(structures=_structures(), ligand="ZIT", **kwargs)


@pytest.mark.parametrize("kwargs, message", [
    ({}, "num_heavy_atoms"),
    ({"num_heavy_atoms": 4}, "num_heavy_atoms"),
    ({"num_heavy_atoms": 81}, "num_heavy_atoms"),
    ({"num_heavy_atoms": 20, "conditioning": "scaffold"}, "needs ligand"),
    ({"num_heavy_atoms": 20, "pocket_cutoff": 8.0}, "pocket_cutoff"),
    ({"num_heavy_atoms": 20, "pocket_residues": "10-20"}, "chain"),
])
def test_bad_parameters_with_a_residue_pocket_are_refused_at_configuration(local_config, isolated_cwd, new_pipeline, kwargs, message):
    from biopipelines.flowr import FLOWR

    with new_pipeline("flowr_bad_sele"):
        with pytest.raises(ValueError, match=message):
            FLOWR(structures=_structures(), **{"pocket_residues": "A10-25", **kwargs})


def test_a_structure_that_is_not_a_pdb_is_refused(local_config, isolated_cwd, new_pipeline):
    from biopipelines.flowr import FLOWR
    from biopipelines.mock import Mock

    with new_pipeline("flowr_cif"):
        cifs = Mock(ids=["s1"], streams={"structures": {"format": "cif", "file": "<id>.cif"}})
        with pytest.raises(ValueError, match="PDB"):
            FLOWR(structures=cifs, ligand="ZIT")


def test_the_size_default_follows_the_task(local_config, isolated_cwd, new_pipeline):
    from biopipelines.flowr import FLOWR

    pipeline = new_pipeline("flowr_sizes")
    with pipeline:
        FLOWR(structures=_structures(), ligand="ZIT")
        FLOWR(structures=_structures(), ligand="ZIT", conditioning="scaffold")
        FLOWR(structures=_structures(), ligand="ZIT", sample_sizes=False)
        FLOWR(structures=_structures(), pocket_residues="A1-12", num_heavy_atoms=20)
        FLOWR(structures=_structures(), pocket_residues="A1-12", num_heavy_atoms=20, sample_sizes=True)
    flowr_tools = [t for t in pipeline.tools if t.TOOL_NAME == "FLOWR"]
    assert [t.sample_sizes for t in flowr_tools] == [True, False, False, False, True]


def test_the_emitted_script_runs_the_pipe_and_writes_the_config(local_config, isolated_cwd, new_pipeline):
    from biopipelines.flowr import FLOWR

    pipeline = new_pipeline("flowr_script")
    with pipeline:
        FLOWR(structures=_structures(), ligand="ZIT", conditioning="substructure", substructure="c1ccccc1", num_molecules=7)
        script_path = pipeline.save()
    text = _step_text(script_path)
    assert "pipe_flowr.py" in text and "--ligand-json" in text and "--container-prefix" in text
    tool = [t for t in pipeline.tools if t.TOOL_NAME == "FLOWR"][0]
    config = yaml.safe_load(open(tool.config_yaml))
    assert config["num_molecules"] == 7 and config["conditioning"] == "substructure" and config["substructure"] == "c1ccccc1"
    assert config["sample_sizes"] is False and config["pocket_residues"] is None


def test_a_residue_pocket_script_passes_no_ligand(local_config, isolated_cwd, new_pipeline):
    from biopipelines.flowr import FLOWR

    pipeline = new_pipeline("flowr_script_sele")
    with pipeline:
        FLOWR(structures=_structures(), pocket_residues="A10-25+B7", num_heavy_atoms=24)
        script_path = pipeline.save()
    text = _step_text(script_path)
    assert "--ligand-json" not in text
    tool = [t for t in pipeline.tools if t.TOOL_NAME == "FLOWR"][0]
    config = yaml.safe_load(open(tool.config_yaml))
    assert config["pocket_residues"] == "A10-25+B7" and config["num_heavy_atoms"] == 24


def test_every_pipe_flag_has_a_matching_argument(local_config, isolated_cwd, new_pipeline):
    from biopipelines.flowr import FLOWR

    pipeline = new_pipeline("flowr_flags")
    with pipeline:
        FLOWR(structures=_structures(), ligand="ZIT")
        script_path = pipeline.save()
    text = _step_text(script_path)
    call = text[text.index("pipe_flowr.py"):]
    call = call[:call.index("\n\n")] if "\n\n" in call else call
    emitted = set(re.findall(r"--[a-z][a-z-]+", call))
    declared = set(re.findall(r'add_argument\("(--[a-z-]+)"', open(SCRIPT, encoding="utf-8").read()))
    assert emitted == declared


def test_the_install_script_pins_the_checkpoint_record_and_checks_the_import(monkeypatch):
    from biopipelines.flowr import FLOWR, CHECKPOINT_MD5, CHECKPOINT_URL

    monkeypatch.setattr(FLOWR, "_install_env", classmethod(lambda cls, env_manager, env_name=None: "FLOWR"))
    script = FLOWR._install_script({"biopipelines": "/bp", "FLOWR": "/data/FLOWR"}, env_manager="mamba")
    assert "https://github.com/jule-c/flowr.git" in script
    assert "records/15737419/" in CHECKPOINT_URL and CHECKPOINT_URL in script and CHECKPOINT_MD5 in script
    assert os.path.join("/data/FLOWR", "checkpoints", "flowr_noHs.ckpt") in script
    assert ", flowr.gen.generate_from_pdb\"" in script and 'touch "$INSTALL_SUCCESS"' in script
    assert script.index("mamba env create") < script.index("md5sum")


def test_the_env_files_keep_upstreams_torch_stack():
    spec = yaml.safe_load(open(os.path.join(ROOT, "environments", "FLOWR.yaml")))
    assert spec["name"] == "FLOWR"
    assert "pytorch::pytorch=2.5.1" in spec["dependencies"] and "pytorch-cuda=11.8" in spec["dependencies"]
    pip = open(os.path.join(ROOT, "environments", "FLOWR.pip.txt")).read()
    assert "biotite==1.0.1" in pip and "numpy==1.26.2" in pip
