"""LigandDdG: wrapper validation, generated script, ligand parameterization and the per-scheme result tables.

Everything here runs without PyRosetta, and the Rosetta side is replaced by recorders. The chemistry tests need RDKit; the last tests need PyRosetta and skip without it. `molfile_to_params.py` itself is never run: what is tested is the SDF it is fed, the command line it is given and the handling of its result.
"""

import glob
import importlib.util
import json
import os
import re
import sys
import types
import xml.etree.ElementTree as ET

import pytest

pytest.importorskip("pandas")
pytest.importorskip("yaml")

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HELPER = os.path.join(REPO_ROOT, "pipe_scripts", "pipe_ligand_ddg.py")


@pytest.fixture(scope="module")
def pipe():
    spec = importlib.util.spec_from_file_location("pipe_ligand_ddg", HELPER)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def step_script(pipeline, tool):
    found = glob.glob(os.path.join(pipeline.folders["output"], "**", f"*_{tool}.sh"), recursive=True)
    assert found, f"no {tool} step script under {pipeline.folders['output']}"
    # the framework writes step scripts in the platform encoding, so an em dash in a comment is not valid UTF-8 on Windows
    return open(found[0], encoding="utf-8", errors="replace").read()


def build(new_pipeline, name, fmt="pdb", ligand="smiles", **kwargs):
    from biopipelines.ligand import Ligand
    from biopipelines.ligand_ddg import LigandDdG
    from biopipelines.mock import Mock

    pipeline = new_pipeline(name)
    with pipeline:
        source = Mock(ids=["s1", "s2"], streams={"structures": {"format": fmt, "file": "<id>.pdb"}})
        if ligand == "smiles":
            ligand = Ligand(smiles="CC(=O)O", codes="LIG")
        LigandDdG(structures=source.streams.structures, ligand=ligand, **kwargs)
        pipeline.save()
    return pipeline, pipeline.tools[-1]


# -- wrapper validation ---------------------------------------------------------------------------

@pytest.mark.parametrize("kwargs,match", [
    ({"scoring": "medium"}, "scoring must be one of"),
    ({"binder_chain": "AB"}, "single-character"),
    ({"ligand_chain": ""}, "single-character"),
    ({"binder_chain": "B", "ligand_chain": "B"}, "different chains"),
    ({"scoring": "hard", "replicates": 5}, "scoring='soft' only"),
    ({"replicates": 0}, "positive integer"),
    ({"ph": 15}, "between 0 and 14"),
    ({"ph": "7.4"}, "between 0 and 14"),
    ({"params_file": "  "}, "params_file must be a path"),
    ({"seed": 0}, "positive integer"),
    ({"seed": 1.5}, "positive integer"),
])
def test_bad_parameters_are_refused_before_a_job_is_queued(local_config, isolated_cwd, new_pipeline, kwargs, match):
    with pytest.raises(ValueError, match=match):
        build(new_pipeline, "lddg_bad", **kwargs)


def test_cif_structures_are_refused(local_config, isolated_cwd, new_pipeline):
    with pytest.raises(ValueError, match="reads PDB files"):
        build(new_pipeline, "lddg_cif", fmt="cif")


def test_a_bare_residue_code_is_shorthand_for_a_ligand(local_config, isolated_cwd, new_pipeline):
    _, tool = build(new_pipeline, "lddg_bare", ligand="LIG")
    assert tool.ligand_stream is not None and len(tool.ligand_stream) == 1


def test_the_default_scheme_is_soft_with_five_replicates(local_config, isolated_cwd, new_pipeline):
    import yaml

    _, tool = build(new_pipeline, "lddg_default")
    config = yaml.safe_load(open(tool.config_yaml))
    assert config == {"scoring": "soft", "binder_chain": "A", "ligand_chain": "B", "replicates": 5,
                      "ph": 7.4, "params_file": "", "seed": 42}


def test_the_hard_scheme_runs_once_and_says_so_in_the_config(local_config, isolated_cwd, new_pipeline):
    import yaml

    _, tool = build(new_pipeline, "lddg_hard", scoring="hard", ph=6.0, params_file="/data/LIG.params")
    config = yaml.safe_load(open(tool.config_yaml))
    assert config["scoring"] == "hard" and config["replicates"] is None
    assert config["ph"] == 6.0 and config["params_file"] == "/data/LIG.params"


# -- generated script and declared outputs ----------------------------------------------------------

def test_script_flags_match_the_pipe_script(local_config, isolated_cwd, new_pipeline):
    pipeline, _ = build(new_pipeline, "lddg_flags")
    script = step_script(pipeline, "LigandDdG")
    sent = set(re.findall(r"--[a-z][a-z-]*", script.split("pipe_ligand_ddg.py", 1)[1].split("Propagating", 1)[0]))
    accepted = set(re.findall(r'add_argument\("(--[a-z-]+)"', open(HELPER, encoding="utf-8").read()))
    assert sent == accepted, (sent, accepted)
    assert "Propagating upstream missing manifest(s)" in script


@pytest.mark.parametrize("scoring,columns", [
    ("hard", ["id", "ddg_norepack"]),
    ("soft", ["id", "ddg_soft", "ddg_soft_std", "n_replicates_used"]),
])
def test_the_result_table_follows_the_scheme(local_config, isolated_cwd, new_pipeline, scoring, columns):
    _, tool = build(new_pipeline, f"lddg_{scoring}_cols", scoring=scoring)
    tables = tool.get_output_files()["tables"]
    assert list(tables["ligand_ddg"].info.columns) == columns
    assert list(tables["missing"].info.columns) == ["id", "removed_by", "kind", "cause"]


def test_install_goes_to_the_shared_env(local_config, monkeypatch):
    from biopipelines.ligand_ddg import LigandDdG

    monkeypatch.setattr(LigandDdG, "_install_env", classmethod(lambda cls, manager, name=None: name or "PyRosetta"))
    script = LigandDdG._install_script({"biopipelines": "/repo"}, env_manager="mamba")
    assert "mamba env create -f" in script and "environments/PyRosetta" in script
    assert "molfile_to_params.py" in script and "sparse-checkout" in script


def test_both_tools_share_one_environment_in_every_shipped_config():
    import yaml

    for variant in ("cluster", "colab", "container"):
        config = yaml.safe_load(open(os.path.join(REPO_ROOT, f"config.{variant}.yaml"), encoding="utf-8"))
        envs = config["environments"]
        assert envs["PPIDdG"] == envs["LigandDdG"] == "PyRosetta", variant


# -- PDB handling -------------------------------------------------------------------------------------

def atom_line(record, serial, name, resname, chain, resseq, xyz, element, bfactor=20.0):
    padded = f" {name:<3}" if len(element) == 1 and len(name) <= 3 else f"{name:<4}"
    return (f"{record:<6}{serial:5d} {padded} {resname:>3s} {chain}{resseq:4d}    "
            f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00{bfactor:6.2f}          {element:>2s}\n")


def write_complex(path, ligand_lines, extra=()):
    protein = [
        atom_line("ATOM", 1, "N", "ALA", "A", 1, (0.0, 0.0, 0.0), "N"),
        atom_line("ATOM", 2, "H", "ALA", "A", 1, (1.0, 0.0, 0.0), "H"),
        atom_line("ATOM", 3, "CA", "ALA", "A", 1, (0.0, 1.4, 0.0), "C"),
        atom_line("HETATM", 4, "O", "HOH", "W", 1, (9.0, 9.0, 9.0), "O"),
    ]
    with open(path, "w") as fh:
        fh.writelines(protein + list(ligand_lines) + list(extra) + ["END\n"])
    return str(path)


def embedded_ligand(smiles, code="LIG", chain="B", resseq=7, with_h=False, names=None):
    """HETATM lines of a posed ligand, with or without hydrogens."""
    Chem = pytest.importorskip("rdkit.Chem")
    from rdkit.Chem import AllChem

    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(mol, randomSeed=11)
    if not with_h:
        mol = Chem.RemoveHs(mol)
    conformer = mol.GetConformer()
    lines = []
    for i, atom in enumerate(mol.GetAtoms()):
        pos = conformer.GetAtomPosition(i)
        name = names[i] if names else f"{atom.GetSymbol().upper()}{i + 1}"
        lines.append(atom_line("HETATM", 100 + i, name, code, chain, resseq, (pos.x, pos.y, pos.z),
                               atom.GetSymbol().upper()))
    return lines


def test_the_binder_loses_its_hydrogens_and_everything_else_on_other_chains_is_dropped(pipe, tmp_path):
    path = write_complex(tmp_path / "c.pdb", embedded_ligand("CC(=O)O"),
                         extra=[atom_line("ATOM", 300, "CA", "GLY", "C", 1, (5.0, 5.0, 5.0), "C")])
    binder, ligand = pipe.split_complex(path, "A", "B", "LIG")
    assert [line[12:16].strip() for line in binder] == ["N", "CA"]
    assert {line[17:20] for line in ligand} == {"LIG"} and len(ligand) == 4


def test_a_hydrogen_without_an_element_column_is_still_a_hydrogen(pipe):
    line = "ATOM      2  HA  ALA A   1       1.000   0.000   0.000  1.00 20.00\n"
    assert pipe.element_of(line) == "H"
    assert pipe.element_of("ATOM      2 1HB  ALA A   1       1.000   0.000   0.000\n") == "H"


@pytest.mark.parametrize("ligand,message", [
    ([], "no residue 'LIG' on chain 'B'"),
    (None, "two residues|2 residues"),
])
def test_a_missing_or_ambiguous_ligand_is_a_clear_failure(pipe, tmp_path, ligand, message):
    if ligand is None:
        ligand = embedded_ligand("CC", resseq=1) + embedded_ligand("CC", resseq=2)
    path = write_complex(tmp_path / "c.pdb", ligand)
    with pytest.raises(ValueError, match=message):
        pipe.split_complex(path, "A", "B", "LIG")


def test_a_missing_binder_chain_is_a_clear_failure(pipe, tmp_path):
    path = write_complex(tmp_path / "c.pdb", embedded_ligand("CC"))
    with pytest.raises(ValueError, match="no atoms on binder chain 'Z'"):
        pipe.split_complex(path, "Z", "B", "LIG")


def test_a_hetatm_record_is_aligned_like_a_pdb_atom_line(pipe):
    line = pipe.hetatm_line(5, "CL1", "Cl", "LIG", "B", 7, (1.5, -2.25, 3.0))
    assert len(line.rstrip("\n")) == 78
    assert (line[:6], line[6:11].strip(), line[12:16], line[17:20], line[21], line[22:26].strip()) == \
        ("HETATM", "5", "CL1 ", "LIG", "B", "7")
    assert (float(line[30:38]), float(line[38:46]), float(line[46:54])) == (1.5, -2.25, 3.0)
    assert line[76:78] == "CL"
    assert pipe.hetatm_line(1, "C1", "C", "LIG", "B", 1, (0, 0, 0))[12:16] == " C1 "


# -- atom naming and the SDF molfile_to_params reads ------------------------------------------------------

def rosetta_element(name):
    """The element Rosetta's mdl_molfile reader infers from an atom name (rosetta_py.io.mdl_molfile.Atom)."""
    if name[0].isalpha():
        element = name[0:2] if len(name) >= 2 and name[1].isalpha() else name[0]
    else:
        element = name[1]
    return element.upper()


def test_names_are_unique_short_and_read_back_as_the_right_element(pipe):
    symbols = ["C", "C", "N", "Cl", "Br", "Si", "O", "H", "H", "S", "P", "F", "I", "Zn"] + ["C"] * 12
    names = pipe.element_first_names(symbols)
    assert len(set(names)) == len(names) and all(len(n) <= 4 for n in names)
    assert [rosetta_element(n) for n in names] == [s.upper() for s in symbols]


def test_a_name_longer_than_four_characters_is_refused(pipe):
    with pytest.raises(ValueError, match="at most 4 characters"):
        pipe.element_first_names(["Cl"] * 100)


def test_the_sdf_carries_the_names_where_rosettas_molfile_reader_looks(pipe):
    Chem = pytest.importorskip("rdkit.Chem")

    mol = Chem.AddHs(Chem.MolFromSmiles("ClCC(=O)O"))
    from rdkit.Chem import AllChem
    AllChem.EmbedMolecule(mol, randomSeed=3)
    names = pipe.element_first_names([a.GetSymbol() for a in mol.GetAtoms()])
    sdf = pipe.sdf_with_atom_names(mol, names)

    lines = sdf.splitlines()
    assert lines[3][34:39] == "V2000" and lines[-1] == "$$$$"
    n_atoms = int(lines[3][0:3])
    assert [lines[4 + i][30:35].strip() for i in range(n_atoms)] == names
    reread = Chem.MolFromMolBlock("\n".join(lines[:-1]), removeHs=False)
    assert reread is not None and reread.GetNumAtoms() == n_atoms
    assert [float(lines[4][c:c + 10]) for c in (0, 10, 20)] == pytest.approx(
        list(mol.GetConformer().GetAtomPosition(0)), abs=1e-3)


def test_a_name_count_that_does_not_match_the_atoms_is_refused(pipe):
    Chem = pytest.importorskip("rdkit.Chem")
    from rdkit.Chem import AllChem

    mol = Chem.AddHs(Chem.MolFromSmiles("CC"))
    AllChem.EmbedMolecule(mol, randomSeed=3)
    with pytest.raises(ValueError, match="names for"):
        pipe.sdf_with_atom_names(mol, ["C1"])


# -- ligand preparation ---------------------------------------------------------------------------------

@pytest.fixture
def dimorphite(monkeypatch):
    """A recorder in place of Dimorphite-DL, which is installed in the Rosetta environment only."""
    calls = []
    states = {"CC(=O)O": ["CC(=O)[O-]"]}
    module = types.ModuleType("dimorphite_dl")

    def protonate_smiles(smiles, **kwargs):
        calls.append((smiles, kwargs))
        return states.get(smiles, [smiles])

    module.protonate_smiles = protonate_smiles
    module.states = states
    module.calls = calls
    monkeypatch.setitem(sys.modules, "dimorphite_dl", module)
    return module


def test_a_hydrogen_free_ligand_is_protonated_at_the_requested_ph_and_given_hydrogens(pipe, tmp_path, dimorphite):
    pytest.importorskip("rdkit")
    lines = embedded_ligand("CC(=O)O")
    sdf, hetatm = pipe.prepare_ligand(lines, "LIG", "B", "CC(=O)O", 6.5, str(tmp_path))

    assert dimorphite.calls == [("CC(=O)O", {"ph_min": 6.5, "ph_max": 6.5, "precision": 0.5})]
    text = open(sdf).read()
    assert "M  CHG  1" in text
    names = [line[12:16].strip() for line in hetatm]
    assert names == ["C1", "C2", "O1", "O2", "H1", "H2", "H3"]
    assert {line[17:20] for line in hetatm} == {"LIG"} and {line[21] for line in hetatm} == {"B"}
    assert {int(line[22:26]) for line in hetatm} == {7}
    assert os.path.basename(sdf) == "LIG.sdf"


def test_a_ligand_that_already_has_hydrogens_is_not_protonated_again(pipe, tmp_path, dimorphite):
    pytest.importorskip("rdkit")
    lines = embedded_ligand("CC(=O)O", with_h=True)
    sdf, hetatm = pipe.prepare_ligand(lines, "LIG", "B", "CC(=O)O", 7.4, str(tmp_path))

    assert dimorphite.calls == []
    assert len(hetatm) == 8 and "M  CHG" not in open(sdf).read()


def test_hydrogens_follow_every_heavy_atom_however_the_structure_orders_them(pipe, tmp_path, dimorphite):
    """molfile_to_params types a hydrogen from the type its carbon already got, so a hydrogen listed first would be typed blind."""
    pytest.importorskip("rdkit")
    lines = list(reversed(embedded_ligand("CC(=O)O", with_h=True)))
    assert lines[0][76:78].strip() == "H"
    sdf, hetatm = pipe.prepare_ligand(lines, "LIG", "B", "CC(=O)O", 7.4, str(tmp_path))

    elements = [line[76:78].strip() for line in hetatm]
    assert elements == sorted(elements, key=lambda e: e == "H")
    assert elements.count("H") == 4
    names = [line[30:35].strip() for line in open(sdf).read().splitlines()[4:4 + len(hetatm)]]
    assert names == [line[12:16].strip() for line in hetatm]


def test_heavy_atom_coordinates_come_from_the_structure_not_from_rdkit(pipe, tmp_path, dimorphite):
    pytest.importorskip("rdkit")
    lines = embedded_ligand("CC(=O)O")
    _, hetatm = pipe.prepare_ligand(lines, "LIG", "B", "CC(=O)O", 7.4, str(tmp_path))
    original = [(float(l[30:38]), float(l[38:46]), float(l[46:54])) for l in lines]
    written = [(float(l[30:38]), float(l[38:46]), float(l[46:54])) for l in hetatm[:4]]
    assert written == pytest.approx(original, abs=1e-3)


def test_no_protonation_state_is_an_error_not_a_silent_fallback(pipe, tmp_path, dimorphite):
    pytest.importorskip("rdkit")
    dimorphite.states["CC(=O)O"] = []
    with pytest.raises(ValueError, match="no protonation state at pH 7.4"):
        pipe.prepare_ligand(embedded_ligand("CC(=O)O"), "LIG", "B", "CC(=O)O", 7.4, str(tmp_path))


def test_a_smiles_that_does_not_match_the_ligand_is_an_error(pipe, tmp_path, dimorphite):
    pytest.importorskip("rdkit")
    with pytest.raises(Exception, match="Check that the SMILES matches|No matching"):
        pipe.prepare_ligand(embedded_ligand("CC(=O)O"), "LIG", "B", "c1ccccc1", 7.4, str(tmp_path))


# -- molfile_to_params ------------------------------------------------------------------------------------

FAKE_MOLFILE_TO_PARAMS = """import json, sys
args = sys.argv[1:]
prefix = args[args.index("-p") + 1]
json.dump(args, open(prefix + ".argv.json", "w"))
if "--fail" in open(args[args.index("-p") + 2]).read():
    sys.stderr.write("atom typing failed for SI1")
    sys.exit(3)
if "--silent" not in open(args[args.index("-p") + 2]).read():
    open(prefix + ".params", "w").write("NAME " + args[args.index("--name") + 1] + chr(10))
"""


@pytest.fixture
def fake_script(pipe, tmp_path, monkeypatch):
    script = tmp_path / "molfile_to_params.py"
    script.write_text(FAKE_MOLFILE_TO_PARAMS)
    monkeypatch.setattr(pipe, "molfile_to_params_path", lambda: str(script))
    return script


def test_molfile_to_params_is_called_with_the_residue_name_the_prefix_and_the_kept_atom_names(pipe, tmp_path, fake_script):
    sdf = tmp_path / "LIG.sdf"
    sdf.write_text("sdf\n")
    params = pipe.generate_params(str(sdf), "LIG", str(tmp_path))

    assert params == str(tmp_path / "LIG.params")
    argv = json.load(open(tmp_path / "LIG.argv.json"))
    assert argv[argv.index("--name") + 1] == "LIG"
    assert argv[argv.index("-p") + 1] == str(tmp_path / "LIG")
    assert argv[argv.index("-p") + 2] == str(sdf)
    assert {"--no-pdb", "--clobber", "--keep-names"} <= set(argv)


def test_a_failing_molfile_to_params_surfaces_its_message(pipe, tmp_path, fake_script):
    sdf = tmp_path / "LIG.sdf"
    sdf.write_text("--fail\n")
    with pytest.raises(RuntimeError, match=r"exit 3\): .*atom typing failed for SI1"):
        pipe.generate_params(str(sdf), "LIG", str(tmp_path))


def test_a_run_that_writes_no_params_is_an_error(pipe, tmp_path, fake_script):
    sdf = tmp_path / "LIG.sdf"
    sdf.write_text("--silent\n")
    with pytest.raises(RuntimeError, match="wrote no LIG.params"):
        pipe.generate_params(str(sdf), "LIG", str(tmp_path))


def test_a_missing_script_tells_the_user_how_to_get_it(pipe, tmp_path, monkeypatch):
    monkeypatch.setattr(pipe, "molfile_to_params_path", lambda: str(tmp_path / "absent.py"))
    with pytest.raises(FileNotFoundError, match=r"LigandDdG\.install\(\).*params_file"):
        pipe.generate_params(str(tmp_path / "LIG.sdf"), "LIG", str(tmp_path))


def test_the_install_location_and_the_lookup_location_agree():
    from biopipelines.rosetta_scoring import MOLFILE_TO_PARAMS, molfile_to_params_path

    assert molfile_to_params_path() == os.path.join(sys.prefix, *MOLFILE_TO_PARAMS.split("/"))
    assert MOLFILE_TO_PARAMS.endswith("share/rosetta_public/source/scripts/python/public/molfile_to_params.py")


# -- scoring orchestration (Rosetta replaced by recorders) -------------------------------------------------

def residue_with(names, hydrogens=()):
    """A stand-in for a Rosetta residue: atom names padded the way Rosetta pads them, numbered from 1."""
    names = list(names) + list(hydrogens)
    return types.SimpleNamespace(natoms=lambda: len(names), atom_name=lambda i: f" {names[i - 1]:<3}",
                                 atom_is_hydrogen=lambda i: names[i - 1] in hydrogens)


@pytest.fixture
def recorders(pipe, monkeypatch):
    log = []

    def load_pose(binder, ligand, params):
        log.append(("pose", len(binder), ligand, params))
        names = [line[12:16].strip() for line in ligand]
        return types.SimpleNamespace(total_residue=lambda: 42, residue=lambda index: residue_with(names))

    monkeypatch.setattr(pipe, "prepare_ligand", lambda *a: (log.append(("prepare", a[1:5])) or ("x.sdf", ["HETATM\n"])))
    monkeypatch.setattr(pipe, "generate_params", lambda sdf, code, work: log.append(("params", sdf, code)) or "gen.params")
    monkeypatch.setattr(pipe, "load_pose", load_pose)
    monkeypatch.setattr(pipe, "seed_rosetta", lambda seed: None)
    monkeypatch.setattr(pipe, "soft_ddg", lambda *a: log.append(("soft", a[1:])) or {"ddg_soft": -4.0, "ddg_soft_std": 0.1,
                                                                                    "n_replicates_used": 3})
    monkeypatch.setattr(pipe, "hard_ddg", lambda p, i: log.append(("hard", i)) or {"ddg_norepack": -6.0})
    return log


def config(**overrides):
    base = {"scoring": "soft", "binder_chain": "A", "ligand_chain": "B", "replicates": 5, "ph": 7.4, "params_file": "",
            "seed": 42}
    base.update(overrides)
    return base


def test_the_soft_scheme_parameterizes_from_the_smiles_and_scores_with_the_chains_and_replicates(
        pipe, tmp_path, recorders):
    path = write_complex(tmp_path / "c.pdb", embedded_ligand("CC(=O)O"))
    result = pipe.score_structure("s1", path, config(replicates=4), "LIG", "CC(=O)O", str(tmp_path / "work"))

    assert result == {"ddg_soft": -4.0, "ddg_soft_std": 0.1, "n_replicates_used": 3}
    assert recorders[0] == ("prepare", ("LIG", "B", "CC(=O)O", 7.4))
    assert recorders[1] == ("params", "x.sdf", "LIG")
    assert recorders[2][3] == "gen.params" and recorders[2][2] == ["HETATM\n"]
    from biopipelines.rosetta_scoring import replicate_seed
    assert recorders[3] == ("soft", ("A", "B", 4, [replicate_seed(42, "s1", r) for r in range(4)]))
    assert os.path.isdir(tmp_path / "work" / "s1")


class PackablePose:
    """A stand-in for a PyRosetta pose (PyRosetta is not installed here): it counts the repacks applied to it, and a clone keeps the count."""

    def __init__(self, repacks=0):
        self.repacks = repacks

    def clone(self):
        return PackablePose(self.repacks)


def test_each_soft_replicate_is_seeded_before_its_repack_and_the_hard_scheme_once(pipe, monkeypatch):
    """The replicate seeds are distinct and fixed, so a rerun repeats exactly; the repack without them drew from /dev/urandom."""
    from biopipelines.rosetta_scoring import replicate_seed

    events = []
    fake = types.SimpleNamespace(create_score_function=lambda name: name)
    monkeypatch.setitem(sys.modules, "pyrosetta", fake)
    monkeypatch.setattr(pipe, "seed_rosetta", lambda seed: events.append(("seed", seed)))
    monkeypatch.setattr(pipe, "repack_around_chain", lambda *a: events.append(("repack",)))
    monkeypatch.setattr(pipe, "unbound_ddg", lambda *a: -20.0 - len(events))
    seeds = [replicate_seed(9, "s1", r) for r in range(3)]
    pipe.soft_ddg(PackablePose(), "A", "B", 3, seeds)
    assert [e for e in events if e[0] == "seed"] == [("seed", s) for s in seeds]
    assert [e[0] for e in events] == ["seed", "repack", "seed", "repack", "seed", "repack"]
    assert len(set(seeds)) == 3


def test_every_soft_replicate_repacks_a_clone_of_the_original_pose(pipe, monkeypatch):
    """Repacking in place made replicate k start from what replicate k-1 left, a chain whose spread says nothing about sampling."""
    started_from, scored = [], []
    monkeypatch.setitem(sys.modules, "pyrosetta", types.SimpleNamespace(create_score_function=lambda name: name))
    monkeypatch.setattr(pipe, "seed_rosetta", lambda seed: None)

    def repack(pose, chain, scorefxn):
        started_from.append(pose.repacks)
        pose.repacks += 1

    def unbound(pose, scorefxn, chain):
        scored.append(pose.repacks)
        return -10.0 - pose.repacks

    monkeypatch.setattr(pipe, "repack_around_chain", repack)
    monkeypatch.setattr(pipe, "unbound_ddg", unbound)
    original = PackablePose()
    result = pipe.soft_ddg(original, "A", "B", 4, [1, 2, 3, 4])

    assert started_from == [0, 0, 0, 0]
    assert scored == [1, 1, 1, 1]
    assert original.repacks == 0
    assert result == {"ddg_soft": -11.0, "ddg_soft_std": 0.0, "n_replicates_used": 2}


def test_the_soft_value_is_the_trimmed_mean_with_a_population_deviation_of_the_replicates(pipe, monkeypatch):
    """The best and the worst replicate are dropped from three on, and the spread is the population deviation (ddof = 0) of what is kept."""
    values = iter([-9.0, -3.0, -1.0, 40.0])
    monkeypatch.setitem(sys.modules, "pyrosetta", types.SimpleNamespace(create_score_function=lambda name: name))
    monkeypatch.setattr(pipe, "seed_rosetta", lambda seed: None)
    monkeypatch.setattr(pipe, "repack_around_chain", lambda *a: None)
    monkeypatch.setattr(pipe, "unbound_ddg", lambda *a: next(values))
    assert pipe.soft_ddg(PackablePose(), "A", "B", 4, [1, 2, 3, 4]) == {
        "ddg_soft": -2.0, "ddg_soft_std": 1.0, "n_replicates_used": 2}


def test_the_hard_scheme_seeds_rosetta_once_from_the_structure_id(pipe, tmp_path, recorders, monkeypatch):
    from biopipelines.rosetta_scoring import replicate_seed

    seeded = []
    monkeypatch.setattr(pipe, "seed_rosetta", lambda seed: seeded.append(seed))
    path = write_complex(tmp_path / "c.pdb", embedded_ligand("CC(=O)O"))
    pipe.score_structure("s1", path, config(scoring="hard", replicates=None, seed=5), "LIG", "CC(=O)O", str(tmp_path / "w"))
    assert seeded == [replicate_seed(5, "s1", 0)]


def test_the_hard_scheme_reads_the_ligand_position_from_the_pose(pipe, tmp_path, recorders):
    path = write_complex(tmp_path / "c.pdb", embedded_ligand("CC(=O)O"))
    result = pipe.score_structure("s1", path, config(scoring="hard", replicates=None), "LIG", "CC(=O)O",
                                  str(tmp_path / "work"))
    assert result == {"ddg_norepack": -6.0}
    assert recorders[-1] == ("hard", 42)


def test_a_premade_params_file_bypasses_the_chemistry_and_keeps_the_ligand_as_written(pipe, tmp_path, recorders):
    params = tmp_path / "LIG.params"
    params.write_text("NAME LIG\n")
    lines = embedded_ligand("CC(=O)O")
    path = write_complex(tmp_path / "c.pdb", lines)
    pipe.score_structure("s1", path, config(params_file=str(params)), "LIG", None, str(tmp_path / "work"))

    assert [entry[0] for entry in recorders] == ["pose", "soft"]
    assert recorders[0][3] == str(params) and recorders[0][2] == lines


def test_a_ligand_with_other_atom_names_than_the_params_is_refused_not_scored(pipe):
    """Rosetta drops an atom it has no name for without a word, and the molecule that is scored is not the one given."""
    lines = [atom_line("HETATM", 1, "C1", "LIG", "B", 7, (1.0, 0.0, 0.0), "C"),
             atom_line("HETATM", 2, "O1", "LIG", "B", 7, (2.0, 0.0, 0.0), "O")]
    with pytest.raises(ValueError, match=r"does not match LIG.params.*\['O1'\].*\['O2'\]"):
        pipe.check_ligand_names(residue_with(["C1", "O2"], hydrogens=["H1"]), lines, "/x/LIG.params")


def test_a_ligand_without_hydrogens_matches_params_that_carry_them(pipe):
    lines = [atom_line("HETATM", 1, "C1", "LIG", "B", 7, (1.0, 0.0, 0.0), "C"),
             atom_line("HETATM", 2, "O2", "LIG", "B", 7, (2.0, 0.0, 0.0), "O")]
    pipe.check_ligand_names(residue_with(["C1", "O2"], hydrogens=["H1"]), lines, "/x/LIG.params")


def test_a_missing_params_file_is_a_failure_of_that_id(pipe, tmp_path, recorders):
    path = write_complex(tmp_path / "c.pdb", embedded_ligand("CC"))
    with pytest.raises(FileNotFoundError, match="params_file not found"):
        pipe.score_structure("s1", path, config(params_file=str(tmp_path / "nope.params")), "LIG", None, str(tmp_path))


@pytest.mark.parametrize("x", [200.0, 8.0])
def test_a_ligand_out_of_range_of_the_binder_chain_is_a_failure_not_a_zero(pipe, tmp_path, recorders, x):
    """The ligand of another copy sits 8 A away in a crystal: no pair of atoms interacts, so the ddG is zero."""
    path = write_complex(tmp_path / "c.pdb", [atom_line("HETATM", 100, "C1", "LIG", "B", 7, (x, 1.4, 0.0), "C")])
    assert pipe.nearest_distance(*pipe.split_complex(path, "A", "B", "LIG")) == pytest.approx(x)
    with pytest.raises(ValueError, match="no interface to score"):
        pipe.score_structure("s1", path, config(), "LIG", "CC", str(tmp_path / "work"))
    assert recorders == []


def test_a_ligand_within_the_contact_range_is_scored(pipe, tmp_path, recorders):
    path = write_complex(tmp_path / "c.pdb", [atom_line("HETATM", 100, "C1", "LIG", "B", 7, (4.0, 1.4, 0.0), "C")])
    assert pipe.score_structure("s1", path, config(), "LIG", "CC", str(tmp_path / "work")) is not None


def test_without_a_smiles_or_a_params_file_the_ligand_cannot_be_parameterized(pipe, tmp_path, recorders):
    path = write_complex(tmp_path / "c.pdb", embedded_ligand("CC"))
    with pytest.raises(ValueError, match="no SMILES"):
        pipe.score_structure("s1", path, config(), "LIG", None, str(tmp_path))
    assert recorders == []


def write_compounds(tmp_path, code="LIG", smiles="CC(=O)O"):
    import pandas as pd

    map_csv = tmp_path / "compounds_map.csv"
    pd.DataFrame([{"id": "lig", "code": code, "smiles": smiles}]).to_csv(map_csv, index=False)
    path = tmp_path / "compounds.json"
    json.dump({"name": "compounds", "ids": ["lig"], "files": [], "map_table": str(map_csv), "format": "csv"},
              open(path, "w"))
    return str(path)


def write_structures(tmp_path, ids):
    import pandas as pd

    files = [str(tmp_path / f"{i}.pdb") for i in ids]
    for f in files:
        open(f, "w").write("END\n")
    map_csv = tmp_path / "structures_map.csv"
    pd.DataFrame({"id": ids, "file": files}).to_csv(map_csv, index=False)
    path = tmp_path / "structures.json"
    json.dump({"name": "structures", "ids": ids, "files": files, "map_table": str(map_csv), "format": "pdb"},
              open(path, "w"))
    return str(path)


@pytest.mark.parametrize("scoring,columns,result", [
    ("hard", ["id", "ddg_norepack"], {"ddg_norepack": -6.0}),
    ("soft", ["id", "ddg_soft", "ddg_soft_std", "n_replicates_used"],
     {"ddg_soft": -4.0, "ddg_soft_std": 0.1, "n_replicates_used": 3}),
])
def test_main_writes_one_row_per_id_in_the_scheme_columns_and_records_the_failure(
        pipe, tmp_path, monkeypatch, scoring, columns, result):
    import pandas as pd
    import yaml

    init = []
    seen = []
    cfg = tmp_path / "config.yaml"
    yaml.safe_dump(config(scoring=scoring, replicates=None if scoring == "hard" else 5), open(cfg, "w"))

    def fake(sid, path, conf, code, smiles, work):
        seen.append((sid, code, smiles, work))
        if sid == "b":
            raise RuntimeError("molfile_to_params failed for LIG (exit 3): atom typing failed")
        return result

    monkeypatch.setattr(pipe, "score_structure", fake)
    monkeypatch.setattr(pipe, "init_pyrosetta", lambda options: init.append(tuple(options)))
    table = tmp_path / "tables" / "ligand_ddg.csv"
    local = tmp_path / "_execution" / "local_missing.csv"
    pipe.main(["--config-yaml", str(cfg), "--structures-json", write_structures(tmp_path, ["a", "b", "c"]),
               "--ligand-json", write_compounds(tmp_path), "--work-dir", str(tmp_path / "work"),
               "--table-csv", str(table), "--local-missing-csv", str(local)])

    from biopipelines.rosetta_scoring import HARD_INIT_OPTIONS, RELAX_INIT_OPTIONS
    assert init == [tuple(HARD_INIT_OPTIONS if scoring == "hard" else RELAX_INIT_OPTIONS)]
    assert [s[:3] for s in seen] == [("a", "LIG", "CC(=O)O"), ("b", "LIG", "CC(=O)O"), ("c", "LIG", "CC(=O)O")]
    frame = pd.read_csv(table)
    assert list(frame.columns) == columns and list(frame["id"]) == ["a", "b", "c"]
    assert frame.loc[frame["id"] == "b", columns[1:]].isna().all().all()
    assert frame.loc[frame["id"] == "c", columns[1]].item() == result[columns[1]]
    failed = pd.read_csv(local)
    assert list(failed["id"]) == ["b"] and "atom typing failed" in failed["cause"].item()


def test_main_does_not_need_a_smiles_when_a_params_file_is_given(pipe, tmp_path, monkeypatch):
    import yaml

    cfg = tmp_path / "config.yaml"
    yaml.safe_dump(config(params_file="/data/LIG.params"), open(cfg, "w"))
    seen = []
    monkeypatch.setattr(pipe, "score_structure", lambda sid, path, conf, code, smiles, work:
                        seen.append(smiles) or {"ddg_soft": 1.0, "ddg_soft_std": 0.0, "n_replicates_used": 5})
    monkeypatch.setattr(pipe, "init_pyrosetta", lambda options: None)
    pipe.main(["--config-yaml", str(cfg), "--structures-json", write_structures(tmp_path, ["a"]),
               "--ligand-json", write_compounds(tmp_path, smiles=""), "--work-dir", str(tmp_path / "w"),
               "--table-csv", str(tmp_path / "t.csv"), "--local-missing-csv", str(tmp_path / "m.csv")])
    assert seen == [None]


def test_main_exits_non_zero_when_every_id_raised(pipe, tmp_path, monkeypatch):
    import yaml

    cfg = tmp_path / "config.yaml"
    yaml.safe_dump(config(), open(cfg, "w"))
    monkeypatch.setattr(pipe, "score_structure", lambda *a: (_ for _ in ()).throw(RuntimeError("boom")))
    monkeypatch.setattr(pipe, "init_pyrosetta", lambda options: None)
    with pytest.raises(SystemExit) as stop:
        pipe.main(["--config-yaml", str(cfg), "--structures-json", write_structures(tmp_path, ["a", "b"]),
                   "--ligand-json", write_compounds(tmp_path), "--work-dir", str(tmp_path / "w"),
                   "--table-csv", str(tmp_path / "t.csv"), "--local-missing-csv", str(tmp_path / "m.csv")])
    assert stop.value.code == 1


# -- the hard protocol --------------------------------------------------------------------------------------

def test_the_hard_protocol_is_well_formed_and_every_reference_resolves(pipe):
    xml = pipe.HARD_PROTOCOL_XML
    assert xml.count("%%ligand_res_number%%") == 1
    root = ET.fromstring(xml.replace("%%ligand_res_number%%", "42"))
    sections = {section.tag: {child.get("name") for child in section if child.get("name")} for section in root}
    scorefxns, selectors = sections["SCOREFXNS"], sections["RESIDUE_SELECTORS"]
    task_ops, movers, filters = sections["TASKOPERATIONS"], sections["MOVERS"], sections["FILTERS"]

    assert set(pipe.HARD_PROTOCOL_MOVERS) <= movers and "ddg_norepack" in filters
    for node in root.iter():
        for attribute in ("scorefxn",):
            if node.get(attribute):
                assert node.get(attribute) in scorefxns, (node.tag, node.get(attribute))
        for attribute in ("selector", "residue_selector"):
            if node.get(attribute):
                assert node.get(attribute) in selectors, (node.tag, node.get(attribute))
        if node.get("selectors"):
            assert set(node.get("selectors").split(",")) <= selectors
        if node.get("task_operations"):
            assert set(node.get("task_operations").split(",")) <= task_ops
        if node.get("relax_mover"):
            assert node.get("relax_mover") in movers
    ddg = next(root.iter("Ddg"))
    assert (ddg.get("repack"), ddg.get("jump"), ddg.get("repeats")) == ("0", "%%ligand_jump%%", "1")


def test_the_hard_protocol_has_the_protocols_section_the_rosetta_schema_requires(pipe):
    """Without it PyRosetta's parser rejects the script ("Missing child element(s)") and every structure fails."""
    root = ET.fromstring(pipe.HARD_PROTOCOL_XML.replace("%%ligand_res_number%%", "42"))
    assert root.find("PROTOCOLS") is not None


def test_the_hard_protocol_uses_the_jump_that_builds_the_ligand_not_jump_1(pipe, monkeypatch):
    """A hetero residue on the binder chain takes jump 1, so across it the Ddg filter would measure that residue."""
    seen = {}

    class Objects:
        def get_mover(self, name):
            return types.SimpleNamespace(apply=lambda pose: None)

        def get_filter(self, name):
            return types.SimpleNamespace(report_sm=lambda pose: -1.5)

    class XmlObjects:
        @staticmethod
        def create_from_string(xml):
            seen["xml"] = xml
            return Objects()

    scripts = types.ModuleType("pyrosetta.rosetta.protocols.rosetta_scripts")
    scripts.XmlObjects = XmlObjects
    for name in ("pyrosetta", "pyrosetta.rosetta", "pyrosetta.rosetta.protocols"):
        monkeypatch.setitem(sys.modules, name, types.ModuleType(name))
    monkeypatch.setitem(sys.modules, "pyrosetta.rosetta.protocols.rosetta_scripts", scripts)

    tree = types.SimpleNamespace(get_jump_that_builds_residue=lambda index: 2 if index == 42 else 1)
    assert pipe.hard_ddg(types.SimpleNamespace(fold_tree=lambda: tree), 42) == {"ddg_norepack": -1.5}
    root = ET.fromstring(seen["xml"])
    jumps = {node.tag: node.get("jump") for node in root.iter() if node.get("jump")}
    # the movers' jump is a boolean ("minimize the rigid-body degrees of freedom"), a number there is refused by the schema
    assert jumps == {"ProteinInterfaceDesign": "2", "TaskAwareMinMover": "1", "Ddg": "2"}
    assert '<Index name="ligand" resnums="42"/>' in seen["xml"]


def test_the_ligand_selection_uses_the_pose_index_and_the_hard_scheme_passes_it(pipe, tmp_path, recorders):
    """Index takes pose numbering for a plain integer; the ligand is the last residue by construction."""
    assert '<Index name="ligand" resnums="%%ligand_res_number%%"/>' in pipe.HARD_PROTOCOL_XML
    path = write_complex(tmp_path / "c.pdb", embedded_ligand("CC(=O)O"))
    pipe.score_structure("s1", path, config(scoring="hard", replicates=None), "LIG", "CC(=O)O", str(tmp_path / "w"))
    assert recorders[-1] == ("hard", 42)


# -- PyRosetta-dependent ------------------------------------------------------------------------------------

def test_the_hard_protocol_parses_in_rosetta_and_exposes_its_movers_and_filter(pipe):
    pyr = pytest.importorskip("pyrosetta", reason="pyrosetta is not installed in this environment (it needs its own licensed one): Rosetta does not parse the hard protocol")
    from biopipelines.rosetta_scoring import HARD_INIT_OPTIONS, init_pyrosetta

    init_pyrosetta(HARD_INIT_OPTIONS)
    objects = pyr.rosetta.protocols.rosetta_scripts.XmlObjects.create_from_string(
        pipe.HARD_PROTOCOL_XML.replace("%%ligand_res_number%%", "2"))
    for name in pipe.HARD_PROTOCOL_MOVERS:
        assert objects.get_mover(name) is not None
    assert hasattr(objects.get_filter("ddg_norepack"), "report_sm")
