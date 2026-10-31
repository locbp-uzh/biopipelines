"""Helpers that more than one test module needs; a test module never imports another test module.

Fixtures live in conftest.py (`pipe_gnina`), plain functions and constants here.
"""

import contextlib
import importlib.util
import json
import os
import pathlib
import shlex
import subprocess
import sys

import numpy as np
import pandas as pd

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPTS = os.path.join(REPO_ROOT, "pipe_scripts")


# ── pipe scripts ─────────────────────────────────────────────────────────────

@contextlib.contextmanager
def sys_path_prepended(path):
    """`path` first on sys.path for the body only: the module-level imports of a test file, without a leak."""
    sys.path.insert(0, str(path))
    try:
        yield
    finally:
        sys.path.remove(str(path))


def load_pipe_script(name, filename):
    """Load a pipe script by path; pipe_scripts is on sys.path while it executes, and only then."""
    with sys_path_prepended(SCRIPTS):
        spec = importlib.util.spec_from_file_location(name, os.path.join(SCRIPTS, filename))
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    return module


def run_pipe_script(script, *args):
    return subprocess.run([sys.executable, os.path.join(SCRIPTS, script), *args], capture_output=True, text=True)


def pipe_script_text(name):
    with open(os.path.join(SCRIPTS, name)) as handle:
        return handle.read()


def propagate_missing(script, upstream_rows):
    """What a generated step script does after its pipe: merge the upstream manifests and the local file into tables/missing.csv."""
    script = script.replace("\r\n", "\n")
    block = script[script.rindex("python", 0, script.index("pipe_propagate_missing.py")):].split("\n\n")[0]
    command = shlex.split(block.replace("\\\n", " "))
    arguments = {}
    for token in command[2:]:
        if token.startswith("--"):
            arguments[token] = []
        else:
            arguments[list(arguments)[-1]].append(token)
    for folder in arguments["--upstream-folders"]:
        if folder:
            pathlib.Path(folder).mkdir(parents=True, exist_ok=True)
            pathlib.Path(folder, "missing.csv").write_text(upstream_rows)
    propagate = load_pipe_script("pipe_propagate_missing_under_test", "pipe_propagate_missing.py")
    propagate.propagate_missing_table(
        upstream_folders=arguments["--upstream-folders"], output_folder=arguments["--output-folder"][0],
        missing_csv_path=arguments["--missing-csv"][0],
        local_missing=(arguments.get("--local-missing") or [None])[0])
    return pd.read_csv(arguments["--missing-csv"][0], dtype=str, keep_default_na=False)


# ── PLACER and the other pair tools ──────────────────────────────────────────

# Ids as upstream tools make them: a protein that is itself a pair (`a+b`), a ligand with `+` and one with `_`.
PROTEINS = ["a+b", "c"]
LIGANDS = ["l+1", "m_2"]
PAIRS = {
    "a+b+l+1": ("a+b", "l+1"), "a+b+m_2": ("a+b", "m_2"), "c+l+1": ("c", "l+1"), "c+m_2": ("c", "m_2"),
}


# ── Boltz2 post-processing ───────────────────────────────────────────────────

BOLTZ_POST = os.path.join(SCRIPTS, "pipe_boltz_postprocessing.py")


def boltz_prediction(root, config_id, samples=1, pae_for=None):
    """boltz_results_<id>/predictions/<id>/ with `samples` models; model k carries a PAE matrix filled with k."""
    folder = root / "execution" / f"boltz_results_{config_id}" / "predictions" / config_id
    folder.mkdir(parents=True)
    for k in range(samples):
        (folder / f"{config_id}_model_{k}.pdb").write_text(f"MODEL {k}\nEND\n")
        scores = {"confidence_score": 0.9 - 0.1 * k, "ptm": 0.8, "iptm": 0.7,
                  "complex_plddt": 0.85, "complex_iplddt": 0.8}
        (folder / f"confidence_{config_id}_model_{k}.json").write_text(json.dumps(scores))
        if pae_for is None or k in pae_for:
            np.savez_compressed(folder / f"pae_{config_id}_model_{k}.npz", pae=np.full((6, 6), float(k)))


# ── Gnina and Vina ───────────────────────────────────────────────────────────

PROTEIN_LINE = "ATOM      1  CA  GLY A   1       0.000   0.000   0.000  1.00  0.00           C\n"
LIGAND_PDB = (
    "HETATM    1  C1  LIG A 100       0.000   0.000   0.000  1.00  0.00           C\n"
    "HETATM    2  C2  LIG A 100       1.500   0.000   0.000  1.00  0.00           C\n"
    "HETATM    3  O1  LIG A 100       2.200   1.200   0.000  1.00  0.00           O\n")


def write_pose_sdf(path, poses):
    """One ethanol pose per (vina, cnn_score, cnn_affinity), carrying the score properties gnina writes."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
    AllChem.EmbedMolecule(mol, randomSeed=7)
    writer = Chem.SDWriter(str(path))
    for vina, cnn_score, cnn_affinity in poses:
        mol.SetProp("minimizedAffinity", str(vina))
        mol.SetProp("CNNscore", str(cnn_score))
        mol.SetProp("CNNaffinity", str(cnn_affinity))
        writer.write(mol)
    writer.close()


def declared_stream(new_pipeline, protein_ids, ligand_ids, cls="Gnina", **kwargs):
    """The structures stream a Gnina or Vina step declares at configuration time."""
    from biopipelines import gnina as gnina_module
    from biopipelines.mock import Mock

    pipeline = new_pipeline("gnina_ids")
    with pipeline:
        struct = Mock(ids=list(protein_ids), streams={"structures": {"format": "pdb", "file": "<id>.pdb"}},
                      map_table_strategy="config")
        comp = Mock(ids=list(ligand_ids), streams={"compounds": {"format": "csv", "file": "<id>.csv"}},
                    map_table_strategy="config")
        out = getattr(gnina_module, cls)(structures=struct.streams.structures,
                                         compounds=comp.streams.compounds, **kwargs)
    return out.streams.structures, struct.streams.structures, comp.streams.compounds


def configured_gnina_step(new_pipeline, cls_name, mode, upstream):
    """Save a pipeline holding one Gnina or Vina step and return (config path, config dict, generated script text, step id)."""
    from biopipelines import gnina as gnina_module
    from biopipelines.mock import Mock

    pipeline = new_pipeline("gnina_missing")
    with pipeline:
        # A table named `missing` is what marks a step as one that drops ids.
        tables = {"missing": {"columns": ["id", "removed_by", "kind", "cause"]}} if upstream else None
        structures = Mock(ids=["good", "bad"], tables=tables,
                          streams={"structures": {"format": "pdb", "file": "<id>.pdb"}}, map_table_strategy="config")
        compounds = Mock(ids=["lig"], streams={"compounds": {"format": "csv", "file": "<id>.csv"}},
                         map_table_strategy="config")
        box = {"center": "0,0,0", "size": 20, "protonate": False} if mode == "docking" else {}
        getattr(gnina_module, cls_name)(structures=structures, compounds=compounds.streams.compounds, mode=mode, **box)
        pipeline.save()
    root = pathlib.Path.cwd()
    (config_path,) = root.rglob("gnina_config.json")
    (script_path,) = root.rglob(f"*_{cls_name}.sh")
    with open(config_path) as handle:
        config = json.load(handle)
    return config_path, config, script_path.read_text(encoding="utf-8", errors="replace"), script_path.stem


def materialize_gnina_inputs(config, upstream):
    """Write the inputs the planned paths point to: a PDB per structure, the compounds with chemistry, the upstream rows."""
    with open(config["structures_json"]) as handle:
        structures = json.load(handle)
    template = structures["files"][0]
    for structure_id in ("good", "bad"):
        path = pathlib.Path(template.replace("<id>", structure_id))
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(PROTEIN_LINE + LIGAND_PDB + "END\n")
    with open(config["compounds_json"]) as handle:
        compounds_map = pathlib.Path(json.load(handle)["map_table"])
    compounds_map.parent.mkdir(parents=True, exist_ok=True)
    compounds_map.write_text("id,smiles,code\nlig,CCO,LIG\n")
    if upstream:
        (path,) = config["upstream_missing"]
        pathlib.Path(path).parent.mkdir(parents=True, exist_ok=True)
        pathlib.Path(path).write_text("id,removed_by,kind,cause\nup,001_Mock,filter,dropped upstream\n")


class GninaBinary:
    """gnina stand-in: exit 1 for the receptor named `bad`, otherwise one pose per call."""

    def __call__(self, cmd, **kwargs):
        receptor = pathlib.Path(cmd[cmd.index("-r") + 1]).stem
        if receptor.startswith("bad"):
            return subprocess.CompletedProcess(cmd, 1, "", "binary failed on this receptor")
        write_pose_sdf(cmd[cmd.index("-o") + 1], [(-7.0, 0.8, 5.0)])
        return subprocess.CompletedProcess(cmd, 0, "", "")


def fake_the_gnina_tools(pipe, monkeypatch, cls_name, mode):
    monkeypatch.setattr(pipe.subprocess, "run", GninaBinary())
    backend = pipe.vina_backend
    if cls_name == "Vina":
        def receptor(src, dst):
            if os.path.basename(src).startswith("bad"):
                raise RuntimeError("receptor has no atoms")
            return dst

        def pdbqt_to_sdf(src, dst, template_sdf=None):
            write_pose_sdf(dst, [(-7.0, 0.8, 5.0)])
            return dst

        monkeypatch.setattr(backend, "receptor_pdbqt", receptor)
        monkeypatch.setattr(backend, "ligand_pdbqt", lambda src, dst: dst)
        monkeypatch.setattr(backend, "box_args", lambda *a, **k: [])
        monkeypatch.setattr(backend, "run",
                            lambda *a, **k: subprocess.CompletedProcess(["vina"], 0, "", ""))
        monkeypatch.setattr(backend, "pdbqt_to_sdf", pdbqt_to_sdf)
    if mode != "docking":
        def extract(complex_pdb, code, out_pdb):
            if os.path.basename(complex_pdb).startswith("bad"):
                return None
            pathlib.Path(out_pdb).write_text(LIGAND_PDB)
            return out_pdb

        monkeypatch.setattr(pipe, "_extract_ligand_hetatm", extract)
        monkeypatch.setattr(pipe, "write_ligand_sdf",
                            lambda lig_pdb, dst, smiles=None: write_pose_sdf(dst, [(-7.0, 0.8, 5.0)]))
