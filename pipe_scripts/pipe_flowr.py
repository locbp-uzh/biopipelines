#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
FLOWR runtime: one upstream run per pocket, then the molecules split into a compounds stream and one SDF pose each.

Upstream's `python -m flowr.gen.generate_from_pdb` takes one protein file and either a reference ligand file (the pocket is cut around it) or a heavy-atom count, and writes `samples_<protein stem>.sdf` (plus `samples_<stem>_protonated.sdf` with `--protonate_generated_ligands`) into `--save_dir`. The molecules are in the coordinate frame of the protein file. It prints its filters' counts but writes no per-molecule score.
"""

import argparse
import math
import os
import subprocess
import sys
from typing import Dict, List, Optional, Tuple

import pandas as pd
import yaml

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from biopipelines.biopipelines_io import (  # noqa: E402
    load_datastream,
    iterate_files,
    step_id_from_table_path,
    container_argv_prefix,
)
from biopipelines.ligand_utils import resolve_ligand_code, resolve_ligand_smiles, write_ligand_sdf  # noqa: E402
from biopipelines.pdb_parser import field_chain, field_i_code, field_res_name, field_res_seq  # noqa: E402
from biopipelines.sele_utils import sele_to_chains_and_list  # noqa: E402

CONDITIONING_FLAGS = {
    "interaction": "--interaction_inpainting",
    "functional_group": "--func_group_inpainting",
    "scaffold": "--scaffold_inpainting",
    "linker": "--linker_inpainting",
    "substructure": "--substructure_inpainting",
}

COMPOUND_COLUMNS = ["id", "structures.id", "sample", "format", "code", "ccd", "smiles"]
LIGAND_CODE = "LIG"


def build_command(config: dict, pdb_file: str, ligand_file: Optional[str], save_dir: str,
                  checkpoint: str, prefix: List[str], num_workers: int) -> List[str]:
    """The upstream command line for one pocket."""
    cmd = prefix + [
        "python", "-m", "flowr.gen.generate_from_pdb",
        "--pdb_file", pdb_file,
        "--arch", "pocket", "--pocket_type", "holo",
        "--ckpt_path", checkpoint,
        "--save_dir", save_dir,
        "--gpus", "1",
        "--batch_cost", str(config["batch_cost"]),
        "--sample_n_molecules_per_target", str(config["num_molecules"]),
        "--integration_steps", str(config["integration_steps"]),
        "--seed", str(config["seed"]),
        "--filter_valid_unique",
    ]
    if ligand_file is not None:
        cmd += ["--ligand_file", ligand_file, "--cut_pocket", "--pocket_cutoff", str(config["pocket_cutoff"])]
    else:
        cmd += ["--num_heavy_atoms", str(config["num_heavy_atoms"])]
    if config["sample_sizes"]:
        cmd.append("--sample_mol_sizes")
    if config["protonate_pocket"]:
        cmd.append("--protonate_pocket")
    conditioning = config["conditioning"]
    if conditioning is not None:
        cmd.append(CONDITIONING_FLAGS[conditioning])
        if conditioning == "interaction":
            cmd.append("--compute_interactions")
        if conditioning == "substructure":
            cmd += ["--substructure", config["substructure"]]
    if config["protonate_ligands"]:
        cmd += ["--protonate_generated_ligands", "--num_workers", str(num_workers)]
    return cmd


def extract_ligand_block(pdb_path: str, code: str) -> str:
    """The HETATM records of the first residue named `code`, as a ligand-only PDB."""
    code = code.strip().upper()
    selected, first_key = [], None
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith("HETATM") or field_res_name(line).upper() != code:
                continue
            key = (field_chain(line), field_res_seq(line))
            if first_key is None:
                first_key = key
            if key == first_key:
                selected.append(line.rstrip("\n"))
    if not selected:
        raise ValueError(f"no HETATM records named {code!r} in {os.path.basename(pdb_path)}")
    return "\n".join(selected) + "\nEND\n"


def stage_reference_ligand(pdb_path: str, code: str, smiles: Optional[str], run_dir: str, sid: str) -> str:
    """The bound ligand as an SDF with bond orders, keeping its coordinates."""
    lig_pdb = os.path.join(run_dir, f"{sid}_ligand.pdb")
    lig_sdf = os.path.join(run_dir, f"{sid}_ligand.sdf")
    with open(lig_pdb, "w") as f:
        f.write(extract_ligand_block(pdb_path, code))
    write_ligand_sdf(lig_pdb, lig_sdf, smiles)
    return lig_sdf


def stage_protein_pdb(pdb_path: str, out_pdb: str) -> None:
    """The structure without its CONECT records, which upstream's reader cannot always resolve (3RY2: IndexError on an atom serial)."""
    with open(pdb_path) as src, open(out_pdb, "w") as dst:
        dst.writelines(line for line in src if not line.startswith("CONECT"))


def stage_pocket_pdb(pdb_path: str, selection: str, out_pdb: str) -> int:
    """A PDB holding only the ATOM records of the selected residues; returns how many residues."""
    whole_chains, residues = sele_to_chains_and_list(selection)
    wanted = set(residues)
    kept, seen = [], set()
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            chain, seq = field_chain(line), field_res_seq(line)
            try:
                key = (chain, int(seq))
            except ValueError:
                continue
            if chain in whole_chains or key in wanted:
                kept.append(line.rstrip("\n"))
                seen.add(key + (field_i_code(line),))
    if not kept:
        raise ValueError(f"no ATOM record of {os.path.basename(pdb_path)} lies in {selection!r}")
    with open(out_pdb, "w") as f:
        f.write("\n".join(kept) + "\nEND\n")
    return len(seen)


def output_sdf(save_dir: str, stem: str, protonated: bool) -> str:
    """Where upstream wrote the molecules: the protonated file when it was asked to protonate."""
    return os.path.join(save_dir, f"samples_{stem}_protonated.sdf" if protonated else f"samples_{stem}.sdf")


def read_records(sdf_path: str) -> list:
    """Every record of an SDF in file order, None where RDKit cannot read one."""
    from rdkit import Chem
    return list(Chem.SDMolSupplier(sdf_path, removeHs=False))


def smiles_of(mol) -> str:
    """Isomeric canonical SMILES of the heavy-atom molecule."""
    from rdkit import Chem
    return Chem.MolToSmiles(Chem.RemoveHs(mol))


def write_pose(mol, mol_id: str, path: str) -> None:
    from rdkit import Chem
    mol = Chem.Mol(mol)
    mol.SetProp("_Name", mol_id)
    writer = Chem.SDWriter(path)
    writer.write(mol)
    writer.close()


def collect_molecules(sid: str, sdf_path: str, num_molecules: int, structures_folder: str, step_id: str) -> Tuple[list, list, list]:
    """Split a pocket's SDF into compounds rows, pose-map rows and missing rows."""
    records = read_records(sdf_path)
    compounds, poses, missing = [], [], []
    for k, mol in enumerate(records, start=1):
        mol_id = f"{sid}_{k}"
        try:
            if mol is None:
                raise ValueError("RDKit cannot read this SDF record")
            smiles = smiles_of(mol)
            pose_path = os.path.join(structures_folder, f"{mol_id}.sdf")
            write_pose(mol, mol_id, pose_path)
        except Exception as e:
            missing.append({"id": mol_id, "removed_by": step_id, "kind": "failure", "cause": str(e)[:200]})
            continue
        compounds.append({"id": mol_id, "structures.id": sid, "sample": k, "format": "smiles",
                          "code": LIGAND_CODE, "ccd": "", "smiles": smiles})
        poses.append({"id": mol_id, "file": pose_path, "structures.id": sid})
    for k in range(len(records) + 1, num_molecules + 1):
        missing.append({"id": f"{sid}_{k}", "removed_by": step_id, "kind": "filter",
                        "cause": "not returned by FLOWR after its validity, uniqueness and diversity filters"})
    return compounds, poses, missing


def subprocess_runner(cmd: List[str], cwd: str, env: Dict[str, str], log_path: str) -> int:
    with open(log_path, "w") as log:
        return subprocess.run(cmd, cwd=cwd, env=env, stdout=log, stderr=subprocess.STDOUT).returncode


def last_line(log_path: str) -> str:
    """The last non-empty line of a log, which is the exception message of a Python traceback."""
    try:
        with open(log_path, errors="replace") as f:
            lines = [l.strip() for l in f if l.strip()]
    except OSError:
        return ""
    return lines[-1] if lines else ""


def run_pocket(sid: str, pdb_path: str, config: dict, reference: Optional[dict], paths: dict,
               prefix: List[str], step_id: str, runner=subprocess_runner) -> Tuple[list, list, list]:
    """Stage one pocket, run upstream once and collect its molecules."""
    run_dir = os.path.join(paths["runs_folder"], sid)
    os.makedirs(run_dir, exist_ok=True)
    if reference is not None:
        ligand_file = stage_reference_ligand(pdb_path, reference["code"], reference["smiles"], run_dir, sid)
        protein = os.path.join(run_dir, f"{sid}.pdb")
        stage_protein_pdb(pdb_path, protein)
    else:
        ligand_file = None
        protein = os.path.join(run_dir, f"{sid}.pdb")
        n_residues = stage_pocket_pdb(pdb_path, config["pocket_residues"], protein)
        print(f"  [flowr] {sid}: pocket of {n_residues} residue(s)")
    stem = os.path.splitext(os.path.basename(protein))[0]
    save_dir = os.path.join(run_dir, "flowr")
    os.makedirs(save_dir, exist_ok=True)
    cmd = build_command(config, protein, ligand_file, save_dir, paths["checkpoint"], prefix, paths["num_workers"])
    print(f"  [flowr] {sid}: {' '.join(cmd)}")
    env = dict(os.environ)
    env["PYTHONPATH"] = paths["flowr_repo"] + (os.pathsep + env["PYTHONPATH"] if env.get("PYTHONPATH") else "")
    log_path = os.path.join(run_dir, "flowr.log")
    returncode = runner(cmd, paths["flowr_repo"], env, log_path)
    if returncode != 0:
        raise RuntimeError(f"FLOWR exited with status {returncode}: {last_line(log_path)}")
    sdf_path = output_sdf(save_dir, stem, config["protonate_ligands"])
    if not os.path.isfile(sdf_path):
        raise RuntimeError(f"FLOWR wrote no {os.path.basename(sdf_path)}: {last_line(log_path)}")
    return collect_molecules(sid, sdf_path, config["num_molecules"], paths["structures_folder"], step_id)


def check_config(config: dict) -> None:
    """Refuse what upstream would only reject after loading its model."""
    if config["substructure"] is not None:
        from rdkit import Chem
        if Chem.MolFromSmiles(config["substructure"]) is None:
            raise ValueError(f"substructure is not a SMILES RDKit can parse: {config['substructure']!r}")


def available_workers() -> int:
    try:
        return max(1, len(os.sched_getaffinity(0)))
    except AttributeError:
        return max(1, os.cpu_count() or 1)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--config-yaml", required=True)
    p.add_argument("--structures-json", required=True)
    p.add_argument("--ligand-json", default=None)
    p.add_argument("--flowr-repo", required=True)
    p.add_argument("--checkpoint", required=True)
    p.add_argument("--runs-folder", required=True)
    p.add_argument("--compounds-csv", required=True)
    p.add_argument("--structures-folder", required=True)
    p.add_argument("--structures-map-csv", required=True)
    p.add_argument("--pockets-csv", required=True)
    p.add_argument("--local-missing-csv", required=True)
    p.add_argument("--container-prefix", default="")
    args = p.parse_args()

    with open(args.config_yaml) as f:
        config = yaml.safe_load(f)
    check_config(config)
    if not os.path.isfile(args.checkpoint):
        print(f"ERROR: checkpoint not found: {args.checkpoint}. Run FLOWR.install().", file=sys.stderr)
        sys.exit(1)

    reference = None
    pocket_label = f"residues:{config['pocket_residues']}"
    if args.ligand_json is not None:
        code = resolve_ligand_code(args.ligand_json)
        reference = {"code": code, "smiles": resolve_ligand_smiles(args.ligand_json)}
        pocket_label = f"ligand:{code}"

    paths = {"runs_folder": args.runs_folder, "structures_folder": args.structures_folder,
             "checkpoint": args.checkpoint, "flowr_repo": args.flowr_repo, "num_workers": available_workers()}
    prefix = container_argv_prefix(args.container_prefix)
    ds = load_datastream(args.structures_json)
    step_id = step_id_from_table_path(args.local_missing_csv)

    for folder in (args.runs_folder, args.structures_folder):
        os.makedirs(folder, exist_ok=True)

    compounds_rows, pose_rows, pocket_rows, missing_rows = [], [], [], []
    attempted = failures = 0
    for sid, pdb_path in iterate_files(ds):
        attempted += 1
        try:
            compounds, poses, missing = run_pocket(sid, pdb_path, config, reference, paths, prefix, step_id)
            compounds_rows += compounds
            pose_rows += poses
            missing_rows += missing
            pocket_rows.append({"id": sid, "pocket": pocket_label, "n_requested": config["num_molecules"],
                                "n_generated": len(compounds)})
            print(f"  [flowr] {sid}: {len(compounds)} of {config['num_molecules']} molecule(s)")
        except Exception as e:
            print(f"WARNING: {sid} failed: {e}", file=sys.stderr)
            failures += 1
            pocket_rows.append({"id": sid, "pocket": pocket_label, "n_requested": config["num_molecules"],
                                "n_generated": math.nan})
            missing_rows.append({"id": sid, "removed_by": step_id, "kind": "failure", "cause": str(e)[:200]})

    for path in (args.compounds_csv, args.structures_map_csv, args.pockets_csv, args.local_missing_csv):
        os.makedirs(os.path.dirname(path), exist_ok=True)
    pd.DataFrame(compounds_rows, columns=COMPOUND_COLUMNS).to_csv(args.compounds_csv, index=False)
    pd.DataFrame(pose_rows, columns=["id", "file", "structures.id"]).to_csv(args.structures_map_csv, index=False)
    pd.DataFrame(pocket_rows, columns=["id", "pocket", "n_requested", "n_generated"]).to_csv(args.pockets_csv, index=False)
    pd.DataFrame(missing_rows, columns=["id", "removed_by", "kind", "cause"]).to_csv(args.local_missing_csv, index=False)

    print(f"Compounds: {args.compounds_csv} ({len(compounds_rows)} molecules)")
    print(f"Pockets: {args.pockets_csv} ({len(pocket_rows)} rows, {failures} failed)")

    if not attempted:
        print("ERROR: no input entities to process (upstream stream is empty or its files are absent)",
              file=sys.stderr)
        sys.exit(1)
    if failures == attempted:
        print(f"ERROR: all {attempted} attempted pocket(s) failed", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
