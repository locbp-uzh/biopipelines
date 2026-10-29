#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Execution half of LigandDdG: protein-ligand binding energy per structure, with the ligand parameterized on the fly.

Follows the Chaperone blocks `ligand_ddg` (scoring="hard") and `ligand_ddg_soft` (scoring="soft"), RosettaCommons/chaperone, MIT license. The soft scheme packs with beta_nov16, scores with beta_nov16_soft and takes a trimmed mean over replicates; the hard scheme runs the minimization protocol below and reads the Ddg filter.

Rosetta needs a .params file for every non-standard residue. It is generated from the ligand in the structure: bond orders come from the compounds stream's SMILES, hydrogens are added at the requested pH when the structure carries none, atoms are renamed to unique element-first names, and Rosetta's molfile_to_params.py (installed beside the environment) writes the params. A pre-made params file bypasses all of that.
"""

import argparse
import math
import os
import subprocess
import sys

import yaml

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from biopipelines.biopipelines_io import (  # noqa: E402
    load_datastream,
    iterate_files,
    step_id_from_table_path,
)
from biopipelines.ligand_utils import resolve_ligand_code, resolve_ligand_smiles  # noqa: E402
from biopipelines.rosetta_scoring import (  # noqa: E402
    RELAX_INIT_OPTIONS,
    HARD_INIT_OPTIONS,
    molfile_to_params_path,
    init_pyrosetta,
    unbound_ddg,
    repack_around_chain,
    trimmed_mean_std,
    replicate_seed,
    seed_rosetta,
    collect_rows,
    write_outputs,
    step_exit_status,
)

COLUMNS = {
    "hard": ["id", "ddg_norepack"],
    "soft": ["id", "ddg_soft", "ddg_soft_std", "n_replicates_used"],
}

# The active steps of Chaperone's ddg_protocol.xml: backbone constraints, soft then hard minimization of the interface, then the no-repack Ddg filter across the ligand's jump, which is not jump 1 when the binder chain carries hetero residues (the movers' `jump` is a boolean, the interface task operation and the filter take a jump number). Its unused selectors, filters and metrics are left out. The schema requires a PROTOCOLS section even though the movers are applied one by one.
HARD_PROTOCOL_XML = """<ROSETTASCRIPTS>
  <SCOREFXNS>
    <ScoreFunction name="sfxn_softish" weights="beta_nov16">
      <Reweight scoretype="fa_rep" weight="0.15"/>
    </ScoreFunction>
    <ScoreFunction name="sfxn" weights="beta_nov16"/>
  </SCOREFXNS>
  <RESIDUE_SELECTORS>
    <Index name="ligand" resnums="%%ligand_res_number%%"/>
    <CloseContact name="interface_by_contact" residue_selector="ligand" contact_threshold="4"/>
    <And name="interface" selectors="interface_by_contact"/>
    <Not name="not_interface" selector="interface"/>
  </RESIDUE_SELECTORS>
  <TASKOPERATIONS>
    <ProteinInterfaceDesign name="pack_long" design_chain1="0" design_chain2="0" jump="%%ligand_jump%%" interface_distance_cutoff="15"/>
    <OperateOnResidueSubset name="restrict_to_interface" selector="not_interface"><PreventRepackingRLT/></OperateOnResidueSubset>
    <OperateOnResidueSubset name="restrict_target2repacking" selector="ligand"><PreventRepackingRLT/></OperateOnResidueSubset>
  </TASKOPERATIONS>
  <MOVERS>
    <AddConstraintsToCurrentConformationMover name="add_bb_cst" use_distance_cst="False" cst_weight="1" bb_only="1" sc_tip_only="0"/>
    <ClearConstraintsMover name="rm_bb_cst"/>
    <TaskAwareMinMover name="softish_min" scorefxn="sfxn_softish" chi="1" bb="1" jump="1" tolerance="0.01" cartesian="false" task_operations="restrict_target2repacking,restrict_to_interface"/>
    <TaskAwareMinMover name="hard_min" scorefxn="sfxn" chi="1" bb="1" jump="1" tolerance="0.01" cartesian="false" task_operations="restrict_target2repacking,restrict_to_interface"/>
    <TaskAwareMinMover name="min" scorefxn="sfxn" bb="0" chi="1" task_operations="pack_long"/>
  </MOVERS>
  <FILTERS>
    <Ddg name="ddg_norepack" threshold="0" jump="%%ligand_jump%%" repeats="1" repack="0" relax_mover="min" confidence="0" scorefxn="sfxn"/>
  </FILTERS>
  <PROTOCOLS/>
</ROSETTASCRIPTS>
"""
HARD_PROTOCOL_MOVERS = ("add_bb_cst", "softish_min", "hard_min", "rm_bb_cst")

# The longest range of a pairwise Rosetta score term: a ligand with no protein atom this close cannot interact, and its ddG is zero.
CONTACT_RANGE_A = 6.0


# ---------------------------------------------------------------------------
# PDB handling (pure Python)
# ---------------------------------------------------------------------------

def element_of(line):
    """The element of a PDB atom line: the element column, else the atom name without digits."""
    element = line[76:78].strip()
    if not element:
        name = line[12:16].strip().lstrip("0123456789")
        element = name[:1]
    return element.upper()


def split_complex(pdb_path, binder_chain, ligand_chain, code):
    """The binder chain without hydrogens, and the ligand residue as written in the structure.

    Everything else (other chains, waters, other hetero residues on the ligand chain) is dropped. Exactly one residue named `code` must sit on `ligand_chain`.
    """
    binder, ligand, residues = [], [], set()
    with open(pdb_path, errors="ignore") as fh:
        for line in fh:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            chain = line[21]
            if chain == binder_chain:
                if element_of(line) != "H":
                    binder.append(line if line.endswith("\n") else line + "\n")
            elif chain == ligand_chain and line[17:20].strip() == code:
                ligand.append(line if line.endswith("\n") else line + "\n")
                residues.add((line[22:26], line[26]))
    if not binder:
        raise ValueError(f"no atoms on binder chain {binder_chain!r}")
    if not ligand:
        raise ValueError(
            f"no residue {code!r} on chain {ligand_chain!r} (a PDB residue name holds at most 3 characters)")
    if len(residues) > 1:
        raise ValueError(f"{len(residues)} residues named {code!r} on chain {ligand_chain!r}; expected one")
    return binder, ligand


def nearest_distance(binder, ligand):
    """Smallest heavy-atom distance in A between the binder chain and the ligand."""
    def position(line):
        return float(line[30:38]), float(line[38:46]), float(line[46:54])

    binder_atoms = [position(line) for line in binder]
    ligand_atoms = [position(line) for line in ligand if element_of(line) != "H"]
    return min(math.dist(a, b) for a in ligand_atoms for b in binder_atoms)


def hetatm_line(serial, name, element, resname, chain, resseq, xyz):
    """A PDB HETATM record with the atom name aligned as molfile_to_params writes it."""
    padded = f" {name:<3}" if len(element) == 1 and len(name) <= 3 else f"{name:<4}"
    return (f"HETATM{serial:5d} {padded} {resname:>3s} {chain:1s}{resseq:4d}    "
            f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00  0.00          {element.upper():>2s}\n")


# ---------------------------------------------------------------------------
# Ligand chemistry and parameters
# ---------------------------------------------------------------------------

def element_first_names(symbols):
    """Unique atom names of the form <ELEMENT><n>: Rosetta's molfile reader infers the element from the leading letters of the name."""
    counts, names = {}, []
    for symbol in symbols:
        symbol = symbol.upper()
        counts[symbol] = counts.get(symbol, 0) + 1
        name = f"{symbol}{counts[symbol]}"
        if len(name) > 4:
            raise ValueError(f"cannot name atom {name!r}: Rosetta atom names hold at most 4 characters")
        names.append(name)
    return names


def sdf_with_atom_names(mol, names):
    """A one-record SDF whose atom lines carry `names` in columns 30-35, where Rosetta's molfile reader looks for the name."""
    from rdkit import Chem

    lines = Chem.MolToMolBlock(mol).splitlines()
    n_atoms = int(lines[3][0:3])
    if n_atoms != len(names):
        raise ValueError(f"{len(names)} names for {n_atoms} atoms")
    for i, name in enumerate(names):
        line = lines[4 + i]
        lines[4 + i] = line[:30] + f"{name:>5}" + line[35:]
    return "\n".join(lines) + "\n$$$$\n"


def protonated_ligand(ligand_pdb, smiles, ph):
    """RDKit mol carrying the structure's coordinates, the SMILES' bond orders and explicit hydrogens.

    A ligand the structure already carries hydrogens for is used as it is. A hydrogen-free one is protonated at `ph` with Dimorphite-DL and given hydrogens placed on its heavy-atom frame.
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from biopipelines.ligand_utils import templated_ligand_mol

    mol = templated_ligand_mol(ligand_pdb, smiles)
    if not any(atom.GetAtomicNum() == 1 for atom in mol.GetAtoms()):
        from dimorphite_dl import protonate_smiles
        states = protonate_smiles(Chem.MolToSmiles(mol), ph_min=ph, ph_max=ph, precision=0.5)
        if not states:
            raise ValueError(f"Dimorphite-DL returned no protonation state at pH {ph}")
        template = Chem.MolFromSmiles(states[0])
        if template is None:
            raise ValueError(f"RDKit cannot parse the protonated SMILES {states[0]!r}")
        mol = AllChem.AssignBondOrdersFromTemplate(template, mol)
        Chem.SanitizeMol(mol)
        mol = Chem.AddHs(mol, addCoords=True)

    # molfile_to_params types a hydrogen from the type its carbon already received, so every heavy atom has to precede every hydrogen.
    heavy = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() != 1]
    hydrogens = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 1]
    return Chem.RenumberAtoms(mol, heavy + hydrogens)


def prepare_ligand(ligand_lines, code, ligand_chain, smiles, ph, work_dir):
    """Write the ligand as an SDF and as renamed PDB records. Returns (sdf_path, hetatm_lines)."""
    resseq = int(ligand_lines[0][22:26])
    raw_pdb = os.path.join(work_dir, "ligand_input.pdb")
    with open(raw_pdb, "w") as fh:
        fh.writelines(ligand_lines)
        fh.write("END".ljust(80) + "\n")

    mol = protonated_ligand(raw_pdb, smiles, ph)
    symbols = [atom.GetSymbol() for atom in mol.GetAtoms()]
    names = element_first_names(symbols)

    sdf_path = os.path.join(work_dir, f"{code}.sdf")
    with open(sdf_path, "w") as fh:
        fh.write(sdf_with_atom_names(mol, names))

    conformer = mol.GetConformer()
    hetatm = []
    for i, (name, symbol) in enumerate(zip(names, symbols)):
        position = conformer.GetAtomPosition(i)
        hetatm.append(hetatm_line(i + 1, name, symbol, code, ligand_chain, resseq,
                                  (position.x, position.y, position.z)))
    return sdf_path, hetatm


def generate_params(sdf_path, code, work_dir):
    """Run Rosetta's molfile_to_params.py on the SDF; returns the .params path."""
    script = molfile_to_params_path()
    if not os.path.isfile(script):
        raise FileNotFoundError(
            f"molfile_to_params.py not found at {script}; run LigandDdG.install() to fetch it, "
            f"or pass params_file=")
    name = code[:3]
    prefix = os.path.join(work_dir, name)
    cmd = [sys.executable, script, "--name", name, "-p", prefix, sdf_path,
           "--no-pdb", "--clobber", "--keep-names"]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        tail = (result.stdout + result.stderr).strip()[-600:]
        raise RuntimeError(f"molfile_to_params failed for {code} (exit {result.returncode}): {tail}")
    params = prefix + ".params"
    if not os.path.isfile(params):
        raise RuntimeError(
            f"molfile_to_params wrote no {os.path.basename(params)}; the ligand may have been split into fragments")
    return params


# ---------------------------------------------------------------------------
# Scoring (PyRosetta)
# ---------------------------------------------------------------------------

def load_pose(binder_lines, ligand_lines, params):
    import pyrosetta as pyr

    text = "".join(binder_lines) + "TER".ljust(80) + "\n" + "".join(ligand_lines) + "END".ljust(80) + "\n"
    pose = pyr.Pose()
    pyr.generate_nonstandard_residue_set(pose, [params])
    pyr.rosetta.core.import_pose.pose_from_pdbstring(pose, text, pose.residue_type_set_for_pose(), "complex.pdb")
    return pose


def check_ligand_names(residue, ligand_lines, params):
    """Refuse a structure whose ligand atom names differ from the ones the params file defines.

    Rosetta drops an atom it has no name for and does not say so, which scores the wrong molecule without an error.
    """
    names = [residue.atom_name(i).strip() for i in range(1, residue.natoms() + 1)]
    heavy = {name for i, name in enumerate(names, start=1) if not residue.atom_is_hydrogen(i)}
    present = {line[12:16].strip() for line in ligand_lines}
    unknown, absent = sorted(present - set(names)), sorted(heavy - present)
    if unknown or absent:
        raise ValueError(
            f"the ligand in the structure does not match {os.path.basename(params)}: atoms {unknown[:6]} are not in the "
            f"params and heavy atoms {absent[:6]} of the params are not in the structure")


def soft_ddg(pose, binder_chain, ligand_chain, replicates, seeds):
    """Pack with beta_nov16, score with beta_nov16_soft, repeat with one seed per replicate; trimmed mean over the replicates."""
    import pyrosetta as pyr

    hard = pyr.create_score_function("beta_nov16")
    soft = pyr.create_score_function("beta_nov16_soft")
    values = []
    for replicate in range(replicates):
        seed_rosetta(seeds[replicate])
        repack_around_chain(pose, ligand_chain, hard)
        values.append(unbound_ddg(pose, soft, binder_chain))
    mean, std, used = trimmed_mean_std(values)
    return {"ddg_soft": mean, "ddg_soft_std": std, "n_replicates_used": used}


def hard_ddg(pose, ligand_index):
    """Minimize the interface with backbone constraints, then read the Ddg filter (no repacking) across the ligand jump."""
    from pyrosetta.rosetta.protocols.rosetta_scripts import XmlObjects

    jump = pose.fold_tree().get_jump_that_builds_residue(ligand_index)
    xml = HARD_PROTOCOL_XML.replace("%%ligand_res_number%%", str(ligand_index)).replace("%%ligand_jump%%", str(jump))
    objects = XmlObjects.create_from_string(xml)
    for mover in HARD_PROTOCOL_MOVERS:
        objects.get_mover(mover).apply(pose)
    return {"ddg_norepack": float(objects.get_filter("ddg_norepack").report_sm(pose))}


def score_structure(sid, path, config, code, smiles, work_root):
    """ddG of one structure; the ligand is parameterized from the structure itself unless a params file was given."""
    scoring = config["scoring"]
    binder_lines, ligand_lines = split_complex(path, config["binder_chain"], config["ligand_chain"], code)
    gap = nearest_distance(binder_lines, ligand_lines)
    if gap > CONTACT_RANGE_A:
        raise ValueError(
            f"ligand {code} on chain {config['ligand_chain']!r} is {gap:.0f} A from chain {config['binder_chain']!r}, "
            f"so there is no interface to score; a crystal with several copies may list another copy's ligand there")

    work_dir = os.path.join(work_root, sid)
    os.makedirs(work_dir, exist_ok=True)
    params = config.get("params_file") or ""
    if params:
        if not os.path.isfile(params):
            raise FileNotFoundError(f"params_file not found: {params}")
    else:
        if not smiles:
            raise ValueError(
                "the ligand compounds stream carries no SMILES to take bond orders from; pass "
                "Ligand(smiles=..., codes=...) or a pre-made params_file")
        sdf_path, ligand_lines = prepare_ligand(ligand_lines, code, config["ligand_chain"], smiles,
                                                config["ph"], work_dir)
        params = generate_params(sdf_path, code, work_dir)

    pose = load_pose(binder_lines, ligand_lines, params)
    if config.get("params_file"):
        check_ligand_names(pose.residue(pose.total_residue()), ligand_lines, params)
    if scoring == "soft":
        seeds = [replicate_seed(config["seed"], sid, r) for r in range(config["replicates"])]
        result = soft_ddg(pose, config["binder_chain"], config["ligand_chain"], config["replicates"], seeds)
    else:
        seed_rosetta(replicate_seed(config["seed"], sid, 0))
        result = hard_ddg(pose, pose.total_residue())
    print("  " + ", ".join(f"{k}={v:.2f}" if isinstance(v, float) else f"{k}={v}" for k, v in result.items()))
    return result


def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument("--config-yaml", required=True)
    p.add_argument("--structures-json", required=True)
    p.add_argument("--ligand-json", required=True)
    p.add_argument("--work-dir", required=True)
    p.add_argument("--table-csv", required=True)
    p.add_argument("--local-missing-csv", required=True)
    args = p.parse_args(argv)

    with open(args.config_yaml) as fh:
        config = yaml.safe_load(fh) or {}
    columns = COLUMNS[config["scoring"]]

    ds = load_datastream(args.structures_json)
    code = resolve_ligand_code(args.ligand_json)
    smiles = None if config.get("params_file") else resolve_ligand_smiles(args.ligand_json)
    step_id = step_id_from_table_path(args.local_missing_csv)
    init_pyrosetta(RELAX_INIT_OPTIONS if config["scoring"] == "soft" else HARD_INIT_OPTIONS)

    rows, missing, attempted, failures = collect_rows(
        iterate_files(ds),
        lambda sid, path: score_structure(sid, path, config, code, smiles, args.work_dir),
        columns, step_id)
    write_outputs(args.table_csv, rows, columns, args.local_missing_csv, missing)

    print(f"ligand_ddg ({config['scoring']}): {args.table_csv} ({len(rows)} rows, {failures} failed)")
    error = step_exit_status(attempted, failures)
    if error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
