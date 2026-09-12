#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""AutoDock Vina backend for the docking pipeline.

Vina differs from GNINA in three ways this module absorbs, so the surrounding
pipeline (conformer prep, pose-consistency clustering, aggregation, complex
writing) runs unchanged:

1. Vina reads and writes PDBQT only. Receptor PDB and ligand SDF are converted
   in, and the result is converted back out to SDF.
2. Vina writes its score in a ``REMARK VINA RESULT`` line rather than an SDF
   property. It is re-attached as ``minimizedAffinity`` so the shared parser
   reads it exactly as it reads GNINA's.
3. ``--autobox`` boxes on the input ligand and takes no reference file, so an
   autobox_ligand reference is turned into an explicit center+size here.

Vina has no CNN, so no CNN properties are ever emitted; the aggregation layer
already treats a missing CNN term as absent rather than zero.
"""

import os
import re
import subprocess
import sys

from rdkit import Chem
from rdkit.Chem import AllChem

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from biopipelines.pdb_parser import parse_pdb_file

_VINA_RESULT_RE = re.compile(
    r"^REMARK\s+VINA\s+RESULT:\s*(-?\d+\.?\d*)", re.MULTILINE
)


def convert(src, dst, extra_args=None):
    """Run an OpenBabel conversion, raising with obabel's own message on failure."""
    cmd = ["obabel", src, "-O", dst] + list(extra_args or [])
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0 or not os.path.exists(dst):
        raise RuntimeError(
            f"obabel failed converting {os.path.basename(src)} -> "
            f"{os.path.basename(dst)}: {result.stderr[:300]}"
        )
    return dst


def receptor_pdbqt(protein_pdb, dst):
    """Convert a prepared receptor PDB to PDBQT.

    ``-xr`` marks the receptor rigid, which suppresses the ROOT/BRANCH torsion
    tree Vina rejects on a receptor.
    """
    return convert(protein_pdb, dst, ["-xr"])


def ligand_pdbqt(ligand_sdf, dst):
    """Convert a ligand SDF to PDBQT, keeping the torsion tree Vina searches over."""
    return convert(ligand_sdf, dst, ["-opdbqt"])


def box_args(box, ligand_sdf=None, receptor_file=None):
    """Build Vina box flags.

    Explicit center+size passes through. An autobox_ligand reference has no Vina
    equivalent (``--autobox`` boxes on the *input* ligand), so its coordinates are
    read here and turned into an explicit box.
    """
    if "center" in box and "size" in box:
        cx, cy, cz = [s.strip() for s in str(box["center"]).split(",")]
        args = ["--center_x", cx, "--center_y", cy, "--center_z", cz]
        size = box["size"]
        if isinstance(size, (int, float)):
            args += ["--size_x", str(size), "--size_y", str(size), "--size_z", str(size)]
        else:
            sx, sy, sz = [s.strip() for s in str(size).split(",")]
            args += ["--size_x", sx, "--size_y", sy, "--size_z", sz]
        return args

    reference = box.get("autobox_ligand")
    from_input_ligand = False
    if not reference and ligand_sdf:
        # Last resort: box on the INPUT ligand -- nothing has been docked yet. Its
        # coordinates mean something only if it was carved from this receptor
        # (Ligand(structures=...)), where they are the crystal pose. A conformer
        # RDKit embedded from SMILES is centred on the origin, so the box would sit
        # tens of angstroms off the protein and Vina would search empty space.
        reference = ligand_sdf
        from_input_ligand = True
    if reference:
        args = _box_from_reference(reference, float(box.get("autobox_add", 4.0)))
        if from_input_ligand and receptor_file:
            _assert_box_overlaps_receptor(args, receptor_file, reference)
        return args
    return []


def _assert_box_overlaps_receptor(box_args_list, receptor_file, reference_file):
    """Refuse a box that does not reach the receptor at all.

    Catches the origin-centred-conformer case decisively rather than warning:
    a docking run whose box misses the protein returns scores, and they look
    like any other scores.
    """
    flags = dict(zip(box_args_list[::2], box_args_list[1::2]))
    try:
        centre = [float(flags[f"--center_{a}"]) for a in "xyz"]
    except (KeyError, ValueError):
        return

    atoms = parse_pdb_file(receptor_file)
    if not atoms:
        return
    lo = [min(getattr(a, c) for a in atoms) for c in ("x", "y", "z")]
    hi = [max(getattr(a, c) for a in atoms) for c in ("x", "y", "z")]

    pad = 10.0
    if all(lo[i] - pad <= centre[i] <= hi[i] + pad for i in range(3)):
        return
    raise RuntimeError(
        f"no binding box was defined and the input ligand "
        f"({os.path.basename(reference_file)}) is not positioned in the receptor: "
        f"its box centre is ({centre[0]:.1f}, {centre[1]:.1f}, {centre[2]:.1f}) while "
        f"{os.path.basename(receptor_file)} spans x[{lo[0]:.1f}, {hi[0]:.1f}] "
        f"y[{lo[1]:.1f}, {hi[1]:.1f}] z[{lo[2]:.1f}, {hi[2]:.1f}]. A generated "
        f"conformer carries no information about where the site is. Pass center= and "
        f"size=, or autobox_ligand=, or use a receptor whose crystal ligand marks the "
        f"pocket.")


def _box_from_reference(reference_file, padding):
    """Center+size covering a reference ligand's heavy atoms, plus padding."""
    coords = _read_coords(reference_file)
    if not coords:
        raise RuntimeError(f"no coordinates found in autobox reference {reference_file}")

    dims = []
    center = []
    for axis in range(3):
        values = [c[axis] for c in coords]
        lo, hi = min(values), max(values)
        center.append((lo + hi) / 2.0)
        dims.append((hi - lo) + 2 * padding)

    args = []
    for flag, value in zip(("--center_x", "--center_y", "--center_z"), center):
        args += [flag, f"{value:.3f}"]
    for flag, value in zip(("--size_x", "--size_y", "--size_z"), dims):
        # Vina refuses a degenerate axis; a planar or single-atom reference
        # would otherwise produce a zero-width box.
        args += [flag, f"{max(value, 2.0):.3f}"]
    return args


def _read_coords(path):
    """Heavy-atom coordinates from a PDB or an RDKit-readable ligand file."""
    if path.lower().endswith((".pdb", ".pdbqt", ".ent")):
        return [(a.x, a.y, a.z) for a in parse_pdb_file(path)]

    mol = Chem.MolFromMolFile(path, removeHs=True)
    if mol is None or mol.GetNumConformers() == 0:
        return []
    conf = mol.GetConformer()
    return [(conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y,
             conf.GetAtomPosition(i).z)
            for i in range(mol.GetNumAtoms())]


def _restore_bond_orders(template, mol):
    """Re-type a PDBQT-derived pose against its input chemistry.

    PDBQT keeps polar hydrogens, so the bare template often fails to match; the
    AddHs form is tried second, matching how the framework templates ligands
    elsewhere.

    When neither matches, the pose is returned with the bond orders PDBQT implied
    -- an all-single-bond skeleton. The coordinates are still valid, so this is
    not a failure worth dropping the pose over, but the chemistry is wrong and
    every downstream consumer that reads the SDF (PoseBusters, RTMScore,
    PoseChange's graph pairing) would take it at face value. The molecule is
    tagged so that is visible in the output rather than only on stderr.
    """
    for tmpl, lig in ((template, mol),
                      (Chem.AddHs(template), Chem.AddHs(mol, addCoords=True))):
        try:
            out = AllChem.AssignBondOrdersFromTemplate(tmpl, lig)
            out.SetProp("bond_orders_restored", "true")
            return out
        except Exception:
            continue
    print("  Warning: could not restore bond orders from the input template; the "
          "pose keeps PDBQT's single-bond perception and is tagged "
          "bond_orders_restored=false", file=sys.stderr)
    mol.SetProp("bond_orders_restored", "false")
    return mol


def pdbqt_to_sdf(pdbqt_path, sdf_path, template_sdf=None):
    """Convert a Vina output PDBQT to SDF, re-attaching each pose's affinity.

    Vina writes one MODEL per pose, each preceded by its own REMARK score. obabel
    emits the models in order, so the Nth affinity belongs to the Nth molecule.
    The score is stored as ``minimizedAffinity`` — the tag GNINA uses — so the
    shared parser needs no backend branch.
    """
    if not os.path.exists(pdbqt_path):
        return None

    with open(pdbqt_path, "r") as handle:
        affinities = [float(m) for m in _VINA_RESULT_RE.findall(handle.read())]

    raw_sdf = sdf_path + ".raw.sdf"
    convert(pdbqt_path, raw_sdf)

    supplier = Chem.SDMolSupplier(raw_sdf, removeHs=False)
    mols = [m for m in supplier if m is not None]
    if not mols:
        return None

    # PDBQT carries no bond orders; recover them from the input chemistry so the
    # emitted pose is a valid molecule rather than an all-single-bond skeleton.
    template = None
    if template_sdf and os.path.exists(template_sdf):
        template = Chem.MolFromMolFile(template_sdf, removeHs=False)

    writer = Chem.SDWriter(sdf_path)
    for index, mol in enumerate(mols):
        if template is not None:
            mol = _restore_bond_orders(template, mol)
        if index < len(affinities):
            mol.SetProp("minimizedAffinity", str(affinities[index]))
        writer.write(mol)
    writer.close()

    os.remove(raw_sdf)
    return sdf_path


def run(vina_binary, receptor_pdbqt_path, ligand_pdbqt_path, out_pdbqt,
        box_flags, scoring="vina", exhaustiveness=8, num_modes=9, seed=0,
        cpu=None, mode="docking", timeout=None, extra_args=None):
    """Invoke the Vina binary. Returns the CompletedProcess.

    ``score_only``/``local_only`` replace the search, but Vina still requires a
    box: without one it exits with "Grid box dimensions must be greater than 0
    Angstrom". This is where it differs from gnina, which autoboxes on the input.
    """
    cmd = [
        vina_binary,
        "--receptor", receptor_pdbqt_path,
        "--ligand", ligand_pdbqt_path,
        "--scoring", scoring,
    ]

    if mode == "score":
        cmd += ["--score_only"] + list(box_flags)
    elif mode == "minimize":
        cmd += ["--local_only", "--out", out_pdbqt] + list(box_flags)
    else:
        cmd += list(box_flags) + [
            "--exhaustiveness", str(exhaustiveness),
            "--num_modes", str(num_modes),
            "--seed", str(seed),
            "--out", out_pdbqt,
        ]

    # Forwarded kwargs as argv tokens; no shell is involved, so nothing re-parses them.
    cmd += list(extra_args or [])

    if cpu:
        cmd += ["--cpu", str(cpu)]

    return subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)


def parse_score_stdout(stdout):
    """Affinity from a ``--score_only`` run, which prints rather than writing a file."""
    match = re.search(r"^\s*Estimated Free Energy of Binding\s*:\s*(-?\d+\.?\d*)",
                      stdout, re.MULTILINE)
    if match:
        return float(match.group(1))
    match = _VINA_RESULT_RE.search(stdout)
    return float(match.group(1)) if match else None
