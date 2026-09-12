#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""OpenMM minimiser. Per input PDB:
  * Load through PDBFixer-equivalent flow (just openmm.app.PDBFile if input is
    already cleanable; otherwise rely on the forcefield to error out and we
    skip the structure).
  * Build the system from the chosen protein forcefield + implicit solvent (or
    vacuum) so no waters/ions are needed. With --ligand-json the bound small
    molecule is parameterised from its SMILES via OpenFF instead of needing a
    residue template.
  * Optionally restrain a residue selection harmonically, and/or freeze atoms
    outright by zeroing their mass.
  * Minimise to the supplied tolerance; write <id>.pdb to output-dir and one
    energy row to the CSV.
"""

import argparse
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from biopipelines.biopipelines_io import load_datastream, iterate_files  # noqa: E402
from biopipelines.sele_utils import sele_to_list  # noqa: E402

import openmm  # noqa: E402
from openmm import LangevinMiddleIntegrator, Platform, unit  # noqa: E402
from openmm.app import (  # noqa: E402
    ForceField, HBonds, Modeller, NoCutoff, PDBFile, Simulation,
)


E_COLS = ["id", "energy_initial_kj_mol", "energy_final_kj_mol", "delta_kj_mol",
          "n_mobile_atoms", "n_frozen_atoms", "charge_method"]

# Map the tool's forcefield/solvent enums onto OpenMM's bundled XML files.
FF_XML = {
    "amber14-all": "amber14-all.xml",
    "amber99sb": "amber99sb.xml",
    "charmm36": "charmm36.xml",
}
SOLVENT_XML = {
    "implicit-gbn2": "implicit/gbn2.xml",
    "implicit-gbn": "implicit/gbn.xml",
    "implicit-obc2": "implicit/obc2.xml",
    "vacuum": None,
}

# Backbone atom names, frozen even for mobile residues so only rotamers relax.
BACKBONE = {"N", "CA", "C", "O", "OXT"}

NAGL_MODEL = "openff-gnn-am1bcc-0.1.0-rc.3.pt"


def _add_restraint(system, modeller, restraint_pairs, k_kjmolnm2):
    """Harmonically restrain heavy atoms of the selected (chain, resnum) pairs."""
    if not restraint_pairs:
        return
    sel = set(restraint_pairs)
    force = openmm.CustomExternalForce("0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    force.addGlobalParameter("k", k_kjmolnm2 * unit.kilojoule_per_mole / unit.nanometer**2)
    for p in ("x0", "y0", "z0"):
        force.addPerParticleParameter(p)
    positions = modeller.positions
    for atom in modeller.topology.atoms():
        if atom.element is not None and atom.element.symbol == "H":
            continue
        chain = atom.residue.chain.id
        try:
            resnum = int(atom.residue.id)
        except (TypeError, ValueError):
            continue
        if (chain, resnum) in sel:
            pos = positions[atom.index].value_in_unit(unit.nanometer)
            force.addParticle(atom.index, [pos[0], pos[1], pos[2]])
    if force.getNumParticles() > 0:
        system.addForce(force)


def _as_reference(value):
    """Return a TableReference when the value is one, else None."""
    from biopipelines.biopipelines_io import TableReference

    if value and str(value).startswith(f"{TableReference.PREFIX}:"):
        return TableReference.from_string(str(value))
    return None


def _matches(atom, sel):
    """True when an atom's residue is in a (chain, resnum) selection set."""
    try:
        resnum = int(atom.residue.id)
    except (TypeError, ValueError):
        return False
    chain = atom.residue.chain.id
    return (chain, resnum) in sel or ("", resnum) in sel


def _apply_freeze(system, topology, mobile_pairs, frozen_pairs, ligand_resnames):
    """Zero the mass of every atom that must not move.

    A zero-mass particle is immovable in OpenMM, unlike a restrained one which
    only pays an energy penalty. Returns (n_mobile, n_frozen).
    """
    if not mobile_pairs and not frozen_pairs:
        return system.getNumParticles(), 0

    freeze = []
    if mobile_pairs:
        sel = set(mobile_pairs)
        for atom in topology.atoms():
            if atom.residue.name in ligand_resnames:
                continue
            # Mobile residues keep their side chains free; backbone stays put so
            # the fold cannot drift and contaminate a downstream pose metric.
            if _matches(atom, sel) and atom.name not in BACKBONE:
                continue
            freeze.append(atom.index)
    else:
        sel = set(frozen_pairs)
        for atom in topology.atoms():
            if atom.residue.name in ligand_resnames:
                continue
            if _matches(atom, sel):
                freeze.append(atom.index)

    frozen = set(freeze)
    for idx in freeze:
        system.setParticleMass(idx, 0.0)

    # A constraint on a massless particle is an error, and HBonds constrains every
    # X-H bond including those inside the frozen region. Drop constraints whose
    # atoms are all frozen -- they are already held by the zero mass.
    for i in reversed(range(system.getNumConstraints())):
        p1, p2, _ = system.getConstraintParameters(i)
        if p1 in frozen or p2 in frozen:
            system.removeConstraint(i)

    return system.getNumParticles() - len(freeze), len(freeze)


def _split_ligand(pdb_path, ligand_resnames, work_dir, sid):
    """Split a complex into (protein-only PDB, ligand-only PDB).

    Returns (protein_path, ligand_path); ligand_path is None when the structure
    carries none of the named residues.
    """
    protein_lines, ligand_lines = [], []
    for line in open(pdb_path):
        if line.startswith(("ATOM", "HETATM")):
            resname = line[17:20].strip()
            if resname in ligand_resnames:
                ligand_lines.append(line)
                continue
            if resname == "HOH":
                continue
            protein_lines.append(line)
        elif line.startswith("TER"):
            protein_lines.append(line)

    protein_path = os.path.join(work_dir, f"{sid}_protein.pdb")
    with open(protein_path, "w") as f:
        f.writelines(protein_lines)
        f.write("TER".ljust(80) + "\n" + "END".ljust(80) + "\n")

    if not ligand_lines:
        return protein_path, None
    ligand_path = os.path.join(work_dir, f"{sid}_ligand.pdb")
    with open(ligand_path, "w") as f:
        f.writelines(ligand_lines)
        f.write("END".ljust(80) + "\n")
    return protein_path, ligand_path


def _fix_protein(protein_path, work_dir, sid):
    """Add terminal OXT and any missing heavy atoms so the forcefield matches.

    Side-chain packers write bare ATOM records with no OXT and no TER, which the
    forcefield reads as a chain broken mid-backbone.
    """
    from pdbfixer import PDBFixer

    fixer = PDBFixer(filename=protein_path)
    fixer.findMissingResidues()
    # Only terminal atoms and hydrogens should be built; inserting missing loop
    # residues would invent coordinates the design never specified.
    fixer.missingResidues = {}
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixed_path = os.path.join(work_dir, f"{sid}_protein_fixed.pdb")
    with open(fixed_path, "w") as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f, keepIds=True)
    return fixed_path


def _name_ligand_residues(topology, n_protein_res, ligand_code):
    """Give every added ligand atom one residue name, number and chain.

    OpenFF carries residue names only for the heavy atoms it read from the PDB; the
    hydrogens it generates land in a separate UNK residue on their own chain. The
    physics is unaffected (bonds come from the molecule graph), but any downstream
    tool selecting the ligand by residue name then sees a fragment of it.
    """
    # Modeller.add appends, so every residue past the protein's count is ligand.
    renamed = 0
    for residue in topology.residues():
        if residue.index < n_protein_res:
            continue
        residue.name = ligand_code
        renamed += sum(1 for _ in residue.atoms())
    return renamed


def _find_covalent_pair(topology, positions, anchor_atom, ligand_resnames, max_distance):
    """Locate the protein-ligand atom pair that carries the covalent link.

    Returns (protein_index, ligand_index, distance_angstrom) for the closest pair
    between `anchor_atom` (e.g. "SG", optionally "SG62" to pin the residue) and any
    ligand heavy atom, or None when nothing is within `max_distance`.
    """
    import numpy as np

    name = "".join(c for c in anchor_atom if not c.isdigit())
    resnum = "".join(c for c in anchor_atom if c.isdigit())

    prot, lig = [], []
    for atom in topology.atoms():
        if atom.residue.name in ligand_resnames:
            if atom.element is not None and atom.element.symbol != "H":
                lig.append(atom.index)
        elif atom.name == name:
            if resnum and str(atom.residue.id) != resnum:
                continue
            prot.append(atom.index)
    if not prot or not lig:
        return None

    coords = np.array(positions.value_in_unit(unit.angstrom))
    best = min(((float(np.linalg.norm(coords[p] - coords[l])), p, l)
                for p in prot for l in lig), key=lambda t: t[0])
    d, p, l = best
    return (p, l, d) if d <= max_distance else None


def _add_covalent_restraint(system, pair, k_kjmolnm2, length_nm):
    """Hold the covalent link with a stiff harmonic bond.

    Not a real bond -- no angle or torsion terms cross the junction, so the ligand
    can still swing about the attachment point -- but the linkage cannot break.

    The 1-2 pair is excluded from the nonbonded terms, as a real bonded pair would
    be. Without that, two heavy atoms held at 0.18 nm still see each other through
    full Lennard-Jones and Coulomb: an enormous spurious repulsion that dominates
    the reported energies and pulls the minimised geometry away from the adduct,
    so the restraint would fight a clash that does not exist in the real molecule.
    """
    p, l, _ = pair
    force = openmm.HarmonicBondForce()
    force.addBond(p, l, length_nm, k_kjmolnm2)
    system.addForce(force)

    excluded = False
    for f in system.getForces():
        if isinstance(f, openmm.NonbondedForce):
            # replace=True: the pair may already carry a 1-4 style exception.
            f.addException(p, l, 0.0, 1.0, 0.0, True)
            excluded = True
        elif isinstance(f, openmm.CustomNonbondedForce):
            existing = {tuple(sorted(f.getExclusionParticles(i)))
                        for i in range(f.getNumExclusions())}
            if tuple(sorted((p, l))) not in existing:
                f.addExclusion(p, l)
            excluded = True
    if not excluded:
        raise RuntimeError(
            'covalent_anchor: the system carries no nonbonded force to exclude the '
            'linked pair from, so the reported energies would be dominated by a '
            'clash between two atoms held a bond length apart')
    return force


def _ligand_positions(mol):
    """The molecule's conformer as an OpenMM-ready quantity."""
    return mol.conformers[0].to_openmm()


_CHARGE_CACHE = {}


def _ligand_molecule(ligand_path, smiles, charge_method):
    """Build an OpenFF Molecule for the bound ligand, with coordinates.

    Returns (molecule, charge_method_used). The second value is not always the
    requested one -- NAGL can decline a molecule outside its element coverage --
    and it reaches the energies table so a mixed-method run is visible there.

    Partial charges depend on the molecular graph, not the conformer, so a run over
    many poses of one ligand computes them once. AM1-BCC on a 44-atom dye takes
    minutes, and a 340-structure run would otherwise spend days recomputing an
    identical answer.
    """
    from openff.toolkit.topology import Molecule
    from rdkit import Chem
    from biopipelines.ligand_utils import posed_ligand_mol

    rdmol = posed_ligand_mol(ligand_path, smiles)
    # Heavy-atom-only input: hydrogens are built onto the existing coordinates so
    # the pose is untouched, and the forcefield needs them present.
    rdmol = Chem.AddHs(rdmol, addCoords=True)
    if rdmol.GetNumConformers() == 0:
        raise ValueError(
            f"{os.path.basename(ligand_path)} yielded no conformer; the bound pose is "
            f"what is being minimised, so it cannot be replaced by a generated one")
    mol = Molecule.from_rdkit(rdmol, allow_undefined_stereo=True)

    # Key on the canonical graph, so a different pose of the same molecule hits.
    cache_key = (mol.to_smiles(isomeric=True, mapped=False), charge_method)
    cached = _CHARGE_CACHE.get(cache_key)
    if cached is not None:
        charges, used = cached
        mol.partial_charges = charges
        return mol, used

    if charge_method == "nagl":
        try:
            mol.assign_partial_charges(NAGL_MODEL)
            _CHARGE_CACHE[cache_key] = (mol.partial_charges, "nagl")
            return mol, "nagl"
        except Exception as e:
            # NAGL's element coverage excludes Si, P and metals, so a silicon
            # rhodamine or a phosphonate lands here rather than being a config error.
            print(f"  NAGL charges unavailable, using AM1-BCC instead: "
                  f"{str(e).splitlines()[0][:160]}", file=sys.stderr)
    try:
        mol.assign_partial_charges("am1bcc")
    except Exception as e:
        raise ValueError(
            f"AM1-BCC charge assignment failed: {str(e).splitlines()[0][:200]}. "
            f"This needs AmberTools on PATH (antechamber/sqm) — check that the "
            f"openmm environment is activated rather than the interpreter "
            f"being called by absolute path.") from e
    _CHARGE_CACHE[cache_key] = (mol.partial_charges, "am1bcc")
    return mol, "am1bcc"


def minimise(pdb_path, out_path, max_iterations, tolerance_kjmolnm,
             forcefield_name, solvent_name, platform_name, restraint_pairs, restraint_k,
             mobile_pairs, frozen_pairs, ligand_resnames, ligand_smiles, charge_method,
             work_dir, sid, ligand_forcefield="gaff-2.11", covalent_anchor="",
             covalent_k=300000.0, covalent_length_nm=0.18,
             covalent_max_distance=3.0):
    xml = [FF_XML[forcefield_name]]
    solvent_xml = SOLVENT_XML[solvent_name]
    if solvent_xml is not None:
        xml.append(solvent_xml)

    charge_method_used = ""

    if ligand_resnames:
        protein_path, ligand_path = _split_ligand(pdb_path, ligand_resnames, work_dir, sid)
        if ligand_path is None:
            raise ValueError(
                f"none of the ligand residues {sorted(ligand_resnames)} are present in "
                f"{os.path.basename(pdb_path)}; the ligand= stream names a molecule this "
                f"structure does not carry")
        from openmmforcefields.generators import SystemGenerator

        mol, charge_method_used = _ligand_molecule(ligand_path, ligand_smiles, charge_method)
        protein_pdb = PDBFile(_fix_protein(protein_path, work_dir, sid))

        forcefield = ForceField(*xml)
        modeller = Modeller(protein_pdb.topology, protein_pdb.positions)
        modeller.addHydrogens(forcefield)
        n_protein_res = sum(1 for _ in modeller.topology.residues())
        # The ligand topology has to come from the OpenFF molecule, not the PDB:
        # template matching is by connectivity, and a PDB HETATM block carries no
        # bonds, so a PDB-derived residue matches nothing.
        modeller.add(mol.to_topology().to_openmm(), _ligand_positions(mol))
        _name_ligand_residues(modeller.topology, n_protein_res,
                              sorted(ligand_resnames)[0])

        covalent_pair = None
        if covalent_anchor:
            covalent_pair = _find_covalent_pair(
                modeller.topology, modeller.positions, covalent_anchor,
                ligand_resnames, covalent_max_distance)
            if covalent_pair is None:
                raise ValueError(
                    f"no {covalent_anchor} atom within {covalent_max_distance} A of the "
                    f"ligand in {os.path.basename(pdb_path)}: there is no covalent link to "
                    f"preserve. Check that the catalytic residue survived upstream design.")

        generator = SystemGenerator(
            forcefields=xml,
            small_molecule_forcefield=ligand_forcefield,
            molecules=[mol],
            forcefield_kwargs={"constraints": HBonds, "rigidWater": False},
        )
        try:
            system = generator.create_system(modeller.topology, molecules=[mol])
        except ValueError as e:
            # GBn2's neck-correction table only covers radii 1-2 A, which excludes
            # silicon; OBC2 has no such table and handles the same ligand.
            if "neck lookup" in str(e) and solvent_name.startswith("implicit-gbn"):
                raise ValueError(
                    f"{solvent_name} cannot solvate this ligand: its neck-correction "
                    f"table has no radius for one of the elements present. Use "
                    f"solvent=\"implicit-obc2\", which parameterises the same molecule.") from e
            raise
        topology = modeller.topology
        if covalent_pair is not None:
            _add_covalent_restraint(system, covalent_pair, covalent_k,
                                    covalent_length_nm)
        if covalent_pair is not None:
            print(f"    covalent restraint: {covalent_anchor}-ligand at "
                  f"{covalent_pair[2]:.2f} A")
            if covalent_pair[2] > 2.5:
                print(f"    WARNING: {covalent_pair[2]:.2f} A is long for a covalent bond; "
                      f"the restraint will pull the ligand toward the anchor",
                      file=sys.stderr)
    else:
        pdb = PDBFile(_fix_protein(pdb_path, work_dir, sid))
        forcefield = ForceField(*xml)
        modeller = Modeller(pdb.topology, pdb.positions)
        modeller.addHydrogens(forcefield)
        system = forcefield.createSystem(
            modeller.topology,
            nonbondedMethod=NoCutoff,
            constraints=HBonds,
        )
        topology = modeller.topology

    _add_restraint(system, modeller, restraint_pairs, restraint_k)
    # Freezing runs after the covalent term is added. A frozen anchor with a mobile
    # ligand is the pivot case and works as intended; the ligand is never frozen, so
    # the restraint can never end up between two immovable atoms.
    n_mobile, n_frozen = _apply_freeze(system, topology, mobile_pairs, frozen_pairs,
                                       ligand_resnames)
    if n_mobile == 0:
        raise ValueError(
            "every atom is frozen — the selection matched nothing mobile. Check that "
            "the selection's chain ids and residue numbers match the structure.")

    integrator = LangevinMiddleIntegrator(
        300 * unit.kelvin, 1.0 / unit.picosecond, 0.002 * unit.picoseconds
    )
    if platform_name and platform_name != "auto":
        sim = Simulation(topology, system, integrator,
                         Platform.getPlatformByName(platform_name))
    else:
        sim = Simulation(topology, system, integrator)
    sim.context.setPositions(modeller.positions)

    e_initial = sim.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    tol = tolerance_kjmolnm * unit.kilojoule_per_mole / unit.nanometer
    sim.minimizeEnergy(tolerance=tol, maxIterations=max_iterations or 0)
    e_final = sim.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)

    positions = sim.context.getState(getPositions=True).getPositions()
    # OpenFF numbers its residues with ints; the PDB writer's keepIds path calls
    # len() on the id and only tolerates strings.
    for res in sim.topology.residues():
        if not isinstance(res.id, str):
            res.id = str(res.id)
    with open(out_path, "w") as f:
        PDBFile.writeFile(sim.topology, positions, f, keepIds=True)
    return e_initial, e_final, n_mobile, n_frozen, charge_method_used


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--structures-json", required=True)
    p.add_argument("--output-dir", required=True)
    p.add_argument("--max-iterations", type=int, default=1000)
    p.add_argument("--tolerance", type=float, default=10.0)
    p.add_argument("--forcefield", default="amber14-all")
    p.add_argument("--solvent", default="implicit-gbn2")
    p.add_argument("--platform", default="auto")
    p.add_argument("--restraint-selection", default="")
    p.add_argument("--restraint-k", type=float, default=1000.0)
    p.add_argument("--mobile-selection", default="")
    p.add_argument("--frozen-selection", default="")
    p.add_argument("--ligand-json", default="")
    p.add_argument("--charge-method", default="am1bcc")
    p.add_argument("--ligand-forcefield", default="gaff-2.11")
    p.add_argument("--covalent-anchor", default="")
    p.add_argument("--covalent-k", type=float, default=300000.0)
    p.add_argument("--covalent-length", type=float, default=0.18)
    p.add_argument("--covalent-max-distance", type=float, default=3.0)
    p.add_argument("--map-csv", required=True)
    p.add_argument("--energies-csv", required=True)
    args = p.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    work_dir = os.path.join(args.output_dir, "_work")
    os.makedirs(work_dir, exist_ok=True)
    ds = load_datastream(args.structures_json)

    restraint_pairs = sele_to_list(args.restraint_selection) if args.restraint_selection else []

    # A TABLE_REFERENCE carries a different selection per structure, so it has to be
    # resolved inside the loop rather than parsed once here.
    mobile_ref = _as_reference(args.mobile_selection)
    frozen_ref = _as_reference(args.frozen_selection)
    mobile_pairs = ([] if mobile_ref else
                    sele_to_list(args.mobile_selection) if args.mobile_selection else [])
    frozen_pairs = ([] if frozen_ref else
                    sele_to_list(args.frozen_selection) if args.frozen_selection else [])

    ligand_resnames, ligand_smiles = set(), None
    if args.ligand_json:
        from biopipelines.ligand_utils import resolve_ligand_codes, resolve_ligand_smiles
        ligand_resnames = set(resolve_ligand_codes(args.ligand_json))
        ligand_smiles = resolve_ligand_smiles(args.ligand_json)
        if not ligand_smiles:
            print("ERROR: ligand stream carries no SMILES; OpenFF cannot parameterise "
                  "a code-only ligand. Pass Ligand(smiles=...).", file=sys.stderr)
            sys.exit(1)
        print(f"Ligand residues: {sorted(ligand_resnames)} charges={args.charge_method}")

    map_rows, energy_rows, failed = [], [], []
    for sid, pdb_path in iterate_files(ds):
        out_path = os.path.join(args.output_dir, f"{sid}.pdb")
        try:
            this_mobile, this_frozen = mobile_pairs, frozen_pairs
            if mobile_ref is not None:
                this_mobile = sele_to_list(mobile_ref.resolve(sid) or "")
                if not this_mobile:
                    # An empty cell would silently minimise the whole structure free,
                    # which is the opposite of what asking for a mobile region means.
                    raise ValueError(
                        f"mobile_selection resolved to nothing for {sid}; nothing would "
                        f"be held and the whole structure would relax")
            if frozen_ref is not None:
                this_frozen = sele_to_list(frozen_ref.resolve(sid) or "")

            e0, e1, n_mob, n_frz, charge_used = minimise(
                pdb_path, out_path, args.max_iterations, args.tolerance,
                args.forcefield, args.solvent, args.platform,
                restraint_pairs, args.restraint_k, this_mobile, this_frozen,
                ligand_resnames, ligand_smiles, args.charge_method, work_dir, sid,
                args.ligand_forcefield, args.covalent_anchor,
                args.covalent_k, args.covalent_length, args.covalent_max_distance)
            map_rows.append({"id": sid, "file": out_path})
            energy_rows.append({
                "id": sid,
                "energy_initial_kj_mol": round(float(e0), 3),
                "energy_final_kj_mol": round(float(e1), 3),
                "delta_kj_mol": round(float(e1 - e0), 3),
                "n_mobile_atoms": n_mob,
                "n_frozen_atoms": n_frz,
                "charge_method": charge_used,
            })
            print(f"  {sid}: E0={e0:.1f} -> E1={e1:.1f} kJ/mol (delta {e1-e0:+.1f}), "
                  f"mobile {n_mob} / frozen {n_frz}")
        except Exception as e:
            print(f"WARNING: {sid} minimisation failed: {e}", file=sys.stderr)
            failed.append(sid)

    os.makedirs(os.path.dirname(args.map_csv), exist_ok=True)
    os.makedirs(os.path.dirname(args.energies_csv), exist_ok=True)
    pd.DataFrame(map_rows, columns=["id", "file"]).to_csv(args.map_csv, index=False)
    pd.DataFrame(energy_rows, columns=E_COLS).to_csv(args.energies_csv, index=False)
    print(f"Map: {args.map_csv} ({len(map_rows)} rows)")
    print(f"Energies: {args.energies_csv} ({len(energy_rows)} rows)")

    if failed:
        print(f"Failed: {len(failed)}/{len(failed)+len(map_rows)}: {failed}", file=sys.stderr)
    if not map_rows:
        sys.exit(1)


if __name__ == "__main__":
    main()
