#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""OpenMM minimiser. Per input PDB:
  * Load through PDBFixer-equivalent flow (just openmm.app.PDBFile if input is
    already cleanable; otherwise rely on the forcefield to error out and we
    skip the structure).
  * Build the system from the chosen protein forcefield + implicit solvent (or
    vacuum) so no waters/ions are needed.
  * Optionally restrain a residue selection harmonically.
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


E_COLS = ["id", "energy_initial_kj_mol", "energy_final_kj_mol", "delta_kj_mol"]

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


def minimise(pdb_path, out_path, max_iterations, tolerance_kjmolnm,
             forcefield_name, solvent_name, platform_name, restraint_pairs, restraint_k):
    pdb = PDBFile(pdb_path)
    xml = [FF_XML[forcefield_name]]
    solvent_xml = SOLVENT_XML[solvent_name]
    if solvent_xml is not None:
        xml.append(solvent_xml)
    forcefield = ForceField(*xml)
    modeller = Modeller(pdb.topology, pdb.positions)
    modeller.addHydrogens(forcefield)
    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=NoCutoff,
        constraints=HBonds,
    )
    _add_restraint(system, modeller, restraint_pairs, restraint_k)
    integrator = LangevinMiddleIntegrator(
        300 * unit.kelvin, 1.0 / unit.picosecond, 0.002 * unit.picoseconds
    )
    if platform_name and platform_name != "auto":
        sim = Simulation(modeller.topology, system, integrator,
                         Platform.getPlatformByName(platform_name))
    else:
        sim = Simulation(modeller.topology, system, integrator)
    sim.context.setPositions(modeller.positions)

    e_initial = sim.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    tol = tolerance_kjmolnm * unit.kilojoule_per_mole / unit.nanometer
    sim.minimizeEnergy(tolerance=tol, maxIterations=max_iterations or 0)
    e_final = sim.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)

    positions = sim.context.getState(getPositions=True).getPositions()
    with open(out_path, "w") as f:
        PDBFile.writeFile(sim.topology, positions, f, keepIds=True)
    return e_initial, e_final


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
    p.add_argument("--map-csv", required=True)
    p.add_argument("--energies-csv", required=True)
    args = p.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    ds = load_datastream(args.structures_json)

    restraint_pairs = sele_to_list(args.restraint_selection) if args.restraint_selection else []

    map_rows, energy_rows, failed = [], [], []
    for sid, pdb_path in iterate_files(ds):
        out_path = os.path.join(args.output_dir, f"{sid}.pdb")
        try:
            e0, e1 = minimise(pdb_path, out_path, args.max_iterations, args.tolerance,
                              args.forcefield, args.solvent, args.platform,
                              restraint_pairs, args.restraint_k)
            map_rows.append({"id": sid, "file": out_path})
            energy_rows.append({
                "id": sid,
                "energy_initial_kj_mol": round(float(e0), 3),
                "energy_final_kj_mol": round(float(e1), 3),
                "delta_kj_mol": round(float(e1 - e0), 3),
            })
            print(f"  {sid}: E0={e0:.1f} -> E1={e1:.1f} kJ/mol (delta {e1-e0:+.1f})")
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
