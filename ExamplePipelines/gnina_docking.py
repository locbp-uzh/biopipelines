# Copyright (c) 2026 Gianluca Quargnali @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Example pipeline for GNINA molecular docking.

Uses PDB 9RTM — the rhodamine-binding protein tag (Rho-tag) from Chlorobaculum
tepidum, bound to a tetramethylrhodamine ligand (CCD: A1EI4).

Three examples:
  1. Basic docking with explicit box (centered on chain C ligand site)
  2. Autobox from a reference ligand extracted from the crystal structure
  3. Conformer generation + multi-run docking with full statistical analysis

Prerequisites:
  - GNINA binary installed at ~/data/gnina/gnina (see config.yaml for details)
  - CUDA modules configured in config.yaml under the 'gnina' section
  - biopipelines environment with RDKit and OpenBabel
"""

from biopipelines.pipeline import *
from biopipelines.entities import *
from biopipelines.gnina import Gnina

# ============================================================================
# Example 1: Basic docking with explicit binding box
# ============================================================================
# The binding site center is taken from the centroid of ligand A1EI4 in chain C
# of PDB 9RTM: approximately (24.9, 21.3, 18.0). The ligand spans ~12 A, so a
# 25 A box provides generous padding.

with Pipeline(project="Examples",
              job="gnina_basic",
              description="GNINA basic docking — 9RTM + tetramethylrhodamine"):

    Resources(gpu="A100",
              time="2:00:00",
              memory="16GB")

    protein = PDB("9RTM", ids="rhotag")

    # Tetramethylrhodamine (A1EI4) via SMILES
    tmr = Ligand(
        smiles="CN(C)c1ccc2c(c1)OC1=CC(=[N+](C)C)C=CC1=C2c1ccccc1C(=O)O",
        ids="TMR",
        codes="LIG"
    )

    docking = Gnina(
        structures=protein,
        compounds=tmr,
        center="24.9,21.3,18.0",
        size=25.0,
        exhaustiveness=32,
        num_runs=5,
        cnn_score_threshold=0.5
    )

# ============================================================================
# Example 2: Autobox from reference ligand
# ============================================================================
# Instead of manually specifying coordinates, use the crystal-structure ligand
# as the autobox reference. GNINA derives the box from the ligand's extent.

with Pipeline(project="Examples",
              job="gnina_autobox",
              description="GNINA autobox docking — 9RTM + reference ligand"):

    Resources(gpu="A100",
              time="2:00:00",
              memory="16GB")

    protein = PDB("9RTM", ids="rhotag")

    # Fetch the crystallographic ligand to use as autobox reference
    ref_ligand = Ligand("A1EI4", ids="A1EI4_ref")

    # Dock a small library against the same binding site
    library = CompoundLibrary({
        'TMR': 'CN(C)c1ccc2c(c1)OC1=CC(=[N+](C)C)C=CC1=C2c1ccccc1C(=O)O',
        'rhodamine_B': 'CCN(CC)c1ccc2c(c1)OC1=CC(=[N+](CC)CC)C=CC1=C2c1ccccc1C(=O)O',
        'fluorescein': 'OC(=O)c1ccccc1-c1c2ccc(=O)cc-2oc2cc(O)ccc12',
    })

    docking = Gnina(
        structures=protein,
        compounds=library,
        autobox_ligand=ref_ligand,
        autobox_add=4.0,
        exhaustiveness=32,
        num_runs=5,
        cnn_score_threshold=0.4
    )

# ============================================================================
# Example 3: Conformer generation with full statistical analysis
# ============================================================================
# Generate RDKit conformers for each ligand, dock each independently across
# multiple runs, and rank by combined Vina + conformer strain energy.

with Pipeline(project="Examples",
              job="gnina_conformers",
              description="GNINA conformer docking — 9RTM with conformer generation"):

    Resources(gpu="A100",
              time="8:00:00",
              memory="16GB")

    protein = PDB("9RTM", ids="rhotag")

    ref_ligand = Ligand("A1EI4", ids="A1EI4_ref")

    tmr = Ligand(
        smiles="CN(C)c1ccc2c(c1)OC1=CC(=[N+](C)C)C=CC1=C2c1ccccc1C(=O)O",
        ids="TMR",
        codes="LIG"
    )

    docking = Gnina(
        structures=protein,
        compounds=tmr,
        autobox_ligand=ref_ligand,
        autobox_add=4.0,
        generate_conformers=True,
        num_conformers=50,
        energy_window=2.0,
        exhaustiveness=64,
        num_runs=10,
        num_modes=9,
        cnn_scoring="rescore",
        cnn_score_threshold=0.5,
        rmsd_threshold=2.0
    )
    # Results:
    #   docking.tables.docking_results   — all accepted poses with Vina & CNN scores
    #   docking.tables.conformer_ranking — per-conformer stats with pseudo_binding_energy
