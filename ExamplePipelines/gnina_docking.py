# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Example pipeline for GNINA molecular docking.

Uses PDB 9RTM — the rhodamine-binding protein tag (Rho-tag) from Chlorobaculum
tepidum, bound to a tetramethylrhodamine ligand (CCD: A1EI4).

Three test cases:
  1. Basic docking — box auto-detected from crystal ligand in the PDB.
  2. Compound library docking — three rhodamine analogs.
  3. Conformer generation + multi-run docking with full statistical analysis.

Prerequisites:
  - GNINA binary installed (see config.yaml gnina section for path and modules).
  - biopipelines environment with RDKit and OpenBabel.
"""

from biopipelines.pipeline import *
from biopipelines.gnina import Gnina

with Pipeline(project="Examples",
              job="Gnina-Docking",
              description="GNINA docking examples — 9RTM + tetramethylrhodamine"):

    Resources(gpu="A100", time="8:00:00", memory="16GB")

    Gnina.install()

    # Shared protein — fetched once, reused across all three tests
    protein = PDB("9RTM", ids="rhotag")

    # =========================================================================
    # Test 1: Basic docking — box auto-detected from crystal ligand
    # =========================================================================
    # No explicit box needed: GNINA auto-detects from the crystal ligand
    # coordinates extracted during protein preparation. The binding box is
    # centered on the ligand bound to the kept chain.
    Suffix("basic")
    tmr_basic = Ligand(
        smiles="CN(C)c1ccc2c(c1)OC1=CC(=[N+](C)C)C=CC1=C2c1ccccc1C(=O)O",
        ids="TMR",
        codes="LIG"
    )
    docking_basic = Gnina(
        structures=protein,
        compounds=tmr_basic,
        exhaustiveness=32,
        num_runs=5,
        cnn_score_threshold=0.5
    )

    # =========================================================================
    # Test 2: Compound library docking — crystal ligand autobox
    # =========================================================================
    # Docks a small library against the same binding site.
    Suffix("library")
    library = CompoundLibrary({
        "TMR":         "CN(C)c1ccc2c(c1)OC1=CC(=[N+](C)C)C=CC1=C2c1ccccc1C(=O)O",
        "rhodamine_B": "CCN(CC)c1ccc2c(c1)OC1=CC(=[N+](CC)CC)C=CC1=C2c1ccccc1C(=O)O",
        "fluorescein": "OC(=O)c1ccccc1-c1c2ccc(=O)cc-2oc2cc(O)ccc12",
    })
    docking_library = Gnina(
        structures=protein,
        compounds=library,
        autobox_add=4.0,
        exhaustiveness=32,
        num_runs=5,
        cnn_score_threshold=0.4
    )

    # =========================================================================
    # Test 3: Conformer generation with full statistical analysis
    # =========================================================================
    # Generates RDKit conformers, docks each independently across multiple runs,
    # and ranks by combined Vina score + conformer strain energy.
    Suffix("conformers")
    tmr_conf = Ligand(
        smiles="CN(C)c1ccc2c(c1)OC1=CC(=[N+](C)C)C=CC1=C2c1ccccc1C(=O)O",
        ids="TMR",
        codes="LIG"
    )
    docking_conformers = Gnina(
        structures=protein,
        compounds=tmr_conf,
        autobox_add=4.0,
        generate_conformers=True,
        num_conformers=50,
        energy_window=2.0,
        conformer_rmsd=1.0,
        exhaustiveness=32,
        num_runs=5,
        num_modes=9,
        cnn_scoring="rescore",
        cnn_score_threshold=0.5,
        rmsd_threshold=2.0
    )
    # Results:
    #   docking_conformers.tables.docking_results   — all poses with Vina & CNN scores
    #   docking_conformers.tables.docking_summary   — per-conformer aggregated stats across runs
