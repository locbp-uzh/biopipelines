# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Minimal Gnina runs exercising every kwarg covered by
tests/tool_parameters/test_gnina_params.py. Reuses the rhodamine-binding
protein 9RTM with TMR (matches example_pipelines/gnina_docking.py).
"""

from biopipelines.pipeline import *
from biopipelines import Gnina

with Pipeline(project="ToolParameters",
              job="Gnina",
              description="Gnina parameter coverage — search/CNN knobs, conformers, pH, box specs"):

    Resources(gpu="A100", time="4:00:00", memory="16GB")

    protein = PDB("9RTM", ids="rhotag", convert="pdb")
    tmr = Ligand(
        smiles="CN(C)c1ccc2c(c1)OC1=CC(=[N+](C)C)C=CC1=C2c1ccccc1C(=O)O",
        ids="TMR",
        codes="LIG",
    )

    # 1: exhaustiveness / num_modes / num_runs / seed
    Suffix("search")
    Gnina(
        structures=protein,
        compounds=tmr,
        exhaustiveness=16,
        num_modes=20,
        num_runs=3,
        seed=7,
    )

    # 2: CNN scoring mode + thresholds
    Suffix("cnn")
    Gnina(
        structures=protein,
        compounds=tmr,
        cnn_scoring="refinement",
        cnn_score_threshold=0.7,
        rmsd_threshold=1.5,
    )

    # 3: conformer generation knobs
    Suffix("conformers")
    Gnina(
        structures=protein,
        compounds=tmr,
        generate_conformers=True,
        num_conformers=12,
        energy_window=3.0,
        conformer_rmsd=0.8,
    )

    # 4: pH-based protonation
    Suffix("ph")
    Gnina(
        structures=protein,
        compounds=tmr,
        protonate=True,
        pH=6.0,
    )

    # 5: explicit autobox padding (larger than default 4Å)
    Suffix("autobox")
    Gnina(
        structures=protein,
        compounds=tmr,
        autobox_add=8.0,
    )
