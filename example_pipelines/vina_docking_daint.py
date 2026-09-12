# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

# tested: Daint (aarch64/GH200), debug partition, 1m22s

"""
AutoDock Vina docking on CSCS Daint.

Vina rather than GNINA because GNINA ships only prebuilt x86-64 binaries, while
conda-forge builds vina, openbabel and rdkit for linux-aarch64. Vina runs on CPU
and has no CNN rescoring, so there are no cnn_* columns in its tables.

Uses PDB 9RTM — the rhodamine-binding protein tag from Chlorobaculum tepidum,
bound to tetramethylrhodamine (CCD: A1EI4).

Submit with:
    BIOPIPELINES_CONFIG_VARIANT=daint ./submit example_pipelines/vina_docking_daint.py
"""

from biopipelines.pipeline import *
from biopipelines import Vina

with Pipeline(project="Examples", job="Vina-Docking",
              description="AutoDock Vina docking on Daint — 9RTM + tetramethylrhodamine"):

    # CPU-only: every Daint node is billed whole, so pack several of these with
    # Parallel(pack=N) rather than giving one docking run a 4-GPU node.
    Resources(gpu="none", partition="debug", time="00:30:00", memory="16GB", cpus=32)

    Vina.install()

    protein = PDB("9RTM", ids="rhotag", convert="pdb")

    tmr = Ligand(
        smiles="CN(C)c1ccc2c(c1)OC1=CC(=[N+](C)C)C=CC1=C2c1ccccc1C(=O)O",
        ids="TMR",
        codes="LIG",
    )

    # No box given: the binding site is recovered from the crystal ligand already
    # in the PDB. Vina's own --autobox has no reference-file form, so the wrapper
    # turns those coordinates into an explicit centre and size.
    Suffix("autobox")
    Vina(structures=protein, compounds=tmr, exhaustiveness=16, num_modes=5)

    # The same pocket stated explicitly, scored with vinardo. Two runs populate
    # the cross-run pose-consistency statistics.
    Suffix("explicit")
    Vina(structures=protein, compounds=tmr, center="22.87,8.65,55.72", size=24.0,
         scoring="vinardo", exhaustiveness=16, num_modes=5, num_runs=2)
