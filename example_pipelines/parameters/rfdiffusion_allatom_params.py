# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Minimal RFdiffusion-AllAtom runs exercising every kwarg covered by
tests/tool_parameters/test_rfdiffusion_allatom_params.py.

Uses 3QRK (ABL1 kinase + 9DP), the same scaffold/ligand pair as
example_pipelines/test_all.py.
"""

from biopipelines.pipeline import *
from biopipelines.rfdiffusion_allatom import RFdiffusionAllAtom

with Pipeline(project="ToolParameters",
              job="RFdiffusionAllAtom",
              description="RFdiffusion-AllAtom parameter coverage — ligand, num_designs, steps, recycles, deterministic"):

    Resources(gpu="A100", time="3:00:00", memory="16GB")

    abl1 = PDB("3QRK", ids="ABL1")
    contigs = "A227-377,10-20"

    # 1: defaults — 9DP ligand context
    Suffix("default")
    RFdiffusionAllAtom(pdb=abl1, ligand="9DP", contigs=contigs, num_designs=2, steps=20)

    # 2: more designs
    Suffix("num_designs")
    RFdiffusionAllAtom(pdb=abl1, ligand="9DP", contigs=contigs, num_designs=4, steps=20)

    # 3: explicit step count
    Suffix("steps")
    RFdiffusionAllAtom(pdb=abl1, ligand="9DP", contigs=contigs, num_designs=2, steps=80)

    # 4: extra recycles
    Suffix("recycles")
    RFdiffusionAllAtom(pdb=abl1, ligand="9DP", contigs=contigs, num_designs=2, steps=20, num_recycles=3)

    # 5: deterministic + design_startnum
    Suffix("deterministic_startnum")
    RFdiffusionAllAtom(
        pdb=abl1,
        ligand="9DP",
        contigs=contigs,
        num_designs=2,
        steps=20,
        deterministic=True,
        design_startnum=5,
    )
