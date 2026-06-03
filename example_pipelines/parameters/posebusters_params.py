# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Minimal PoseBusters runs exercising every kwarg covered by
tests/tool_parameters/test_posebusters_params.py. PDB 9RTM ships with a
co-crystal LIG (TMR), so PoseBusters can validate it directly without
needing a docking pre-step.
"""

from biopipelines.pipeline import *
from biopipelines import PoseBusters

with Pipeline(project="ToolParameters",
              job="PoseBusters",
              description="PoseBusters parameter coverage — ligand, mode"):

    Resources(gpu=None, time="1:00:00", memory="8GB")

    complex_struct = PDB("9RTM", ids="rhotag", convert="pdb")

    # 1: defaults — dock mode, ligand code from the compounds stream
    Suffix("dock_lig")
    PoseBusters(structures=complex_struct, ligand="LIG", mode="dock")

    # 2: alternative ligand code
    Suffix("hem")
    hbar = PDB("1A3N", ids="HBA")
    PoseBusters(structures=hbar, ligand="HEM", mode="dock")
