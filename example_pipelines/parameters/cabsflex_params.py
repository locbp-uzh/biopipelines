# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Minimal CABS-Flex runs exercising every kwarg covered by
tests/tool_parameters/test_cabsflex_params.py. Lysozyme (1AKI) is small
enough that even with reduced MC settings each run finishes in minutes.
"""

from biopipelines.pipeline import *
from biopipelines.cabsflex import CABSflex

with Pipeline(project="ToolParameters",
              job="CABSflex",
              description="CABS-Flex parameter coverage — MC knobs, temperature, restraints, output format"):

    Resources(gpu=None, time="6:00:00", memory="8GB")

    target = PDB("1AKI", ids="LYZ", convert="pdb")

    # 1: smaller ensemble + tighter MC
    Suffix("small_mc")
    CABSflex(structures=target, num_models=5, mc_cycles=20, mc_steps=30)

    # 2: annealing + temperature
    Suffix("annealing_temp")
    CABSflex(structures=target, num_models=5, mc_annealing=10, temperature="1.5")

    # 3: filtering count
    Suffix("filtering")
    CABSflex(structures=target, num_models=5, filtering_count=200)

    # 4: all-atom rebuild + medoid output
    Suffix("aa_rebuild")
    CABSflex(structures=target, num_models=5, aa_rebuild=True, pdb_output="M")

    # 5: secondary-structure restraints
    Suffix("restraints")
    CABSflex(
        structures=target,
        num_models=5,
        restraints="ss1",
        restraints_gap=4,
        restraints_min=4.0,
        restraints_max=9.0,
    )
