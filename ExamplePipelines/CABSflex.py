# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
This pipeline shows how to determine the flexibility of a protein (RMSF) using a coarse grained simulation with cabs flex.
"""

from biopipelines.pipeline import *
from biopipelines.cabsflex import CABSflex

with Pipeline(project="Examples",
              job="PerResidueFlexibility",
              description="determination of per residue flexibility of a few proteins"):
    Resources(gpu=None,
              time="24:00:00",
              memory="16GB")
    proteins = PDB(["1a2j",
                    "1gfl",
                    "1pii",
                    "2b3p",
                    "2b3q",
                    "5mbn",
                    "1rx4"])
    CABSflex(proteins)