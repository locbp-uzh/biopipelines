# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Install the 8 new CPU tools added on the package-two branch.

UniProt, RDKit, and OpenBabel reuse the shared `biopipelines` env and need
no install step; the 8 tools below each create their own conda env from
environments/<tool>.yaml.
"""

from biopipelines.pipeline import *
from biopipelines import (
    PLIP,
    Prodigy,
    ProLIF,
    FPocket,
    DSSP,
    APBS,
    OpenMM,
    Reduce,
)


with Pipeline(project="Setup", job="InstallNewCPUTools"):
    Resources(time="2:00:00", memory="8GB")
    Reduce.install()
    DSSP.install()
    FPocket.install()
    OpenMM.install()
    APBS.install()
    Prodigy.install()
    ProLIF.install()
    PLIP.install()
