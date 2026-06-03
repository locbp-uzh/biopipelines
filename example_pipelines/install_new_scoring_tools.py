# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Install the new ligand-scoring tools added on package-two:
XTB (GFN2-xTB), RTMScore, GEMS.

XTB is a single conda binary; RTMScore/GEMS clone repos + create envs.
"""

from biopipelines.pipeline import *
from biopipelines import XTB, RTMScore, GEMS


with Pipeline(project="Setup", job="InstallNewScoringTools"):
    Resources(time="4:00:00", memory="16GB")
    XTB.install()
    RTMScore.install()
    GEMS.install()
