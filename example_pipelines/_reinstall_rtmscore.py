# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Force-reinstall RTMScore with the corrected torchdata conda pin."""

from biopipelines.pipeline import *
from biopipelines import RTMScore


with Pipeline(project="Setup", job="ReinstallRTMScore"):
    Resources(time="3:00:00", memory="16GB")
    RTMScore.install(force_reinstall=True)
