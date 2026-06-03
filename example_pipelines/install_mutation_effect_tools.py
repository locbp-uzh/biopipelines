# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Install the mutation-effect prediction tools added on package-two:
ThermoMPNN (fold-stability ddG, structure-based) and VespaG (zero-shot
single-substitution fitness, sequence-based).

ThermoMPNN clones the repo (model weights ship with it) + creates a torch env.
VespaG installs from GitHub source (not on PyPI); its ESM-2 3B weights (~11 GB)
download lazily on the first prediction, so run VespaG itself on a GPU node/runtime.
"""

from biopipelines.pipeline import *
from biopipelines import ThermoMPNN, VespaG


with Pipeline(project="Setup", job="InstallMutationEffectTools"):
    Resources(time="4:00:00", memory="16GB")
    ThermoMPNN.install()
    VespaG.install()
