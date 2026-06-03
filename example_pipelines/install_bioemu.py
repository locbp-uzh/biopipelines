# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Install BioEmu (Microsoft's generative equilibrium-ensemble emulator).

Creates the bioemu conda env and pip-installs bioemu[cuda,md]. AlphaFold2
weights (~3.5 GB) download lazily to ~/.cache/colabfold on first sampling
run, not at install time.
"""

from biopipelines.pipeline import *
from biopipelines import BioEmu


with Pipeline(project="Setup", job="InstallBioEmu"):
    Resources(time="2:00:00", memory="16GB")
    BioEmu.install()
