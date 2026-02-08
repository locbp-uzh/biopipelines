# Copyright (c) 2026 Gianluca Quargnali @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# This software is freely available for use, modification, and redistribution.
# If you use this software or any derivative of it in published work,
# you must cite the original author and this repository.
"""
This pipeline should not be modified as it is used by the gpu server
"""

from PipelineScripts.pipeline import *
from PipelineScripts.mmseqs2 import MMseqs2Server

with Pipeline(project="MMseqs2Server",
              job="GPU",
              description="Runs MMseqs2 server with GPU accelaration for 1 against many sequence alignment"):

    Resources(gpu="80GB|96GB",
              time="4:00:00",
              memory="32GB",
              cpus=4)

    MMseqs2Server("gpu")

