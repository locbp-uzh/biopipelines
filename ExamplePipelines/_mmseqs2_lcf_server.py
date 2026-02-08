# Copyright (c) 2026 Gianluca Quargnali @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# This software is freely available for use, modification, and redistribution.
# If you use this software or any derivative of it in published work,
# you must cite the original author and this repository.
"""
This pipeline should not be modified as it is used by the LCF gpu server.

MMseqs2 LCF (LocalColabFold) server uses colabfold_search for MSA generation,
searching both UniRef30 and ColabFoldDB for richer alignments.
"""

from PipelineScripts.pipeline import *
from PipelineScripts.mmseqs2_lcf import MMseqs2ServerLCF

with Pipeline(project="MMseqs2LCFServer",
              job="GPU",
              description="Runs MMseqs2 LCF server with GPU acceleration using colabfold_search for UniRef30 + ColabFoldDB"):

    Resources(gpu="80GB|96GB",
              time="4:00:00",
              memory="32GB",
              cpus=16)

    MMseqs2ServerLCF()
