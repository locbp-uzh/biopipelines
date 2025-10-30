"""
This pipeline shows how to extract metrics from a past run
"""

#!/usr/bin/env python3
import os, sys
sys.path.insert(0, os.getcwd()) #to see scripts in current folder

from PipelineScripts.pipeline import *
from PipelineScripts.rfdiffusion import RFdiffusion
from PipelineScripts.protein_mpnn import ProteinMPNN

with Pipeline("MyProject", 
              "Test", 
              "Random modelling"):

    Resources(time="24:00:00",memory="16GB")

    rfd = RFdiffusion(contigs="50-100", num_designs=5)
    pmd = ProteinMPNN(structures=rfd,num_sequences=2)
