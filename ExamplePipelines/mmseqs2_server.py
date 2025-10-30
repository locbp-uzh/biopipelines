"""
This pipeline should not be modified as it is used by the gpu server
"""

import os, sys
sys.path.insert(0, os.getcwd()) #to see scripts in current folder

from PipelineScripts.pipeline import *
from PipelineScripts.mmseqs2 import MMseqs2Server

with Pipeline(pipeline_name="MMseqs2Server",
              job_name="GPU",
              job_description="Runs MMseqs2 server with GPU accelaration for 1 against many sequence alignment"):

    Resources(gpu="80GB",
              time="4:00:00",
              memory="32GB",
              cpus=4)

    MMseqs2Server("gpu")

