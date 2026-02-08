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

