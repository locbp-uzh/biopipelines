"""
This pipeline shows how to extract metrics from a past run
"""

#!/usr/bin/env python3
import os, sys
sys.path.insert(0, os.getcwd()) #to see scripts in current folder

from PipelineScripts.pipeline import Pipeline
from PipelineScripts.rfdiffusion import RFdiffusion
from PipelineScripts.load_output import LoadOutput
from PipelineScripts.extract_metrics import ExtractMetrics

pipeline = Pipeline(
    pipeline_name="Pipeline", 
    job_name="Test", 
    job_description="Average by datasheets and extract metrics",
    debug=True)

pipeline.resources(
    time="24:00:00",
    memory="16GB"
)

from PipelineScripts.protein_mpnn import ProteinMPNN
rfd = pipeline.add(RFdiffusion(contigs="50-100", num_designs=5))
pmd = pipeline.add(ProteinMPNN(structures=rfd,num_sequences=2))
print(pmd)

pipeline.slurm()