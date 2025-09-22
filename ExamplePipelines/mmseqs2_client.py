"""
This pipeline shows how to improve the difference in predicted binding affinity between open and close form of a carbocyanine 7 chloride and halotag7 starting from a Boltz model of the open form.
"""

#!/usr/bin/env python3
import os, sys
sys.path.insert(0, os.getcwd()) #to see scripts in current folder

from PipelineScripts.pipeline import Pipeline
from PipelineScripts.mmseqs2 import MMseqs2

pipeline = Pipeline(
    pipeline_name="MMseqs2Server", #Will create a folder in /shares/USER/<pipeline_name>
    job_name="CPU", #Unique job folder in /shares/USER/<pipeline_name>/job_name_NNN
    job_description="Test of MMseqs2 local")

pipeline.resources(
    time="24:00:00",
    memory="2GB",
    cpus=32,
    nodes=1
)

pipeline.add(MMseqs2(sequences="MAEIGTGFPFDPHYVEVLGERMHYVDVGPRDGTPVLFLHGNPTSSYVWRNIIPHVAPTHRCIAPDLIGMGKSDKPDLGYFFDDHVRFMDAFIEALGLEEVVLVIHDWGSALGFHWAKRNPERVKGIAFMEFIRPIPTWDEWPEFARETFQAFRTTDVGRKLIIDQNVFIEGTLPMGVVRPLTEVEMDHYREPFLNPVDREPLWRFPNELPIAGEPANIVALVEEYMDWLHQSPVPKLLFWGTPGVLIPPAEAARLAKSLPNCKAVDIGPGLNLLQEDNPDLIGSEIARWLSTLEISG"))

#Prints
pipeline.save()
pipeline.slurm() 
