"""
This pipeline shows how predict to generate mutations at a specific position using the site directed mutagenesis tool
"""

import os, sys
sys.path.insert(0, os.getcwd()) #to see scripts in current folder

from PipelineScripts.pipeline import Pipeline
from PipelineScripts.load_output import LoadOutput
from PipelineScripts.site_directed_mutagenesis import SDM
from PipelineScripts.mmseqs2 import MMseqs2
from PipelineScripts.boltz2 import Boltz2

pp = Pipeline(
    pipeline_name="Examples",
    job_name="SDM-MMseqs-Boltz",
    job_description="Saturation mutagenesis of HaloTag7 at position 175 with Cy7 methyl amide open enantiomer R")

pp.resources(
    gpu="32GB",
    time="24:00:00",
    memory="16GB"
)

#we load this
best_R = pp.add(LoadOutput('/shares/locbp.chem.uzh/public/BioPipelines/Boltz/HT7_Cy7_C_R_001/ToolOutputs/1_Boltz2_output.json'))


single_point_mutants = pp.add(SDM(original=best_R,
                                  position=175,
                                  mode="saturation"))

msas = pp.add(MMseqs2(sequences=single_point_mutants))

boltz2 = pp.add(Boltz2(proteins=single_point_mutants,
                        ligands=best_R,
                        msas=msas))

pp.save()
pp.slurm()