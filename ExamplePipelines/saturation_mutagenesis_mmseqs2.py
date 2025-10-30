"""
This pipeline shows how predict to generate mutations at a specific position using the site directed mutagenesis tool
"""

import os, sys
sys.path.insert(0, os.getcwd()) #to see scripts in current folder

from PipelineScripts.pipeline import *
from PipelineScripts.load_output import LoadOutput
from PipelineScripts.site_directed_mutagenesis import SDM
from PipelineScripts.mmseqs2 import MMseqs2
from PipelineScripts.boltz2 import Boltz2

with Pipeline(project="Examples",
              job="SDM-MMseqs-Boltz",
              description="Saturation mutagenesis of HaloTag7 at position 175 with Cy7 methyl amide open enantiomer R"):

    Resources(gpu="32GB",
              time="24:00:00",
              memory="16GB")

    best_R = LoadOutput('/shares/locbp.chem.uzh/public/BioPipelines/Boltz/HT7_Cy7_C_R_001/ToolOutputs/1_Boltz2_output.json')

    single_point_mutants = SDM(original=best_R,
                               position=175,
                               mode="saturation")

    msas = MMseqs2(sequences=single_point_mutants)

    boltz2 = Boltz2(proteins=single_point_mutants,
                    ligands=best_R,
                    msas=msas)
