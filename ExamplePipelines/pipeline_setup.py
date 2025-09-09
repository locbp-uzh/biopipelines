#!/usr/bin/env python3
import os, sys
sys.path.insert(0, os.getcwd()) #to see scripts in current folder

from PipelineScripts.pipeline import Pipeline
from PipelineScripts.ligand_mpnn import LigandMPNN
from PipelineScripts.boltz2 import Boltz2
from PipelineScripts.filter import Filter
from PipelineScripts.residue_atom_distance import ResidueAtomDistance
from PipelineScripts.confidence import Confidence
from PipelineScripts.merge_datasheets import MergeDatasheets
from PipelineScripts.select_best import SelectBest

pipeline = Pipeline(
    pipeline_name="LigandMPNN-Boltz-Cycle", #Will create a folder in /shares/USER/<pipeline_name>
    job_name="HT7_Cy7_ChlorineFilter", #Unique job folder in /shares/USER/<pipeline_name>/job_name_NNN
    job_description="Test on filter based on distance between chlorine and aspartate")

pipeline.resources(
    gpu="V100",
    time="24:00:00",
    memory="16GB"
)

lmpnn = pipeline.add(
    LigandMPNN(structures="HT7_Cy7_CHF2_S_noncov.pdb",
        ligand="HAL",
        num_sequences=3, 
        redesigned="145-180", 
        design_within=4))

boltz_apo = pipeline.add(Boltz2(proteins=lmpnn.output))
boltz_holo = pipeline.add(Boltz2(proteins=lmpnn.output,
                                ligands="CN(C(/C=C/C=C/C=C/C=C(C1(C)C)/N(C)C2=C1C=CC=C2)(N3CC(F)F)C4(CC3=O)CC5=CN(CCOCCOCCCCCCCl)N=N5)C6=C4C=CC=C6",
                                msas=boltz_apo.output.datasheets.msas,
                                affinity=True))
holo_chlorine_aspartate_distance = pipeline.add(ResidueAtomDistance(input=boltz_holo.output,
                                                                    atom='LIG.Cl',
                                                                    residue='protein.D in IGDWG',
                                                                    metric_name='chlorine_distance'))
analysis = pipeline.add(MergeDatasheets(datasheets=[boltz_apo.output.datasheets.affinity,
                                                      boltz_holo.output.datasheets.affinity,
                                                      holo_chlorine_aspartate_distance.output.datasheets.analysis],
                                        prefixes=["apo_","holo_",""],
                                        calculate = {"affinity_difference":"holo_affinity_pred_value-apo_affinity_pred_value"} ))
filtered  = pipeline.add(Filter(input=analysis,expression="chlorine_distance < 5.0"))
best = pipeline.add(SelectBest(filtered,metric="affinity_difference",mode="max"))

lmpnn2 = pipeline.add(
    LigandMPNN(input=best.output,
        ligand="HAL",
        num_sequences=3, 
        redesigned="145-180", 
        design_within=4))

#Prints
pipeline.save()
pipeline.slurm(email="") 


