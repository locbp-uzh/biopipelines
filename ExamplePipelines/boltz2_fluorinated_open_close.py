"""
This pipeline shows how improve the difference in predicted binding affinity between open and close form of a carbocyanine 7 chloride and halotag7 starting from a Boltz model of the open form.
"""

#!/usr/bin/env python3
import os, sys
sys.path.insert(0, os.getcwd()) #to see scripts in current folder

from PipelineScripts.pipeline import Pipeline
from PipelineScripts.load_output import LoadOutput
from PipelineScripts.ligand_mpnn import LigandMPNN
from PipelineScripts.mutation_profiler import MutationProfiler
from PipelineScripts.mutation_composer import MutationComposer
from PipelineScripts.mmseqs2 import MMseqs2
from PipelineScripts.boltz2 import Boltz2
from PipelineScripts.residue_atom_distance import ResidueAtomDistance
from PipelineScripts.merge_datasheets import MergeDatasheets
from PipelineScripts.concatenate_datasheets import ConcatenateDatasheets
from PipelineScripts.remove_duplicates import RemoveDuplicates
from PipelineScripts.filter import Filter
from PipelineScripts.select_best import SelectBest
from PipelineScripts.average_by_datasheet import AverageByDatasheet
from PipelineScripts.extract_metrics import ExtractMetrics

pipeline = Pipeline(
    pipeline_name="Boltz", #Will create a folder in /shares/USER/<pipeline_name>
    job_name="LigandMPNN-MutationComposer-MMseqs-Cycle_HT7_Cy7_C_R_Fluorinated", #Unique job folder in /shares/USER/<pipeline_name>/job_name_NNN
    job_description="Folding of best mutant from a cycle with only methy amide, with the fluoroethyl analogues")

pipeline.resources(
    gpu="A100",
    time="24:00:00",
    memory="16GB"
)

"""
We load both open and close form so that we calculate the delta in affinity and use it as benchmark
"""
best = pipeline.add(LoadOutput('/shares/locbp.chem.uzh/gquarg/BioPipelines/LigandMPNN-MutationComposer-MMseqs-Cycle/HT7_Cy7_C_R_003/ToolOutputs/42_SelectBest_output.json'))

msa = pipeline.add(MMseqs2(sequences=best))

#all_analyses=[]

for substitution in ['CCF','CCF2','CCF3']:
    pipeline.set_suffix(substitution)
    open = pipeline.add(LoadOutput(f"/shares/locbp.chem.uzh/public/BioPipelines/Boltz/HT7_Cy7_{substitution}_R_001/ToolOutputs/001_Boltz2_output.json"))
    close = pipeline.add(LoadOutput(f"/shares/locbp.chem.uzh/public/BioPipelines/Boltz/HT7_Cy7_{substitution}_RR_001/ToolOutputs/001_Boltz2_output.json"))

    boltz_holo_open = pipeline.add(Boltz2(proteins=best,
                                    ligands=open,
                                    msas=msa))
    boltz_holo_close = pipeline.add(Boltz2(proteins=best,
                                    ligands=close,
                                    msas=msa))

    open_chlorine_aspartate_distance = pipeline.add(ResidueAtomDistance(input=boltz_holo_open,
                                                                        residue='D in IHDWG',
                                                                        atom='LIG.Cl',
                                                                        metric_name='open_chlorine_distance'))
    open_cap_aspartate_distance = pipeline.add(ResidueAtomDistance(input=boltz_holo_open,
                                                                   residue='D in IHDWG',
                                                                   atom='LIG.N88',
                                                                   metric_name='open_cap_distance'))
    current_analysis = pipeline.add(MergeDatasheets(datasheets=[boltz_holo_open.datasheets.affinity,
                                                        boltz_holo_close.datasheets.affinity,
                                                        open_chlorine_aspartate_distance.datasheets.analysis,
                                                        open_cap_aspartate_distance.datasheets.analysis],
                                            prefixes=["open_","close_","",""],
                                            calculate = {"affinity_delta":"open_affinity_pred_value-close_affinity_pred_value"} ))
    
    #all_analyses.append(current_analysis)
    
#Prints
pipeline.save()
pipeline.slurm() 
