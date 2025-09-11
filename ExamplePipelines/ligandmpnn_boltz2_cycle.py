###IN PROGRESS
"""
This pipeline shows how improve the difference in predicted binding affinity between open and close form of a carbocyanine 7 chloride and halotag7 starting from a Boltz model of the open form.
"""

#!/usr/bin/env python3
import os, sys
sys.path.insert(0, os.getcwd()) #to see scripts in current folder

from PipelineScripts.pipeline import Pipeline
from PipelineScripts.load_output import LoadOutput
from PipelineScripts.ligand_mpnn import LigandMPNN
from PipelineScripts.boltz2 import Boltz2
from PipelineScripts.residue_atom_distance import ResidueAtomDistance
from PipelineScripts.merge_datasheets import MergeDatasheets
from PipelineScripts.concatenate_datasheets import ConcatenateDatasheets
from PipelineScripts.remove_duplicates import RemoveDuplicates
from PipelineScripts.filter import Filter
from PipelineScripts.select_best import SelectBest

pipeline = Pipeline(
    pipeline_name="LigandMPNN-Boltz-Cycle", #Will create a folder in /shares/USER/<pipeline_name>
    job_name="HT7_Cy7_C_R", #Unique job folder in /shares/USER/<pipeline_name>/job_name_NNN
    job_description="Test on filter based on distance between chlorine and aspartate")

pipeline.resources(
    gpu="V100",
    time="24:00:00",
    memory="16GB"
)

"""
It is best practise to start from a Boltz2 output with the open form, to have a benchmark affinity.
One can then load it with the LoadOutput tool, which will contain the same structures (pdbs), ids, and datasheets as the Boltz2 tool of the past pipeline.  
"""
best_open = pipeline.add(LoadOutput(
    '/shares/locbp.chem.uzh/gquarg/BioPipelines/Boltz/HT_Cy7_C_R_001/ToolOutputs/1_Boltz2_output.json'
    #'path/to/job/ToolOutputs/<Job>_Boltz2_output.json'
))
best_closed = pipeline.add(LoadOutput(
    '/shares/locbp.chem.uzh/gquarg/BioPipelines/Boltz/HT_Cy7_C_RR_001/ToolOutputs/1_Boltz2_output.json'
    #'path/to/job/ToolOutputs/<Job>_Boltz2_output.json'
))

"""
Calculate baseline affinity delta for original structures
"""
original_analysis = pipeline.add(MergeDatasheets(
    datasheets=[best_open.output.datasheets.affinity,
               best_closed.output.datasheets.affinity],
    prefixes=["open_", "close_"],
    calculate={"affinity_delta": "open_affinity_pred_value - close_affinity_pred_value"}
))

NUM_CYCLES = 3
all_sequences_seen = None  # Track all sequences across cycles
previous_analysis = original_analysis  # Start with original baseline for comparison

for CYCLE in range(NUM_CYCLES):
    """
    Diversify with LigandMPNN
    """
    lmpnn = pipeline.add(LigandMPNN(structures=best_open.output, #this is equivalent to boltz2.output
                                    ligand="LIG", #in ligand mpnn you should always specify the ligand name, which is LIG if from Boltz
                                    num_sequences=3, 
                                    redesigned="145-180", #similarly you can specify fixed=...
                                    design_within=4))
    
    """
    Filter NEW sequences only (remove duplicates against historical)
    """
    unique_new_sequences = pipeline.add(RemoveDuplicates(
        pool=lmpnn.output,           # Current cycle sequences  
        history=all_sequences_seen.output if all_sequences_seen else None,  # History from previous cycles (None for first cycle)
        compare="sequence"           # Compare protein sequences
    ))
    
    """
    Update history with unique sequences AFTER deduplication
    """
    if all_sequences_seen is None:
        # First cycle - initialize history with ConcatenateDatasheets (single input)
        all_sequences_seen = pipeline.add(ConcatenateDatasheets(
            datasheets=[unique_new_sequences.output.datasheets.sequences]
        ))
    else:
        # Subsequent cycles - concatenate unique sequences with existing history
        all_sequences_seen = pipeline.add(ConcatenateDatasheets(
            datasheets=[unique_new_sequences.output.datasheets.sequences, 
                       all_sequences_seen.output.datasheets.concatenated]
        ))
    
    """
    Fold only the NEW unique structures
    """
    boltz_apo = pipeline.add(Boltz2(proteins=unique_new_sequences.output))
    boltz_holo_open = pipeline.add(Boltz2(proteins=unique_new_sequences.output,
                                    ligands=r"C[N+]1=C(/C=C/C=C/C=C/C=C(C2(C)C)/N(C)C3=C2C=CC=C3)[C@](CC4=CN(CCOCCOCCCCCCCl)N=N4)(CC(NC)=O)C5=C1C=CC=C5",
                                    msas=boltz_apo.output,
                                    affinity=True))
    boltz_holo_close = pipeline.add(Boltz2(proteins=unique_new_sequences.output,
                                    ligands=r"CN([C@@]1(/C=C/C=C/C=C/C=C(C2(C)C)/N(C)C3=C2C=CC=C3)[C@@]4(CC(N1C)=O)CC5=CN(CCOCCOCCCCCCCl)N=N5)C6=C4C=CC=C6",
                                    msas=boltz_apo.output,
                                    affinity=True))
    
    """
    Calculate distances and analysis (unchanged)
    """
    open_chlorine_aspartate_distance = pipeline.add(ResidueAtomDistance(input=boltz_holo_open.output,
                                                                        residue='D in IHDWG',
                                                                        atom='LIG.Cl',
                                                                        metric_name='open_chlorine_distance'))
    close_chlorine_aspartate_distance = pipeline.add(ResidueAtomDistance(input=boltz_holo_close.output,
                                                                        residue='D in IHDWG',
                                                                        atom='LIG.Cl',
                                                                        metric_name='close_chlorine_distance'))
    current_analysis = pipeline.add(MergeDatasheets(datasheets=[boltz_holo_open.output.datasheets.affinity,
                                                        boltz_holo_close.output.datasheets.affinity,
                                                        open_chlorine_aspartate_distance.output.datasheets.analysis,
                                                        close_chlorine_aspartate_distance.output.datasheets.analysis],
                                            prefixes=["open_","close_","",""],
                                            calculate = {"affinity_delta":"open_affinity_pred_value-close_affinity_pred_value"} ))
    current_filtered = pipeline.add(Filter(pool=boltz_holo_open.output,
                                    data=current_analysis.output,
                                    expression="open_chlorine_distance < 5.0"))
    
    """
    Combine current results with previous best (if not first cycle)
    """
    if CYCLE > 0:
        combined_candidates = pipeline.add(ConcatenateDatasheets(
            datasheets=[current_analysis.output.datasheets.merged,
                       previous_analysis.output.datasheets.merged]
        ))
        selection_pool = combined_candidates.output
        selection_data = combined_candidates.output
    else:
        selection_pool = current_filtered.output
        selection_data = current_analysis.output
    
    """
    Select best for next cycle (FIXED)
    """
    best_open = pipeline.add(SelectBest(
        pool=selection_pool,
        data=selection_data,
        metric='affinity_delta',
        mode='max',
        name=f'{CYCLE+1}_best'
    ))
    
    # Store analysis for next cycle comparison
    previous_analysis = current_analysis

#Prints
pipeline.save()
pipeline.slurm(email="") 
