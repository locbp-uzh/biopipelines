"""
This pipeline shows how to run calculations on a past cycle run containing lots of outputs in case you forgot or want to add something.
The idea is to load tool outputs based on their name, order them based on the execution order, and feed the datasheets into analysis tools.
"""

import os, sys
sys.path.insert(0, os.getcwd()) #to see scripts in current folder

from PipelineScripts.pipeline import Pipeline
from PipelineScripts.load_output import LoadOutput
from PipelineScripts.average_by_datasheet import AverageByDatasheet
from PipelineScripts.extract_metrics import ExtractMetrics

pipeline = Pipeline(
    pipeline_name="PostProcessing", 
    job_name="Cycle_HT7_Cy7_C_R_003", 
    job_description="Average by datasheets and extract metrics")

pipeline.resources(
    time="24:00:00",
    memory="16GB"
)

tool_outputs = '/shares/locbp.chem.uzh/gquarg/BioPipelines/LigandMPNN-MutationComposer-MMseqs-Cycle/HT7_Cy7_C_R_003/ToolOutputs'

original = pipeline.add(LoadOutput(tool_outputs+'/3_MergeDatasheets_output.json'))
# a filter was applied at each cycle. we can load all the Filter outputs in the folder .../ToolOutputs
# filtered now contains tuples (str,ToolOutput) e.g. (008_Filter.json, <output of 008_Filter>)
filtered=[(x,pipeline.add(LoadOutput(os.path.join(tool_outputs,x)))) for x in os.listdir(tool_outputs) if "Filter" in x]
# we order based on NNN extracted from the first element of the tuple
filtered_sorted = sorted(filtered,key=lambda x: int(x[0].split('_')[0])) # we sort them base on NNN_... (they are not loaded in order)
# we analyse based on 2nd element of the tuple
all_merged = [original.datasheets.merged] + [x[1].datasheets.merged for x in filtered_sorted]

pipeline.add(AverageByDatasheet(all_merged))
metrics = ["affinity_delta",
           "open_affinity_pred_value",
           "close_affinity_pred_value"]
pipeline.add(ExtractMetrics(datasheets=all_merged,
                            metrics=metrics))

pipeline.save()
pipeline.slurm()