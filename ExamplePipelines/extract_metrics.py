"""
This pipeline shows how to extract metrics from a past run
"""

#!/usr/bin/env python3
import os, sys
sys.path.insert(0, os.getcwd()) #to see scripts in current folder

from PipelineScripts.pipeline import Pipeline
from PipelineScripts.load_output import LoadOutput
from PipelineScripts.average_by_datasheet import AverageByDatasheet
from PipelineScripts.extract_metrics import ExtractMetrics

pipeline = Pipeline(
    pipeline_name="Calculations", 
    job_name="LigandMPNN-MutationComposer-MMseqs-Cycle_HT7_Cy7_C_R_003", 
    job_description="Average by datasheets and extract metrics")

pipeline.resources(
    time="24:00:00",
    memory="16GB"
)

tool_outputs = '/shares/locbp.chem.uzh/gquarg/BioPipelines/LigandMPNN-MutationComposer-MMseqs-Cycle/HT7_Cy7_C_R_003/ToolOutputs'

original = pipeline.add(LoadOutput('/shares/locbp.chem.uzh/gquarg/BioPipelines/LigandMPNN-MutationComposer-MMseqs-Cycle/HT7_Cy7_C_R_003/ToolOutputs/3_MergeDatasheets_output.json'))
filtered=[(x,pipeline.add(LoadOutput(os.path.join(tool_outputs,x)))) for x in os.listdir(tool_outputs) if "Filter" in x]
filtered_sorted = sorted(filtered,key=lambda x: int(x[0].split('_')[0]))
all_merged = [original.datasheets.merged] + [x[1].datasheets.merged for x in filtered_sorted]


pipeline.add(AverageByDatasheet(all_merged))
metrics = ["affinity_delta",
           "open_affinity_pred_value",
           "close_affinity_pred_value"]
pipeline.add(ExtractMetrics(datasheets=all_merged,
                            metrics=metrics))


pipeline.save()
pipeline.slurm()