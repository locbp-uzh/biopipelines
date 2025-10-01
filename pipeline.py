"""
This pipeline shows how to extract metrics from a past run
"""

#!/usr/bin/env python3
import os, sys
sys.path.insert(0, os.getcwd()) #to see scripts in current folder

from PipelineScripts.pipeline import Pipeline
from PipelineScripts.load_output import LoadOutput
from PipelineScripts.extract_metrics import ExtractMetrics

pipeline = Pipeline(
    pipeline_name="Calculations", 
    job_name="LigandMPNN-MutationComposer-MMseqs-Cycles", 
    job_description="Average by datasheets and extract metrics")

pipeline.resources(
    time="24:00:00",
    memory="16GB"
)

tool_outputs = [f'/shares/locbp.chem.uzh/gquarg/BioPipelines/LigandMPNN-MutationComposer-MMseqs-Cycle/HT7_Cy7_C_R_{x}/ToolOutputs' for x in ['003','008','009','010','011']] + [f'//shares/locbp.chem.uzh/gquarg/BioPipelines/LigandMPNN-MutationComposer-MMseqs-Cycle/HT7_Cy7_C_R_W2_00{i}/ToolOutputs' for i in [1,2,3,4,5]]

original = pipeline.add(LoadOutput('/shares/locbp.chem.uzh/gquarg/BioPipelines/LigandMPNN-MutationComposer-MMseqs-Cycle/HT7_Cy7_C_R_003/ToolOutputs/3_MergeDatasheets_output.json'))

for folder in tool_outputs:
    pipeline.set_suffix(os.path.basename(folder.replace('/ToolOutputs','')))
    selected=[(x,pipeline.add(LoadOutput(os.path.join(folder,x)))) for x in os.listdir(folder) if "Select" in x]
    selected_sorted = sorted(selected,key=lambda x: int(x[0].split('_')[0]))
    all_merged = [original.datasheets.merged] + [x[1].datasheets.selected for x in selected_sorted]
    metrics = ["affinity_delta",
            "open_affinity_pred_value",
            "close_affinity_pred_value"]
    pipeline.add(ExtractMetrics(datasheets=all_merged,
                                metrics=metrics))


pipeline.save()
pipeline.slurm()