#!/usr/bin/env python3

import sys
import os
sys.path.insert(0, os.getcwd())

from PipelineScripts.pipeline import Pipeline
from PipelineScripts.load_output import LoadOutput
from PipelineScripts.merge_datasheets import MergeDatasheets

# Create minimal test to see what datasheets.merged returns
pipeline = Pipeline(
    pipeline_name="Debug-Test",
    job_name="debug", 
    job_description="Debug datasheet types"
)

# Load some data
best_open = pipeline.add(LoadOutput(
    '/shares/locbp.chem.uzh/gquarg/BioPipelines/Boltz/HT_Cy7_C_R_001/ToolOutputs/1_Boltz2_output.json'
))
best_close = pipeline.add(LoadOutput(
    '/shares/locbp.chem.uzh/gquarg/BioPipelines/Boltz/HT_Cy7_C_RR_001/ToolOutputs/1_Boltz2_output.json'
))

# Create analysis
analysis = pipeline.add(MergeDatasheets(
    datasheets=[best_open.output.datasheets.affinity,
               best_close.output.datasheets.affinity],
    prefixes=["open_", "close_"],
    calculate={"affinity_delta": "open_affinity_pred_value - close_affinity_pred_value"}
))

# Check what type the merged datasheet is
merged_ds = analysis.output.datasheets.merged
print(f"Type of analysis.output.datasheets.merged: {type(merged_ds)}")
print(f"merged_ds: {merged_ds}")

if hasattr(merged_ds, '__dict__'):
    print(f"merged_ds.__dict__: {merged_ds.__dict__}")

# Import the types to check
from PipelineScripts.base_config import DatasheetInfo, ToolOutput, StandardizedOutput

print(f"Is ToolOutput? {isinstance(merged_ds, ToolOutput)}")
print(f"Is StandardizedOutput? {isinstance(merged_ds, StandardizedOutput)}")
print(f"Is DatasheetInfo? {isinstance(merged_ds, DatasheetInfo)}")

# Check if it has expected attributes
print(f"Has path? {hasattr(merged_ds, 'path')}")
if hasattr(merged_ds, 'path'):
    print(f"Path: {merged_ds.path}")