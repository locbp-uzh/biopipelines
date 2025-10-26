"""
Example usage of SelectionEditor tool.

This example shows how to modify PyMOL selection strings using SelectionEditor
with various operations like expand, shrink, shift, and invert.
"""

import os, sys
sys.path.insert(0, os.getcwd())

from PipelineScripts.pipeline import Pipeline
from PipelineScripts.load_output import LoadOutput
from PipelineScripts.rfdiffusion_allatom import RFdiffusionAllAtom
from PipelineScripts.distance_selector import DistanceSelector
from PipelineScripts.selection_editor import SelectionEditor
from PipelineScripts.ligand_mpnn import LigandMPNN

pipeline = Pipeline(
    pipeline_name="Examples",
    job_name="SelectionEditor-Example",
    job_description="Example showing SelectionEditor tool usage"
)

pipeline.resources(
    gpu="80GB",
    time="12:00:00",
    memory="16GB"
)

# Load a previous structure with ligand
rifampicin = pipeline.add(LoadOutput(
    "/shares/locbp.chem.uzh/public/BioPipelines/Boltz/rifampicin_001/ToolOutputs/1_Boltz2_output.json"
))

# Run RFdiffusion to generate new designs
rfdaa = pipeline.add(RFdiffusionAllAtom(
    pdb=rifampicin,
    ligand='LIG',
    contigs='10-20,A6-140',
    num_designs=2,
    steps=20
))

# Get residues within 5Å of ligand
distances = pipeline.add(DistanceSelector(
    structures=rfdaa,
    ligand="LIG",
    distance=5
))

# Example 1: Expand the selection by 2 residues on each side
# If original selection is "3-45+58-60", this becomes "1-47+56-62"
# (assuming residues 1,2,56,57,61,62 exist in PDB)
expanded = pipeline.add(SelectionEditor(
    selection=distances.datasheets.selections.within,
    expand=2
))

# Example 2: Shrink a selection by 1 residue from each side
# This makes the binding site tighter
shrunk = pipeline.add(SelectionEditor(
    selection=distances.datasheets.selections.within,
    shrink=1
))

# Example 3: Invert a selection to get everything EXCEPT the binding site
# Useful for fixing distant residues while redesigning the binding site
inverted = pipeline.add(SelectionEditor(
    selection=distances.datasheets.selections.within,
    invert=True
))

# Example 4: Shift a selection by +3 residues
# Useful for offsetting selections
shifted = pipeline.add(SelectionEditor(
    selection=distances.datasheets.selections.within,
    shift=3
))

# Example 5: Combine operations - expand then use with LigandMPNN
# Expand binding site selection and use for redesign
binding_site_expanded = pipeline.add(SelectionEditor(
    selection=distances.datasheets.selections.within,
    expand=2
))

# Use the expanded selection for LigandMPNN redesign
lmpnn = pipeline.add(LigandMPNN(
    structures=rfdaa,
    ligand="LIG",
    num_sequences=5,
    redesigned=binding_site_expanded.datasheets.selections.within
))

# Example 6: Auto-detect structures from selection source
# You don't need to provide structures explicitly
auto_detect = pipeline.add(SelectionEditor(
    selection=distances.datasheets.selections.beyond,
    expand=1
    # structures parameter is optional - will be auto-detected from distances
))

print("\n" + "="*60)
print("SelectionEditor Examples")
print("="*60)
print("\nExamples created:")
print("1. Expanded selection (expand=2)")
print("2. Shrunk selection (shrink=1)")
print("3. Inverted selection (invert=True)")
print("4. Shifted selection (shift=3)")
print("5. Expanded binding site for LigandMPNN")
print("6. Auto-detected structures from selection source")
print("\nUsage in LigandMPNN:")
print("  redesigned=expanded.datasheets.selections.within")
print("\nNote: Operations are structure-aware and validate against")
print("      actual PDB residue numbering, merging overlapping ranges.")
print("="*60 + "\n")

# Uncomment to generate and submit SLURM jobs
# pipeline.slurm()
