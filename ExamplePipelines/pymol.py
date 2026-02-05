"""
This pipeline shows how to run RFdiffusion-AllAtom, generate sequences far from the ligand with ProteinMPNN and close with LigandMPNN, obtain alignments with MMseqs2 and fold the sequences with Boltz2.
"""

import os, sys
sys.path.insert(0, os.getcwd()) #to see scripts in current folder

from PipelineScripts.pipeline import *
from PipelineScripts.entities import * # PDB, Ligand, Sequence
from PipelineScripts.load import LoadOutput
from PipelineScripts.pymol import PyMOL

with Pipeline(project="Examples",
              job="PyMOL",
              description="redesign of N terminus of rifampicin binding protein"):

    Resources(gpu="any",
              time="4:00:00",
              memory="16GB")
   
    boltz_rifampicin = LoadOutput("/shares/locbp.chem.uzh/gquarg/BioPipelines/Examples/RFDAA-ProteinMPNN-LigandMPNN-MMseqs-Boltz_006/ToolOutputs/011_Boltz2.json")
    # PyMOL visualization: render holo structures with pLDDT coloring and ligand
    PyMOL(
        PyMOL.RenderEach(
            structures=boltz_rifampicin,
            orient_selection="resn LIG",  # Orient towards ligand
            color_protein="plddt",
            color_ligand="byatom",
            ligand_selection="resn LIG",
            plddt_upper=1,  # Boltz2 uses 0-1 confidence scores
            title="Rifampicin Binder",
            caption="Probability: {affinity_probability_binary:.2f} | Affinity: {affinity_pred_value:.2f}",
            data_table=boltz_rifampicin.tables.affinity,
            width=1920,
            height=1080,
            background="white"
        ),
        session="rifampicin_holo_renders"
    )
