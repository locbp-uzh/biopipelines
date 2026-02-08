# Copyright (c) 2026 Gianluca Quargnali @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# This software is freely available for use, modification, and redistribution.
# If you use this software or any derivative of it in published work,
# you must cite the original author and this repository.
"""
This pipeline shows how to run RFdiffusion-AllAtom, generate sequences far from the ligand with ProteinMPNN and close with LigandMPNN, obtain alignments with MMseqs2 and fold the sequences with Boltz2.
"""

from PipelineScripts.pipeline import *
from PipelineScripts.entities import * # PDB, Ligand, Sequence
from PipelineScripts.rfdiffusion import RFdiffusion
from PipelineScripts.protein_mpnn import ProteinMPNN
from PipelineScripts.alphafold import AlphaFold
from PipelineScripts.panda import Panda
from PipelineScripts.plot import Plot
from PipelineScripts.pymol import PyMOL

with Pipeline(project="Examples",
              job="RFD-ProteinMPNN-AlphaFold2",
              description="Redesign of N terminus domain of lysozyme"):

    Resources(gpu="any", 
              time="4:00:00",
              memory="16GB")

    lysozyme = PDB("168L") # loads from PDB

    rfd = RFdiffusion(pdb=lysozyme,
                        contigs='50-70/A81-140',
                        num_designs=3)

    pmpnn = ProteinMPNN(structures=rfd,
                        num_sequences=2)
    
    af = AlphaFold(proteins=pmpnn)

    # Generate plots of AlphaFold confidence metrics
    Plot(
        # pLDDT vs pTM scatter
        Plot.Scatter(
            data=af.tables.confidence,
            x="plddt",
            y="ptm",
            title="pLDDT vs pTM",
            xlabel="pLDDT",
            ylabel="pTM",
            grid=True
        ),
        # pLDDT distribution
        Plot.Histogram(
            data=af.tables.confidence,
            x="plddt",
            bins=20,
            title="pLDDT Distribution",
            xlabel="pLDDT",
            ylabel="Count"
        ),
        # Max PAE distribution
        Plot.Histogram(
            data=af.tables.confidence,
            x="max_pae",
            bins=20,
            title="Max PAE Distribution",
            xlabel="Max PAE",
            ylabel="Count"
        )
    )

    # PyMOL visualization
    PyMOL(
        PyMOL.Load(af),
        PyMOL.ColorAF(af),
        PyMOL.Align(),
        session="Final results"
    )
