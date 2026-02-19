# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
This pipeline shows how to run RFdiffusion, generate sequences with ProteinMPNN, and fold the sequences with Boltz2.
"""

from biopipelines.pipeline import *
from biopipelines.rfdiffusion import RFdiffusion
from biopipelines.protein_mpnn import ProteinMPNN
from biopipelines.alphafold import AlphaFold
from biopipelines.panda import Panda
from biopipelines.plot import Plot
from biopipelines.pymol import PyMOL

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
                        num_sequences=2,
                        redesigned=rfd.tables.structures.designed)
    
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
