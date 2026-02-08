# Copyright (c) 2026 Gianluca Quargnali @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# This software is freely available for use, modification, and redistribution.
# If you use this software or any derivative of it in published work,
# you must cite the original author and this repository.
"""
This pipeline shows how to run RFdiffusion3, generate sequences far from the ligand with ProteinMPNN and close with LigandMPNN, and fold the sequences with Boltz2.
"""

from PipelineScripts.pipeline import *
from PipelineScripts.entities import * # PDB, Ligand, Sequence
from PipelineScripts.rfdiffusion3 import RFdiffusion3
from PipelineScripts.distance_selector import DistanceSelector
from PipelineScripts.protein_mpnn import ProteinMPNN
from PipelineScripts.ligand_mpnn import LigandMPNN
from PipelineScripts.stitch_sequences import StitchSequences
from PipelineScripts.boltz2 import Boltz2
from PipelineScripts.panda import Panda
from PipelineScripts.plot import Plot
from PipelineScripts.pymol import PyMOL

with Pipeline(project="Examples",
              job="RFD3-ProteinMPNN-LigandMPNN-Boltz",
              description="Redesign of a portion of adenilate kinase"):

    Resources(gpu="any", 
              time="4:00:00",
              memory="16GB")

    adenylate_kinase = PDB("3BE4") # ATP + AMP <-> 2ADP
    atp = Ligand("ATP")
    adp = Ligand("ADP")
    amp = Ligand("AMP")
    ap5 = Ligand("AP5")

    adenylate_kinase_boltz = Boltz2(proteins=adenylate_kinase,
                                    ligands=ap5) # the sequence is also extracted with PDB
    
    adenylate_kinase_boltz_renamed = PDB(adenylate_kinase_boltz,
                                         PDB.Rename("LIG",":L:")) # RFdiffusion3 cannot handle ccd-like ligand codes

    rfd3 = RFdiffusion3(pdb=adenylate_kinase_boltz_renamed, #RFdiffusion3 often needs some PDB cleanup. The easiest solution is to start from a Boltz prediction
                        ligand_code=':L:', 
                        contig='A1-121,1-10,A170-214', #They have renamed contigs -> contig
                        num_designs=3)

    #this generates a table showing for each structure id a pymol selection for residues within and beyond the distance from the ligand
    distances = DistanceSelector(structures=rfd3,
                                  ligand=":L:",
                                  distance=5,
                                  restrict_to=rfd3.tables.structures.designed)
    pmpnn = ProteinMPNN(structures=rfd3,
                        num_sequences=2,
                        redesigned=distances.tables.selections.beyond)
    lmpnn = LigandMPNN(structures=rfd3,
                       ligand=":L:", #in ligand mpnn you should always specify the ligand name.
                       num_sequences=2,
                       redesigned=distances.tables.selections.within)

    #stitch sequences by replacing the rfd3 sequences with lmpnn sequences at the "within" positions, and pmpnn at the "beyond" positions
    #we'll get 2x2 = 4 sequences per design
    sequences = StitchSequences(template=rfd3,
                                substitutions={
                                    distances.tables.selections.beyond: pmpnn,
                                    distances.tables.selections.within: lmpnn
                                })

    boltz_apo = Boltz2(proteins=sequences)
    boltz_atp = Boltz2(proteins=sequences,
                       ligands=atp,
                       msas=boltz_apo)
    boltz_amp = Boltz2(proteins=sequences,
                       ligands=amp,
                       msas=boltz_apo)
    boltz_adp = Boltz2(proteins=sequences,
                       ligands=Bundle(adp,adp), #put together two copies of adp, get affinity against the first one
                       msas=boltz_apo)


    # Merge affinity predictions across ligand conditions
    merged = Panda(
        tables=[boltz_atp.tables.affinity,
                boltz_amp.tables.affinity,
                boltz_adp.tables.affinity],
        operations=[
            Panda.merge(on="id", prefixes=["atp_", "amp_", "adp_"]),
            Panda.sort("atp_affinity_pred_value", ascending=True)
        ]
    )

    # Generate plots comparing Boltz2 metrics
    Plot(
        # Scatter: ATP vs AMP affinity
        Plot.Scatter(
            data=merged.tables.result,
            x="atp_affinity_pred_value",
            y="amp_affinity_pred_value",
            title="ATP vs AMP Binding Affinity",
            xlabel="ATP Affinity",
            ylabel="AMP Affinity",
            grid=True
        ),
        # Scatter: ATP vs ADP affinity
        Plot.Scatter(
            data=merged.tables.result,
            x="atp_affinity_pred_value",
            y="adp_affinity_pred_value",
            title="ATP vs ADP Binding Affinity",
            xlabel="ATP Affinity",
            ylabel="ADP Affinity",
            grid=True
        ),
        # Column: confidence across conditions
        Plot.Column(
            data=[boltz_apo.tables.confidence, boltz_atp.tables.confidence,
                  boltz_amp.tables.confidence, boltz_adp.tables.confidence],
            y="complex_plddt",
            labels=["Apo", "ATP", "AMP", "ADP"],
            title="Complex pLDDT Across Conditions",
            ylabel="complex_pLDDT",
            style="column"
        ),
        # Histogram: apo confidence score distribution
        Plot.Histogram(
            data=boltz_apo.tables.confidence,
            x="confidence_score",
            bins=20,
            title="Apo Confidence Score Distribution",
            xlabel="Confidence Score",
            ylabel="Count"
        )
    )

    PyMOL(
        PyMOL.Load(boltz_adp),
        PyMOL.ColorAF(boltz_adp,upper=1),
        session="adenylate_kinase_adp"
    )
