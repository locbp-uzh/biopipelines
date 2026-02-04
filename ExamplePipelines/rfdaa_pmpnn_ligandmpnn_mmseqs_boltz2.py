"""
This pipeline shows how to run RFdiffusion-AllAtom, generate sequences far from the ligand with ProteinMPNN and close with LigandMPNN, obtain alignments with MMseqs2 and fold the sequences with Boltz2.
"""

import os, sys
sys.path.insert(0, os.getcwd()) #to see scripts in current folder

from PipelineScripts.pipeline import *
from PipelineScripts.entities import * # PDB, Ligand, Sequence
from PipelineScripts.rfdiffusion_allatom import RFdiffusionAllAtom
from PipelineScripts.distance_selector import DistanceSelector
from PipelineScripts.protein_mpnn import ProteinMPNN
from PipelineScripts.ligand_mpnn import LigandMPNN
from PipelineScripts.stitch_sequences import StitchSequences
from PipelineScripts.mmseqs2 import MMseqs2
from PipelineScripts.boltz2 import Boltz2
from PipelineScripts.panda import Panda
from PipelineScripts.plot import Plot
from PipelineScripts.pymol import PyMOL

with Pipeline(project="Examples",
              job="RFDAA-ProteinMPNN-LigandMPNN-MMseqs-Boltz",
              description="redesign of N terminus of rifampicin binding protein"):

    Resources(gpu="any", 
              time="4:00:00",
              memory="16GB")

    rifampicin_binding_protein = PDB("9IAF") # loads from PDB
    rifampicin = Ligand("rifampicin") # loads from pubchem. Also possible: CAS, CCD
    quinolinone = Ligand("A1I1V") # Loads CCF from https://www.rcsb.org/ligand/A1I1V

    rfdaa = RFdiffusionAllAtom(pdb=rifampicin_binding_protein, #can also be a path, preferentially to PDBs folder inside biopipelines folder
                               ligand='RFP', #in rfdaa always specify the ligand name
                               contigs='10-20,A6-140',
                               num_designs=3)

    #this generates a table showing for each structure id a pymol selection for residues within and beyond the distance from the ligand
    distances = DistanceSelector(structures=rfdaa,
                                  ligand="RFP",
                                  distance=5,
                                  restrict_to=rfdaa.tables.structures.designed)

    pmpnn = ProteinMPNN(structures=rfdaa,
                        num_sequences=2,
                        redesigned=distances.tables.selections.beyond)

    lmpnn = LigandMPNN(structures=rfdaa,
                       ligand="RFP", #in ligand mpnn you should always specify the ligand name.
                       num_sequences=2,
                       redesigned=distances.tables.selections.within)

    #stitch sequences by replacing the pmpnn sequences with lmpnn sequences at the "within" positions
    #we'll get 2x2 = 4 sequences per design
    sequences = StitchSequences(template=pmpnn,
                                substitutions={
                                    distances.tables.selections.within: lmpnn
                                })
    """
    rfdiffusion3 gives also sequences in output, so we can write:
    sequences = StitchSequences(template=rfd3,
                                substitutions={
                                    distances.tables.selections.beyond: pmpnn,
                                    distances.tables.selections.within: lmpnn
                                })
    """

    # MMseqs2 has to wait for the server to initialize and provide the msas.
    # This is often not immediate. It is better to switch to non-GPU mode when waiting.
    Resources(gpu=None,
              time="4:00:00",
              memory="16GB")
    msas = MMseqs2(sequences=sequences) 

    # Now we need again a GPU
    # This will be a separate job, dependent on the msa job, dependent on the design job.
    Resources(gpu="any", 
            time="4:00:00",
            memory="16GB")
    boltz_quino = Boltz2(proteins=sequences,
                       ligands=quinolinone,
                        msas=msas) 
    boltz_rifampicin = Boltz2(proteins=sequences,
                        ligands=rifampicin,
                        msas=msas) 
    # we also take the affinity of quinoline in presence of rifampicin and viceversa
    boltz_quino_w_rifampicin = Boltz2(proteins=sequences,
                                      ligands=Bundle(quinolinone,rifampicin), #affinity: quinoline
                                      msas=msas) 
    boltz_rifampicin_w_quino = Boltz2(proteins=sequences,
                                      ligands=Bundle(rifampicin,quinolinone), #affinity: rifampicin
                                      msas=msas) 
    """
    For up to 250-500 sequences (depending on length) one can simply use the public MMseqs2 server while removing the MMseqs2 tool usage:
    boltz_quino = Boltz2(proteins=sequences,
                       ligands=quinolinone) 
    boltz_rifampicin = Boltz2(proteins=sequences,
                        ligands=rifampicin, #ligand smiles taken from original
                        msas=boltz_quino) #recycle the msas
    ...
    """

    # Merge metrics to calculate affinity delta
    merged = Panda(
        tables=[boltz_quino.tables.affinity, 
                boltz_rifampicin.tables.affinity,
                boltz_quino_w_rifampicin.tables.affinity,
                boltz_rifampicin_w_quino.tables.affinity],
        operations=[
            Panda.merge(on="id", prefixes=["quino_", "rif_","quino_w_r_","rif_w_q_"]), # without <metric> -> <prefix><metric>. We don't want overlaps
            Panda.calculate({"affinity_delta": "rif_affinity_pred_value - quino_affinity_pred_value",
                             "affinity_delta_both": "rif_w_affinity_pred_value - quino_w_affinity_pred_value"}),
            Panda.sort("affinity_delta", ascending=True)
        ]
    )

    # Generate plots comparing  metrics
    Plot(
        # Scatter plot:  affinity
        Plot.Scatter(
            data=merged.tables.result,
            x="rif_affinity_pred_value",
            y="quino_affinity_pred_value",
            title="Rif vs Quino Binding Affinity",
            xlabel="Rif Affinity",
            ylabel="Quino Affinity",
            grid=True
        ),
        # Scatter plot:  affinity in presence of other analyte
        Plot.Scatter(
            data=merged.tables.result,
            x="rif_w_affinity_pred_value",
            y="quino_w_affinity_pred_value",
            title="Rif vs Quino Binding Affinity in presence of the other",
            xlabel="Rif Affinity",
            ylabel="Quino Affinity",
            grid=True
        ),
        # Histogram of affinity delta
        Plot.Histogram(
            data=merged.tables.result,
            x="affinity_delta",
            bins=15,
            title="Binding Affinity Delta Distribution",
            xlabel="Affinity Delta (Rif -  Quino)",
            ylabel="Count"
        ),
        # Histogram of affinity delta
        Plot.Histogram(
            data=merged.tables.result,
            x="affinity_delta_both",
            bins=15,
            title="Binding Affinity Delta Distribution",
            xlabel="Affinity Delta (Rif -  Quino) with both present",
            ylabel="Count"
        ),
        # Confidence score distributions
        Plot.Histogram(
            data=boltz_rifampicin.tables.confidence,
            x="confidence_score",
            bins=20,
            title="Quino Structure Confidence Distribution",
            xlabel="Confidence Score",
            ylabel="Count"
        )
    )

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
