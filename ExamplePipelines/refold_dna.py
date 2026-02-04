"""
This shows how to refold the filtered results from a BoltzGen run with Boltz (full MSA from MMseqs2 public server)

"""

import os, sys
sys.path.insert(0, os.getcwd()) #to see scripts in current folder

from PipelineScripts.pipeline import *
from PipelineScripts.entities import *
from PipelineScripts.boltz2 import Boltz2
from PipelineScripts.panda import Panda
from PipelineScripts.plot import Plot
from PipelineScripts.dna_encoder import DNAEncoder
from PipelineScripts.conformational_change import ConformationalChange

with Pipeline(project="Examples",
              job="BoltzGen-Refold-DNA",
              description="Refold of filtered designs from a BoltzGen run and DNA output"):

    Resources(gpu="any",
              time="4:00:00",
              memory="16GB")

    # Load structures from previous BoltzGen run
    cif_folder_path = "/shares/locbp.chem.uzh/gquarg/BioPipelines/DeNovo-Gentamicin-Sensor/Gentamicin_BoltzGen_20x500designs_Filtering_001/041_BoltzGenMerge/final_ranked_designs/final_100_designs"
    proteins = PDB(cif_folder_path) #PDB will also extract the sequences
    DNAEncoder(proteins, organism="EC") # Generate DNA sequences

    # Load previous results table (BoltzGen metrics with partial MSA)
    table_path = "/shares/locbp.chem.uzh/gquarg/BioPipelines/DeNovo-Gentamicin-Sensor/Gentamicin_BoltzGen_20x500designs_Filtering_001/041_BoltzGenMerge/final_ranked_designs/final_designs_metrics_100.csv"
    boltzgen_data = Table(table_path)
    ligand = Ligand("gentamicin")

    # Refold with full MSA: apo and holo
    boltz_apo  = Boltz2(proteins=proteins)
    boltz_holo = Boltz2(proteins=proteins,
                        ligands=ligand,
                        msas=boltz_apo)

    # Calculate conformational change between apo and holo
    conf_change = ConformationalChange(reference_structures=boltz_apo,
                                       target_structures=boltz_holo)

    # Merge all metrics: BoltzGen (partial MSA) vs Boltz apo/holo (full MSA)
    merged = Panda(
        tables=[boltzgen_data.tables.data,
                boltz_apo.tables.confidence,
                boltz_holo.tables.confidence,
                boltz_holo.tables.affinity,
                conf_change.tables.conformational_analysis],
        operations=[
            Panda.merge(on="id", prefixes=["", "apo_", "holo_", "holo_", ""]),
            Panda.rank("holo_affinity_pred_value", ascending=True)
        ]
    )

    # Generate comparison plots
    Plot(
        # Scatter: BoltzGen vs Boltz pTM
        Plot.Scatter(
            data=merged.tables.result,
            x="ptm",
            y="holo_ptm",
            title="pTM: BoltzGen (partial MSA) vs Boltz (full MSA)",
            xlabel="BoltzGen pTM",
            ylabel="Boltz Holo pTM",
            grid=True
        ),
        # Scatter: BoltzGen vs Boltz affinity
        Plot.Scatter(
            data=merged.tables.result,
            x="affinity_pred_value",
            y="holo_affinity_pred_value",
            title="Affinity: BoltzGen vs Boltz",
            xlabel="BoltzGen Affinity",
            ylabel="Boltz Holo Affinity",
            grid=True
        ),
        # Scatter: BoltzGen vs Boltz affinity probability
        Plot.Scatter(
            data=merged.tables.result,
            x="affinity_probability_binary",
            y="holo_affinity_probability_binary",
            title="Affinity Probability: BoltzGen vs Boltz",
            xlabel="BoltzGen Affinity Probability",
            ylabel="Boltz Holo Affinity Probability",
            grid=True
        ),
        # Column plot: complex_pLDDT comparison (apo vs holo)
        Plot.Column(
            data=[boltz_apo.tables.confidence, boltz_holo.tables.confidence],
            y="complex_plddt",
            labels=["Apo", "Holo"],
            title="Complex pLDDT: Apo vs Holo",
            ylabel="complex_pLDDT",
            style="column"
        ),
        # Histogram: RMSD distribution
        Plot.Histogram(
            data=conf_change.tables.conformational_analysis,
            x="rmsd",
            bins=20,
            title="Apo-Holo RMSD Distribution",
            xlabel="RMSD (Å)",
            ylabel="Count"
        )
    )

