# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
This shows how to refold the filtered results from a BoltzGen run with Boltz (full MSA from MMseqs2 public server),
plot a metrics comparison between the two, add some analysis steps (apo vs holo),
and obtain codon-optimized sequences for e coli
"""

from biopipelines.pipeline import *
from biopipelines.boltz2 import Boltz2
from biopipelines.panda import Panda
from biopipelines.plot import Plot
from biopipelines.dna_encoder import DNAEncoder
from biopipelines.conformational_change import ConformationalChange

with Pipeline(project="GentamicinBinder", job="BoltzGen-Refold-DNA"):
    Resources(gpu="any", time="24:00:00", memory="16GB")
    final_metrics_csv = "/path/to/final_designs_metrics_100.csv"
    final_sequences = Sequence(final_metrics_csv)
    DNAEncoder(final_sequences, 
               organism="EC")
    
    ligand = Ligand("gentamicin")
    boltz_apo  = Boltz2(proteins=final_sequences)
    boltz_holo = Boltz2(proteins=final_sequences,
                        ligands=ligand,
                        msas=boltz_apo)
    conf_change = ConformationalChange(reference_structures=boltz_apo,
                                       target_structures=boltz_holo)
    boltzgen_data = Table(final_metrics_csv) 
    merged = Panda(tables=[boltzgen_data.tables.data,
                           boltz_apo.tables.confidence,
                           boltz_holo.tables.confidence,
                           boltz_holo.tables.affinity,
                           conf_change.tables.changes],
                   operations=[Panda.merge(on="id", 
                                           prefixes=["boltzgen_", 
                                                     "apo_", 
                                                     "holo_", 
                                                     "holo_", 
                                                     ""]),
                               Panda.rank("holo_affinity_pred_value", 
                                          ascending=True)])
    Plot(Plot.Scatter(data=merged.tables.result,
                      x="boltzgen_affinity_pred_value",
                      y="holo_affinity_pred_value",
                      title="Affinity: BoltzGen vs Boltz",
                      xlabel="BoltzGen Affinity",
                      ylabel="Boltz Holo Affinity",
                      grid=True),
        Plot.Scatter(data=merged.tables.result,
                     x="boltzgen_affinity_probability_binary",
                     y="holo_affinity_probability_binary",
                     title="Affinity Probability: BoltzGen vs Boltz",
                     xlabel="BoltzGen Affinity Probability",
                     ylabel="Boltz Holo Affinity Probability",
                     grid=True),
        Plot.Histogram(data=conf_change.tables.changes,
                       x="RMSD",
                       bins=20,
                       title="Apo-Holo RMSD Distribution",
                       xlabel="RMSD (Å)",
                       ylabel="Count"))

