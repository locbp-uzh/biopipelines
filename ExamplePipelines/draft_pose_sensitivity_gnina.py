# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

# Pose sensitivity pipeline comparing Boltz2 affinity predictions with GNINA
# docking scores across saturation mutagenesis at the ABL1 gatekeeper residue.
#
# Workflow:
#   1. Mutagenesis at gatekeeper position (AUTH 315, internal 93)
#   2. Boltz2 predicts each mutant complex with affinity
#   3. GNINA re-docks imatinib into each Boltz2-predicted structure
#   4. Panda merges all metrics into a single table
#   5. Plots: Boltz2 affinity, GNINA scores, cross-method scatter

from biopipelines.pipeline import *
from biopipelines.mutagenesis import Mutagenesis
from biopipelines.boltz2 import Boltz2
from biopipelines.gnina import Gnina
from biopipelines.panda import Panda
from biopipelines.plot import Plot

with Pipeline(project="Imatinib", job="PoseSensitivityGnina"):
    Resources(gpu="A100", time="12:00:00", memory="32GB")
    abl1 = PDB("3QRK")
    imatinib = Ligand("STI")
    original = Boltz2(proteins=abl1,
                      ligands=imatinib)
    mutant_sequences = Mutagenesis(original=abl1,
                                   position=93, #AUTH 315
                                   mode="saturation")
    mutants = Boltz2(proteins=mutant_sequences,
                     ligands=imatinib)
    docking = Gnina(structures=mutants,
                    compounds=imatinib)

    mutants_docking = Panda(tables=[mutants.tables.affinity,
                                    mutants.tables.confidence,
                                    docking.tables.docking_summary],
                            operations=[Panda.merge(on=["id",
                                                        "id",
                                                        "structures.id"])])

    analysis = Panda(tables=[mutant_sequences.tables.sequences,
                             mutants_docking.tables.result],
                     operations=[Panda.merge(on=["id",
                                                 "proteins.id"])])

    # --- Plots ---

    # 1. Boltz2 predicted affinity per mutation
    Plot(Plot.Bar(data=analysis.tables.result,
                  x="mutations",
                  y="affinity_pred_value",
                  title="Boltz2 Predicted Affinity per Gatekeeper Mutant",
                  xlabel="Mutation",
                  ylabel="Boltz2 Affinity (pred)",
                  x_tick_rotation=45,
                  grid=True))

    # 2. GNINA Vina score and CNN score side-by-side per mutation
    Plot(Plot.Bar(data=analysis.tables.result,
                  x="mutations",
                  y="best_vina",
                  y_right="best_cnn_score",
                  title="GNINA Docking Scores per Gatekeeper Mutant",
                  xlabel="Mutation",
                  ylabel="Best Vina Score (kcal/mol)",
                  ylabel_right="CNN Score",
                  x_tick_rotation=45,
                  grid=True))

    # 3. Cross-method scatter: Boltz2 affinity vs GNINA Vina score
    Plot(Plot.Scatter(data=analysis.tables.result,
                      x="affinity_pred_value",
                      y="best_vina",
                      color="mutations",
                      title="Boltz2 Affinity vs GNINA Vina Score",
                      xlabel="Boltz2 Affinity (pred)",
                      ylabel="GNINA Best Vina (kcal/mol)",
                      grid=True,
                      color_legend_outside=True))
