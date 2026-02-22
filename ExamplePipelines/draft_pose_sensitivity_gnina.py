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
#   4. PoseChange compares ligand poses between mutants and wild-type
#   5. Panda merges all metrics into a single table
#   6. Plots: Boltz2 affinity, GNINA scores, cross-method scatter, pose RMSD

from biopipelines.pipeline import *
from biopipelines.mutagenesis import Mutagenesis
from biopipelines.boltz2 import Boltz2
from biopipelines.gnina import Gnina
from biopipelines.pose_change import PoseChange
from biopipelines.panda import Panda
from biopipelines.plot import Plot

with Pipeline(project="Imatinib", job="PoseSensitivityGnina"):
    Resources(gpu="A100", time="12:00:00", memory="32GB")

    abl1 = PDB("3QRK")
    imatinib = Ligand("STI")

    # Predict wild-type complex as pose reference
    original = Boltz2(proteins=abl1,
                      ligands=imatinib,
                      affinity=True)

    # Saturation mutagenesis at gatekeeper position (internal 93, AUTH 315)
    mutant_sequences = Mutagenesis(original=abl1,
                                   position=93,
                                   mode="saturation")

    # Boltz2: predict all mutant complexes with affinity
    mutants = Boltz2(proteins=mutant_sequences,
                     ligands=imatinib,
                     affinity=True)

    # GNINA: re-dock imatinib into each Boltz2-predicted mutant structure.
    # Autobox is derived from the co-predicted ligand coordinates in the PDB.
    docking = Gnina(structures=mutants,
                    compounds=imatinib,
                    exhaustiveness=32,
                    num_runs=5,
                    cnn_scoring="rescore")

    # PoseChange: measure ligand RMSD between each mutant and wild-type pose
    pose = PoseChange(reference_structure=original,
                      sample_structures=mutants,
                      reference_ligand="LIG")

    # Merge Boltz2 affinity + confidence + GNINA best scores + pose RMSD + mutation labels
    boltz_metrics = Panda(tables=[mutants.tables.affinity,
                                  mutants.tables.confidence],
                          operations=[Panda.merge()])

    gnina_best = Panda(tables=docking.tables.docking_summary,
                       operations=[Panda.sort("best_vina", ascending=True),
                                   Panda.drop_duplicates(subset="protein_id", keep="first"),
                                   Panda.select_columns(["protein_id", "best_vina",
                                                         "best_cnn_score", "pose_consistency"])])

    analysis = Panda(tables=[boltz_metrics.tables.result,
                              gnina_best.tables.result,
                              pose.tables.changes,
                              mutant_sequences.tables.sequences],
                     operations=[Panda.merge()])

    # --- Plots ---

    # 1. Boltz2 predicted affinity per mutation
    Plot(Plot.Bar(data=analysis.tables.result,
                  x="mutation",
                  y="affinity_pred_value",
                  title="Boltz2 Predicted Affinity per Gatekeeper Mutant",
                  xlabel="Mutation",
                  ylabel="Boltz2 Affinity (pred)",
                  x_tick_rotation=45,
                  grid=True))

    # 2. GNINA Vina score and CNN score side-by-side per mutation
    Plot(Plot.Bar(data=analysis.tables.result,
                  x="mutation",
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
                      color="mutation",
                      title="Boltz2 Affinity vs GNINA Vina Score",
                      xlabel="Boltz2 Affinity (pred)",
                      ylabel="GNINA Best Vina (kcal/mol)",
                      grid=True,
                      color_legend_outside=True))

    # 4. Ligand RMSD per mutation (pose displacement from wild-type)
    Plot(Plot.Bar(data=analysis.tables.result,
                  x="mutation",
                  y="ligand_rmsd",
                  y_right="pose_consistency",
                  title="Ligand Pose Change per Gatekeeper Mutant",
                  xlabel="Mutation",
                  ylabel="Ligand RMSD vs WT (Å)",
                  ylabel_right="GNINA Pose Consistency",
                  x_tick_rotation=45,
                  grid=True))

    # 5. Boltz2 confidence score per mutation
    Plot(Plot.Bar(data=analysis.tables.result,
                  x="mutation",
                  y="confidence_score",
                  title="Boltz2 Confidence Score per Gatekeeper Mutant",
                  xlabel="Mutation",
                  ylabel="Confidence Score",
                  x_tick_rotation=45,
                  grid=True))
