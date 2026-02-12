from biopipelines.pipeline import *
from biopipelines.entities import *
from biopipelines.mutagenesis import Mutagenesis
from biopipelines.boltz2 import Boltz2
from biopipelines.pose_change import PoseChange
from biopipelines.contacts import Contacts
from biopipelines.panda import Panda
from biopipelines.plot import Plot
from biopipelines.pymol import PyMOL

with Pipeline(project="Imatinib", job="PoseSensitivity"):
    Resources(gpu="A100", time="8:00:00", memory="32GB")

    # ABL1 kinase with imatinib — study how gatekeeper mutations affect drug binding
    abl1 = PDB("1IEP")
    imatinib = Ligand("STI")  # STI is the PDB ligand code for imatinib

    # Saturation mutagenesis at the gatekeeper position (Thr315 in ABL1)
    mutants = Mutagenesis(original=abl1,
                          position=315,
                          mode="saturation",
                          include_original=True)

    # Predict all mutant–imatinib complexes
    predictions = Boltz2(proteins=mutants,
                         ligands=imatinib,
                         affinity=True)

    # Measure ligand pose deviation relative to the wild-type crystal structure
    pose = PoseChange(reference_structure=abl1,
                      sample_structures=predictions,
                      reference_ligand="STI",
                      calculate_centroid=True,
                      calculate_orientation=True)

    # Count binding contacts for each mutant
    contacts = Contacts(structures=predictions,
                        ligand="STI",
                        contact_threshold=4.0)

    # Merge pose change, contacts, and affinity into a single analysis table
    analysis = Panda(
        tables=[
            pose.tables.changes,
            contacts.tables.contacts,
            predictions.tables.affinity,
        ],
        operations=[
            Panda.merge(on="id"),
            Panda.sort("ligand_rmsd", ascending=True),
        ],
        pool=predictions
    )

    # Visualize: which mutations preserve the binding pose?
    Plot(
        Plot.Scatter(
            data=analysis.tables.result,
            x="ligand_rmsd", y="contacts",
            title="Pose Deviation vs Binding Contacts",
            xlabel="Ligand RMSD (A)", ylabel="Contacts (<4 A)",
            grid=True
        ),
        Plot.Scatter(
            data=analysis.tables.result,
            x="ligand_rmsd", y="affinity_pred_value",
            title="Pose Deviation vs Predicted Affinity",
            xlabel="Ligand RMSD (A)", ylabel="Predicted Affinity",
            grid=True
        ),
        Plot.Histogram(
            data=analysis.tables.result,
            x="ligand_rmsd", bins=20,
            title="Imatinib Pose RMSD Distribution Across Mutants"
        ),
    )

    # PyMOL session of all predicted complexes
    PyMOL(
        PyMOL.Load(predictions),
        PyMOL.ColorAF(predictions),
        PyMOL.Align(),
        session="abl1_imatinib_mutants"
    )
