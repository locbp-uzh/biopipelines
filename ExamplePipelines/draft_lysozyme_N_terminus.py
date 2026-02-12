from biopipelines.pipeline import *
from biopipelines.entities import *
from biopipelines.rfdiffusion import RFdiffusion
from biopipelines.protein_mpnn import ProteinMPNN
from biopipelines.alphafold import AlphaFold
from biopipelines.conformational_change import ConformationalChange
from biopipelines.plot import Plot
from biopipelines.pymol import PyMOL

with Pipeline(project="Lysozyme", job="Redesign"):
    Resources(gpu="A100", time="4:00:00", memory="16GB")

    lysozyme = PDB("168L")

    backbones = RFdiffusion(pdb=lysozyme,
                            contigs='50-70/A81-140',
                            num_designs=10)

    sequences = ProteinMPNN(structures=backbones,
                            num_sequences=5,
                            redesigned=backbones.tables.structures.designed)

    refolded = AlphaFold(proteins=sequences)

    # Validate refolding: compare designed region RMSD between backbone and prediction
    conf_change = ConformationalChange(reference_structures=backbones,
                                       target_structures=refolded,
                                       selection=backbones.tables.structures.designed)

    # Plot confidence metrics and conformational change
    Plot(
        Plot.Scatter(
            data=refolded.tables.confidence,
            x="plddt", y="ptm",
            title="AlphaFold Confidence",
            xlabel="pLDDT", ylabel="pTM", grid=True
        ),
        Plot.Histogram(
            data=refolded.tables.confidence,
            x="plddt", bins=20,
            title="pLDDT Distribution"
        ),
        Plot.Histogram(
            data=conf_change.tables.changes,
            x="RMSD", bins=20,
            title="Designed Region RMSD (backbone vs refolded)"
        ),
    )

    # PyMOL session: non-designed region in white, designed region colored by pLDDT
    PyMOL(
        PyMOL.Load(refolded),
        PyMOL.ColorAF(refolded),
        PyMOL.Color(refolded, selection=backbones.tables.structures.fixed, color="white"),
        PyMOL.Align(),
        session="lysozyme_redesign"
    )
