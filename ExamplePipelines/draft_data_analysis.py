from biopipelines.pipeline import *
from biopipelines.entities import *
from biopipelines.rfdiffusion import RFdiffusion
from biopipelines.protein_mpnn import ProteinMPNN
from biopipelines.alphafold import AlphaFold
from biopipelines.panda import Panda
from biopipelines.plot import Plot
from biopipelines.pymol import PyMOL

with Pipeline(project="Lysozyme", job="Analysis"):
    Resources(gpu="A100", time="4:00:00", memory="16GB")

    lysozyme = PDB("168L")

    backbones = RFdiffusion(pdb=lysozyme,
                            contigs='50-70/A81-140',
                            num_designs=10)

    sequences = ProteinMPNN(structures=backbones,
                            num_sequences=5,
                            redesigned=backbones.tables.structures.designed)

    predictions = AlphaFold(proteins=sequences)

    # Filter by confidence, sort, and keep top designs with their structure files
    best_designs = Panda(
        table=predictions.tables.confidence,
        operations=[
            Panda.filter("plddt > 80"),
            Panda.sort("ptm", ascending=False),
            Panda.head(20)
        ],
        pool=predictions,
        rename="top"
    )

    # Generate publication-ready plots
    Plot(
        Plot.Scatter(data=predictions.tables.confidence,
                     x="plddt", y="ptm",
                     title="AlphaFold Confidence Metrics"),
        Plot.Histogram(data=predictions.tables.confidence,
                       x="plddt", bins=20,
                       title="pLDDT Distribution")
    )

    # Create a PyMOL session with aligned and colored structures
    PyMOL(
        PyMOL.Load(best_designs),
        PyMOL.ColorAF(best_designs),
        PyMOL.Align(),
        session="Top designs"
    )
