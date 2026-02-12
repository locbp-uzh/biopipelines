from biopipelines.pipeline import *
from biopipelines.entities import *
from biopipelines.fuse import Fuse
from biopipelines.alphafold import AlphaFold
from biopipelines.distance import Distance
from biopipelines.angle import Angle
from biopipelines.panda import Panda
from biopipelines.plot import Plot
from biopipelines.pymol import PyMOL

with Pipeline(project="Biosensor", job="LinkerScreen"):
    Resources(gpu="A100", time="8:00:00", memory="16GB")

    # FRET biosensor: mClover3 donor and mRuby3 acceptor, linked by flexible GS linkers
    donor = PDB("5WJ2")    # mClover3
    acceptor = PDB("5WJ4") # mRuby3

    # Screen linker lengths from 5 to 25 residues (GS repeats)
    fusions = Fuse(
        proteins=[donor, acceptor],
        name="FRET_sensor",
        linker="GGGGSGGGGSGGGGSGGGGSGGGGSGGGGS",
        linker_lengths=["5-25"]
    )

    # Predict structure of each fusion variant
    folded = AlphaFold(proteins=fusions)

    # Measure donor-acceptor distance (chromophore-bearing residues)
    # mClover3 Tyr66 to mRuby3 Tyr66 (approximate chromophore positions)
    chromophore_dist = Distance(
        structures=folded,
        residue=["66", "-66"],
        method="min",
        metric_name="FRET_distance"
    )

    # Measure inter-domain orientation angle
    orientation = Angle(
        structures=folded,
        atoms=["66.CA", "1.CA", "-1.CA", "-66.CA"],
        metric_name="domain_orientation"
    )

    # Merge distance and angle data, rank by FRET distance
    combined = Panda(
        tables=[chromophore_dist.tables.distances, orientation.tables.angles],
        operations=[
            Panda.merge(on="id"),
            Panda.sort("FRET_distance", ascending=True),
        ],
        pool=folded
    )

    # Visualize how linker length affects FRET distance and domain orientation
    Plot(
        Plot.Scatter(
            data=combined.tables.result,
            x="FRET_distance", y="domain_orientation",
            title="Linker Length vs FRET Geometry",
            xlabel="Chromophore Distance (A)", ylabel="Domain Orientation (deg)",
            grid=True
        ),
        Plot.Histogram(
            data=folded.tables.confidence,
            x="plddt", bins=20,
            title="AlphaFold Confidence Across Linker Variants"
        ),
    )

    # PyMOL session: color domains, highlight linker
    PyMOL(
        PyMOL.Load(folded),
        PyMOL.ColorAF(folded),
        PyMOL.Color(folded, selection=fusions.tables.sequences.L1, color="yellow"),
        PyMOL.Align(),
        session="FRET_linker_screen"
    )
