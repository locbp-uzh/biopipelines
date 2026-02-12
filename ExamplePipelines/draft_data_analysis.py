from biopipelines.pipeline import *
from biopipelines.entities import *
from biopipelines.panda import Panda
from biopipelines.plot import Plot
from biopipelines.pymol import PyMOL

best_designs = Panda(
    tables=predictions.tables.confidence,
    operations=[
        Panda.filter("plddt > 80"),
        Panda.sort("ptm", ascending=False),
        Panda.head(20)
    ],
    pool=predictions,
    rename="top"
)

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
