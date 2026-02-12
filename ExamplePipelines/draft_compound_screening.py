from biopipelines.pipeline import *
from biopipelines.entities import *
from biopipelines.boltz2 import Boltz2
from biopipelines.panda import Panda
from biopipelines.plot import Plot
from biopipelines.pymol import PyMOL

with Pipeline(project="Screening", job="FragmentScreen"):
    Resources(gpu="A100", time="8:00:00", memory="32GB")

    target = Sequence("MDPLNLS...", ids="DRD2")
    cofactor = Ligand("ATP")

    library = CompoundLibrary(
        library={
            "candidate": "<aryl><linker><amide>",
            "aryl": ["<methoxyphenyl>", "<fluorophenyl>"],
            "methoxyphenyl": r"c1cc(OC)cc(c1)",
            "fluorophenyl": r"c1cc(F)ccc1",
            "linker": ["CC", "CCC", "CCCC"],
            "amide": ["C(=O)N", "C(=O)NC"],
        },
        primary_key="candidate"
    )
    # Generates 2 x 3 x 2 = 12 compounds with fragment tracking

    boltz = Boltz2(
        proteins=target,
        ligands=Bundle(Each(library), cofactor),
        affinity=True
    )
    # 12 predictions, each compound co-predicted with the ATP cofactor

    # Select top hits by predicted affinity
    top_hits = Panda(
        table=boltz.tables.affinity,
        operations=[
            Panda.sort("affinity_pred_value", ascending=True),
            Panda.head(5)
        ],
        pool=boltz,
        rename="top"
    )

    # Plot affinity and confidence distributions
    Plot(
        Plot.Histogram(
            data=boltz.tables.affinity,
            x="affinity_pred_value", bins=12,
            title="Predicted Binding Affinity Distribution"
        ),
        Plot.Scatter(
            data=boltz.tables.confidence,
            x="complex_plddt", y="ptm",
            title="Boltz-2 Confidence",
            xlabel="Complex pLDDT", ylabel="pTM", grid=True
        ),
    )

    # PyMOL session with top hits
    PyMOL(
        PyMOL.Load(top_hits),
        PyMOL.ColorAF(top_hits),
        PyMOL.Align(),
        session="fragment_screen_top_hits"
    )
