# Copyright (c) 2026 Gianluca Quargnali @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

# tested:

from biopipelines.pipeline import *
from biopipelines.fuse import Fuse
from biopipelines.boltz2 import Boltz2
from biopipelines.distance import Distance
from biopipelines.panda import Panda
from biopipelines.plot import Plot
from biopipelines.pymol import PyMOL

with Pipeline(project="Biosensor", job="CaFRET"):
    Resources(gpu="A100", time="8:00:00", memory="16GB")
    donor = Sequence("VSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTWGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYISHNVYITADKQKNGIKANFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITLGMDELYK")     # ECFP
    cam = PDB("1CFD") # Calmodulin
    acceptor = Sequence("VSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFGYGLQCFARYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSYQSALSKDPNEKRDHMVLLEFVTAAGITLGMDELYK")   # EYFP
    fusions = Fuse(sequences=[donor, cam, acceptor],
                   name="CaFRET",
                   linker="GSGAG",
                   linker_lengths=["3-5", "3-5"])
    apo = Boltz2(proteins=fusions)
    calcium = Ligand("CA")
    holo = Boltz2(proteins=fusions,
                  ligands=Bundle(calcium,calcium,calcium,calcium),
                  msas=apo)
    dist_apo = Distance(structures=apo,
                        residue=["66", "-173"], 
                        metric_name="FRET_distance_apo")
    dist_holo = Distance(structures=holo,
                         residue=["66", "-173"], 
                         metric_name="FRET_distance_holo")
    R0 = 49.0  # Forster radius for CFP-YFP pair (Angstrom), assumes kappa2 = 2/3
    derived_metrics = {"FRET_E_apo": f"1 / (1 + (FRET_distance_apo / {R0}) ** 6)",
                       "FRET_E_holo": f"1 / (1 + (FRET_distance_holo / {R0}) ** 6)",
                       "delta_FRET": "FRET_E_holo - FRET_E_apo"}
    analysis = Panda(tables=[fusions.tables.sequences,
                             dist_apo.tables.distances,
                             dist_holo.tables.distances],
                     operations=[Panda.merge(on="id"),
                                 Panda.calculate(derived_metrics),
                                 Panda.sort(by="delta_FRET", ascending=False)])
    Plot(Plot.Bar(data=analysis.tables.result,
                  title="Calcium-Induced FRET Change by Linker Length",
                  x="lengths",
                  y="FRET_E_apo",
                  y_right="FRET_E_holo",
                  xlabel="Linker Lengths",
                  ylabel="FRET apo",
                  ylabel_right="FRET holo"))
    best_apo = Panda(tables=[analysis.tables.result],
                     operations=[Panda.head(1)],
                     pool=apo)
    best_holo = Panda(tables=[analysis.tables.result],
                      operations=[Panda.head(1)],
                      pool=holo)
    PyMOL(PyMOL.Load(best_apo),
          PyMOL.Load(best_holo),
          PyMOL.Align(selection=fusions.tables.sequences.S2),
          PyMOL.Color("white"),
          PyMOL.Color("cyan", selection=fusions.tables.sequences.S1),
          PyMOL.Color("pink", selection=fusions.tables.sequences.S2),
          PyMOL.Color("yellow", selection=fusions.tables.sequences.S3))
    