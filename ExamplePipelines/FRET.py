# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
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

donor = Sequence("VSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTHGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNFNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAA",
                     ids="EBFP")  
acceptor = Sequence("VSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFGYGLQCFARYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSYQSALSKDPNEKRDHMVLLEFVTAA",
                        ids="EYFP") 


with Pipeline(project="Biosensor", job="CaFRET"):
    Resources(gpu="A100", time="8:00:00", memory="16GB")
    donor = Sequence("VSKGEELFTG...", ids="EBFP")     
    cam = PDB("1CFD", ids="CaM")
    acceptor = Sequence("VSKGEELFTG...", ids="EYFP")   
    fusions = Fuse(sequences=[donor, cam, acceptor],
                   name="CaFRET",
                   linker="GSG",
                   linker_lengths=["0-3", "0-3"])
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
    R0 = 35.4  # Forster radius for CFP-YFP pair (Angstrom), assumes kappa2 = 2/3
    derived_metrics = {"FRET_E_apo": f"1 / (1 + (FRET_distance_apo / {R0}) ** 6)",
                       "FRET_E_holo": f"1 / (1 + (FRET_distance_holo / {R0}) ** 6)",
                       "delta_FRET": "abs(FRET_E_holo - FRET_E_apo)"}
    analysis = Panda(tables=[fusions.tables.sequences,
                             dist_apo.tables.distances,
                             dist_holo.tables.distances],
                     operations=[Panda.merge(),
                                 Panda.calculate(derived_metrics)])
    Plot(...)
    Plot(Plot.Bar(data=analysis.tables.result,
                  title="FRET efficiency by Linker Length",
                  x="lengths",
                  y="FRET_E_apo",
                  y_right="FRET_E_holo",
                  xlabel="Linker Lengths",
                  ylabel="FRET apo",
                  ylabel_right="FRET holo"),
         Plot.Bar(data=analysis.tables.result,
                  title="Calcium-Induced FRET Change by Linker Length",
                  x="lengths",
                  y="delta_FRET",
                  xlabel="Linker Lengths",
                  ylabel="FRET difference"))
    
    best = Panda(tables=[analysis.tables.result],
                 operations=[Panda.sort("delta_FRET",ascending=False)])
    
    best_apo = Panda(tables=[best.tables.result],
                     operations=[Panda.head(1)],
                     pool=apo)
    best_holo = Panda(tables=[best.tables.result],
                      operations=[Panda.head(1)],
                      pool=holo)
    PyMOL(PyMOL.Load(best_apo),
          PyMOL.Load(best_holo),
          PyMOL.Align(selection=fusions.tables.sequences.S2),
          PyMOL.Color("white"),
          PyMOL.Color("blue", selection=fusions.tables.sequences.S1),
          PyMOL.Color("pink", selection=fusions.tables.sequences.S2),
          PyMOL.Color("yellow", selection=fusions.tables.sequences.S3))
    

"""
with Pipeline(project="Biosensor", job="CaFRET"):
    Resources(gpu="A100", time="8:00:00", memory="16GB")
    donor = Sequence("VSKGEELFTG...", ids="EBFP")     
    cam = PDB("1CFD", ids="CaM")
    acceptor = Sequence("VSKGEELFTG...", ids="EYFP")   
    fusions = Fuse(sequences=[donor, cam, acceptor],
                   name="CaFRET",
                   linker="GS",
                   linker_lengths=["0-2", "0-2"])
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
    R0 = 35.4  # Forster radius for CFP-YFP pair (Angstrom), assumes kappa2 = 2/3
    derived_metrics = {"FRET_E_apo": f"1 / (1 + (FRET_distance_apo / {R0}) ** 6)",
                       "FRET_E_holo": f"1 / (1 + (FRET_distance_holo / {R0}) ** 6)",
                       "delta_FRET": "FRET_E_holo - FRET_E_apo"}
    analysis = Panda(tables=[fusions.tables.sequences,
                             dist_apo.tables.distances,
                             dist_holo.tables.distances],
                     operations=[Panda.merge(on="id"),
                                 Panda.calculate(derived_metrics)])
    Plot(Plot.Bar(data=analysis.tables.result,
                  title="Calcium-Induced FRET Change by Linker Length",
                  x="lengths",
                  y="FRET_E_apo",
                  y_right="FRET_E_holo",
                  xlabel="Linker Lengths",
                  ylabel="FRET apo",
                  ylabel_right="FRET holo"))


    best = Panda(tables=[analysis.tables.result],
                 operations=[Panda.sort("delta_FRET",ascending=False)])
    best_apo = Panda(tables=[best.tables.result],
                     operations=[Panda.head(1)],
                     pool=apo)
    best_holo = Panda(tables=[best.tables.result],
                      operations=[Panda.head(1)],
                      pool=holo)
    PyMOL(PyMOL.Load(best_apo),
          PyMOL.Load(best_holo),
          PyMOL.Align(selection=fusions.tables.sequences.S2),
          PyMOL.Color("white"),
          PyMOL.Color("blue", selection=fusions.tables.sequences.S1),
          PyMOL.Color("pink", selection=fusions.tables.sequences.S2),
          PyMOL.Color("yellow", selection=fusions.tables.sequences.S3))
"""