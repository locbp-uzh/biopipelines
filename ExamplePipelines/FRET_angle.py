# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

# tested:

from biopipelines.pipeline import *
from biopipelines.fuse import Fuse
from biopipelines.boltz2 import Boltz2
from biopipelines.distance import Distance
from biopipelines.angle import Angle
from biopipelines.panda import Panda
from biopipelines.plot import Plot
from biopipelines.pymol import PyMOL 

with Pipeline(project="Biosensor", job="CaFRET"):
    Resources(gpu="A100", time="8:00:00", memory="16GB")
    donor = Sequence("VSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTHGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNFNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAA",
                     ids="EBFP") 
    cam = PDB("1CFD", 
              ids="CaM")
    acceptor = Sequence("VSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFGYGLQCFARYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSYQSALSKDPNEKRDHMVLLEFVTAA",
                        ids="EYFP") 
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
    orientation_apo = Angle(structures=apo,
                            atoms=(('66.NE1', '66.CA'), ('-173.OH', '-173.CA')),
                            metric_name="orientation_apo",
                            unit="radians")
    orientation_holo = Angle(structures=holo,
                            atoms=(('66.NE1', '66.CA'), ('-173.OH', '-173.CA')),
                            metric_name="orientation_holo",
                            unit="radians")
    R0 = 35.4  # Forster radius for BFP-YFP pair (Angstrom), assumes kappa2 = 2/3
    # kappa2 from the inter-chromophore transition dipole angle (orientation_apo/holo, in rad):
    #   kappa2 = (2/3) * (1 + cos^2(theta))  -- approximation for random azimuthal averaging
    # R0_eff = R0 * (kappa2 / (2/3))^(1/6) rescales R0 by the orientation factor
    derived_metrics = {
        "kappa2_apo":  "(2/3) * (1 + cos(orientation_apo) ** 2)",
        "kappa2_holo": "(2/3) * (1 + cos(orientation_holo) ** 2)",
        "R0_eff_apo":  f"{R0} * (kappa2_apo / (2/3)) ** (1/6)",
        "R0_eff_holo": f"{R0} * (kappa2_holo / (2/3)) ** (1/6)",
        "FRET_E_apo":  "1 / (1 + (FRET_distance_apo / R0_eff_apo) ** 6)",
        "FRET_E_holo": "1 / (1 + (FRET_distance_holo / R0_eff_holo) ** 6)",
        "delta_FRET":  "abs(FRET_E_holo - FRET_E_apo)",
    }
    analysis = Panda(tables=[fusions.tables.sequences,
                             dist_apo.tables.distances,
                             dist_holo.tables.distances,
                             orientation_apo.tables.angles,
                             orientation_holo.tables.angles],
                     operations=[Panda.merge(),
                                 Panda.calculate(derived_metrics)])
    #Plot(...)
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
    PyMOL(PyMOL.Load(best_holo),
          PyMOL.Color("white"),
          PyMOL.Color("blue", selection=fusions.tables.sequences.S1),
          PyMOL.Color("pink", selection=fusions.tables.sequences.S2),
          PyMOL.Color("yellow", selection=fusions.tables.sequences.S3))