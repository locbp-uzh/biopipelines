from biopipelines.pipeline import *
from biopipelines.entities import *
from biopipelines.fuse import Fuse
from biopipelines.boltz2 import Boltz2
from biopipelines.distance import Distance
from biopipelines.angle import Angle
from biopipelines.panda import Panda
from biopipelines.plot import Plot

with Pipeline(project="Biosensor", job="CaFRET"):
    Resources(gpu="A100", time="8:00:00", memory="16GB")
    donor = Sequence("MVSKGEELFTGV...")     # ECFP
    cam = PDB("1CFD")                       # Calmodulin
    acceptor = Sequence("MVSKGEELFTG...")   # EYFP
    fusions = Fuse(proteins=[donor, cam, acceptor],
                   name="CaFRET",
                   linker="GSGAG",
                   linker_lengths=["3-5", "3-5"])
    apo = Boltz2(proteins=fusions)
    calcium = Ligand("CA")
    holo = Boltz2(proteins=fusions,
                  ligands=calcium,
                  msas=apo)
    dist_apo = Distance(structures=apo,
                        residue=["66", "-173"], 
                        metric_name="FRET_distance_apo")
    dist_holo = Distance(structures=holo,
                         residue=["66", "-173"], 
                         metric_name="FRET_distance_holo")
    dihedral_apo = Angle(structures=apo,
                         atoms=["64.NE1", "66.CA", "-173.OH", "-173.CA"],
                         metric_name="domain_orientation_apo")
    dihedral_holo = Angle(structures=holo,
                          atoms=["64.NE1", "66.CA", "-173.OH", "-173.CA"],
                          metric_name="domain_orientation_holo")
    derived_metrics = {"delta_distance": "FRET_distance_holo - FRET_distance_apo",
                       "delta_orientation": "domain_orientation_holo - domain_orientation_apo"}
    analysis = Panda(tables=[dist_apo.tables.distances, 
                             dist_holo.tables.distances,
                             dihedral_apo.tables.angles, 
                             dihedral_holo.tables.angles],
                     operations=[Panda.merge(on="id"),
                                 Panda.calculate(derived_metrics)])
    Plot(Plot.Scatter(data=analysis.tables.result,
                      x="delta_distance", 
                      y="delta_orientation",
                      title="Calcium-Induced FRET Geometry Change",
                      xlabel="Distance Change, Apo to Holo (A)",
                      ylabel="Orientation Change (deg)", 
                      grid=True))
