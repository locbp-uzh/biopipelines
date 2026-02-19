# 10.1111/mmi.12043 Atu4243 9 ± 2 μM

from biopipelines.pipeline import *
from biopipelines.boltz2 import Boltz2
from biopipelines.ligand_mpnn import LigandMPNN
from biopipelines.distance_selector import DistanceSelector
from biopipelines.mutation_profiler import MutationProfiler
from biopipelines.mutation_composer import MutationComposer
from biopipelines.panda import Panda

with Pipeline(project="Optimization", job="Non-IterativeBinding"):
    Resources(gpu="A100", time="24:00:00", memory="16GB")
    Atu4243 = PDB("4EQ7")
    GABA = Ligand("GABA")
    original = Boltz2(proteins=Atu4243, 
                      ligands=GABA)
    pocket = DistanceSelector(structures=original,
                                ligand="LIG", 
                                distance=5)
    variants = LigandMPNN(structures=original,
                            ligand="LIG",
                            num_sequences=1000,
                            redesigned=pocket.tables.selections.within)
    Suffix("50_composer")
    profile = MutationProfiler(original=original, 
                                mutants=variants)
    candidates = MutationComposer(frequencies=profile.tables.absolute_frequencies,
                                    num_sequences=50, 
                                    mode="weighted_random", 
                                    max_mutations=3)
    predicted = Boltz2(proteins=candidates, 
                        ligands=GABA)
    Suffix("50_LigandMPNN")
    sequences = Panda(tables=variants.tables.sequences,
                      operations=[Panda.head(50)],
                      pool=variants)
    predicted = Boltz2(proteins=candidates, 
                        ligands=GABA)


