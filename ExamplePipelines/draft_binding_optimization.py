from biopipelines.pipeline import *
from biopipelines.entities import *
from biopipelines.boltz2 import Boltz2
from biopipelines.ligand_mpnn import LigandMPNN
from biopipelines.distance_selector import DistanceSelector
from biopipelines.mutation_profiler import MutationProfiler
from biopipelines.mutation_composer import MutationComposer
from biopipelines.panda import Panda

with Pipeline(project="Optimization", job="IterativeBinding"):
    Resources(gpu="A100", time="24:00:00", memory="16GB")
    halotag = PDB("6U32")
    ligand = Ligand("123626542", ids="SiR-CA")
    original = Boltz2(proteins=halotag, 
                      ligands=ligand)
    
    current_best = Panda(tables=original.tables.affinity,
                         operations=[], 
                         pool=original)
    for cycle in range(3):
        Suffix(f"Cycle{cycle+1}")
        pocket = DistanceSelector(structures=current_best,
                                  ligand="LIG", 
                                  distance=5)
        variants = LigandMPNN(structures=current_best,
                              ligand="LIG",
                              num_sequences=1000,
                              redesigned=pocket.tables.selections.within)
        profile = MutationProfiler(original=current_best, 
                                   mutants=variants)
        candidates = MutationComposer(frequencies=profile.tables.absolute_frequencies,
                                      num_sequences=5, 
                                      mode="weighted_random", 
                                      max_mutations=3)
        predicted = Boltz2(proteins=candidates, 
                           ligands=ligand)
        current_best = Panda(tables=[current_best.tables.result, 
                                     predicted.tables.affinity],
                             operations=[Panda.concat(add_source=True),
                                         Panda.sort("affinity_pred_value", ascending=True),
                                         Panda.head(1)],
                             pool=[current_best, predicted])

