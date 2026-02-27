# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

from biopipelines.pipeline import *
from biopipelines.boltz2 import Boltz2
from biopipelines.ligand_mpnn import LigandMPNN
from biopipelines.distance_selector import DistanceSelector
from biopipelines.mutation_profiler import MutationProfiler
from biopipelines.mutation_composer import MutationComposer
from biopipelines.panda import Panda
from biopipelines.gnina import Gnina

def drop_duplicates_history(new_sequences, all_sequences_seen):
    if all_sequences_seen is None: # First cycle
        unique_new_sequences = Panda(table=new_sequences.tables.sequences,
                                     operations=[Panda.drop_duplicates(subset="sequence", 
                                                                       keep="first")],
                                     pool=new_sequences)
        all_sequences_updated = Panda(table=unique_new_sequences.tables.result)
    else: # Subsequent cycles - concatenate with history and remove duplicates
        combined = Panda(tables=[new_sequences.tables.sequences, 
                                 all_sequences_seen.tables.result],
                         operations=[Panda.concat(add_source=True)])
        unique_new_sequences = Panda(table=combined.tables.result,
                                     operations=[Panda.drop_duplicates(subset="sequence", 
                                                                       keep="first"),
                                                 Panda.filter("source_table == 0")], # keep only new
                                     pool=new_sequences)
        all_sequences_updated = Panda(table=combined.tables.result,
                                      operations=[Panda.drop_duplicates(subset="sequence", 
                                                                        keep="first")])
    return unique_new_sequences,all_sequences_updated

with Pipeline(project="NocT", job=f"IterativeBindingOptimization"):
    Resources(gpu="A100", time="24:00:00", memory="16GB")
    protein = PDB("5OT9")
    ligand = Ligand("Histopine")
    original = Boltz2(proteins=protein, 
                    ligands=ligand)
    
    current_best = Panda(tables=original.tables.affinity,
                        pool=original)
    all_candidates = None
    for cycle in range(5):
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
                                    num_sequences=10, 
                                    mode="weighted_random", 
                                    max_mutations=3)
        unique_new_candidates,all_candidates = drop_duplicates_history(candidates,all_candidates)
        predicted = Boltz2(proteins=unique_new_candidates, 
                        ligands=ligand)
        current_best = Panda(tables=[current_best.tables.result, 
                                    predicted.tables.affinity],
                            operations=[Panda.concat(add_source=True),
                                        Panda.sort("affinity_pred_value"),
                                        Panda.head(1)],
                            pool=[current_best, 
                                    predicted])


