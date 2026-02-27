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
from biopipelines.remap import ReMap


with Pipeline(project="NocT", job=f"IterativeBindingOptimization_Gnina"):
    Resources(gpu="A100", time="24:00:00", memory="16GB")
    protein = PDB("5OT9")
    ligand = Ligand("Histopine")
    original = Boltz2(proteins=protein,
                    ligands=ligand)
    docked_original = Gnina(structures=original,
                            compounds=ligand)
    
    current_best_sequence = original
    current_best_structure = Panda(tables=docked_original.tables.docking_summary,
                        pool=docked_original)
    for cycle in range(5):
        Suffix(f"Cycle{cycle+1}")
        pocket = DistanceSelector(structures=current_best_structure,
                                ligand="LIG",
                                distance=5)
        variants = LigandMPNN(structures=current_best_structure,
                            ligand="LIG",
                            num_sequences=1000,
                            redesigned=pocket.tables.selections.within)
        profile = MutationProfiler(original=current_best_sequence,
                                    mutants=variants)
        candidates = MutationComposer(frequencies=profile.tables.absolute_frequencies,
                                    num_sequences=3,
                                    mode="weighted_random",
                                    max_mutations=3)
        predicted = Boltz2(proteins=candidates,
                        ligands=ligand)
        docked = Gnina(structures=predicted,
                        compounds=ligand)
        best_structure = Panda(tables=[current_best_structure.tables.result,
                                    docked.tables.docking_summary],
                            operations=[Panda.concat(add_source=True),
                                        Panda.sort("mean_cnn_affinity", ascending=True),
                                        Panda.head(1)],
                            pool=[current_best_structure,
                                    docked])
        docked_sequences = ReMap(candidates, docked)
        best_sequence = Panda(tables=[current_best_structure.tables.result,
                                      docked.tables.docking_summary],
                              operations=[Panda.concat(add_source=True),
                                          Panda.sort("mean_cnn_affinity", ascending=True),
                                          Panda.head(1)],
                              pool=[current_best_sequence,
                                    docked_sequences])
        current_best_sequence,current_best_structure=best_sequence,best_structure
