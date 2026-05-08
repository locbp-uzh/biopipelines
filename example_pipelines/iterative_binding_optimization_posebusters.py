# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

# tested:


from biopipelines.pipeline import *
from biopipelines.boltz2 import Boltz2
from biopipelines.ligand_mpnn import LigandMPNN
from biopipelines.distance_selector import DistanceSelector
from biopipelines.mutation_profiler import MutationProfiler
from biopipelines.mutation_composer import MutationComposer
from biopipelines.panda import Panda
from biopipelines.posebusters import PoseBusters


# Variant of iterative_binding_optimization.py that gates each cycle's
# Boltz2-predicted pose through PoseBusters before the
# affinity-based ranking. Poses that fail physical-validity checks
# (bond lengths, bond angles, clashes, volume overlap with protein, etc.)
# are dropped via Panda.filter("all_pass == True") so they cannot be
# carried forward as the next cycle's best.

with Pipeline(project="NocT", job=f"IterativeBindingOptimization_PoseBusters"):
    Resources(gpu="A100", time="24:00:00", memory="16GB")
    protein = Sequence("5OT9")
    ligand = Ligand("Histopine")
    original = Boltz2(proteins=protein,
                      ligands=ligand)

    current_best = Panda(tables=original.tables.affinity,
                         pool=original)
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
                                      num_sequences=3,
                                      mode="weighted_random",
                                      max_mutations=3)
        predicted = Boltz2(proteins=candidates,
                           ligands=ligand)
        busted = PoseBusters(structures=predicted,
                             ligand="LIG")
        validated = Panda(tables=[predicted.tables.affinity,
                                  busted.tables.posebusters],
                          operations=[Panda.merge(),
                                      Panda.filter("all_pass == True")],
                          pool=predicted)
        current_best = Panda(tables=[current_best.tables.result,
                                     validated.tables.result],
                             operations=[Panda.concat(),
                                         Panda.sort("affinity_probability_binary",
                                                    ascending=False),
                                         Panda.head(1)],
                             pool=[current_best,
                                   validated])
