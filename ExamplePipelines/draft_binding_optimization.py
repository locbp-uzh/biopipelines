from biopipelines.pipeline import *
from biopipelines.entities import *
from biopipelines.boltz2 import Boltz2
from biopipelines.ligand_mpnn import LigandMPNN
from biopipelines.distance_selector import DistanceSelector
from biopipelines.mutation_profiler import MutationProfiler
from biopipelines.mutation_composer import MutationComposer
from biopipelines.panda import Panda

with Pipeline(project="Optimization", job="IterativeBinding"):
    Resources(gpu="A100", time="4:00:00", memory="16GB")

    protein = PDB("6U32")
    ligand = Ligand("JF646-HaloTag ligand", ids="JF646")

    current_best = Boltz2(proteins=protein, ligands=ligand) # problem: id of sequence is 6U32 instead of actual id
    current_best = Panda(tables=current_best.tables.affinity,
                         operations=[], pool=current_best)

    for cycle in range(3):
        Suffix(f"Cycle{cycle+1}")

        # Identify residues near the ligand binding site
        pocket = DistanceSelector(structures=current_best,
                                  ligand="LIG", distance=5)

        # Generate 1000 sequence variants near the binding site
        variants = LigandMPNN(structures=current_best,
                              ligand="LIG",
                              num_sequences=1000,
                              redesigned=pocket.tables.selections.within)

        # Analyze mutation frequencies and compose new candidates
        profile = MutationProfiler(original=current_best, 
                                   mutants=variants)
        candidates = MutationComposer(
            frequencies=profile.tables.absolute_frequencies,
            num_sequences=3, mode="weighted_random", max_mutations=3)

        # Predict structures and select the best across all cycles
        predicted = Boltz2(proteins=candidates, ligands=ligand)
        current_best = Panda(
            tables=[current_best.tables.result, predicted.tables.affinity],
            operations=[
                Panda.concat(add_source=True),
                Panda.sort("affinity_pred_value", ascending=True),
                Panda.head(1)
            ],
            pool=[current_best, predicted]
        )
