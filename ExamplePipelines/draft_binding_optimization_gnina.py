# 10.1111/mmi.12043 Atu4243 9 ± 2 μM

from biopipelines.pipeline import *
from biopipelines.boltz2 import Boltz2
from biopipelines.ligand_mpnn import LigandMPNN
from biopipelines.distance_selector import DistanceSelector
from biopipelines.mutation_profiler import MutationProfiler
from biopipelines.mutation_composer import MutationComposer
from biopipelines.panda import Panda
from biopipelines.gnina import Gnina
from biopipelines.remap import ReMap

# provenence propagation

examples = {
    "Atu4243-GABA": ("4EQ7","GABA"), # no gnina: (0.86,0.88)->(0.60,0.86)
    "Atu4243-TACA": ("4EQ7","TACA"), # no gnina: (1.26,0.79)->(1.08,0.77)
    "NocT-Histopine": ("5OT9","Histopine"), #no gnina: (-0.22,0.88)->(-0.48,0.93)
    "Atu4661-Galactinol": ("6EQ8","galactinol"), #no gnina: (1.33,0.56)->
    "TeaA-hydroxyectoine": ("2VPO","hydroxyectoine"),
}
for example, comp in examples.items():
    if example.startswith("Atu"): continue
    with Pipeline(project="Optimization", job=f"IterativeBinding_{example}_Gnina_ReMap"):
        Resources(gpu="A100", time="24:00:00", memory="16GB")
        protein = PDB(comp[0])
        ligand = Ligand(comp[1])
        original = Boltz2(proteins=protein,
                        ligands=ligand)
        current_best_sequence = protein
        docked_original = Gnina(structures=original,
                                compounds=ligand)
        current_best = Panda(tables=docked_original.tables.docking_summary,
                            operations=[Panda.sort("mean_cnn_affinity", ascending=True),
                                        Panda.head(1)],
                            pool=docked_original)
        for cycle in range(5):
            Suffix(f"Cycle{cycle+1}")
            pocket = DistanceSelector(structures=current_best,
                                    ligand="LIG",
                                    distance=5)
            variants = LigandMPNN(structures=current_best,
                                ligand="LIG",
                                num_sequences=1000,
                                redesigned=pocket.tables.selections.within)
            profile = MutationProfiler(original=current_best_sequence,
                                       mutants=variants)
            candidates = MutationComposer(frequencies=profile.tables.absolute_frequencies,
                                        num_sequences=10,
                                        mode="weighted_random",
                                        max_mutations=3)
            predicted = Boltz2(proteins=candidates,
                            ligands=ligand)
            docked = Gnina(structures=predicted,
                           compounds=ligand)
            # Select best structure for next cycle (pool from docked for structures)
            current_best = Panda(tables=[current_best.tables.result,
                                        docked.tables.docking_summary],
                                operations=[Panda.concat(add_source=True),
                                            Panda.sort("mean_cnn_affinity", ascending=True),
                                            Panda.head(1)],
                                pool=[current_best,
                                      docked])
            docked_sequences = ReMap(candidates, docked)
            current_best_sequence = Panda(
                                tables=[current_best.tables.result,
                                        docked.tables.docking_summary],
                                operations=[Panda.concat(add_source=True),
                                            Panda.sort("mean_cnn_affinity", ascending=True),
                                            Panda.head(1)],
                                pool=[current_best_sequence,
                                      docked_sequences])
            
