from biopipelines.pipeline import *
from biopipelines.entities import *
from biopipelines.protein_mpnn import ProteinMPNN
from biopipelines.alphafold import AlphaFold
from biopipelines.conformational_change import ConformationalChange
from biopipelines.panda import Panda
from biopipelines.dna_encoder import DNAEncoder

with Pipeline(project="GFPSensor", job="GeneSynthesis"):
    Resources(gpu="A100", time="4:00:00", memory="16GB")

    caII = PDB("3KS3")

    sequences = ProteinMPNN(structures=caII,
                            fixed="94+96+119",
                            num_sequences=10,
                            soluble_model=True)

    folded = AlphaFold(proteins=sequences)

    dna = DNAEncoder(sequences=sequences, organism="EC")

    conf_change = ConformationalChange(reference_structures = caII,
                                       target_structures = folded)

    filtered = Panda(
        tables=[folded.tables.confidence,conf_change.tables.changes],
        operations=[
            Panda.merge(),
            Panda.filter("RMSD < 1.5 and plddt > 80")
        ],
        pool=folded
    )

    dna = DNAEncoder(sequences=filtered, organism="EC")