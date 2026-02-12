from biopipelines.pipeline import *
from biopipelines.entities import *
from biopipelines.rfdiffusion import RFdiffusion
from biopipelines.protein_mpnn import ProteinMPNN
from biopipelines.alphafold import AlphaFold
from biopipelines.panda import Panda
from biopipelines.dna_encoder import DNAEncoder

with Pipeline(project="GFPSensor", job="GeneSynthesis"):
    Resources(gpu="A100", time="4:00:00", memory="16GB")

    # Split-GFP complementation sensor: redesign the loop connecting strands 7-8
    gfp = PDB("1GFL")

    backbones = RFdiffusion(pdb=gfp,
                            contigs='A1-142/10-20/A157-230',
                            num_designs=10)

    sequences = ProteinMPNN(structures=backbones,
                            num_sequences=10,
                            redesigned=backbones.tables.structures.designed)

    predictions = AlphaFold(proteins=sequences)

    # Select high-confidence designs
    best = Panda(
        table=predictions.tables.confidence,
        operations=[
            Panda.filter("plddt > 85"),
            Panda.sort("ptm", ascending=False),
            Panda.head(5)
        ],
        pool=predictions,
        rename="top"
    )

    # Reverse-translate to E. coli-optimized DNA for gene synthesis ordering
    dna = DNAEncoder(sequences=best, organism="EC")
