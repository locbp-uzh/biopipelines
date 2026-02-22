# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

# tested: 17.2.2026
# pubs.acs.org/doi/10.1021/jacs.5c19875
# 4LCD: add a keep parameter e.g. chain A and resi 100-150 keep chain E
# Conformational change alignment add option CA for resi pymol selection (is it standardized?)

from biopipelines.pipeline import *
from biopipelines.protein_mpnn import ProteinMPNN
from biopipelines.alphafold import AlphaFold
from biopipelines.conformational_change import ConformationalChange
from biopipelines.panda import Panda
from biopipelines.dna_encoder import DNAEncoder


with Pipeline(project="Ubiquitin", job="InverseFolding"):
    Resources(gpu="A100", time="4:00:00", memory="16GB")
    ubiquitin = PDB("4LCD",
                    chain="E")
    sequences = ProteinMPNN(structures=ubiquitin,
                            num_sequences=50,
                            soluble_model=True)
    folded = AlphaFold(proteins=sequences)
    dna = DNAEncoder(sequences=sequences, 
                     organism="EC")
    conf_change = ConformationalChange(reference_structures = ubiquitin,
                                       target_structures = folded)
    filtered_sequences = Panda(tables=[folded.tables.confidence,
                                       conf_change.tables.changes],
                               operations=[Panda.merge(),
                                           Panda.filter("RMSD < 1.5 and plddt > 80")],
                               pool=sequences)
    dna = DNAEncoder(sequences=filtered_sequences, 
                     organism="EC") 









