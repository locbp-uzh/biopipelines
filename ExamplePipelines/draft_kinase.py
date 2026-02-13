# Copyright (c) 2026 Gianluca Quargnali @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

# debugged: 13.2.2026

from biopipelines.pipeline import *
from biopipelines.entities import *
from biopipelines.rfdiffusion import RFdiffusion
from biopipelines.protein_mpnn import ProteinMPNN
from biopipelines.alphafold import AlphaFold
from biopipelines.conformational_change import ConformationalChange
from biopipelines.panda import Panda

with Pipeline(project="AdenylateKinase", job="LID_Redesign"):
    Resources(gpu="A100", time="4:00:00", memory="16GB")
    kinase = PDB("4AKE")
    backbones = RFdiffusion(pdb=kinase,
                            contigs='A1-117/50-70/A161-214', 
                            num_designs=10)
    sequences = ProteinMPNN(structures=backbones,
                            num_sequences=2,
                            redesigned=backbones.tables.structures.designed)
    refolded = AlphaFold(proteins=sequences)    
    conf_change = ConformationalChange(reference_structures = kinase,
                                       target_structures = refolded)
    top3 = Panda(tables=[refolded.tables.confidence,
                         conf_change.tables.changes],
                 operations=[Panda.merge(),
                             Panda.filter("RMSD < 1.5 and plddt > 80"),
                             Panda.sort("plddt"),
                             Panda.head(3)],
                 pool=refolded)
    
    

