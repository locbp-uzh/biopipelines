# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

# tested: 

from biopipelines.pipeline import *
from biopipelines.rfdiffusion import RFdiffusion
from biopipelines.protein_mpnn import ProteinMPNN
from biopipelines.alphafold import AlphaFold

with Pipeline(project="Examples", job="ContextAutoregistrationDemonstration"):
    Resources(gpu="A100", time="4:00:00", memory="16GB")
    AlphaFold(ProteinMPNN(RFdiffusion(contigs='100-150',num_designs=10), num_sequences=2))    


    


    

