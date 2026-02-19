# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
This pipeline shows how to run RFdiffusion-AllAtom, generate sequences far from the ligand with ProteinMPNN and close with LigandMPNN, obtain alignments with MMseqs2 and fold the sequences with Boltz2.
"""

from biopipelines.pipeline import *
from biopipelines.mmseqs2 import MMseqs2
from biopipelines.boltz2 import Boltz2

with Pipeline(project="Examples",
              job="MMseqs-Boltz",
              description="redesign of N terminus of rifampicin binding protein"):
    # MMseqs2 has to wait for the server to initialize and provide the msas.
    # This is often not immediate. It is better to switch to non-GPU mode when waiting.
    Resources(gpu=None,
              time="4:00:00",
              memory="16GB")
    
    sequences = Sequence(["GNSKKHNLILIGAPGSGKGTQCEFIKKEYGLAHLSTGDMLREAIKNGTKIGLEAKSIIESGNFVGDEIVLGLVKEK",
                          "FDLGVCVNGFVLDGFPRTIPQAEGLAKILSEIGDSLTSVIYFEIDDSEIIERISGRCTHPASGRIYHVKYNPPKQPGIDDVTGEPLVWRDDDNAEAVKVRLDVFHKQTAPLVKFYEDLGILKRVNAKLPPKEVTEQIKKIL"])
    msas = MMseqs2(sequences=sequences) # as of now provides only reduced MSA, not properly implemented

    # Now we need a GPU. This will be a separate job, dependent on the msa job
    Resources(gpu="any", 
            time="4:00:00",
            memory="16GB")
    boltz = Boltz2(proteins=sequences,
                        msas=msas) 
    