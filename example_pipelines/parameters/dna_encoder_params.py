# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Minimal DNAEncoder runs exercising every kwarg covered by
tests/tool_parameters/test_dna_encoder_params.py.
"""

from biopipelines.pipeline import *
from biopipelines.dna_encoder import DNAEncoder

with Pipeline(project="ToolParameters",
              job="DNAEncoder",
              description="DNAEncoder parameter coverage — organism codon tables"):

    Resources(gpu=None, time="0:30:00", memory="4GB")

    seq = Sequence("MKTAYIAKQRQISFVKSHFSRQLEERLGL", type="protein", ids="probe")

    # 1: E. coli
    Suffix("EC")
    DNAEncoder(sequences=seq, organism="EC")

    # 2: H. sapiens
    Suffix("HS")
    DNAEncoder(sequences=seq, organism="HS")

    # 3: Combined EC + HS
    Suffix("EC_HS")
    DNAEncoder(sequences=seq, organism="EC&HS")

    # 4: S. cerevisiae
    Suffix("SC")
    DNAEncoder(sequences=seq, organism="SC")
