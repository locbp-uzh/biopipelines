# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Minimal RBSDesigner runs exercising every kwarg covered by
tests/tool_parameters/test_rbs_designer_params.py.
"""

from biopipelines.pipeline import *
from biopipelines.rbs_designer import RBSDesigner

with Pipeline(project="ToolParameters",
              job="RBSDesigner",
              description="RBSDesigner parameter coverage — TIR levels, pre_sequence, start codon"):

    Resources(gpu=None, time="1:00:00", memory="4GB")

    gene = Sequence("ATGAAAGCATCAGCTGCAGCT", type="dna", ids="g1")

    # 1: TIR preset low
    Suffix("low")
    RBSDesigner(sequences=gene, tir="low")

    # 2: TIR preset high
    Suffix("high")
    RBSDesigner(sequences=gene, tir="high")

    # 3: numeric TIR
    Suffix("tir_2500")
    RBSDesigner(sequences=gene, tir=2500)

    # 4: explicit pre_sequence + start codon
    Suffix("pre_seq_start")
    RBSDesigner(sequences=gene, pre_sequence="AAGGAG", add_start_codon=True)
