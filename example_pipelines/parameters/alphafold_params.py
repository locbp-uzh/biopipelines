# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Minimal AlphaFold runs exercising every kwarg covered by
tests/tool_parameters/test_alphafold_params.py.
"""

from biopipelines.pipeline import *
from biopipelines.alphafold import AlphaFold

with Pipeline(project="ToolParameters",
              job="AlphaFold",
              description="AlphaFold parameter coverage — relax, recycle, seed"):

    Resources(gpu="A100", time="4:00:00", memory="32GB")

    target = Sequence("MKTAYIAKQRQISFVKSHFSRQLEERLGL", ids="short")

    # 1: defaults
    Suffix("defaults")
    AlphaFold(proteins=target)

    # 2: amber relax (num_relax > 0 implies --amber)
    Suffix("relax")
    AlphaFold(proteins=target, num_relax=2)

    # 3: extra recycles
    Suffix("recycle")
    AlphaFold(proteins=target, num_recycle=5)

    # 4: deterministic random seed
    Suffix("seed")
    AlphaFold(proteins=target, rand_seed=2026)

    # 5: full smoke
    Suffix("all")
    AlphaFold(proteins=target, num_relax=1, num_recycle=4, rand_seed=99)
