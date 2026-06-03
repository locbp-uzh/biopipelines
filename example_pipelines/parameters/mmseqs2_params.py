# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Minimal MMseqs2 runs exercising every kwarg covered by
tests/tool_parameters/test_mmseqs2_params.py.
"""

from biopipelines.pipeline import *
from biopipelines import MMseqs2

with Pipeline(project="ToolParameters",
              job="MMseqs2",
              description="MMseqs2 parameter coverage — output_format, mask"):

    Resources(gpu=None, time="2:00:00", memory="16GB")

    seq = Sequence("MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEK", ids="probe")

    # 1: CSV output
    Suffix("csv")
    MMseqs2(sequences=seq, output_format="csv")

    # 2: A3M MSA output
    Suffix("a3m")
    MMseqs2(sequences=seq, output_format="a3m")

    # 3: Masked region
    Suffix("mask")
    MMseqs2(sequences=seq, output_format="a3m", mask="10-20+30-40")
