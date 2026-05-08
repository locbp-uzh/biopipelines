# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Minimal RFdiffusion3 runs exercising every kwarg covered by
tests/tool_parameters/test_rfdiffusion3_params.py.
"""

from biopipelines.pipeline import *
from biopipelines.rfdiffusion3 import RFdiffusion3

with Pipeline(project="ToolParameters",
              job="RFdiffusion3",
              description="RFdiffusion3 parameter coverage — length, num_designs, num_models, design_startnum"):

    Resources(gpu="A100", time="2:00:00", memory="16GB")

    # 1: de novo design with length range
    Suffix("length")
    RFdiffusion3(length="80-100", num_designs=2, num_models=1)

    # 2: more designs
    Suffix("num_designs")
    RFdiffusion3(length=80, num_designs=4, num_models=1)

    # 3: multiple models per design
    Suffix("num_models")
    RFdiffusion3(length=80, num_designs=2, num_models=2)

    # 4: explicit start number
    Suffix("startnum")
    RFdiffusion3(length=80, num_designs=2, num_models=1, design_startnum=99)
