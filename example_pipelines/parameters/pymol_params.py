# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Minimal PyMOL runs exercising every op-type covered by
tests/tool_parameters/test_pymol_params.py.
"""

from biopipelines.pipeline import *
from biopipelines import PyMOL

with Pipeline(project="ToolParameters",
              job="PyMOL",
              description="PyMOL parameter coverage — Load, Color, Show, Align, session"):

    Resources(gpu=None, time="1:00:00", memory="8GB")

    a = PDB("1AKI", ids="LYZ_A")
    b = PDB("2LZM", ids="LYZ_B")

    # 1: bare Load
    Suffix("load")
    PyMOL(PyMOL.Load(structures=a))

    # 2: Color
    Suffix("color")
    PyMOL(
        PyMOL.Load(structures=a),
        PyMOL.Color(color="cyan", structures=a),
    )

    # 3: Show representation
    Suffix("show")
    PyMOL(
        PyMOL.Load(structures=a),
        PyMOL.Show(representation="cartoon"),
    )

    # 4: Align (cealign)
    Suffix("align")
    PyMOL(
        PyMOL.Load(structures=a),
        PyMOL.Load(structures=b),
        PyMOL.Align(method="cealign"),
    )

    # 5: Custom session name + multi-op
    Suffix("multi")
    PyMOL(
        PyMOL.Load(structures=a),
        PyMOL.Load(structures=b),
        PyMOL.Color(color="magenta", structures=a),
        PyMOL.Show(representation="cartoon"),
        PyMOL.Align(method="align"),
        session="custom_view",
    )
