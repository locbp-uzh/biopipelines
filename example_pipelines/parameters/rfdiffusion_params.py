# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Minimal RFdiffusion runs exercising every kwarg covered by
tests/tool_parameters/test_rfdiffusion_params.py.
"""

from biopipelines.pipeline import *
from biopipelines.rfdiffusion import RFdiffusion

with Pipeline(project="ToolParameters",
              job="RFdiffusion",
              description="RFdiffusion parameter coverage — contigs, num_designs, steps, deterministic, startnum, active_site"):

    Resources(gpu="A100", time="3:00:00", memory="16GB")

    # 1: short contig + small num_designs
    Suffix("short_contig")
    RFdiffusion(contigs="50-100", num_designs=2)

    # 2: explicit step count
    Suffix("steps")
    RFdiffusion(contigs="50-100", num_designs=2, steps=25)

    # 3: deterministic seed via reproducible flag
    Suffix("reproducible")
    RFdiffusion(contigs="50-100", num_designs=2, reproducible=True)

    # 4: design_startnum (output indexing offset)
    Suffix("startnum")
    RFdiffusion(contigs="50-100", num_designs=2, design_startnum=42)

    # 5: smoke combining everything (no input PDB → length-only contig)
    Suffix("all")
    RFdiffusion(
        contigs="80-100",
        num_designs=2,
        steps=30,
        reproducible=True,
        design_startnum=10,
    )

    # 6: inpaint sequence + secondary structure on a motif scaffold
    # Uses 1AKI as the scaffold (lysozyme) so the chain+range tokens resolve.
    Suffix("inpaint")
    lyz = PDB("1AKI", ids="LYZ", convert="pdb")
    RFdiffusion(
        pdb=lyz,
        contigs="A1-50/30-50",
        inpaint="A10-20",
        inpaint_str="A30-40",
        num_designs=2,
        steps=30,
    )
