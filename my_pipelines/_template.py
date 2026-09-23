# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Template pipeline.

A minimal starting point: fetch a structure, then run a custom Scripting step
on it. Copy this file, rename it, and replace the body with your protocol.

Submission:
    bp_submit(script="my_pipelines/my_run.py")   # from an agent with bp-mcp registered
    bp-submit my_pipelines/my_run.py             # or from a shell in the repo
"""

from biopipelines.pipeline import *
from biopipelines import PDB, Scripting

with Pipeline(project="Template",
              job="example",
              description="Template pipeline — fetch a structure and run a Scripting step"):

    # Compute resources for the batch. Drop or adjust as needed; CPU-only here.
    Resources(time="1:00:00", memory="8GB", cpus=4)

    # Entities (PDB / Sequence / Ligand / ...) provide typed inputs.
    protein = PDB("168L")

    # A custom step. Scripting("_template.py") resolves against the configured
    # scripts folder (folders.infrastructure.scripts, default my_scripts/), so
    # no absolute path is needed. The inputs dict keys are the names the script
    # reads back via inputs[...].
    custom = Scripting("_template.py",
                       inputs={"structures": protein})

    # custom.streams.structures / custom.tables.scores feed downstream tools
    # exactly like any other tool's output.
