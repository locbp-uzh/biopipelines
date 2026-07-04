# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Scripting example — table-only output.

Counts atoms per input structure into one standalone table. Uses only the inputs/outputs proxies (no biopipelines import in the script itself). See ``biopipelines/scripting_api.py`` for the full contract.

    Scripting("count_atoms.py", inputs={"structures": some_tool})
"""


def configuration(inputs):
    from biopipelines.scripting_api import Table
    if inputs["structures"].format != "pdb":
        raise ValueError("count_atoms expects a pdb structures stream")
    return {"count": Table(columns=["id", "n_atoms"])}

def execution(inputs, outputs):
    for structure_id, path in inputs["structures"].iterate():
        n_atoms = sum(1 for line in open(path) if line.startswith("ATOM"))
        outputs["count"].row({"id": structure_id, "n_atoms": n_atoms})
