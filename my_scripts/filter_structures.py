# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Scripting example — file-stream output.

Copies through only the structures that clear an atom-count threshold, emitting a filtered ``structures`` stream (per-id PDB files). Dropped ids are reported via ``outputs.drop`` so the completion check excuses their absent files instead of failing the step. Uses only the inputs/outputs proxies (no biopipelines import in the script itself).

    Scripting("filter_structures.py", inputs={"structures": some_tool})
"""

MIN_ATOMS = 500

def configuration(inputs):
    from biopipelines.scripting_api import Stream
    return {"structures": Stream("pdb", inputs["structures"].ids)}

def execution(inputs, outputs):
    for structure_id, path in inputs["structures"].iterate():
        lines = open(path).readlines()
        n_atoms = sum(1 for line in lines if line.startswith("ATOM"))
        if n_atoms < MIN_ATOMS:
            outputs.drop(structure_id, cause=f"n_atoms {n_atoms} < {MIN_ATOMS}")
            continue
        out_path = outputs["structures"].file(structure_id, f"{structure_id}.pdb")
        with open(out_path, "w") as dst:
            dst.writelines(lines)
