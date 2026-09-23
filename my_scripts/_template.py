# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Scripting template — copy this to start a step of your own.

Passes every input structure through unchanged and records one number per structure, so it runs against any pdb stream and shows both output kinds at once: a file stream downstream tools can consume, and a table you can filter or plot on. Replace the body of ``execution`` with your protocol and adjust what ``configuration`` declares. Uses only the inputs/outputs proxies; see ``biopipelines/scripting_api.py`` for the full contract.

    Scripting("_template.py", inputs={"structures": some_tool})

``configuration`` runs at config time in the pipeline's own Python and must predict shapes without touching disk. ``execution`` runs later, in the configured env, against real files.
"""


def configuration(inputs):
    from biopipelines.scripting_api import Stream, Table
    ids = inputs["structures"].ids
    return {"structures": Stream("pdb", ids),
            "metrics": Table(columns=["id", "n_residues"])}


def execution(inputs, outputs):
    for structure_id, path in inputs["structures"].iterate():
        lines = open(path).readlines()
        # One CA per residue, summed over every chain in the file — 168L's asymmetric unit has
        # five, so this reports 820 rather than the 164 of one chain.
        n_residues = sum(1 for line in lines if line.startswith("ATOM") and line[12:16].strip() == "CA")

        out_path = outputs["structures"].file(structure_id, f"{structure_id}.pdb")
        with open(out_path, "w") as dst:
            dst.writelines(lines)

        outputs["metrics"].row({"id": structure_id, "n_residues": n_residues})

        # Dropping an id here instead would need outputs.drop(structure_id, cause="...")
        # so the completion check excuses its missing file rather than failing the step.
